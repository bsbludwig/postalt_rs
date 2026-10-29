// SPDX-License-Identifier: GPL-3.0-only
// Part of the BWA/bwakit postalt Rust port; see NOTICE.md for provenance.

//! ALT projection, overlap grouping and read-end-local MAPQ adjustment.
//!
//! Most coordinates here are zero-based, half-open intervals; `Hit::parse` takes
//! a one-based SAM/XA position. Lifted coordinates deliberately include r985's
//! extra subtraction (`Alt::offset = SAM POS - 2`) and must not be "fixed" as a
//! routine coordinate conversion. Comments marked "r985" identify preserved
//! reference behavior; comments marked "Deviation" describe intentional changes.
//!
//! The reference repeatedly scans ALT CIGARs and uses a binned interval index.
//! This implementation precomputes query segments and prefix-maximum interval
//! indexes. Stable ordering is retained because the first lift and tied groups
//! can determine output order and MAPQ.
//!
//! Most reads in a whole-genome stream never touch an ALT contig. Those records
//! bypass hit parsing entirely (see [`Processor::push`]) and, for BAM input, are
//! not even decoded; only records that can change are converted to SAM text.
use crate::io::Record;
use crate::sam_text::{SamRecord, parse_number};
use anyhow::{Context, Result, bail};
use rustc_hash::FxHashMap;
use seqair::bam::{
    BamHeader, CigarOp, RawBamRecord,
    aux::{AuxValue, find_tag},
};
use seqair::sam::{format, parse_cigar};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    sync::Arc,
};

/// An interval `[start, end)` carrying an ALT mapping or primary-contig label.
#[derive(Clone, Debug)]
struct Interval<T> {
    start: i64,
    end: i64,
    value: T,
}
/// Intervals sorted stably by start, with a monotone maximum-end prefix.
///
/// Replaces r985's `intv_ovlp` bin table. Queries skip impossible prefixes and
/// suffixes with binary searches, then filter the remaining candidates; the
/// candidate scan can still be linear for deeply nested intervals.
#[derive(Clone, Debug)]
struct Index<T> {
    entries: Vec<Interval<T>>,
    max_end: Vec<i64>,
}
impl<T> Index<T> {
    /// Build once, retaining insertion order among equal-start intervals.
    fn new(mut entries: Vec<Interval<T>>) -> Self {
        entries.sort_by_key(|x| x.start);
        let mut max = i64::MIN;
        let max_end = entries
            .iter()
            .map(|x| {
                max = max.max(x.end);
                max
            })
            .collect();
        Self { entries, max_end }
    }
    /// Return overlaps in stable start order; touching boundaries do not overlap.
    fn overlaps(&self, start: i64, end: i64) -> impl Iterator<Item = &Interval<T>> {
        // Starts at/after the query end cannot overlap. Before `lo`, even the
        // largest end seen so far is at/before the query start.
        let hi = self.entries.partition_point(|x| x.start < end);
        let lo = self.max_end[..hi].partition_point(|&x| x <= start);
        self.entries[lo..hi].iter().filter(move |x| x.end > start)
    }
}
/// Query-consuming piece of an ALT CIGAR, prepared for `cigar2pos` lookup.
///
/// `[query, end)` includes clips; `reference` is the reference offset before the
/// operation. Deletions advance that offset without needing a query segment.
#[derive(Clone, Debug)]
struct Segment {
    end: i64,
    query: i64,
    reference: i64,
    op: u8,
}
/// One ALT-to-primary alignment; an ALT contig can have several such mappings.
#[derive(Clone, Debug)]
struct Alt {
    ctg: Vec<u8>,
    reverse: bool,
    /// Primary SAM POS minus two, preserving r985's extra subtraction.
    offset: i64,
    /// Aligned query length plus soft/hard clips, used to reflect reverse hits.
    query_len: i64,
    segments: Vec<Segment>,
}
impl Alt {
    /// Map an ALT query base to a relative primary offset, or `-1` in a clip.
    ///
    /// Equivalent to r985's `cigar2pos`: insertions map to their reference anchor;
    /// out-of-range positions and clips return `-1`. Segment construction omits
    /// reference skips (`N`) from this projection, exactly as the script does.
    fn project(&self, pos: i64) -> i64 {
        let i = self.segments.partition_point(|x| x.end <= pos);
        match self.segments.get(i) {
            Some(s) if pos >= s.query => match s.op {
                b'M' => s.reference + pos - s.query,
                b'I' => s.reference,
                _ => -1,
            },
            _ => -1,
        }
    }
}
/// Projected interval/strand used for grouping and emitted verbatim in `lt:Z`.
///
/// Coordinates retain the reference's offset and clipping quirks; an end that
/// projects into a clip can precede the start. Do not normalize these intervals.
#[derive(Clone, Copy, Debug)]
struct Lift<'a> {
    ctg: &'a [u8],
    reverse: bool,
    start: i64,
    end: i64,
}
/// Reported alignment or XA candidate, with its score and zero or more lifts.
///
/// Borrows contig and CIGAR text from the record and lifted contigs from the
/// ALT index, so evaluating hits allocates nothing per hit.
#[derive(Debug)]
struct Hit<'a> {
    ctg: &'a [u8],
    start: i64,
    end: i64,
    reverse: bool,
    cigar: &'a [u8],
    nm: i64,
    score: i64,
    hard: bool,
    /// This hit's lifts within the shared lift list.
    lifts: std::ops::Range<usize>,
    /// Grouping key `(contig, start, end)`: the first lift, else the hit itself.
    key: (&'a [u8], i64, i64),
    /// Pre-sort index: zero identifies the SAM-reported hit, others are XA hits.
    original: usize,
    /// Overlap-component ID assigned after sorting by projected contig/start.
    group: usize,
}
impl<'a> Hit<'a> {
    /// Decode a hit using r985's fixed scoring parameters, not its SAM AS tag.
    ///
    /// `pos` is one-based and positive for valid mapped hits; strand is supplied
    /// separately (the sign of an XA position is not part of the coordinate).
    /// Returns errors from CIGAR decoding rather than JS numeric coercions.
    /// `ops` is scratch space for the decoded CIGAR.
    fn parse(
        ctg: &'a [u8],
        pos: i64,
        reverse: bool,
        cigar: &'a [u8],
        nm: i64,
        original: usize,
        ops: &mut Vec<CigarOp>,
    ) -> Result<Self> {
        let mut mat = 0;
        let mut ins = 0;
        let mut del = 0;
        let mut gaps = 0;
        let mut skip = 0;
        let mut hard = false;
        decode_cigar(cigar, ops)?;
        for (len, op) in ops.iter().map(op_parts) {
            match op {
                b'M' => mat += len,
                b'I' => {
                    ins += len;
                    gaps += 1;
                }
                b'D' => {
                    del += len;
                    gaps += 1;
                }
                b'N' => skip += len,
                b'H' => hard = true,
                _ => (),
            }
        }
        // r985 clamps NM to the indel length before counting mismatches.
        // With a=1, b=4, gap-open=6 and gap-extension=1, its score is integral;
        // floor(score / a + .499) therefore simplifies to the expression below.
        let nm = nm.max(del + ins);
        let (start, end) = (pos - 1, pos - 1 + mat + del + skip);
        Ok(Self {
            ctg,
            start,
            end,
            reverse,
            cigar,
            nm,
            score: mat - 5 * (nm - del - ins) - 6 * gaps - del - ins,
            hard,
            lifts: 0..0,
            key: (ctg, start, end),
            original,
            group: 0,
        })
    }
}
/// Decode CIGAR with seqair; callers apply r985's operation-specific semantics.
///
/// Unlike the reference regex, the parser recognizes `=`, `X` and `P`, but the
/// scoring/projection callers intentionally ignore them. Invalid operations
/// rejected by seqair produce errors instead of being skipped by the JS regex.
fn decode_cigar(cigar: &[u8], ops: &mut Vec<CigarOp>) -> Result<()> {
    Ok(parse_cigar(cigar, ops)?)
}
/// Decode a BAM record into `out` as SAM text, using `text` as scratch.
fn decode_bam(
    bam: &RawBamRecord,
    header: &BamHeader,
    text: &mut Vec<u8>,
    out: &mut SamRecord,
) -> Result<()> {
    text.clear();
    format::format_sam_record(bam, header, text)?;
    out.parse_into(text)
}
/// Split a decoded operation into its length and SAM operation character.
fn op_parts(op: &CigarOp) -> (i64, u8) {
    let packed = op.to_bam_u32();
    // Op codes past X (B and the reserved ones) are not scored by r985.
    let op = b"MIDNSHP=X"
        .get((packed & 15) as usize)
        .copied()
        .unwrap_or(b'?');
    (i64::from(packed >> 4), op)
}
/// Parse an integer field, returning `None` for anything `i64::from_str` rejects.
fn parse_i64(field: &[u8]) -> Option<i64> {
    parse_number(field).ok()
}
/// Recognize the gene prefix in r985-style names such as `HLA-A*01:01`.
/// The allele suffix only needs to begin with a digit, matching the JS regex.
fn hla_gene(name: &str) -> Option<&str> {
    let (gene, allele) = name.split_once('*')?;
    (gene.starts_with("HLA-")
        && gene.len() > 4
        && !gene.bytes().any(|b| b.is_ascii_whitespace())
        && allele.starts_with(|c: char| c.is_ascii_digit()))
    .then_some(gene)
}

/// Immutable ALT-to-primary mappings and primary intervals for HLA collection.
///
/// Contigs mentioned by unmapped ALT records remain classified as ALT even
/// though they have no usable projection. This distinction affects MAPQ rules.
/// Contig names are compared as bytes, exactly as they appear in SAM text.
pub struct AltIndex {
    /// Every ALT contig; the value holds its usable mappings, if any.
    contigs: FxHashMap<Vec<u8>, Option<Index<Alt>>>,
    /// Primary intervals labelled with the index of their HLA gene, if any.
    primary: FxHashMap<Vec<u8>, Index<Option<u32>>>,
    hla_genes: Vec<String>,
    /// r985 tracks only the chromosome of the last mapped HLA allele it reads.
    hla_chr: Option<Vec<u8>>,
}
impl AltIndex {
    /// Load a plain or gzip-compressed ALT SAM file, detected by its magic bytes.
    ///
    /// # Errors
    /// Returns file/decompression errors or the parsing errors of [`Self::read`].
    pub fn from_path(path: &std::path::Path) -> Result<Self> {
        let mut reader = BufReader::new(
            File::open(path).with_context(|| format!("open ALT file {}", path.display()))?,
        );
        if reader.fill_buf()?.starts_with(&[31, 139]) {
            Self::read(BufReader::new(flate2::read::MultiGzDecoder::new(reader)))
        } else {
            Self::read(reader)
        }
    }

    /// Build both interval indexes from ALT SAM text; no `@SQ` header is required.
    ///
    /// Headers and lines with fewer than eleven fields are skipped like r985.
    /// Unlike the alignment-input loop, this accepts an unterminated final line.
    /// ALT hard clips contribute to full query length and behave like soft clips
    /// during projection, as in the reference's explicit H-to-S conversion.
    ///
    /// # Errors
    /// Returns read/field/CIGAR errors. A mapped ALT record with no recognized
    /// MIDSHN operations is rejected explicitly instead of failing in JS later.
    pub fn read(reader: impl BufRead) -> Result<Self> {
        let mut alt: FxHashMap<Vec<u8>, Vec<Interval<Alt>>> = FxHashMap::default();
        let mut primary: FxHashMap<Vec<u8>, Vec<Interval<Option<u32>>>> = FxHashMap::default();
        let mut contigs = FxHashMap::default();
        let mut hla_genes: Vec<String> = Vec::new();
        let mut hla_chr = None;
        let mut ops = Vec::new();
        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            if line.starts_with('@') || line.split('\t').count() < 11 {
                continue;
            }
            let r =
                SamRecord::parse(line.as_bytes()).with_context(|| format!("ALT line {}", i + 1))?;
            let f: Vec<&str> = line.split('\t').collect();
            // Classification happens before the unmapped check in r985.
            contigs.insert(f[0].as_bytes().to_vec(), None);
            if r.flag() & 4 != 0 || f[2] == "*" {
                continue;
            }
            let mut gene_id = None;
            if let Some(gene) = hla_gene(f[0]) {
                let id = match hla_genes.iter().position(|x| x == gene) {
                    Some(id) => id,
                    None => {
                        hla_genes.push(gene.to_owned());
                        hla_genes.len() - 1
                    }
                };
                gene_id = Some(id as u32);
                hla_chr = Some(f[2].as_bytes().to_vec());
            }
            decode_cigar(f[5].as_bytes(), &mut ops)?;
            let ops: Vec<_> = ops
                .iter()
                .map(op_parts)
                .filter(|(_, op)| b"MIDSHN".contains(op))
                .collect();
            if ops.is_empty() {
                bail!(
                    "ALT line {}: mapped alignment has no MIDSHN CIGAR operations",
                    i + 1
                );
            }
            // Interval spans and projection offsets intentionally use separate
            // counters: N consumes reference span but r985 cigar2pos ignores it.
            let mut qaln = 0;
            let mut clip = 0;
            let mut tlen = 0;
            let mut x = 0;
            let mut y = 0;
            let mut segments = Vec::new();
            for &(len, op) in &ops {
                match op {
                    b'M' => {
                        qaln += len;
                        tlen += len;
                    }
                    b'I' => qaln += len,
                    b'S' | b'H' => clip += len,
                    b'D' | b'N' => tlen += len,
                    _ => (),
                }
                match op {
                    b'M' | b'I' | b'S' | b'H' => {
                        segments.push(Segment {
                            query: y,
                            end: y + len,
                            reference: x,
                            op,
                        });
                        y += len;
                        if op == b'M' {
                            x += len;
                        }
                    }
                    b'D' => x += len,
                    // r985 cigar2pos ignores N, although the interval includes it.
                    _ => (),
                }
            }
            let reverse = r.flag() & 16 != 0;
            let &(len, op) = if reverse {
                ops.last().unwrap()
            } else {
                &ops[0]
            };
            // Only the outermost operation on the relevant strand determines
            // the clip offset; do not sum consecutive clips here.
            let start = if op == b'S' || op == b'H' { len } else { 0 };
            alt.entry(f[0].as_bytes().to_vec())
                .or_default()
                .push(Interval {
                    start,
                    end: start + qaln,
                    value: Alt {
                        ctg: f[2].as_bytes().to_vec(),
                        reverse,
                        // r985 first computes pos = SAM POS - 1, then stores pos - 1.
                        // Example: 1001 + ALT query offset 100 becomes 1099, not 1100.
                        offset: r.pos() - 2,
                        query_len: qaln + clip,
                        segments,
                    },
                });
            primary
                .entry(f[2].as_bytes().to_vec())
                .or_default()
                .push(Interval {
                    start: r.pos() - 1,
                    end: r.pos() - 1 + tlen,
                    value: gene_id,
                });
        }
        for (k, v) in alt {
            contigs.insert(k, Some(Index::new(v)));
        }
        Ok(Self {
            contigs,
            primary: primary
                .into_iter()
                .map(|(k, v)| (k, Index::new(v)))
                .collect(),
            hla_genes,
            hla_chr,
        })
    }
    /// HLA gene names in first-seen order; FASTQ gene buffers use these indexes.
    pub fn hla_genes(&self) -> &[String] {
        &self.hla_genes
    }
    /// Whether the ALT file names this contig (mapped or not).
    fn is_alt(&self, ctg: &[u8]) -> bool {
        self.contigs.contains_key(ctg)
    }
    /// Usable ALT-to-primary mappings of this contig.
    fn alt_maps(&self, ctg: &[u8]) -> Option<&Index<Alt>> {
        self.contigs.get(ctg)?.as_ref()
    }
    /// Whether this is the single chromosome r985 searches for HLA genes.
    fn is_hla_chr(&self, ctg: &[u8]) -> bool {
        self.hla_chr.as_deref() == Some(ctg)
    }
    /// Union the HLA genes whose primary intervals overlap this query span.
    fn collect_hla(&self, ctg: &[u8], start: i64, end: i64, genes: &mut Vec<u32>) {
        if let Some(idx) = self.primary.get(ctg) {
            for gene in idx.overlaps(start, end).filter_map(|a| a.value) {
                if !genes.contains(&gene) {
                    genes.push(gene);
                }
            }
        }
    }
}

/// Stateful postprocessor for one contiguous read-name/read-end group at a time.
///
/// Records cannot be emitted immediately: later ALT hits can lower or restore
/// earlier MAPQs. Call [`Self::starts_new_group`] before each [`Self::push`] and
/// [`Self::flush`] at boundaries and EOF. `push` does not enforce this contract.
/// Memory scales with the ALT index plus the current group and its generated hits;
/// emitted records are kept for reuse, so steady-state processing does not allocate.
pub struct Processor {
    index: Arc<AltIndex>,
    ratio: f64,
    k8_compat: bool,
    /// Whether HLA genes are collected and FASTQ text produced.
    fastq: bool,
    /// Dictionary of BAM input records, needed to decode them.
    bam_header: Option<Arc<BamHeader>>,
    /// Per BAM reference: whether records there need hit processing even
    /// without XA; the last entry covers unplaced (`*`) records.
    tid_needed: Vec<bool>,
    buffer: Vec<Record>,
    /// HLA gene indexes collected for the current group.
    hla: Vec<u32>,
    /// Records generated by the current `push`.
    generated: Vec<Record>,
    spare_sam: Vec<SamRecord>,
    spare_bam: Vec<RawBamRecord>,
    ops: Vec<CigarOp>,
    reversed: (Vec<u8>, Vec<u8>),
    scratch: SamRecord,
    /// SAM text of a decoded BAM record, before it is parsed.
    scratch_text: Vec<u8>,
}
impl Processor {
    /// Configure MAPQ adjustment for SAM input, without HLA FASTQ extraction.
    ///
    /// Supply a finite `ratio`; values above one are capped, and negative values
    /// remain allowed as in r985. Correct sequence/quality reversal is enabled
    /// by default.
    pub fn new(index: Arc<AltIndex>, ratio: f64) -> Self {
        let mut processor = Self {
            index,
            ratio: ratio.min(1.0),
            k8_compat: false,
            fastq: false,
            bam_header: None,
            tid_needed: Vec::new(),
            buffer: Vec::new(),
            hla: Vec::new(),
            generated: Vec::new(),
            spare_sam: Vec::new(),
            spare_bam: Vec::new(),
            ops: Vec::new(),
            reversed: (Vec::new(), Vec::new()),
            scratch: SamRecord::default(),
            scratch_text: Vec::new(),
        };
        processor.classify_references();
        processor
    }
    /// Select installed-k8 compatibility for subsequently generated reverse hits.
    ///
    /// Deviation: the default actually reverse-complements sequence and reverses
    /// quality. k8 1.2 no longer implements `Bytes[i]`, so r985 leaves both strings
    /// unchanged while toggling the strand flag. `enabled = true` reproduces that
    /// output without reproducing the runtime's memory corruption.
    pub fn set_k8_compat(&mut self, enabled: bool) {
        self.k8_compat = enabled;
    }
    /// Collect HLA genes and produce FASTQ text in [`Self::flush`].
    ///
    /// Collection only feeds FASTQ output, so it is skipped when disabled.
    pub fn set_hla_fastq(&mut self, enabled: bool) {
        self.fastq = enabled;
        self.classify_references();
    }
    /// Set the dictionary that [`Record::Bam`] inputs refer to.
    pub fn set_bam_header(&mut self, header: Option<Arc<BamHeader>>) {
        self.bam_header = header;
        self.classify_references();
    }
    /// Whether records on this contig need hit processing even without XA.
    fn contig_needed(&self, ctg: &[u8]) -> bool {
        self.index.is_alt(ctg) || (self.fastq && self.index.is_hla_chr(ctg))
    }
    /// Precompute [`Self::contig_needed`] for each BAM reference index.
    fn classify_references(&mut self) {
        let names: Vec<&str> = match &self.bam_header {
            Some(header) => header.target_names().collect(),
            None => Vec::new(),
        };
        let mut needed: Vec<bool> = names
            .iter()
            .map(|n| self.contig_needed(n.as_bytes()))
            .collect();
        needed.push(self.contig_needed(b"*"));
        self.tid_needed = needed;
    }
    /// An empty SAM record for the caller to fill, reusing earlier allocations.
    pub fn spare_sam(&mut self) -> SamRecord {
        self.spare_sam.pop().unwrap_or_default()
    }
    /// An empty BAM record for the caller to fill, reusing earlier allocations.
    pub fn spare_bam(&mut self) -> RawBamRecord {
        self.spare_bam.pop().unwrap_or_default()
    }
    /// Whether this record changes QNAME or read-end bits (`FLAG & 0xc0`).
    /// Returns false for an empty buffer; strand and secondary flags do not split
    /// groups. Flush the old group before pushing a record for which this is true.
    pub fn starts_new_group(&self, record: &Record) -> bool {
        self.buffer.first().is_some_and(|s| {
            s.qname() != record.qname() || (s.flag() & 0xc0) != (record.flag() & 0xc0)
        })
    }
    /// Emit staged alignments in order, append HLA FASTQ text, then clear the group.
    ///
    /// r985 writes the first buffered record's sequence/quality to every HLA gene
    /// collected across the group, keeping that record's orientation and `/end`
    /// suffix. It does not select the best HLA alignment or reorient FASTQ reads.
    /// `fastq[i]` receives text for gene `i` of [`AltIndex::hla_genes`].
    ///
    /// # Errors
    /// Propagates callback or decoding failures. Output may be partially
    /// emitted; the group is cleared regardless, so abort rather than retry.
    pub fn flush(
        &mut self,
        mut emit: impl FnMut(&Record) -> Result<()>,
        fastq: &mut [Vec<u8>],
    ) -> Result<()> {
        let result = self.emit_group(&mut emit, fastq);
        self.clear();
        result
    }
    fn emit_group(
        &mut self,
        emit: &mut impl FnMut(&Record) -> Result<()>,
        fastq: &mut [Vec<u8>],
    ) -> Result<()> {
        for record in &self.buffer {
            emit(record)?;
        }
        let Some(first) = self.buffer.first() else {
            return Ok(());
        };
        if self.hla.is_empty() {
            return Ok(());
        }
        let first = match first {
            Record::Sam(text) => text,
            Record::Bam(bam) => {
                let header = self
                    .bam_header
                    .as_deref()
                    .context("BAM record without header")?;
                decode_bam(bam, header, &mut self.scratch_text, &mut self.scratch)?;
                &self.scratch
            }
        };
        for &gene in &self.hla {
            if let Some(out) = fastq.get_mut(gene as usize) {
                out.push(b'@');
                out.extend_from_slice(first.qname());
                out.push(b'/');
                out.push(b'0' + ((first.flag() >> 6) & 3) as u8);
                out.push(if first.flag() & 16 != 0 { b'-' } else { b'+' });
                out.push(b'\n');
                out.extend_from_slice(first.field(9));
                out.extend_from_slice(b"\n+\n");
                out.extend_from_slice(first.field(10));
                out.push(b'\n');
            }
        }
        Ok(())
    }
    /// Discard the current group, keeping its records' allocations for reuse.
    pub fn clear(&mut self) {
        for record in self.buffer.drain(..) {
            match record {
                Record::Sam(r) => self.spare_sam.push(r),
                Record::Bam(r) => self.spare_bam.push(r),
            }
        }
        self.hla.clear();
    }

    /// Process one alignment, updating earlier records and staging new XA hits.
    ///
    /// The stages follow `bwa_postalt`: parse the reported/XA hits, lift ALT
    /// positions, group overlapping loci, estimate MAPQ, adjust buffered records,
    /// and append the reported record followed by generated supplementary hits.
    /// Only ALT candidates in the reported hit's group produce new SAM records.
    ///
    /// Records that are unmapped, or have no XA tag and lie on a contig that is
    /// neither ALT nor (with FASTQ enabled) the HLA chromosome, are staged
    /// unchanged without parsing their CIGAR: r985 would also pass them through.
    /// Deviation: their malformed CIGARs are therefore not reported, matching the
    /// script, which never rejects CIGAR text.
    ///
    /// # Errors
    /// Returns malformed XA/CIGAR errors, with the read name as context. This
    /// method is not transactional (HLA membership may already have changed);
    /// abort processing on an error.
    ///
    /// # Caller requirements
    /// Flush an existing group if [`Self::starts_new_group`] is true. Records must
    /// follow aligner order so later ALT records can adjust earlier primary hits.
    pub fn push(&mut self, record: Record) -> Result<()> {
        let mut record = match record {
            Record::Bam(bam) => {
                let needed = bam.tid().map_or(self.tid_needed.last(), |t| {
                    self.tid_needed.get(t.as_usize())
                });
                let has_xa = matches!(
                    find_tag(bam.aux().as_bytes(), *b"XA"),
                    Some(AuxValue::String(_))
                );
                if bam.flags().is_unmapped() || (needed == Some(&false) && !has_xa) {
                    self.buffer.push(Record::Bam(bam));
                    return Ok(());
                }
                let mut text = self.spare_sam();
                let header = self
                    .bam_header
                    .as_deref()
                    .context("BAM record without header")?;
                let decoded = decode_bam(&bam, header, &mut self.scratch_text, &mut text);
                let name = String::from_utf8_lossy(bam.qname()).into_owned();
                self.spare_bam.push(bam);
                decoded.with_context(|| format!("read {name}"))?;
                text
            }
            Record::Sam(text) => text,
        };
        if record.flag() & 4 != 0
            || (!self.contig_needed(record.field(2)) && record.tag(b"XA:Z:").is_none())
        {
            self.buffer.push(Record::Sam(record));
            return Ok(());
        }
        match self.process(&mut record) {
            Ok(()) => {
                self.buffer.push(Record::Sam(record));
                self.buffer.append(&mut self.generated);
                Ok(())
            }
            Err(e) => {
                let name = String::from_utf8_lossy(record.qname()).into_owned();
                self.spare_sam.push(record);
                Err(e.context(format!("read {name}")))
            }
        }
    }

    /// The r985 per-record algorithm for a mapped record that may involve ALTs.
    /// Leaves generated records in `self.generated`; the caller stages `record`.
    fn process(&mut self, record: &mut SamRecord) -> Result<()> {
        let index = Arc::clone(&self.index);
        let index: &AltIndex = &index;
        let rec: &SamRecord = record;
        // NM is bounded like a BAM `i` tag so scoring cannot overflow.
        let nm = rec
            .tag(b"NM:i:")
            .and_then(|v| parse_number::<i32>(v).ok())
            .map_or(0, i64::from);
        let reported = Hit::parse(
            rec.field(2),
            rec.pos(),
            rec.flag() & 16 != 0,
            rec.field(5),
            nm,
            0,
            &mut self.ops,
        )?;
        if self.fastq && index.is_hla_chr(reported.ctg) {
            index.collect_hla(reported.ctg, reported.start, reported.end, &mut self.hla);
        }
        // HLA collection precedes the hard-clip bypass in r985. Such a record
        // can contribute to FASTQ extraction despite skipping all XA processing.
        if reported.hard {
            return Ok(());
        }
        let span = reported.end - reported.start;
        let mut hits = vec![reported];
        if let Some(xa) = rec.tag(b"XA:Z:") {
            parse_xa(xa, &mut hits, &mut self.ops)?;
        }
        if !hits.iter().any(|h| index.is_alt(h.ctg)) {
            return Ok(());
        }
        let lifts = lift(index, &mut hits);
        // Save reported-hit information before sorting; original index zero
        // identifies it afterward. Ambiguous lifts suppress buffered adjustment.
        let n_lifted = hits[0].lifts.len();
        let reported_lift = (n_lifted == 1).then(|| lifts[hits[0].lifts.start]);
        let groups = group_hits(&mut hits);
        let reported = hits.iter().position(|h| h.original == 0).unwrap();
        let reported_group = hits[reported].group;
        let n_group = hits.iter().filter(|h| h.group == reported_group).count();
        let original_mapq = rec.mapq();
        let mut mapq = original_mapq;
        // Recompute only when multiple hits collapse into the reported locus.
        if n_group > 1 {
            mapq = group_mapq(&hits, groups, reported_group);
            // A usable ALT map permits MAPQ to rise, while a primary hit can
            // only fall. Membership in the ALT list alone is insufficient here.
            mapq = if index.alt_maps(rec.field(2)).is_some() {
                mapq.max(original_mapq)
            } else {
                mapq.min(original_mapq)
            };
        }
        // HLA membership uses the bounding span of the whole reported group,
        // including gaps between individual hits, and r985's single hla_chr.
        if self.fastq && index.is_hla_chr(hits[reported].key.0) {
            let mut start = 1 << 30;
            let mut end = 0;
            for h in hits.iter().filter(|h| h.group == reported_group) {
                start = start.min(h.key.1);
                end = end.max(h.key.2);
            }
            index.collect_hla(hits[reported].key.0, start, end, &mut self.hla);
        }
        // Preserve the reference's directional update: only already-buffered
        // records are adjusted. With no lift, all are treated as nonoverlapping;
        // with multiple lifts, none are adjusted by this record.
        if n_lifted <= 1 {
            self.adjust_buffered(reported_lift, span, mapq)?;
        }
        self.generate(rec, &hits, &lifts, reported, mapq, index)?;
        let mut lt = Vec::new();
        write_lifts(&mut lt, &lifts[hits[reported].lifts.clone()]);
        // Append, rather than replace, tags to retain r985's duplicate-tag and
        // optional-field ordering behavior on previously processed input.
        record.set_mapq(mapq);
        if n_group > 1 {
            record.push_field_fmt(format_args!("om:i:{original_mapq}"));
        }
        if !lt.is_empty() {
            record.push_field(&lt);
        }
        Ok(())
    }

    /// Restore or lower MAPQs of records already staged for this read end.
    ///
    /// BAM records staged unchanged are decoded first so they can be edited.
    fn adjust_buffered(&mut self, lift: Option<Lift<'_>>, span: i64, mapq: u8) -> Result<()> {
        for i in 0..self.buffer.len() {
            if let Record::Bam(bam) = &self.buffer[i] {
                let mut text = self.spare_sam.pop().unwrap_or_default();
                let header = self
                    .bam_header
                    .as_deref()
                    .context("BAM record without header")?;
                decode_bam(bam, header, &mut self.scratch_text, &mut text)?;
                if let Record::Bam(bam) = std::mem::replace(&mut self.buffer[i], Record::Sam(text))
                {
                    self.spare_bam.push(bam);
                }
            }
            let Record::Sam(s) = &mut self.buffer[i] else {
                unreachable!("decoded above")
            };
            // r985 uses t[5] (the CURRENT record's CIGAR), not s[5] (the
            // buffered record's). Thus `span` below intentionally belongs to
            // the current hit, even when the earlier record is much shorter.
            let overlap = lift.is_some_and(|l| {
                l.ctg == s.field(2)
                    && (s.flag() & 16 != 0) == l.reverse
                    && s.pos() - 1 < l.end
                    && l.start < s.pos() - 1 + span
            });
            // Last matching duplicate wins in this scan, unlike the first
            // NM/XA/RG matches. Absent pa defaults to 10, so it is not reduced.
            let mut om = -1;
            let mut pa = 10.0;
            for tag in s.tags() {
                if let Some(v) = tag
                    .strip_prefix(b"om:i:")
                    .and_then(parse_i64)
                    .filter(|&v| v >= 0)
                {
                    om = v;
                } else if let Some(v) = tag.strip_prefix(b"pa:f:") {
                    pa = parse_number(v).unwrap_or(f64::NAN);
                }
            }
            if overlap {
                // Strictly positive om is restored; om=0 is not restored. The
                // minimum is taken before narrowing, so om above 255 cannot wrap.
                let restored = if om > 0 { om } else { i64::from(s.mapq()) };
                s.set_mapq(restored.min(i64::from(mapq)) as u8);
            } else if pa < self.ratio {
                if om < 0 {
                    let old = s.mapq();
                    s.push_field_fmt(format_args!("om:i:{old}"));
                }
                s.set_mapq(0);
            }
        }
        Ok(())
    }

    /// Stage supplementary records for ALT hits in the reported group.
    fn generate(
        &mut self,
        rec: &SamRecord,
        hits: &[Hit<'_>],
        lifts: &[Lift<'_>],
        reported: usize,
        mapq: u8,
        index: &AltIndex,
    ) -> Result<()> {
        let group = hits[reported].group;
        let rg = rec.tag(b"RG:Z:");
        let mut reversed = false;
        for (i, h) in hits.iter().enumerate() {
            if h.group != group || i == reported || index.alt_maps(h.ctg).is_none() {
                continue;
            }
            // Keep all original flags and add supplementary (0x800), including
            // any secondary flag already present. Reversal is cached per record.
            let mut flag = rec.flag() | 0x800;
            let opposite = h.reverse != hits[reported].reverse;
            if opposite {
                flag ^= 0x10;
                if !self.k8_compat && !reversed {
                    let (seq, qual) = &mut self.reversed;
                    reverse_read(rec.field(9), rec.field(10), seq, qual);
                    reversed = true;
                }
            }
            let (seq, qual) = if opposite && !self.k8_compat {
                (&self.reversed.0[..], &self.reversed.1[..])
            } else {
                (rec.field(9), rec.field(10))
            };
            // RNEXT '=' is relative to RNAME; a new contig needs the original
            // reference name spelled out. PNEXT/TLEN remain unchanged as in r985.
            let rnext = if rec.field(6) == b"=" && h.ctg != rec.field(2) {
                rec.field(2)
            } else {
                rec.field(6)
            };
            let pos = h.start + 1;
            let mut out = self.spare_sam.pop().unwrap_or_default();
            out.assemble(flag, pos, mapq, |t| {
                let mut n = itoa::Buffer::new();
                // Generated hits inherit only the core fields and RG; do not copy
                // XA, AS, SA or other optional tags from the reported alignment.
                for (i, field) in [
                    rec.qname(),
                    n.format(flag).as_bytes(),
                    h.ctg,
                    itoa::Buffer::new().format(pos).as_bytes(),
                    itoa::Buffer::new().format(mapq).as_bytes(),
                    h.cigar,
                    rnext,
                    rec.field(7),
                    rec.field(8),
                    seq,
                    qual,
                ]
                .into_iter()
                .enumerate()
                {
                    if i > 0 {
                        t.push(b'\t');
                    }
                    t.extend_from_slice(field);
                }
                t.extend_from_slice(b"\tNM:i:");
                t.extend_from_slice(itoa::Buffer::new().format(h.nm).as_bytes());
                if !h.lifts.is_empty() {
                    t.push(b'\t');
                    write_lifts(t, &lifts[h.lifts.clone()]);
                }
                if let Some(rg) = rg {
                    t.extend_from_slice(b"\tRG:Z:");
                    t.extend_from_slice(rg);
                }
            })?;
            self.generated.push(Record::Sam(out));
        }
        Ok(())
    }
}

/// Append XA hits (`ctg,±pos,CIGAR,NM;...`), skipping empty entries.
fn parse_xa<'a>(xa: &'a [u8], hits: &mut Vec<Hit<'a>>, ops: &mut Vec<CigarOp>) -> Result<()> {
    for entry in xa.split(|&b| b == b';').filter(|x| !x.is_empty()) {
        let text = || String::from_utf8_lossy(entry).into_owned();
        let mut fields = [&b""[..]; 4];
        let mut n = 0;
        for field in entry.split(|&b| b == b',') {
            if let Some(slot) = fields.get_mut(n) {
                *slot = field;
            }
            n += 1;
        }
        let [ctg, pos, cigar, nm] = fields;
        if n != 4 || pos.len() < 2 {
            bail!("invalid XA hit: {}", text());
        }
        // Bounded so that scoring arithmetic cannot overflow; r985's floats
        // never wrapped, they just lost precision far beyond any real value.
        let (sign, digits) = pos.split_first().unwrap_or((&b'+', b""));
        if !matches!(sign, b'+' | b'-') || !digits.iter().all(u8::is_ascii_digit) {
            bail!("invalid XA position: {}", text());
        }
        let start = parse_number::<i32>(digits)
            .ok()
            .filter(|&p| p > 0)
            .with_context(|| format!("invalid XA position: {}", text()))?;
        let nm = parse_number::<u32>(nm).with_context(|| format!("invalid XA NM: {}", text()))?;
        let (start, nm) = (i64::from(start), i64::from(nm));
        let hit = Hit::parse(ctg, start, *sign == b'-', cigar, nm, hits.len(), ops)?;
        hits.push(hit);
    }
    Ok(())
}

/// Lift each ALT hit through every overlapping ALT-to-primary mapping and set
/// its grouping key. Reflect coordinates across the full ALT query length for
/// reverse maps. Returns the lifts that `Hit::lifts` ranges refer to.
fn lift<'a>(index: &'a AltIndex, hits: &mut [Hit<'a>]) -> Vec<Lift<'a>> {
    let mut lifts = Vec::new();
    for h in hits.iter_mut() {
        let first = lifts.len();
        if let Some(idx) = index.alt_maps(h.ctg) {
            for a in idx.overlaps(h.start, h.end) {
                let a = &a.value;
                let (s, e) = if a.reverse {
                    (a.query_len - h.end, a.query_len - h.start - 1)
                } else {
                    (h.start, h.end - 1)
                };
                let s = a.project(s);
                let e = a.project(e) + 1;
                // r985 adds one BEFORE checking the projected end. A clipped
                // end returns -1, becomes 0, and is accepted. Keep that quirk;
                // checking the unincremented result would alter lt tags/groups.
                if s < 0 || e < 0 {
                    continue;
                }
                lifts.push(Lift {
                    ctg: &a.ctg,
                    reverse: h.reverse != a.reverse,
                    start: s + a.offset,
                    end: e + a.offset,
                });
            }
        }
        h.lifts = first..lifts.len();
        // Grouping key: only the first lift is used, matching the reference's
        // TODO. Unlifted hits keep their own contig and interval, even on an ALT.
        if let Some(l) = lifts.get(first) {
            h.key = (l.ctg, l.start, l.end);
        }
    }
    lifts
}

/// Sort hits by grouping key and number connected overlap components.
/// Returns the number of groups.
fn group_hits(hits: &mut [Hit<'_>]) -> usize {
    // Stable sorting matters for first lifts, equal-start hits and output
    // order. Grouping deliberately ignores strand, as the reference does.
    hits.sort_by(|a, b| a.key.0.cmp(b.key.0).then(a.key.1.cmp(&b.key.1)));
    // Sweep connected overlap components: extend the running maximum end
    // for transitive overlaps, but split intervals that only touch.
    let mut group = 0;
    let mut end = 0;
    for i in 0..hits.len() {
        let (chr, start, stop) = hits[i].key;
        if i == 0 || hits[i - 1].key.0 != chr {
            if i > 0 {
                group += 1;
            }
            end = 0;
        } else if start >= end {
            group += 1;
        }
        end = end.max(stop);
        hits[i].group = group;
    }
    group + 1
}

/// MAPQ from the reported group's best score versus the best other group.
/// Each group's best score competes; stable score ties retain group order.
fn group_mapq(hits: &[Hit<'_>], groups: usize, reported_group: usize) -> u8 {
    let mut maxima = vec![i64::MIN; groups];
    for h in hits {
        maxima[h.group] = maxima[h.group].max(h.score);
    }
    let mut order: Vec<_> = maxima.into_iter().enumerate().collect();
    order.sort_by_key(|a| std::cmp::Reverse(a.1));
    if order[0].0 != reported_group {
        0
    } else if order.len() == 1 {
        60
    } else {
        // Scores are bounded i64s; saturate rather than wrap on extreme input.
        let diff = order[0].1.saturating_sub(order[1].1).saturating_mul(6);
        u8::try_from(diff.clamp(0, 60)).unwrap_or(60)
    }
}

/// Append an `lt:Z` tag for these lifts, in index order, with the reference's
/// trailing semicolons. Appends nothing for no lifts.
fn write_lifts(out: &mut Vec<u8>, lifts: &[Lift<'_>]) {
    if lifts.is_empty() {
        return;
    }
    out.extend_from_slice(b"lt:Z:");
    let mut n = itoa::Buffer::new();
    for l in lifts {
        out.extend_from_slice(l.ctg);
        out.push(b',');
        out.extend_from_slice(n.format(l.start).as_bytes());
        out.push(b',');
        out.extend_from_slice(n.format(l.end).as_bytes());
        out.push(b',');
        out.push(if l.reverse { b'-' } else { b'+' });
        out.push(b';');
    }
}

/// IUPAC complement for every byte; other bytes map to themselves.
static COMPLEMENT: [u8; 256] = {
    let from = b"WSATUGCYRKMBDHVNwsatugcyrkmbdhvn";
    let to = b"WSTAACGRYMKVHDBNwstaacgrymkvhdbn";
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = i as u8;
        i += 1;
    }
    let mut j = 0;
    while j < from.len() {
        table[from[j] as usize] = to[j];
        j += 1;
    }
    table
};

/// Perform the intended r985 IUPAC reversal on ASCII SAM sequence and quality.
///
/// Deviation from the broken k8 1.2 runtime: actually reverse both strings.
/// Deviation from the old Bytes table's zero fallback: preserve unknown bytes
/// (including SAM '=') and the missing-value sentinel '*', rather than emitting
/// NULs. Quality is reversed independently, so missing quality stays missing.
fn reverse_read(seq: &[u8], qual: &[u8], rev_seq: &mut Vec<u8>, rev_qual: &mut Vec<u8>) {
    rev_seq.clear();
    if seq == b"*" {
        rev_seq.push(b'*');
    } else {
        rev_seq.extend(seq.iter().rev().map(|&b| COMPLEMENT[b as usize]));
    }
    rev_qual.clear();
    rev_qual.extend(qual.iter().rev());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interval_queries_match_exhaustive_search() {
        let mut seed = 985u64;
        let entries = (0..500)
            .map(|id| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                let start = (seed >> 32) as i64 % 10000;
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                Interval {
                    start,
                    end: start + 1 + (seed >> 32) as i64 % 2000,
                    value: id,
                }
            })
            .collect();
        let idx = Index::new(entries);
        for start in (-10..12000).step_by(13) {
            for len in [1, 20, 1000, 20000] {
                let expected: Vec<_> = idx
                    .entries
                    .iter()
                    .filter(|x| x.start < start + len && x.end > start)
                    .map(|x| x.value)
                    .collect();
                let actual: Vec<_> = idx.overlaps(start, start + len).map(|x| x.value).collect();
                assert_eq!(actual, expected);
            }
        }
    }

    #[test]
    fn projection_boundaries_include_insertions_and_ignore_reference_skips() {
        let text = b"alt\t0\tchr1\t1001\t60\t3S4M2I3M2D4M2N2M2H\t*\t0\t0\t*\t*\n";
        let idx = AltIndex::read(&text[..]).unwrap();
        let alt = &idx.alt_maps(b"alt").unwrap().entries[0].value;
        let expected = [
            -1, -1, -1, 0, 1, 2, 3, 4, 4, 4, 5, 6, 9, 10, 11, 12, 13, 14, -1, -1,
        ];
        for (pos, reference) in expected.into_iter().enumerate() {
            assert_eq!(alt.project(pos as i64), reference, "query {pos}");
        }
        assert_eq!(alt.project(-1), -1);
        assert_eq!(alt.project(20), -1);
        assert_eq!(alt.offset, 999);
    }

    #[test]
    fn corrected_reverse_preserves_missing_fields_and_ambiguity_codes() {
        let reverse = |seq: &str, qual: &str| {
            let (mut s, mut q) = (Vec::new(), Vec::new());
            reverse_read(seq.as_bytes(), qual.as_bytes(), &mut s, &mut q);
            (String::from_utf8(s).unwrap(), String::from_utf8(q).unwrap())
        };
        assert_eq!(
            reverse("ACGTRYMKBDHVN", "ABCDEFGHIJKLM"),
            ("NBDHVMKRYACGT".into(), "MLKJIHGFEDCBA".into())
        );
        assert_eq!(reverse("*", "*"), ("*".into(), "*".into()));
        assert_eq!(reverse("acgt", "*"), ("acgt".into(), "*".into()));
    }

    #[test]
    fn restored_om_above_255_does_not_wrap() {
        // Overlapping earlier record with om:i:300: r985 restores 300, then
        // caps it at the new MAPQ (60 here), rather than wrapping 300 to 44.
        let alt = b"altA\t0\tchr1\t1001\t60\t1000M\t*\t0\t0\t*\t*\n";
        let mut p = Processor::new(Arc::new(AltIndex::read(&alt[..]).unwrap()), 1.0);
        let earlier = b"q\t0\tchr1\t1100\t0\t100M\t*\t0\t0\t*\t*\tom:i:300\tpa:f:0.9";
        let current = b"q\t2048\taltA\t101\t60\t100M\t*\t0\t0\t*\t*\tNM:i:0";
        for line in [&earlier[..], &current[..]] {
            p.push(Record::Sam(SamRecord::parse(line).unwrap()))
                .unwrap();
        }
        let mut mapqs = Vec::new();
        p.flush(
            |r| {
                mapqs.push(match r {
                    Record::Sam(s) => s.mapq(),
                    Record::Bam(b) => b.mapq().as_byte(),
                });
                Ok(())
            },
            &mut [],
        )
        .unwrap();
        assert_eq!(mapqs, [60, 60]);
    }
}
