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
use anyhow::{Context, Result, bail};
use seqair::sam::stream::SamRecord;
use std::fmt::Write as _;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
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
    ctg: String,
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
#[derive(Clone, Debug)]
struct Lift {
    ctg: String,
    reverse: bool,
    start: i64,
    end: i64,
}
/// Reported alignment or XA candidate, with its score and zero or more lifts.
#[derive(Debug)]
struct Hit {
    ctg: String,
    start: i64,
    end: i64,
    reverse: bool,
    cigar: String,
    nm: i64,
    score: i64,
    hard: bool,
    lifted: Vec<Lift>,
    /// Pre-sort index: zero identifies the SAM-reported hit, others are XA hits.
    original: usize,
    /// Overlap-component ID assigned after sorting by projected contig/start.
    group: usize,
}
impl Hit {
    /// Decode a hit using r985's fixed scoring parameters, not its SAM AS tag.
    ///
    /// `pos` is one-based and positive for valid mapped hits; strand is supplied
    /// separately (the sign of an XA position is not part of the coordinate).
    /// Returns errors from CIGAR decoding rather than JS numeric coercions.
    fn parse(
        ctg: &str,
        pos: i64,
        reverse: bool,
        cigar: &str,
        nm: i64,
        original: usize,
    ) -> Result<Self> {
        let mut mat = 0;
        let mut ins = 0;
        let mut del = 0;
        let mut gaps = 0;
        let mut skip = 0;
        let mut hard = false;
        for (len, op) in cigar_ops(cigar)? {
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
        Ok(Self {
            ctg: ctg.into(),
            start: pos - 1,
            end: pos - 1 + mat + del + skip,
            reverse,
            cigar: cigar.into(),
            nm,
            score: mat - 5 * (nm - del - ins) - 6 * gaps - del - ins,
            hard,
            lifted: Vec::new(),
            original,
            group: 0,
        })
    }
    /// Grouping key: only the first lift is used, matching the reference's TODO.
    /// Unlifted hits keep their own contig and interval, even on an ALT contig.
    fn primary(&self) -> (&str, i64, i64) {
        self.lifted
            .first()
            .map(|x| (x.ctg.as_str(), x.start, x.end))
            .unwrap_or((&self.ctg, self.start, self.end))
    }
    /// Format every lift, in index order, with the reference's trailing semicolon.
    fn lift_tag(&self) -> Option<String> {
        if self.lifted.is_empty() {
            return None;
        }
        let mut s = String::from("lt:Z:");
        for l in &self.lifted {
            write!(
                s,
                "{},{},{},{};",
                l.ctg,
                l.start,
                l.end,
                if l.reverse { '-' } else { '+' }
            )
            .unwrap();
        }
        Some(s)
    }
}
/// Decode CIGAR with seqair; callers apply r985's operation-specific semantics.
///
/// Unlike the reference regex, the parser recognizes `=`, `X` and `P`, but the
/// scoring/projection callers intentionally ignore them. Invalid operations
/// rejected by seqair produce errors instead of being skipped by the JS regex.
fn cigar_ops(cigar: &str) -> Result<Vec<(i64, u8)>> {
    let ops = seqair::sam::stream::parse_cigar(cigar).map_err(anyhow::Error::from_boxed)?;
    Ok(ops
        .iter()
        .map(|op| {
            let packed = op.to_bam_u32();
            ((packed >> 4) as i64, b"MIDNSHP=X"[(packed & 15) as usize])
        })
        .collect())
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
pub struct AltIndex {
    alt: HashMap<String, Index<Alt>>,
    primary: HashMap<String, Index<String>>,
    is_alt: HashSet<String>,
    hla_genes: Vec<String>,
    /// r985 tracks only the chromosome of the last mapped HLA allele it reads.
    hla_chr: Option<String>,
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
        let mut alt: HashMap<String, Vec<Interval<Alt>>> = HashMap::new();
        let mut primary: HashMap<String, Vec<Interval<String>>> = HashMap::new();
        let mut is_alt = HashSet::new();
        let mut hla_genes = Vec::new();
        let mut hla_chr = None;
        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            if line.starts_with('@') || line.split('\t').count() < 11 {
                continue;
            }
            let r = SamRecord::parse(&line)
                .map_err(anyhow::Error::from_boxed)
                .with_context(|| format!("ALT line {}", i + 1))?;
            let f = &r.fields;
            // Classification happens before the unmapped check in r985.
            is_alt.insert(f[0].clone());
            if r.flag & 4 != 0 || f[2] == "*" {
                continue;
            }
            if let Some(gene) = hla_gene(&f[0]) {
                if !hla_genes.iter().any(|x| x == gene) {
                    hla_genes.push(gene.to_owned());
                }
                hla_chr = Some(f[2].clone());
            }
            let ops: Vec<_> = cigar_ops(&f[5])?
                .into_iter()
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
            let reverse = r.flag & 16 != 0;
            let &(len, op) = if reverse {
                ops.last().unwrap()
            } else {
                &ops[0]
            };
            // Only the outermost operation on the relevant strand determines
            // the clip offset; do not sum consecutive clips here.
            let start = if op == b'S' || op == b'H' { len } else { 0 };
            alt.entry(f[0].clone()).or_default().push(Interval {
                start,
                end: start + qaln,
                value: Alt {
                    ctg: f[2].clone(),
                    reverse,
                    // r985 first computes pos = SAM POS - 1, then stores pos - 1.
                    // Example: 1001 + ALT query offset 100 becomes 1099, not 1100.
                    offset: r.pos - 2,
                    query_len: qaln + clip,
                    segments,
                },
            });
            primary.entry(f[2].clone()).or_default().push(Interval {
                start: r.pos - 1,
                end: r.pos - 1 + tlen,
                value: f[0].clone(),
            });
        }
        Ok(Self {
            alt: alt.into_iter().map(|(k, v)| (k, Index::new(v))).collect(),
            primary: primary
                .into_iter()
                .map(|(k, v)| (k, Index::new(v)))
                .collect(),
            is_alt,
            hla_genes,
            hla_chr,
        })
    }
    /// Union the HLA genes whose primary intervals overlap this query span.
    fn collect_hla(&self, ctg: &str, start: i64, end: i64, genes: &mut HashSet<String>) {
        if let Some(idx) = self.primary.get(ctg) {
            for a in idx.overlaps(start, end) {
                if let Some(gene) = hla_gene(&a.value) {
                    genes.insert(gene.to_owned());
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
/// Memory scales with the ALT index plus the current group and its generated hits.
pub struct Processor {
    index: AltIndex,
    ratio: f64,
    buffer: Vec<SamRecord>,
    hla: HashSet<String>,
    fastq: HashMap<String, BufWriter<File>>,
    k8_compat: bool,
}
impl Processor {
    /// Configure MAPQ adjustment and optionally create one FASTQ per HLA gene.
    ///
    /// Supply a finite `ratio`; values above one are capped, and negative values
    /// remain allowed as in r985. A nonempty prefix creates/truncates
    /// `<prefix>.<gene>.fq` even if no reads eventually match that gene.
    /// Correct sequence/quality reversal is enabled by default.
    ///
    /// # Errors
    /// Returns HLA file-creation errors. Files already created are not rolled back.
    pub fn new(index: AltIndex, ratio: f64, prefix: Option<&str>) -> Result<Self> {
        let mut fastq = HashMap::new();
        if let Some(prefix) = prefix.filter(|x| !x.is_empty()) {
            for gene in &index.hla_genes {
                let path = format!("{prefix}.{gene}.fq");
                fastq.insert(
                    gene.clone(),
                    BufWriter::new(File::create(&path).with_context(|| format!("create {path}"))?),
                );
            }
        }
        Ok(Self {
            index,
            ratio: ratio.min(1.0),
            buffer: Vec::new(),
            hla: HashSet::new(),
            fastq,
            k8_compat: false,
        })
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
    /// Whether this record changes QNAME or read-end bits (`FLAG & 0xc0`).
    /// Returns false for an empty buffer; strand and secondary flags do not split
    /// groups. Flush the old group before pushing a record for which this is true.
    pub fn starts_new_group(&self, record: &SamRecord) -> bool {
        self.buffer.first().is_some_and(|s| {
            s.fields[0] != record.fields[0] || (s.flag & 0xc0) != (record.flag & 0xc0)
        })
    }
    /// Emit staged alignments in order, write HLA FASTQs, then clear group state.
    ///
    /// r985 writes the first buffered record's sequence/quality to every HLA gene
    /// collected across the group, keeping that record's orientation and `/end`
    /// suffix. It does not select the best HLA alignment or reorient FASTQ reads.
    ///
    /// # Errors
    /// Propagates callback or FASTQ write failures. Output may be partially
    /// written and the buffer remains populated on failure; abort rather than
    /// retrying, which could duplicate already emitted records.
    pub fn flush(&mut self, mut emit: impl FnMut(&SamRecord) -> Result<()>) -> Result<()> {
        if let Some(first) = self.buffer.first() {
            for record in &self.buffer {
                emit(record)?;
            }
            for gene in &self.hla {
                if let Some(out) = self.fastq.get_mut(gene) {
                    writeln!(
                        out,
                        "@{}/{}{}\n{}\n+\n{}",
                        first.fields[0],
                        (first.flag >> 6) & 3,
                        if first.flag & 16 != 0 { '-' } else { '+' },
                        first.fields[9],
                        first.fields[10]
                    )?;
                }
            }
        }
        self.buffer.clear();
        self.hla.clear();
        Ok(())
    }
    /// Flush HLA file buffers so final I/O failures are observable.
    ///
    /// Call after [`Self::flush`]; this does not emit staged alignments or prevent
    /// future writes. Returns any FASTQ writer's flush error.
    pub fn finish_fastq(&mut self) -> Result<()> {
        for writer in self.fastq.values_mut() {
            writer.flush()?;
        }
        Ok(())
    }
    /// Process one alignment, updating earlier records and staging new XA hits.
    ///
    /// The stages follow `bwa_postalt`: parse the reported/XA hits, lift ALT
    /// positions, group overlapping loci, estimate MAPQ, adjust buffered records,
    /// and append the reported record followed by generated supplementary hits.
    /// Only ALT candidates in the reported hit's group produce new SAM records.
    ///
    /// # Errors
    /// Returns malformed XA/CIGAR errors. This method is not transactional (HLA
    /// membership may already have changed); abort processing on an error.
    ///
    /// # Caller requirements
    /// Flush an existing group if [`Self::starts_new_group`] is true. Records must
    /// follow aligner order so later ALT records can adjust earlier primary hits.
    pub fn push(&mut self, mut record: SamRecord) -> Result<()> {
        if record.flag & 4 != 0 {
            self.buffer.push(record);
            return Ok(());
        }
        let nm = record
            .tag("NM:i:")
            .and_then(|x| x.parse().ok())
            .unwrap_or(0);
        let hit = Hit::parse(
            &record.fields[2],
            record.pos,
            record.flag & 16 != 0,
            &record.fields[5],
            nm,
            0,
        )?;
        if self.index.hla_chr.as_deref() == Some(&hit.ctg) {
            self.index
                .collect_hla(&hit.ctg, hit.start, hit.end, &mut self.hla);
        }
        // HLA collection precedes the hard-clip bypass in r985. Such a record
        // can contribute to FASTQ extraction despite skipping all XA processing.
        if hit.hard {
            self.buffer.push(record);
            return Ok(());
        }
        let span = hit.end - hit.start;
        let mut hits = vec![hit];
        if let Some(xa) = record.tag("XA:Z:") {
            for xa in xa.split(';').filter(|x| !x.is_empty()) {
                let fields: Vec<_> = xa.split(',').collect();
                if fields.len() != 4 || fields[1].len() < 2 {
                    bail!("invalid XA hit: {xa}");
                }
                hits.push(Hit::parse(
                    fields[0],
                    fields[1][1..]
                        .parse()
                        .with_context(|| format!("invalid XA position: {xa}"))?,
                    fields[1].starts_with('-'),
                    fields[2],
                    fields[3].parse()?,
                    hits.len(),
                )?);
            }
        }
        if !hits.iter().any(|h| self.index.is_alt.contains(&h.ctg)) {
            self.buffer.push(record);
            return Ok(());
        }
        // Lift each ALT hit through every overlapping ALT-to-primary mapping.
        // Reflect coordinates across the full ALT query length for reverse maps.
        for h in &mut hits {
            if let Some(idx) = self.index.alt.get(&h.ctg) {
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
                    h.lifted.push(Lift {
                        ctg: a.ctg.clone(),
                        reverse: h.reverse != a.reverse,
                        start: s + a.offset,
                        end: e + a.offset,
                    });
                }
            }
        }
        // Save reported-hit information before sorting; original index zero
        // identifies it afterward. Ambiguous lifts suppress buffered adjustment.
        let n_lifted = hits[0].lifted.len();
        let reported_lift = (n_lifted == 1).then(|| hits[0].lifted[0].clone());
        // Stable sorting matters for first lifts, equal-start hits and output
        // order. Grouping deliberately ignores strand, as the reference does.
        hits.sort_by(|a, b| {
            let a = a.primary();
            let b = b.primary();
            a.0.cmp(b.0).then(a.1.cmp(&b.1))
        });
        // Sweep connected overlap components: extend the running maximum end
        // for transitive overlaps, but split intervals that only touch.
        let mut group = 0;
        let mut end = 0;
        for i in 0..hits.len() {
            let (chr, start, stop) = hits[i].primary();
            if i == 0 || hits[i - 1].primary().0 != chr {
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
        let reported = hits.iter().position(|h| h.original == 0).unwrap();
        let reported_group = hits[reported].group;
        let n_group = hits.iter().filter(|h| h.group == reported_group).count();
        let original_mapq = record.mapq;
        let mut mapq = original_mapq;
        // Recompute only when multiple hits collapse into the reported locus.
        // Each group's best score competes; stable score ties retain group order.
        if n_group > 1 {
            let mut maxima = vec![i64::MIN; group + 1];
            for h in &hits {
                maxima[h.group] = maxima[h.group].max(h.score);
            }
            let mut order: Vec<_> = maxima.into_iter().enumerate().collect();
            order.sort_by_key(|a| std::cmp::Reverse(a.1));
            mapq = if order[0].0 != reported_group {
                0
            } else if order.len() == 1 {
                60
            } else {
                (6 * (order[0].1 - order[1].1)).min(60) as u8
            };
            // A usable ALT map permits MAPQ to rise, while a primary hit can
            // only fall. Membership in is_alt alone is insufficient here.
            mapq = if self.index.alt.contains_key(&record.fields[2]) {
                mapq.max(original_mapq)
            } else {
                mapq.min(original_mapq)
            };
        }
        // HLA membership uses the bounding span of the whole reported group,
        // including gaps between individual hits, and r985's single hla_chr.
        if self.index.hla_chr.as_deref() == Some(hits[reported].primary().0) {
            let mut start = 1 << 30;
            let mut end = 0;
            for h in hits.iter().filter(|h| h.group == reported_group) {
                start = start.min(h.primary().1);
                end = end.max(h.primary().2);
            }
            self.index
                .collect_hla(hits[reported].primary().0, start, end, &mut self.hla);
        }
        // Preserve the reference's directional update: only already-buffered
        // records are adjusted. With no lift, all are treated as nonoverlapping;
        // with multiple lifts, none are adjusted by this record.
        if n_lifted <= 1 {
            for s in &mut self.buffer {
                // r985 uses t[5] (the CURRENT record's CIGAR), not s[5] (the
                // buffered record's). Thus `span` below intentionally belongs to
                // the current hit, even when the earlier record is much shorter.
                let overlap = reported_lift.as_ref().is_some_and(|l| {
                    l.ctg == s.fields[2]
                        && (s.flag & 16 != 0) == l.reverse
                        && s.pos - 1 < l.end
                        && l.start < s.pos - 1 + span
                });
                // Last matching duplicate wins in this scan, unlike the first
                // NM/XA/RG matches. Absent pa defaults to 10, so it is not reduced.
                let mut om = -1;
                let mut pa = 10.0;
                for tag in &s.fields[11..] {
                    if let Some(v) = tag
                        .strip_prefix("om:i:")
                        .and_then(|x| x.parse::<i64>().ok())
                        .filter(|&v| v >= 0)
                    {
                        om = v;
                    } else if let Some(v) = tag.strip_prefix("pa:f:") {
                        pa = v.parse().unwrap_or(f64::NAN);
                    }
                }
                if overlap {
                    // Strictly positive om is restored; om=0 is not restored.
                    if om > 0 {
                        s.mapq = om as u8;
                    }
                    s.mapq = s.mapq.min(mapq);
                } else if pa < self.ratio {
                    if om < 0 {
                        s.fields.push(format!("om:i:{}", s.mapq));
                    }
                    s.mapq = 0;
                }
            }
        }
        // Append, rather than replace, tags to retain r985's duplicate-tag and
        // optional-field ordering behavior on previously processed input.
        record.mapq = mapq;
        if n_group > 1 {
            record.fields.push(format!("om:i:{original_mapq}"));
        }
        if let Some(tag) = hits[reported].lift_tag() {
            record.fields.push(tag);
        }
        let rg = record.tag("RG:Z:").map(|s| format!("RG:Z:{s}"));
        let mut reversed = None;
        let mut generated = Vec::new();
        for (i, h) in hits.iter().enumerate() {
            if h.group != reported_group || i == reported || !self.index.alt.contains_key(&h.ctg) {
                continue;
            }
            let mut fields = record.fields[..11].to_vec();
            fields[2] = h.ctg.clone();
            fields[5] = h.cigar.clone();
            // RNEXT '=' is relative to RNAME; a new contig needs the original
            // reference name spelled out. PNEXT/TLEN remain unchanged as in r985.
            if fields[6] == "=" && h.ctg != record.fields[2] {
                fields[6] = record.fields[2].clone();
            }
            // Keep all original flags and add supplementary (0x800), including
            // any secondary flag already present. Reversal is cached per record.
            let mut flag = record.flag | 0x800;
            if h.reverse != hits[reported].reverse {
                let (seq, qual) = reversed.get_or_insert_with(|| {
                    if self.k8_compat {
                        (record.fields[9].clone(), record.fields[10].clone())
                    } else {
                        reverse_read(&record.fields[9], &record.fields[10])
                    }
                });
                fields[9] = seq.clone();
                fields[10] = qual.clone();
                flag ^= 0x10;
            }
            // Generated hits inherit only the core fields and RG; do not copy
            // XA, AS, SA or other optional tags from the reported alignment.
            fields.push(format!("NM:i:{}", h.nm));
            if let Some(tag) = h.lift_tag() {
                fields.push(tag);
            }
            if let Some(rg) = &rg {
                fields.push(rg.clone());
            }
            generated.push(SamRecord {
                fields,
                flag,
                pos: h.start + 1,
                mapq,
            });
        }
        self.buffer.push(record);
        self.buffer.extend(generated);
        Ok(())
    }
}
/// Perform the intended r985 IUPAC reversal on ASCII SAM sequence and quality.
///
/// Deviation from the broken k8 1.2 runtime: actually reverse both strings.
/// Deviation from the old Bytes table's zero fallback: preserve unknown bytes
/// (including SAM '=') and the missing-value sentinel '*', rather than emitting
/// NULs. Quality is reversed independently, so missing quality stays missing.
fn reverse_read(seq: &str, qual: &str) -> (String, String) {
    let from = b"WSATUGCYRKMBDHVNwsatugcyrkmbdhvn";
    let to = b"WSTAACGRYMKVHDBNwstaacgrymkvhdbn";
    let rs: Vec<u8> = seq
        .bytes()
        .rev()
        .map(|b| {
            from.iter()
                .position(|&x| x == b)
                .map(|i| to[i])
                .unwrap_or(b)
        })
        .collect();
    let rs = if seq == "*" {
        "*".into()
    } else {
        String::from_utf8(rs).unwrap()
    };
    let rq = qual.chars().rev().collect();
    (rs, rq)
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
        let alt = &idx.alt["alt"].entries[0].value;
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
        assert_eq!(
            reverse_read("ACGTRYMKBDHVN", "ABCDEFGHIJKLM"),
            ("NBDHVMKRYACGT".into(), "MLKJIHGFEDCBA".into())
        );
        assert_eq!(reverse_read("*", "*"), ("*".into(), "*".into()));
        assert_eq!(reverse_read("acgt", "*"), ("acgt".into(), "*".into()));
    }
}
