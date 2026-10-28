// SPDX-License-Identifier: GPL-3.0-only
// Part of the BWA/bwakit postalt Rust port; see NOTICE.md for provenance.

//! Sequential SAM/BAM transport with bounded buffers and explicit finalization.
//!
//! Deviation from `bwa-postalt.js`: BAM input/output, compression workers and
//! explicit output files are supported in addition to plain/gzip SAM streams.
//! SAM text is preserved by the seqair adapter; BAM serialization can normalize
//! base case, reference aliases and optional-field numeric spelling.
//!
//! Input is read in [`Batch`]es of whole read-end groups so that groups can be
//! processed independently; [`Serializer`] turns processed records back into
//! output bytes. BAM records stay in wire form until something needs their text,
//! and unchanged BAM records are copied to BAM output verbatim.
//!
//! Format detection replays bytes consumed while sniffing, so pipes work without
//! seeking. Neither input nor output sorts records or adds program/header lines.
//! Call [`Output::finish`] to observe compression and buffered-write errors.
use crate::sam_text::SamRecord;
use anyhow::{Context, Result, bail};
use seqair::{
    bam::{
        BamHeader, RawBamRecord,
        bgzf::BgzfStream,
        stream::BamStreamReader,
        writer::{BamWriter, BamWriterBuilder},
    },
    sam::{SamRecordEncoder, format},
};
use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    num::NonZeroU8,
    path::Path,
    sync::Arc,
};

/// One alignment: SAM text, or a BAM record still in wire form.
#[derive(Clone, Debug)]
pub enum Record {
    /// Parsed SAM text (SAM input, generated hits, or decoded BAM).
    Sam(SamRecord),
    /// Validated BAM record that has not needed decoding.
    Bam(RawBamRecord),
}
impl Record {
    /// QNAME bytes.
    pub fn qname(&self) -> &[u8] {
        match self {
            Self::Sam(r) => r.qname(),
            Self::Bam(r) => r.qname(),
        }
    }
    /// FLAG bits.
    pub fn flag(&self) -> u16 {
        match self {
            Self::Sam(r) => r.flag(),
            Self::Bam(r) => r.flags().raw(),
        }
    }
}

/// Serialized form of a [`Batch`] or of output bytes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    /// Newline-terminated SAM lines.
    Sam,
    /// BAM records, each a little-endian `block_size` and its payload.
    Bam,
}

/// Consecutive input records forming whole read-end groups.
///
/// A group (same QNAME and `FLAG & 0xc0`) never spans two batches, so batches
/// can be processed independently. Records have been framed but not parsed.
#[derive(Debug)]
pub struct Batch {
    /// Zero-based position in the input; output must follow this order.
    pub index: u64,
    /// One-based SAM line number (headers included) or BAM record number of
    /// the first record, for error messages.
    pub first: u64,
    /// How `data` is framed.
    pub format: Format,
    /// The records' bytes.
    pub data: Vec<u8>,
}
impl Batch {
    /// Iterate the framed records' bytes: SAM lines without their newline, or
    /// BAM payloads without their size prefix.
    pub fn records(&self) -> impl Iterator<Item = &[u8]> {
        let mut rest = &self.data[..];
        let format = self.format;
        std::iter::from_fn(move || {
            let (record, tail) = match format {
                Format::Sam => {
                    if rest.is_empty() {
                        return None;
                    }
                    let end = memchr::memchr(b'\n', rest).unwrap_or(rest.len());
                    (&rest[..end], rest.get(end + 1..).unwrap_or_default())
                }
                // Input framed every record with a size it checked.
                Format::Bam => {
                    let (size, rest) = rest.split_first_chunk::<4>()?;
                    rest.split_at_checked(u32::from_le_bytes(*size) as usize)?
                }
            };
            rest = tail;
            Some(record)
        })
    }
}

/// Grouping key of a framed record, or `None` if it cannot be read.
fn group_key(format: Format, record: &[u8]) -> Option<(&[u8], u16)> {
    match format {
        Format::Sam => {
            let mut fields = record.splitn(3, |&b| b == b'\t');
            let name = fields.next()?;
            let flag = std::str::from_utf8(fields.next()?)
                .ok()?
                .parse::<u16>()
                .ok()?;
            Some((name, flag & 0xc0))
        }
        Format::Bam => {
            let payload = record.get(4..)?;
            let name_len = usize::from(*payload.get(8)?);
            let flag = u16::from_le_bytes(payload.get(14..16)?.try_into().ok()?);
            Some((payload.get(32..32 + name_len.checked_sub(1)?)?, flag & 0xc0))
        }
    }
}

/// Buffered file/stdin/replayed stream shared by both input formats.
type InputStream = BufReader<Box<dyn Read + Send>>;

/// Alignment input whose header has already been consumed.
///
/// Use the constructors, then read either whole groups with
/// [`Self::next_batch`] or single records with [`Self::next_record`]; stop at
/// the first `None`. No alignment index or coordinate sorting is required.
pub struct Input {
    source: Source,
    /// Input record count so far: SAM line number (headers included) or BAM
    /// record number of the most recently framed record.
    number: u64,
    /// Set once a terminating condition has been reached.
    done: bool,
    /// Record read past the end of the previous batch, and its number.
    carry: Vec<u8>,
    carry_number: u64,
    /// Read error deferred until the records before it have been returned.
    error: Option<anyhow::Error>,
    batches: u64,
}
/// Format-specific input state.
enum Source {
    /// Plain or decompressed SAM text, with original leading header lines.
    Sam {
        /// Buffered text after any format-detection replay/decompression.
        reader: InputStream,
        /// First nonheader line read during header detection, if any.
        pending: Option<Vec<u8>>,
        /// Leading header text with line endings normalized to LF.
        header: Vec<u8>,
    },
    /// Sequential BGZF BAM decoded through seqair.
    Bam {
        /// Record framing over the decompressed stream, possibly with
        /// read-ahead workers.
        reader: Box<BamStreamReader<InputStream>>,
        /// Reference dictionary of the BAM header.
        header: Arc<BamHeader>,
        /// The header as SAM text, as `samtools view -h` prints it: cut at its
        /// first NUL and newline-terminated.
        text: Vec<u8>,
    },
}
impl Input {
    /// Open a path, or stdin when `path` is `None` or `-`, and read its header.
    ///
    /// BAM input is inflated on rayon's global pool, up to `threads` batches of
    /// blocks at once (zero inflates on the reading thread); SAM input ignores
    /// it.
    ///
    /// # Errors
    /// Returns open errors or the detection/header errors of [`Self::from_reader`].
    pub fn open(path: Option<&Path>, threads: usize) -> Result<Self> {
        let reader: Box<dyn Read + Send> = match path.filter(|p| *p != Path::new("-")) {
            Some(p) => Box::new(File::open(p).with_context(|| format!("open {}", p.display()))?),
            None => Box::new(io::stdin()),
        };
        Self::from_reader(reader, threads)
    }
    /// Detect SAM, gzip/BGZF SAM or BGZF BAM by content and consume its header.
    ///
    /// Partial reads, including one byte at a time, are supported. A gzip sniff
    /// is replayed to the actual decoder; extensions are never consulted. Empty
    /// and headerless SAM streams are accepted. CRAM and raw, non-BGZF BAM are
    /// not supported.
    ///
    /// # Errors
    /// Returns stream, decompression, UTF-8 header text or BAM-header decoding
    /// errors, or a failure to start decompression workers. Reading can block
    /// until sufficient bytes arrive to identify the format.
    pub fn from_reader(mut reader: Box<dyn Read + Send>, threads: usize) -> Result<Self> {
        // Read the magic bytes explicitly: fill_buf() may return only one byte
        // on a pipe. Treating that as plain SAM would misclassify a valid BAM.
        let prefix = read_prefix(&mut reader, 4)?;
        let gzip = prefix.starts_with(&[31, 139]);
        if prefix == b"CRAM" {
            bail!("CRAM input is not supported; convert it with `samtools view -b`");
        }
        if prefix == b"BAM\x01" {
            bail!("uncompressed BAM input is not supported; compress it with `samtools view -b`");
        }
        let reader: Box<dyn Read + Send> = Box::new(io::Cursor::new(prefix).chain(reader));
        let mut reader = BufReader::with_capacity(256 * 1024, reader);
        let mut is_bam = false;
        if gzip {
            // Inspect decompressed magic without assuming a seekable stream or
            // relying on the filename. Replay consumed compressed bytes to seqair.
            let mut consumed = Vec::new();
            let mut tee = Tee {
                reader: &mut reader,
                copy: &mut consumed,
            };
            let magic = read_prefix(&mut flate2::read::MultiGzDecoder::new(&mut tee), 4)?;
            is_bam = magic == b"BAM\x01";
            let replay = io::Cursor::new(consumed).chain(reader);
            let source: Box<dyn Read + Send> = if is_bam {
                Box::new(replay)
            } else {
                Box::new(flate2::read::MultiGzDecoder::new(replay))
            };
            reader = BufReader::with_capacity(256 * 1024, source);
        }
        let mut number = 0;
        let source = if is_bam {
            // Decompression shares rayon's global pool; `threads` bounds how
            // many batches it inflates at once (255 at most).
            let builder = BgzfStream::builder(reader);
            let stream = match NonZeroU8::new(u8::try_from(threads).unwrap_or(u8::MAX)) {
                Some(n) => builder.threads(n),
                None => builder,
            }
            .build();
            let reader = BamStreamReader::from_bgzf(stream).context("read BAM header")?;
            let header = Arc::new(reader.header().clone());
            let mut text = Vec::new();
            format::format_sam_header(&header, &mut text).context("BAM header text")?;
            Source::Bam {
                reader: Box::new(reader),
                header,
                text,
            }
        } else {
            let mut header = Vec::new();
            let mut line = Vec::new();
            loop {
                line.clear();
                if reader.read_until(b'\n', &mut line)? == 0 {
                    break;
                }
                number += 1;
                // Preserve the supplied k8 readline()>0 loop: a final line
                // without a delimiter is omitted, even if it is a header.
                if !line.ends_with(b"\n") {
                    warn_stop(number, "its last line has no newline");
                    line.clear();
                    break;
                }
                if !line.starts_with(b"@") {
                    break;
                }
                // Any bytes are kept (k8 does not decode); BAM output checks
                // UTF-8 when it needs a dictionary.
                let end = line
                    .iter()
                    .rposition(|&b| b != b'\r' && b != b'\n')
                    .map_or(0, |i| i + 1);
                header.extend_from_slice(&line[..end]);
                header.push(b'\n');
            }
            // Header detection has already read the first alignment; retain it
            // rather than requiring a seek back to its start.
            Source::Sam {
                reader,
                pending: (!line.is_empty()).then_some(line),
                header,
            }
        };
        Ok(Self {
            source,
            number,
            done: false,
            carry: Vec::new(),
            carry_number: 0,
            error: None,
            batches: 0,
        })
    }
    /// Return the leading header without adding @PG records or changing sort order.
    /// SAM CRLF endings are normalized to LF; headerless input returns no bytes.
    /// A BAM header is returned as `samtools view -h` prints it: up to its first
    /// NUL, newline-terminated, with `@SQ` lines generated if the text has none.
    pub fn header_text(&self) -> &[u8] {
        match &self.source {
            Source::Sam { header, .. } => header,
            Source::Bam { text, .. } => text,
        }
    }
    /// The BAM reference dictionary that [`Record::Bam`] values refer to.
    pub fn bam_header(&self) -> Option<&Arc<BamHeader>> {
        match &self.source {
            Source::Sam { .. } => None,
            Source::Bam { header, .. } => Some(header),
        }
    }
    /// Framing of the records this input yields.
    pub fn format(&self) -> Format {
        match self.source {
            Source::Sam { .. } => Format::Sam,
            Source::Bam { .. } => Format::Bam,
        }
    }

    /// Append the next framed record to `buf` and return whether there was one.
    ///
    /// r985 compatibility: a blank SAM line terminates processing and an
    /// unterminated final SAM line is discarded; nothing is read afterwards.
    /// SAM lines keep their newline; BAM records keep their size prefix.
    ///
    /// # Errors
    /// Returns read/decompression errors and BAM truncation or size errors.
    /// BAM payloads must be 32 bytes to 2 MiB, matching seqair's record-size limit.
    fn read_framed(&mut self, buf: &mut Vec<u8>) -> Result<bool> {
        if self.done {
            return Ok(false);
        }
        let start = buf.len();
        match &mut self.source {
            Source::Sam {
                reader, pending, ..
            } => {
                if let Some(line) = pending.take() {
                    buf.extend_from_slice(&line);
                } else if reader.read_until(b'\n', buf)? == 0 {
                    self.done = true;
                    return Ok(false);
                } else {
                    self.number += 1;
                }
                // Keep both of the reference runtime's stopping conditions;
                // do not silently turn a blank line into a skipped record.
                let line = &buf[start..];
                let unterminated = !line.ends_with(b"\n");
                if unterminated || line.iter().all(|&b| b == b'\r' || b == b'\n') {
                    let more = !unterminated && !reader.fill_buf()?.is_empty();
                    if unterminated {
                        warn_stop(self.number, "its last line has no newline");
                    } else if more {
                        warn_stop(self.number, "it is blank");
                    }
                    buf.truncate(start);
                    self.done = true;
                    return Ok(false);
                }
            }
            Source::Bam { reader, .. } => {
                // The size prefix is filled in once the payload is read.
                buf.extend_from_slice(&[0; 4]);
                if !reader.read_payload_into(buf)? {
                    buf.truncate(start);
                    self.done = true;
                    return Ok(false);
                }
                // Payloads are at most seqair's 2 MiB record limit.
                let size = u32::try_from(buf.len() - start - 4).context("BAM record size")?;
                buf[start..start + 4].copy_from_slice(&size.to_le_bytes());
                self.number += 1;
            }
        }
        Ok(true)
    }

    /// Read whole read-end groups totalling at least `target` bytes, if available.
    ///
    /// Batches end only where QNAME or `FLAG & 0xc0` changes, so a single large
    /// group can exceed `target`. Records are framed but not validated here;
    /// unparsable grouping fields simply allow a batch boundary. After a read
    /// error, the records before it are returned first and the error on the next
    /// call.
    ///
    /// # Errors
    /// Returns read errors as described for single records.
    pub fn next_batch(&mut self, target: usize) -> Result<Option<Batch>> {
        self.next_batch_into(target, Vec::new())
    }
    /// As [`Self::next_batch`], reusing `buffer`'s allocation for the batch data.
    ///
    /// # Errors
    /// As [`Self::next_batch`].
    pub fn next_batch_into(&mut self, target: usize, mut buffer: Vec<u8>) -> Result<Option<Batch>> {
        if let Some(e) = self.error.take() {
            return Err(e);
        }
        let format = self.format();
        let mut data = {
            buffer.clear();
            buffer.extend_from_slice(&self.carry);
            self.carry.clear();
            buffer
        };
        let mut first = self.carry_number;
        let mut previous = (!data.is_empty()).then_some(0);
        loop {
            let start = data.len();
            match self.read_framed(&mut data) {
                Ok(true) => {}
                Ok(false) => break,
                Err(e) => {
                    data.truncate(start);
                    self.done = true;
                    if data.is_empty() {
                        return Err(e);
                    }
                    self.error = Some(e);
                    break;
                }
            }
            let Some(prev) = previous else {
                first = self.number;
                previous = Some(start);
                continue;
            };
            // A mid-stream header line belongs with the group it interrupts:
            // r985 prints it while that group is still buffered.
            if start >= target && !is_header_line(format, &data[start..]) {
                let key = group_key(format, &data[start..]);
                if key.is_none() || key != group_key(format, &data[prev..start]) {
                    self.carry.extend_from_slice(&data[start..]);
                    self.carry_number = self.number;
                    data.truncate(start);
                    break;
                }
            }
            previous = Some(start);
        }
        if data.is_empty() {
            return Ok(None);
        }
        self.batches += 1;
        Ok(Some(Batch {
            index: self.batches - 1,
            first,
            format,
            data,
        }))
    }

    /// Read and decode one record as SAM text, for simple sequential use.
    ///
    /// Do not interleave with [`Self::next_batch`].
    ///
    /// # Errors
    /// Returns read errors and record parsing/validation errors.
    pub fn next_record(&mut self) -> Result<Option<SamRecord>> {
        let mut buf = Vec::new();
        if !self.read_framed(&mut buf)? {
            return Ok(None);
        }
        let number = self.number;
        match &self.source {
            Source::Sam { .. } => SamRecord::parse(&buf)
                .map(Some)
                .with_context(|| format!("SAM line {number}")),
            Source::Bam { header, .. } => {
                let record = RawBamRecord::from_payload(buf.split_off(4), header)
                    .with_context(|| format!("BAM record {number}"))?;
                let mut text = Vec::new();
                format::format_sam_record(&record, header, &mut text)
                    .with_context(|| format!("BAM record {number}"))?;
                SamRecord::parse(&text)
                    .map(Some)
                    .with_context(|| format!("BAM record {number}"))
            }
        }
    }
}

/// Whether a framed record is a SAM header line (`@...`).
pub fn is_header_line(format: Format, record: &[u8]) -> bool {
    format == Format::Sam && record.first() == Some(&b'@')
}

/// Read up to `n` bytes, fewer only at the end of input.
fn read_prefix(reader: &mut impl Read, n: usize) -> io::Result<Vec<u8>> {
    let mut prefix = Vec::with_capacity(n);
    reader.by_ref().take(n as u64).read_to_end(&mut prefix)?;
    Ok(prefix)
}

/// Report that SAM input stops early for an r985 compatibility reason.
fn warn_stop(line: u64, why: &str) {
    eprintln!(
        "postalt_rs: warning: input stops at line {line} because {why} (as bwa-postalt.js does)"
    );
}

/// Converts records to output bytes, keeping unchanged BAM records verbatim.
///
/// BAM-to-BAM copying needs the output dictionary to list the same reference
/// names in the same order as the input's, which holds when the output header is
/// the input's own text. Otherwise records go through SAM text and are re-encoded.
pub struct Serializer {
    input: Option<Arc<BamHeader>>,
    output: Option<Arc<BamHeader>>,
    verbatim: bool,
    text: Vec<u8>,
    encoder: SamRecordEncoder,
    encoded: RawBamRecord,
}
impl Serializer {
    /// Serialize SAM, or BAM for the `output` dictionary (see
    /// [`Output::bam_header`]), decoding BAM input records against `input`.
    pub fn new(input: Option<Arc<BamHeader>>, output: Option<Arc<BamHeader>>) -> Self {
        let verbatim = match (&input, &output) {
            (Some(i), Some(o)) => i.target_names().eq(o.target_names()),
            _ => false,
        };
        Self {
            input,
            output,
            verbatim,
            text: Vec::new(),
            encoder: SamRecordEncoder::default(),
            encoded: RawBamRecord::default(),
        }
    }
    /// Append `record` to `out` in the output format: a SAM line, or a BAM
    /// `block_size` and payload.
    ///
    /// # Errors
    /// Returns BAM encoding errors (for example references missing from the
    /// output dictionary) and BAM-to-SAM decoding errors (malformed aux data).
    /// `out` is unchanged on error.
    pub fn write(&mut self, record: &Record, out: &mut Vec<u8>) -> Result<()> {
        self.write_inner(record, out)
            .with_context(|| format!("read {}", String::from_utf8_lossy(record.qname())))
    }
    fn write_inner(&mut self, record: &Record, out: &mut Vec<u8>) -> Result<()> {
        let text = match record {
            Record::Bam(bam) if self.verbatim => {
                append_payload(bam, out);
                return Ok(());
            }
            Record::Bam(bam) => {
                let header = self
                    .input
                    .as_deref()
                    .context("BAM record without a header")?;
                self.text.clear();
                format::format_sam_record(bam, header, &mut self.text)?;
                &self.text[..]
            }
            Record::Sam(text) => text.as_bytes(),
        };
        match &self.output {
            Some(header) => {
                self.encoder.encode(text, header, &mut self.encoded)?;
                append_payload(&self.encoded, out);
            }
            None => {
                out.extend_from_slice(text);
                out.push(b'\n');
            }
        }
        Ok(())
    }
}

impl Serializer {
    /// Append a SAM header line found among the alignments (without its line
    /// ending) to SAM output.
    ///
    /// # Errors
    /// BAM output has nowhere to put it.
    pub fn write_header_line(&self, line: &[u8], out: &mut Vec<u8>) -> Result<()> {
        if self.output.is_some() {
            bail!("a header line after the first alignment cannot be written to BAM");
        }
        let end = line
            .iter()
            .rposition(|&b| b != b'\r' && b != b'\n')
            .map_or(0, |i| i + 1);
        out.extend_from_slice(&line[..end]);
        out.push(b'\n');
        Ok(())
    }
}

/// Append a record's `block_size` and payload.
fn append_payload(record: &RawBamRecord, out: &mut Vec<u8>) {
    let payload = record.as_payload();
    // A valid record is at most seqair's 2 MiB limit.
    out.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    out.extend_from_slice(payload);
}

/// Buffer small field writes before sending them to a file or stdout.
type OutputStream = BufWriter<Box<dyn Write>>;
/// SAM or BAM output with the header written eagerly at construction.
///
/// This sink preserves emission order and does not build an index: supplementary
/// alignments generally make the output unsuitable for coordinate indexing until
/// it is sorted. Always call [`Self::finish`] to complete the stream.
pub enum Output {
    /// Text SAM with original optional-field spelling and order.
    Sam(OutputStream),
    /// Compressed BAM with a dictionary used to resolve emitted reference names.
    Bam {
        /// Seqair serializer and single- or multiple-worker BGZF compressor.
        writer: Box<BamWriter<OutputStream>>,
        /// Same dictionary used to initialize `writer`; do not replace midstream.
        header: Arc<BamHeader>,
        /// A record being written, reused.
        record: RawBamRecord,
    },
}
impl Output {
    /// Create/truncate a file or use stdout (`None`/`-`), then emit the header.
    ///
    /// `bam` selects format independently of the filename. SAM accepts empty
    /// header text; BAM needs a valid UTF-8 `@SQ` dictionary. `threads` (at
    /// least one) and `compression` (0..=12) configure BAM compression only.
    /// Callers must ensure the output path does not alias an input they are
    /// still reading.
    ///
    /// # Errors
    /// Returns invalid-header, create, compression-setup or header-write errors.
    /// BAM header validation happens before creating/truncating the destination.
    pub fn create(
        path: Option<&Path>,
        bam: bool,
        text: &[u8],
        threads: usize,
        compression: i32,
    ) -> Result<Self> {
        // Validate header before creating/truncating any output.
        let header = if bam {
            let text = std::str::from_utf8(text).context("BAM output requires a UTF-8 header")?;
            Some(
                BamHeader::from_sam_text(text)
                    .context("BAM output requires a SAM @SQ dictionary")?,
            )
        } else {
            None
        };
        let out: Box<dyn Write> = match path.filter(|p| *p != Path::new("-")) {
            Some(p) => {
                Box::new(File::create(p).with_context(|| format!("create {}", p.display()))?)
            }
            None => Box::new(io::stdout()),
        };
        let mut out = BufWriter::with_capacity(256 * 1024, out);
        if let Some(header) = header {
            let writer = BamWriterBuilder::to_writer(out, &header)
                .compression_level(compression)
                .compression_threads(threads)
                .build()?;
            Ok(Self::Bam {
                writer: Box::new(writer),
                header: Arc::new(header),
                record: RawBamRecord::default(),
            })
        } else {
            out.write_all(text)?;
            Ok(Self::Sam(out))
        }
    }
    /// The BAM dictionary, or `None` for SAM output.
    pub fn bam_header(&self) -> Option<&Arc<BamHeader>> {
        match self {
            Self::Sam(_) => None,
            Self::Bam { header, .. } => Some(header),
        }
    }
    /// Write bytes produced by a [`Serializer`] for this output.
    ///
    /// # Errors
    /// Propagates I/O, record validation and compression errors. Writes can be
    /// partial; callers should abort on failure.
    pub fn write_serialized(&mut self, data: &[u8]) -> Result<()> {
        match self {
            Self::Sam(out) => out.write_all(data)?,
            Self::Bam {
                writer,
                header,
                record,
            } => {
                let mut rest = data;
                while let Some((size, tail)) = rest.split_first_chunk::<4>() {
                    let (payload, tail) = tail
                        .split_at_checked(u32::from_le_bytes(*size) as usize)
                        .context("serialized BAM record is truncated")?;
                    // Reuse the record's allocation for each payload.
                    let mut buffer = std::mem::take(record).into_payload();
                    buffer.clear();
                    buffer.extend_from_slice(payload);
                    *record = RawBamRecord::from_payload(buffer, header)?;
                    writer.write_raw(record)?;
                    rest = tail;
                }
                if !rest.is_empty() {
                    bail!("serialized BAM record size is truncated");
                }
            }
        }
        Ok(())
    }
    /// Serialize and write one SAM record.
    ///
    /// # Errors
    /// Propagates I/O and serialization errors, including references absent from
    /// the BAM dictionary. Writes can be partial; callers should abort on failure.
    pub fn write(&mut self, record: &SamRecord) -> Result<()> {
        match self {
            Self::Sam(out) => record.write_sam(out)?,
            Self::Bam {
                writer,
                header,
                record: raw,
            } => {
                SamRecordEncoder::default().encode(record.as_bytes(), header, raw)?;
                writer.write_raw(raw)?;
            }
        }
        Ok(())
    }
    /// Finish compression, emit the BAM EOF block when applicable, and flush.
    ///
    /// # Errors
    /// Returns pending compressor or buffered I/O failures that dropping a writer
    /// would hide. Consumes the sink to prevent writing after finalization.
    pub fn finish(self) -> Result<()> {
        match self {
            Self::Sam(mut out) => out.flush()?,
            Self::Bam { writer, .. } => {
                let (mut out, _) = writer.finish()?;
                out.flush()?;
            }
        }
        Ok(())
    }
}

/// Copy compressed bytes read during format detection so stdin can be replayed.
struct Tee<'a, R> {
    reader: &'a mut R,
    copy: &'a mut Vec<u8>,
}
impl<R: Read> Read for Tee<'_, R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.reader.read(buf)?;
        self.copy.extend_from_slice(&buf[..n]);
        Ok(n)
    }
}
