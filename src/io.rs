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
use anyhow::{Context, Result, bail};
use seqair::{
    bam::{
        BamHeader,
        bgzf::BgzfReader,
        writer::{BamWriter, BamWriterBuilder},
    },
    sam::stream::{BamEncoder, BamRecord, BgzfStream, MAX_BAM_RECORD, SamRecord, read_bam_header},
};
use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    path::Path,
    sync::Arc,
};

/// One alignment: SAM text, or a BAM record still in wire form.
#[derive(Clone, Debug)]
pub enum Record {
    /// Parsed SAM text (SAM input, generated hits, or decoded BAM).
    Sam(SamRecord),
    /// Validated BAM payload that has not needed decoding.
    Bam(BamRecord),
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
            Self::Bam(r) => r.flag(),
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
            if rest.is_empty() {
                return None;
            }
            let (record, tail) = match format {
                Format::Sam => {
                    let end = memchr::memchr(b'\n', rest).unwrap_or(rest.len());
                    (&rest[..end], rest.get(end + 1..).unwrap_or_default())
                }
                Format::Bam => {
                    // Input framed every record with a validated size.
                    let size = i32::from_le_bytes(rest[..4].try_into().unwrap()) as usize;
                    (&rest[4..4 + size], &rest[4 + size..])
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

/// Adapt a forward-only reader to seqair's offset-tracking `Read + Seek` API.
///
/// Only position queries and seeks to the current offset succeed. Actual
/// repositioning returns `Unsupported`; this is sufficient for sequential BGZF
/// decoding and avoids materializing stdin merely to satisfy a `Seek` bound.
pub struct ForwardReader<R> {
    inner: R,
    offset: u64,
}
impl<R: Read> Read for ForwardReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.offset += n as u64;
        Ok(n)
    }
}
impl<R> Seek for ForwardReader<R> {
    fn seek(&mut self, from: SeekFrom) -> io::Result<u64> {
        match from {
            SeekFrom::Current(0) => Ok(self.offset),
            SeekFrom::Start(n) if n == self.offset => Ok(n),
            _ => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "input stream cannot seek",
            )),
        }
    }
}
/// Buffered file/stdin/replayed stream shared by both input formats.
type InputStream = BufReader<Box<dyn Read + Send>>;
/// Decompressed BAM records after the header.
type BamStream = BgzfStream<BufReader<ForwardReader<InputStream>>>;

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
        header: String,
    },
    /// Sequential BGZF BAM decoded through seqair.
    Bam {
        /// Decompressed record stream, possibly with read-ahead workers.
        reader: Box<BamStream>,
        /// Reference dictionary and original header text from the BAM header.
        header: Arc<BamHeader>,
    },
}
impl Input {
    /// Open a path, or stdin when `path` is `None` or `-`, and read its header.
    ///
    /// `threads` read-ahead workers decompress BAM input (zero decompresses on
    /// the reading thread); SAM input ignores it.
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
        // Read both magic bytes explicitly: fill_buf() may return only one byte
        // on a pipe. Treating that as plain SAM would misclassify a valid BAM.
        let mut prefix = [0u8; 2];
        let mut used = 0;
        while used < prefix.len() {
            let n = reader.read(&mut prefix[used..])?;
            if n == 0 {
                break;
            }
            used += n;
        }
        let gzip = used == 2 && prefix == [31, 139];
        let reader: Box<dyn Read + Send> =
            Box::new(std::io::Cursor::new(prefix[..used].to_vec()).chain(reader));
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
            let mut decoder = flate2::read::MultiGzDecoder::new(&mut tee);
            let mut magic = [0u8; 4];
            let mut n = 0;
            while n < 4 {
                let got = decoder.read(&mut magic[n..])?;
                if got == 0 {
                    break;
                }
                n += got;
            }
            if n == 4 && magic == *b"BAM\x01" {
                is_bam = true;
            }
            let replay = std::io::Cursor::new(consumed).chain(reader);
            let source: Box<dyn Read + Send> = if is_bam {
                Box::new(replay)
            } else {
                Box::new(flate2::read::MultiGzDecoder::new(replay))
            };
            reader = BufReader::with_capacity(256 * 1024, source);
        }
        let mut number = 0;
        let source = if is_bam {
            let mut reader = BgzfReader::from_reader(ForwardReader {
                inner: reader,
                offset: 0,
            });
            let header = read_bam_header(&mut reader)
                .map_err(anyhow::Error::from_boxed)
                .context("read BAM header")?;
            // Hand the rest of the stream, including any already-decompressed
            // remainder of the header's last block, to the read-ahead decoder.
            let (rest, inner, ended) = reader.into_parts();
            Source::Bam {
                reader: Box::new(BgzfStream::new(inner, rest, ended, threads)?),
                header: Arc::new(header),
            }
        } else {
            let mut header = String::new();
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
                    line.clear();
                    break;
                }
                if !line.starts_with(b"@") {
                    break;
                }
                let text = std::str::from_utf8(&line).context("SAM header is not UTF-8")?;
                header.push_str(text.trim_end_matches(['\r', '\n']));
                header.push('\n');
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
    /// SAM CRLF endings are normalized to LF; headerless input returns an empty str.
    pub fn header_text(&self) -> &str {
        match &self.source {
            Source::Sam { header, .. } => header,
            Source::Bam { header, .. } => header.header_text(),
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
                if !line.ends_with(b"\n") || line.iter().all(|&b| b == b'\r' || b == b'\n') {
                    buf.truncate(start);
                    self.done = true;
                    return Ok(false);
                }
            }
            Source::Bam { reader, .. } => {
                let mut size = [0u8; 4];
                // EOF before the first size byte is normal. EOF partway through
                // a size prefix or payload must instead report truncation.
                let n = loop {
                    match reader.read(&mut size) {
                        Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                        result => break result?,
                    }
                };
                if n == 0 {
                    self.done = true;
                    return Ok(false);
                }
                reader
                    .read_exact(&mut size[n..])
                    .context("truncated BAM record size")?;
                let len = i32::from_le_bytes(size);
                if !(32..=MAX_BAM_RECORD as i32).contains(&len) {
                    bail!("invalid BAM record size {len}");
                }
                buf.extend_from_slice(&size);
                buf.resize(start + 4 + len as usize, 0);
                reader
                    .read_exact(&mut buf[start + 4..])
                    .context("truncated BAM record")?;
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
            if start >= target {
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
                .map_err(anyhow::Error::from_boxed)
                .with_context(|| format!("SAM line {number}")),
            Source::Bam { header, .. } => {
                let mut record = BamRecord::default();
                *record.buffer_mut() = buf.split_off(4);
                let mut text = SamRecord::default();
                record
                    .validate(header)
                    .and_then(|()| record.to_sam(header, &mut text))
                    .map_err(anyhow::Error::from_boxed)
                    .with_context(|| format!("BAM record {number}"))?;
                Ok(Some(text))
            }
        }
    }
}

/// Converts records to output bytes, keeping unchanged BAM records verbatim.
///
/// BAM-to-BAM copying needs the output dictionary to list the same reference
/// names in the same order as the input's, which holds when the output header is
/// the input's own text. Otherwise records go through SAM text and are re-encoded.
pub struct Serializer {
    format: Format,
    input: Option<Arc<BamHeader>>,
    output: Option<Arc<BamHeader>>,
    verbatim: bool,
    text: SamRecord,
    encoder: BamEncoder,
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
            format: if output.is_some() {
                Format::Bam
            } else {
                Format::Sam
            },
            input,
            output,
            verbatim,
            text: SamRecord::default(),
            encoder: BamEncoder::default(),
        }
    }
    /// Append `record` to `out` in the output format.
    ///
    /// # Errors
    /// Returns BAM encoding errors (for example references missing from the
    /// output dictionary). `out` is unchanged on error.
    pub fn write(&mut self, record: &Record, out: &mut Vec<u8>) -> Result<()> {
        let text = match record {
            Record::Bam(bam) if self.verbatim => {
                bam.append_bam(out);
                return Ok(());
            }
            Record::Bam(bam) => {
                let header = self
                    .input
                    .as_deref()
                    .context("BAM record without a header")?;
                bam.to_sam(header, &mut self.text)
                    .map_err(anyhow::Error::from_boxed)?;
                &self.text
            }
            Record::Sam(text) => text,
        };
        match (self.format, &self.output) {
            (Format::Bam, Some(header)) => text
                .append_bam(header, &mut self.encoder, out)
                .map_err(anyhow::Error::from_boxed),
            _ => {
                text.append_sam(out);
                Ok(())
            }
        }
    }
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
    },
}
impl Output {
    /// Create/truncate a file or use stdout (`None`/`-`), then emit the header.
    ///
    /// `bam` selects format independently of the filename. SAM accepts empty
    /// header text; BAM needs a valid `@SQ` dictionary. `threads` (at least one)
    /// and `compression` (0..=12) configure BAM compression only. Callers must
    /// ensure the output path does not alias an input they are still reading.
    ///
    /// # Errors
    /// Returns invalid-header, create, compression-setup or header-write errors.
    /// BAM header validation happens before creating/truncating the destination.
    pub fn create(
        path: Option<&Path>,
        bam: bool,
        text: &str,
        threads: usize,
        compression: i32,
    ) -> Result<Self> {
        // Validate header before creating/truncating any output.
        let header = if bam {
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
            })
        } else {
            out.write_all(text.as_bytes())?;
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
    /// Propagates I/O, framing and compression errors. Writes can be partial;
    /// callers should abort on failure.
    pub fn write_serialized(&mut self, data: &[u8]) -> Result<()> {
        match self {
            Self::Sam(out) => out.write_all(data)?,
            Self::Bam { writer, header } => writer
                .write_bam_records(data, header)
                .map_err(anyhow::Error::from_boxed)?,
        }
        Ok(())
    }
    /// Serialize and write one record using its authoritative flag/position/MAPQ.
    ///
    /// # Errors
    /// Propagates I/O and serialization errors, including references absent from
    /// the BAM dictionary. Writes can be partial; callers should abort on failure.
    pub fn write(&mut self, record: &SamRecord) -> Result<()> {
        let mut data = Vec::new();
        match self {
            Self::Sam(out) => record.write_sam(out)?,
            Self::Bam { header, .. } => {
                record
                    .append_bam(header, &mut BamEncoder::default(), &mut data)
                    .map_err(anyhow::Error::from_boxed)?;
                self.write_serialized(&data)?;
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
