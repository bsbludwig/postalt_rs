// SPDX-License-Identifier: GPL-3.0-only
// Part of the BWA/bwakit postalt Rust port; see NOTICE.md for provenance.

//! Sequential SAM/BAM transport with bounded buffers and explicit finalization.
//!
//! Deviation from `bwa-postalt.js`: BAM input/output, compression workers and
//! explicit output files are supported in addition to plain/gzip SAM streams.
//! SAM text is preserved by the seqair adapter; BAM serialization can normalize
//! base case, reference aliases and optional-field numeric spelling.
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
    sam::stream::{SamRecord, read_bam_header},
};
use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    path::Path,
};

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
type InputStream = BufReader<Box<dyn Read>>;
/// Alignment input whose header has already been consumed.
///
/// Use the constructors and stop on the first `None` from [`Self::next_record`].
/// SAM keeps one pending line after header detection; BAM reuses a raw-record
/// buffer. No alignment index or coordinate sorting is required.
pub enum Input {
    /// Plain or decompressed SAM text, with original leading header lines.
    Sam {
        /// Buffered text after any format-detection replay/decompression.
        reader: InputStream,
        /// Reused line storage; initially the first nonheader line, if present.
        line: String,
        /// Whether `line` must be processed before reading another line.
        pending: bool,
        /// One-based source line count for error messages, including headers.
        number: u64,
        /// Leading header text with line endings normalized to LF.
        header: String,
    },
    /// Sequential BGZF BAM decoded through seqair.
    Bam {
        /// BGZF reader over a possibly nonseekable compressed byte stream.
        reader: Box<BgzfReader<ForwardReader<InputStream>>>,
        /// Reused BAM record payload, excluding its four-byte size prefix.
        raw: Vec<u8>,
        /// Reference dictionary and original header text from the BAM header.
        header: BamHeader,
    },
}
impl Input {
    /// Open a path, or stdin when `path` is `None` or `-`, and read its header.
    ///
    /// # Errors
    /// Returns open errors or the detection/header errors of [`Self::from_reader`].
    pub fn open(path: Option<&Path>) -> Result<Self> {
        let reader: Box<dyn Read> = match path.filter(|p| *p != Path::new("-")) {
            Some(p) => Box::new(File::open(p).with_context(|| format!("open {}", p.display()))?),
            None => Box::new(io::stdin()),
        };
        Self::from_reader(reader)
    }
    /// Detect SAM, gzip/BGZF SAM or BGZF BAM by content and consume its header.
    ///
    /// Partial reads, including one byte at a time, are supported. A gzip sniff
    /// is replayed to the actual decoder; extensions are never consulted. Empty
    /// and headerless SAM streams are accepted. CRAM and raw, non-BGZF BAM are
    /// not supported.
    ///
    /// # Errors
    /// Returns stream, decompression, UTF-8 text or BAM-header decoding errors.
    /// Reading can block until sufficient bytes arrive to identify the format.
    pub fn from_reader(mut reader: Box<dyn Read>) -> Result<Self> {
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
        let reader: Box<dyn Read> =
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
            let source: Box<dyn Read> = if is_bam {
                Box::new(replay)
            } else {
                Box::new(flate2::read::MultiGzDecoder::new(replay))
            };
            reader = BufReader::with_capacity(256 * 1024, source);
        }
        if is_bam {
            let mut reader = BgzfReader::from_reader(ForwardReader {
                inner: reader,
                offset: 0,
            });
            let header = read_bam_header(&mut reader)
                .map_err(anyhow::Error::from_boxed)
                .context("read BAM header")?;
            Ok(Self::Bam {
                reader: Box::new(reader),
                raw: Vec::new(),
                header,
            })
        } else {
            let mut header = String::new();
            let mut line = String::new();
            let mut number = 0;
            loop {
                line.clear();
                if reader.read_line(&mut line)? == 0 {
                    break;
                }
                number += 1;
                // Preserve the supplied k8 readline()>0 loop: a final line
                // without a delimiter is omitted, even if it is a header.
                if !line.ends_with('\n') {
                    line.clear();
                    break;
                }
                if !line.starts_with('@') {
                    break;
                }
                header.push_str(line.trim_end_matches(['\r', '\n']));
                header.push('\n');
            }
            // Header detection has already read the first alignment; retain it
            // rather than requiring a seek back to its start.
            let pending = !line.is_empty();
            Ok(Self::Sam {
                reader,
                line,
                pending,
                number,
                header,
            })
        }
    }
    /// Return the leading header without adding @PG records or changing sort order.
    /// SAM CRLF endings are normalized to LF; headerless input returns an empty str.
    pub fn header_text(&self) -> &str {
        match self {
            Self::Sam { header, .. } => header,
            Self::Bam { header, .. } => header.header_text(),
        }
    }
    /// Decode the next alignment, returning `None` when processing should stop.
    ///
    /// r985 compatibility: a blank SAM line terminates processing and an
    /// unterminated final SAM line is discarded. Stop at the first `None`; the
    /// reader does not latch a separate terminal state after a blank line.
    /// Headers must precede records; interspersed headers are not supported.
    ///
    /// # Errors
    /// Returns record parsing, read/decompression or truncation errors. BAM
    /// payloads must be 32 bytes to 2 MiB, matching seqair's record-size limit.
    pub fn next_record(&mut self) -> Result<Option<SamRecord>> {
        match self {
            Self::Sam {
                reader,
                line,
                pending,
                number,
                ..
            } => {
                if !*pending {
                    line.clear();
                    if reader.read_line(line)? == 0 {
                        return Ok(None);
                    }
                    *number += 1;
                }
                *pending = false;
                // Keep both of the reference runtime's stopping conditions;
                // do not silently turn a blank line into a skipped record.
                if !line.ends_with('\n') || line.trim_end_matches(['\r', '\n']).is_empty() {
                    return Ok(None);
                }
                SamRecord::parse(line)
                    .map(Some)
                    .map_err(anyhow::Error::from_boxed)
                    .with_context(|| format!("SAM line {number}"))
            }
            Self::Bam {
                reader,
                raw,
                header,
            } => {
                let mut size = [0u8; 4];
                // EOF before the first size byte is normal. EOF partway through
                // a size prefix or payload must instead report truncation.
                let n = reader.read_up_to(&mut size)?;
                if n == 0 {
                    return Ok(None);
                }
                reader
                    .read_exact_into(&mut size[n..])
                    .context("truncated BAM record size")?;
                let size = i32::from_le_bytes(size);
                if !(32..=2 * 1024 * 1024).contains(&size) {
                    bail!("invalid BAM record size {size}");
                }
                raw.resize(size as usize, 0);
                reader
                    .read_exact_into(raw)
                    .context("truncated BAM record")?;
                SamRecord::from_bam(raw, header)
                    .map(Some)
                    .map_err(anyhow::Error::from_boxed)
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
        header: BamHeader,
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
                header,
            })
        } else {
            out.write_all(text.as_bytes())?;
            Ok(Self::Sam(out))
        }
    }
    /// Serialize one staged record using its authoritative flag/position/MAPQ.
    ///
    /// # Errors
    /// Propagates I/O and serialization errors, including references absent from
    /// the BAM dictionary. Writes can be partial; callers should abort on failure.
    pub fn write(&mut self, record: &SamRecord) -> Result<()> {
        match self {
            Self::Sam(out) => record.write_sam(out)?,
            Self::Bam { writer, header } => writer
                .write_sam_record(record, header)
                .map_err(anyhow::Error::from_boxed)?,
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
