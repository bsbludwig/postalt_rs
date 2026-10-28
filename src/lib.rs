// SPDX-License-Identifier: GPL-3.0-only
// Part of the BWA/bwakit postalt Rust port; see NOTICE.md for provenance.

//! Streaming ALT-alignment postprocessing compatible with `bwa-postalt.js` r985.
//!
//! [`postalt::AltIndex`] loads ALT-to-primary mappings once. [`postalt::Processor`]
//! then stages alignments for one read name and read end, allowing subsequent ALT
//! hits to adjust earlier records before output. [`io`] handles sequential SAM,
//! gzip SAM and BGZF BAM without requiring an index or seekable input.
//!
//! # Compatibility contract
//!
//! The grouping, scoring, tag order and coordinate quirks of r985 are retained.
//! The default intentionally repairs opposite-strand sequence/quality reversal:
//! k8 1.2 removed the indexed `Bytes` access that the script relies on. Use
//! [`postalt::Processor::set_k8_compat`] to reproduce that runtime's unchanged
//! sequence/quality instead. Both modes preserve the script's extra coordinate
//! subtraction and its use of the current CIGAR to test earlier-record overlap.
//!
//! # Streaming lifecycle
//!
//! Keep records contiguous by read name **and** read end, in aligner order.
//! [`pipeline::run`] reads whole groups in batches, processes them on worker
//! threads and writes results in input order; [`postalt::Processor`] can also be
//! driven directly. Finalize the HLA writers and alignment output explicitly so
//! buffered I/O errors propagate.
//!
//! ```no_run
//! use postalt_rs::{io::{Input, Output}, pipeline::{self, Settings, Worker}, postalt::AltIndex};
//! use std::{path::Path, sync::Arc};
//!
//! # fn main() -> anyhow::Result<()> {
//! let index = Arc::new(AltIndex::from_path(Path::new("reference.fa.alt"))?);
//! let input = Input::open(None, 1)?; // stdin, one BGZF decompression worker
//! let mut output = Output::create(None, false, input.header_text(), 1, 6)?;
//! let settings = Settings { ratio: 1.0, k8_compat: false, fastq: false };
//! let headers = (input.bam_header().cloned(), output.bam_header().cloned());
//! let make = || Worker::new(Arc::clone(&index), &settings, headers.0.clone(), headers.1.clone());
//! pipeline::run(input, &mut output, &mut [], 2, make)?;
//! output.finish()?;
//! # Ok(())
//! # }
//! ```
pub mod io;
pub mod pipeline;
pub mod postalt;
pub mod sam_text;
