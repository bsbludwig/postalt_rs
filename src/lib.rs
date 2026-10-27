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
//! Flush before pushing a different group, and again at EOF. Finalize both the
//! HLA writers and alignment output explicitly so buffered I/O errors propagate.
//!
//! ```no_run
//! use postalt_rs::{io::{Input, Output}, postalt::{AltIndex, Processor}};
//! use std::path::Path;
//!
//! # fn main() -> anyhow::Result<()> {
//! let index = AltIndex::from_path(Path::new("reference.fa.alt"))?;
//! let mut input = Input::open(None)?; // stdin
//! let mut output = Output::create(None, false, input.header_text(), 1, 6)?;
//! let mut processor = Processor::new(index, 1.0, None)?;
//! while let Some(record) = input.next_record()? {
//!     if processor.starts_new_group(&record) {
//!         processor.flush(|r| output.write(r))?;
//!     }
//!     processor.push(record)?;
//! }
//! processor.flush(|r| output.write(r))?;
//! processor.finish_fastq()?;
//! output.finish()?;
//! # Ok(())
//! # }
//! ```
pub mod io;
pub mod postalt;
