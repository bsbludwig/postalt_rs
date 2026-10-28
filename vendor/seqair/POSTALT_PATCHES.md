# Pinned seqair source

Source: https://github.com/Softleif/seqair
Revision: `94f9317bed689e0a6977eb3470bf1ec8382a253f`
Packages: seqair 0.8.0 and seqair-types 0.4.0.

This snapshot includes production source, upstream licenses, README files and the two upstream test-support modules referenced from source (needed by rustfmt). Unrelated crates, fixtures, benchmarks and dev-dependencies are omitted. Dependency versions are locked by the root Cargo.lock. The original workspace dependency declarations are retained; package manifests retain only production dependencies and the cfg lint declarations.

Upstream currently exposes indexed SAM readers, but no public sequential SAM parser/writer. Its pileup-oriented owned base type only represents A/C/G/T/N; generic RecordStore decoding also requires placed positions. Those APIs alone would lose unmapped records, IUPAC ambiguity bases, or SAM lexical details in this application.

The production source differs in exactly six files:

- `Cargo.toml`: adds `memchr` (already in the dependency graph) for SAM field splitting.
- `src/sam.rs`: exports the added `stream` module.
- `src/sam/reader.rs`: gives sibling-module access to the existing CIGAR and auxiliary parsers (`pub(super)`). No parser logic changes.
- `src/sam/stream.rs`: added streaming adapter.
  - `SamRecord` keeps one line in a single reusable buffer with field offsets. FLAG/POS/MAPQ are decoded and kept in canonical decimal text, so SAM output writes the buffer verbatim. `append_bam` encodes directly to BAM wire format with seqair's CIGAR/aux parsers, fixed-header encoder, bin computation and sequence codec, preserving full IUPAC sequence.
  - `BamRecord` keeps a validated BAM payload in wire form (validation handles POS=-1 by temporarily supplying zero to the indexed parser, restoring the sentinel afterward). It exposes the fields postalt needs without decoding, re-emits the payload unchanged, and decodes to `SamRecord` with seqair's full IUPAC codec on demand.
  - `BgzfStream` continues a stream after header parsing and inflates blocks ahead of the reader on optional worker threads. Framing, CRC and size checks mirror `BgzfReader::read_block`, including its end-of-stream rules.
  - Also: BAM header decoding for sequential streams and a CIGAR parse helper.
- `src/bam/bgzf.rs`: adds `BgzfReader::into_parts`, returning the unread decompressed remainder, the inner reader and the EOF flag. No reader logic changes.
- `src/bam/writer.rs`: adds `write_bam_records`, which writes size-prefixed BAM payloads through the existing size validation, BGZF writer and optional index path. Existing upstream APIs are unchanged.

The application owns stream orchestration and read-end grouping; the adapter owns alignment serialization. No alternate BAM/SAM library is used. Tests in the parent project exercise these additions, including all BAM auxiliary types supported by the adapter, IUPAC, missing sequence/quality, unmapped records, BAM pipes, gzip SAM, one/two-worker writers, read-ahead decompression, verbatim BAM copying and multi-batch streams.

Upstream MIT and Apache-2.0 licenses remain in this directory. The parent project does not run the omitted upstream integration suite; it tests the used I/O paths and cross-checks real outputs with samtools.
