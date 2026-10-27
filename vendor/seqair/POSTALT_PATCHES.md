# Pinned seqair source

Source: https://github.com/Softleif/seqair
Revision: `94f9317bed689e0a6977eb3470bf1ec8382a253f`
Packages: seqair 0.8.0 and seqair-types 0.4.0.

This snapshot includes production source, upstream licenses, README files and the two upstream test-support modules referenced from source (needed by rustfmt). Unrelated crates, fixtures, benchmarks and dev-dependencies are omitted. Dependency versions are locked by the root Cargo.lock. The original workspace dependency declarations are retained; package manifests retain only production dependencies and the cfg lint declarations.

Upstream currently exposes indexed SAM readers, but no public sequential SAM parser/writer. Its pileup-oriented owned base type only represents A/C/G/T/N; generic RecordStore decoding also requires placed positions. Those APIs alone would lose unmapped records, IUPAC ambiguity bases, or SAM lexical details in this application.

The production source differs in exactly four files:

- `src/sam.rs`: exports the added `stream` module.
- `src/sam/reader.rs`: gives sibling-module access to the existing CIGAR and auxiliary parsers (`pub(super)`). No parser logic changes.
- `src/sam/stream.rs`: added text-preserving `SamRecord` parsing/writing, conversions using the upstream CIGAR/aux parsers, BAM header decoder, wire validation, and sequence codecs. Handles BAM POS=-1 by temporarily supplying zero to the indexed parser, restoring the original sentinel afterward. Decodes raw sequence with seqair's full IUPAC codec, bypassing the five-base pileup representation.
- `src/bam/writer.rs`: adds `write_sam_record`. Uses the upstream owned encoder and BGZF writer, then preserves the full IUPAC sequence in the packed record before serialization. Existing upstream APIs are unchanged.

The application owns stream orchestration and read-end grouping; the adapter owns alignment serialization. No alternate BAM/SAM library is used. Tests in the parent project exercise these additions, including all BAM auxiliary types supported by the adapter, IUPAC, missing sequence/quality, unmapped records, BAM pipes, gzip SAM, and one/two-worker writers.

Upstream MIT and Apache-2.0 licenses remain in this directory. The parent project does not run the omitted upstream integration suite; it tests the used I/O paths and cross-checks real outputs with samtools.
