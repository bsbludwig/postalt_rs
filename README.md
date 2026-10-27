# postalt_rs

A streaming Rust replacement for `bwa-postalt.js` r985, using [seqair](https://github.com/Softleif/seqair) for alignment I/O. It lifts ALT hits, groups competing alignments, recalculates MAPQ, emits supplementary ALT alignments, and optionally extracts HLA FASTQ files.

Default behavior fixes reverse-complementing of generated opposite-strand records. Use `--k8-1.2-compat` to reproduce the supplied k8 1.2 runtime's output, including its failure to reverse sequence and quality. All other postalt algorithm quirks are retained.

## Build and run

The pinned toolchain is Rust 1.98.1. A C compiler, pkg-config, and liblzma development headers are needed for seqair's dependencies. No k8 or samtools installation is needed to run the replacement.

```bash
cargo build --release --locked

# SAM input/output; input defaults to stdin
./target/release/postalt_rs reference.fa.alt alignments.sam > adjusted.sam
bwa-mem2 mem reference.fa reads_1.fq.gz reads_2.fq.gz \
  | ./target/release/postalt_rs reference.fa.alt \
  | samtools sort -o adjusted.bam

# BAM input and output, with two BAM compression workers
./target/release/postalt_rs -O bam -@ 2 -o adjusted.bam reference.fa.alt alignments.bam

# Reproduce the installed k8 1.2 behavior
./target/release/postalt_rs --k8-1.2-compat reference.fa.alt alignments.sam > compatible.sam

# Optional HLA extraction and pa threshold
./target/release/postalt_rs -p hla_reads -r 0.8 reference.fa.alt alignments.sam > adjusted.sam
```

Input is detected by content: SAM, gzip/BGZF SAM, or BGZF BAM, from a file or pipe. The ALT mapping file may also be gzip-compressed. Output is SAM by default; `-O bam` selects BAM independently of the filename. BAM output needs an `@SQ` dictionary containing every emitted reference. Headerless SAM remains supported for SAM output. CRAM is not supported.

**Keep all records for each read name and read end contiguous, in aligner order.** Sorting only by read name can interleave mates. Primary records must precede the ALT records that adjust them. Coordinate-sorted input is unsuitable. The program holds only the current read-end group and the ALT index in memory. Output order and headers are preserved as in the script; output is not coordinate-sorted, even if an input header claims otherwise. Sort the output before indexing it.

| Option | Meaning |
| --- | --- |
| `-p STR` | Write `STR.HLA-<gene>.fq` files when the ALT file contains `HLA-<gene>*<allele>` names |
| `-r FLOAT` | Reduce nonoverlapping buffered records with `pa < FLOAT` to MAPQ 0; default 1, capped at 1 |
| `-v` | Print reference compatibility version `r985` |
| `--version` | Print replacement version |
| `-o FILE` | Output file; default stdout, `-` explicitly selects stdout |
| `-O sam\|bam` | Output format; default SAM |
| `-@ INT` | BAM compression workers; default 1; does not parallelize postalt processing |
| `--compression INT` | BAM compression level 0–12; default 6 |
| `--k8-1.2-compat` | Leave generated opposite-strand sequence and quality unchanged |

## Compatibility choices

The reference is the installed `/opt/apps/manual/bwa-mem2/2.3.0/noarch/bin/bwa-postalt.js`, SHA-256 `50660424f540ad38e08c181f628289754a39b38e2092f3010cf7478e1d0d83b7`.

Preserved behavior includes:

- The extra subtraction in lift coordinates: the primary offset is SAM POS minus **two**.
- The use of the current ALT record's CIGAR span when testing overlap of earlier buffered records.
- First-lift-only grouping, stable tie ordering, transitive overlap grouping, MAPQ floors/caps, `om` restoration, and duplicate optional tags.
- Hard-clipped reported hits bypassing adjustment; `=`, `X`, and `P` ignored by postalt scoring/projection; `N` ignored by ALT CIGAR projection but counted in reference spans.
- The reference's end-clipping boundary behavior and its use of the first buffered record for HLA FASTQ sequence/quality.
- Blank lines terminating SAM input, and an unterminated final SAM line being omitted. Ensure every input record ends in a newline.

[k8 1.x removed indexed access to `Bytes`](https://github.com/attractivechaos/k8/releases/tag/v1.0), which the reference uses for reversal. The default corrects that operation using the reference's IUPAC complement table. Missing sequence/quality (`*`) remain missing. Compatibility mode reproduces the unchanged bytes from the supplied k8 1.2 runtime. Neither mode reproduces runtime memory corruption.

SAM output retains original optional-field order, duplicate tags, and decimal spelling. BAM conversion uses the BAM types, so case, reference aliases, integer storage widths, and decimal spelling may normalize; comparisons should use decoded values or samtools normalization. Invalid records and unsupported CLI options produce explicit errors rather than reproducing JavaScript coercions or crashes.

## Tests

```bash
cargo fmt --all --check
cargo clippy --all-targets -- -D warnings
cargo test --locked

# Live differential checks: requires k8 and the reference script
module load bwa-mem2
python3 scripts/compatibility.py --random 10000 --seed 20261027
```

The checked-in synthetic fixtures were generated independently by k8, not by the Rust implementation. `cargo test` needs no external reference runtime or patient data. It covers both modes, projection boundaries, interval indexing, grouping and score ties, mate/name boundaries, MAPQ restoration, HLA output, reversal and missing fields, SAM/BAM streams, gzip, auxiliary types, multithreaded BAM output, truncation, and CLI errors.

The differential runner compares SAM and HLA files byte-for-byte. For default mode it makes a temporary reference copy replacing only the reversal code. For HLA cases that crash k8's string-writing path, a documented test-only wrapper writes the same ASCII bytes through `Uint8Array`/`ArrayBuffer`. The ALT/grouping/MAPQ/HLA algorithms are unchanged. `--write-golden` explicitly regenerates fixtures; `tests/fixtures/manifest.json` records which cases use the safe write wrapper. `--reference`, `--k8`, `POSTALT_REFERENCE`, and `K8` allow other installations.

## Real data and performance

[VALIDATION.md](VALIDATION.md) records the supplied NA12878 evaluation, including exact comparison results and timing methodology. To reproduce with an indexed BAM:

```bash
module load samtools bwa-mem2
python3 scripts/prepare_real_data.py sample.bam target/validation/real
python3 scripts/validate_real.py reference.fa.alt target/validation/real/input.sam \
  --bam target/validation/real/input.bam --json target/validation/real-validation.json
python3 scripts/benchmark.py reference.fa.alt target/validation/real/input.sam \
  --bam target/validation/real/input.bam --repeats 3 \
  --json target/validation/benchmark.json
```

Extraction includes all `_alt` references and `chr6:30000000-31000000`. The script sorts by name, then groups read ends and puts original primary records first. This simulates streaming input; it cannot reconstruct the aligner's original order or restore records/tags already changed by earlier processing. It selects alignments in these regions, not every distant alignment of the same template.

The implementation precomputes ALT CIGAR segments for binary-search projection and uses sorted interval indexes with prefix maxima. Buffered I/O and read-end-local processing avoid whole-file storage. Seqair supplies BGZF decompression and multithreaded compression. See [vendor/seqair/POSTALT_PATCHES.md](vendor/seqair/POSTALT_PATCHES.md) for the pinned upstream revision and the streaming adapter.

## License

GPLv3 (`GPL-3.0-only`), with the full license in [LICENSE](LICENSE). This project
ports `bwa-postalt.js` from BWA/bwakit, whose [upstream README declares GPLv3](https://github.com/lh3/bwa/blob/27dd1da7020151be9358fce67cbe2c10801c348d/README.md#availability).
The script's installation alongside BWA-MEM2 does not supply an Intel/MIT license
for this port. See [NOTICE.md](NOTICE.md) for authorship, the pinned upstream
source, and the comparison with the installed reference. Vendored seqair's
upstream code retains its MIT/Apache-2.0 notices.
