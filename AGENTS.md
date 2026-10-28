# Working on postalt_rs

This repository implements a streaming Rust replacement for BWA/bwakit's
`bwa-postalt.js` r985. Read [README.md](README.md) for behavior and CLI usage,
[VALIDATION.md](VALIDATION.md) for validation evidence, and
[NOTICE.md](NOTICE.md) for source provenance.

## Source map

- `src/main.rs`: CLI validation, file setup and the global allocator (mimalloc).
- `src/pipeline.rs`: batches of whole read-end groups, worker threads, ordered
  output, and sequential error semantics (batch "tails").
- `src/postalt.rs`: ALT index, CIGAR projection, alignment grouping, MAPQ
  adjustment, supplementary records, HLA FASTQ text, pass-through fast path.
- `src/io.rs`: SAM/gzip/BAM stream detection, batch framing, `Record`
  (SAM text or BAM wire form), serialization, and output.
- `vendor/seqair/`: pinned alignment parser/writer and local streaming adapter.
  Read `vendor/seqair/POSTALT_PATCHES.md` before modifying or updating it.
- `tests/compatibility.rs`: reference-derived golden cases and I/O regression tests.
- `scripts/compatibility.py`: live k8 differential tests and fixture generation.
- `scripts/prepare_real_data.py`, `validate_real.py`, `benchmark.py`: real-data
  preparation, output comparison, and performance measurement.
- `validation/`: saved output hashes and benchmark results. Large data and
  disposable experiments belong under ignored `target/`.

## Development and validation

Use the pinned Rust toolchain from `rust-toolchain.toml` (currently 1.98.1).
Build dependencies include a C compiler, pkg-config, and liblzma headers.
Run from the repository root:

```bash
cargo fmt --all --check
cargo clippy --all-targets -- -D warnings
cargo test --locked
cargo build --release --locked
```

For algorithm or serialization changes, also run differential validation:

```bash
module load bwa-mem2
python3 scripts/compatibility.py --random 10000 --seed 20261027
```

The differential runner defaults to the release binary; rebuild it before
testing a change. `--reference` / `POSTALT_REFERENCE` and `--k8` / `K8` allow
alternate installations. Expected outputs must come from the reference,
never from the Rust implementation being tested. Explicit regeneration is:

```bash
module load bwa-mem2
python3 scripts/compatibility.py --write-golden
```

Check fixture availability on a fresh checkout: as observed on 2026-10-28,
`.gitignore` excludes `tests/fixtures`, and only a few fixture files are tracked.
Most golden files therefore need local regeneration before `cargo test`.
Do not mistake missing fixture files for an algorithm regression, or silently
rewrite existing expected output to make a failing test pass.

Test changes at the relevant level; documentation-only edits need a content
and diff review rather than rerunning expensive data processing. Benchmark
release builds with matched input/output formats, compression settings, and
worker counts. Record repetitions, cache state, hardware, and output sink;
check output equivalence before interpreting timing differences.

## Memory

Maintain this section across sessions. Save explicit user preferences,
accepted behavior decisions, and important verified findings here when they
change. Date measurements and environment-specific observations; include
evidence paths and limitations. Replace stale entries instead of accumulating
contradictory notes. Keep transient progress and large logs out of this file.
Current user instructions take precedence over remembered instructions.

### Standing user instructions and decisions
- Do not commit to git on the user's behalf. Do not run any mutating git commands, unless explicitly told to do so.
- Use **seqair** for production BAM/SAM parsing and writing. Performance and
  comprehensive compatibility testing against `bwa-postalt.js` are priorities.
- Add proper Rust doc comments, Python docstrings, and explanatory comments
  throughout the source, especially wherever behavior deviates from the script.
  Explain compatibility quirks and their reasons at the relevant code.
- Correct reverse-complementing of generated opposite-strand sequence and
  reversal of quality **by default**. Keep explicit `--k8-1.2-compat` mode to
  reproduce the installed runtime's unchanged sequence/quality output.
- The script is from **BWA/bwakit**, despite being installed with BWA-MEM2.
  The port uses **GPL-3.0-only**, not the Intel/BWA-MEM2 MIT license. Preserve
  `LICENSE`, `NOTICE.md`, and vendored seqair's upstream MIT/Apache-2.0 notices.
- The htslib comparison was an experiment, not a request to change the
  production backend.

### Compatibility findings

- Reference: `/opt/apps/manual/bwa-mem2/2.3.0/noarch/bin/bwa-postalt.js`;
  load it and k8 with `module load bwa-mem2`. Reference SHA-256:
  `50660424f540ad38e08c181f628289754a39b38e2092f3010cf7478e1d0d83b7`.
- Preserve the extra coordinate subtraction: the lift's primary offset is
  SAM POS minus **two**, not the usual POS minus one. Also preserve use of the
  **current ALT record's CIGAR span** when checking overlap of earlier buffered
  records. These are intentional compatibility choices, not cleanup targets.
- Preserve first-lift grouping, stable ties, duplicate SAM tags, lexical SAM
  fields, and the reference's clipping/scoring quirks. Blank SAM lines end input;
  an unterminated final SAM line is omitted. See README for the complete list.
- Records must be contiguous by **read name and read end**, with primaries
  before the ALT records that can adjust them. Name sorting alone can interleave
  mates. Flush the final group at EOF. Keep memory bounded to the ALT index and
  the batches in flight (`pipeline.rs`); avoid whole-file buffering.
- k8 1.2 lacks the indexed `Bytes` access used by the script's reversal code.
  It also has string-writing heap corruption on some HLA cases. The differential
  runner uses documented, temporary reference repairs for reversal and safe
  byte writing where appropriate; do not modify the installed reference or
  attempt to reproduce memory corruption.
- SAM output comparisons are byte-for-byte. BAM can normalize lexical details;
  compare decoded values or normalize with `samtools view --no-PG -h`.
- Seqair is pinned at `94f9317bed689e0a6977eb3470bf1ec8382a253f` (0.8.0,
  seqair-types 0.4.0). Its adapter preserves unmapped records, IUPAC bases, and
  SAM lexical details that the upstream indexed/pileup APIs alone would lose.
- Seqair's `compression_threads(0)` is synchronous; positive values create
  that many compression workers. The application CLI currently requires at
  least one worker. `-@` sets BAM compression workers only; `-j` sets
  processing and BAM decompression workers (`-j 0`: single thread).
- Output must be byte-identical for every `-j`. Within a batch, a group's
  output is only written once sequential processing would have flushed it
  (next record parsed, or clean EOF); keep that rule when changing errors.

### Real-data validation and performance

The supplied NA12878 inputs are:

- BAM: `/users/bschuster/sharedscratch/WM_GIAB/Results/2.0/Alignments/ce3a42619fd344bbb0f4b134e14569a6.bam`
- ALT: `/mnt/lustre/shared/genomes/hg38+HPV+HBV+HVC+EBV+TAPS_SI_v2/hg38_full_gatk_HPV_HBV_HCV_spike-ins_v2.fa.gz.alt`

Prepared data, when present, is under `target/validation/real/`. It includes all
261 ALT references plus `chr6:30000000-31000000`: 1,042,930 input records and
3,628,954 output records. The source is coordinate-sorted and previously
processed; regrouping cannot recover original aligner order or missing distant
alignments of the same template. It is not a pristine bwa-mem2 stream.

The 2026-10-27 validation passed both reversal modes, 49 golden cases, and
10,000 randomized templates. Full-data SAM comparisons matched the respective
references, and normalized BAM output matched corrected SAM. Evidence lives in
`VALIDATION.md` and `validation/real-validation.json`. Historical test counts
and results are not substitutes for checks after a code change.

The separate htslib experiment used a **custom Rust–C bridge to htslib 1.21**,
**not rust-htslib**, and the unchanged production postalt core. The harness and
raw results are under `target/htslib-comparison/` if retained; `cargo clean`
can remove them. Median elapsed seconds over three warm-cache runs, writing to
`/dev/null`, with BAM compression level 6 and synchronous input decompression:

| Path | seqair | htslib |
| --- | ---: | ---: |
| SAM to SAM | 4.89 | 7.40 |
| BAM to SAM | 6.13 | 7.62 |
| BAM to BAM, synchronous compression | 15.51 | 16.47 |
| BAM to BAM, two compression workers | 9.88 | 7.39 |

Pass-through benchmark (2026-10-28): the ALT/MHC dataset above is ~85% XA
records and unrepresentative of WGS, where most reads are pass-through.
`target/validation/passthrough/input.{sam,bam}` holds `chr1:100000000-115000000`
from the same BAM (`prepare_real_data.py --regions ...`; 2,893,003 records, 5%
with XA). Profiling before optimization showed ~40% glibc malloc/free from
`Vec<String>` records and a full BAM->text->BAM round trip per record.

Optimization (2026-10-28, user-requested): single-buffer `SamRecord`, raw
`BamRecord` pass-through, record/buffer pools, exact fast path for records
without XA off ALT/HLA contigs, mimalloc, read-ahead BGZF decompression and a
batch pipeline (`-j`). Outputs are byte-identical to the previous binary on
both datasets in all formats and `-j` values, including error paths; k8
differential (10,000 templates, `-j 0/1/4`) and `validate_real.py` passed.
Warm-cache medians, `/dev/null` sink (old -> new): pass-through SAM->SAM
2.91 -> 0.83 s (`-j 1`) / 0.45 s (`-j 8`); BAM->BAM `-@ 4` 8.52 -> 2.36 s
(samtools view 4.29 s); ALT-rich SAM->SAM 4.97 -> 0.30 s (`-j 8`). Peak RSS
rose from ~7 MiB to ~30 MiB (`-j 0/1`) and up to ~140 MiB (`-j 8`): in-flight
batches plus per-thread allocator caches. Evidence:
`validation/optimization-benchmark.json`, `VALIDATION.md`.

Normalized outputs matched on this dataset. These measurements include the
bridge's extra conversions between htslib records and the core's SAM text
fields. They do not measure rust-htslib or a core operating directly on htslib
records, and do not establish a universal htslib slowdown. The two-worker BAM
advantage is elapsed time, not lower total CPU time. Production still uses seqair.
