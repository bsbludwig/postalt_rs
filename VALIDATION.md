# Validation record — 2026-10-27, updated 2026-10-28

Validated on the supplied BAM and ALT mapping file. Raw sample reads and large outputs stay under ignored `target/validation/`; checked-in fixtures are synthetic. Machine-readable timing and output hashes are in [validation/](validation/).

## Reference and inputs

- Reference: BWA/bwakit's `bwa-postalt.js`, version r985, installed alongside BWA-MEM2 2.3.0, SHA-256 `50660424f540ad38e08c181f628289754a39b38e2092f3010cf7478e1d0d83b7`. Its source and GPLv3 provenance are documented in [NOTICE.md](NOTICE.md).
- Installed runtime: k8 1.2-r137, V8 10.2.154.26-node.28.
- Original BAM: `/users/bschuster/sharedscratch/WM_GIAB/Results/2.0/Alignments/ce3a42619fd344bbb0f4b134e14569a6.bam`.
- ALT file: `/mnt/lustre/shared/genomes/hg38+HPV+HBV+HVC+EBV+TAPS_SI_v2/hg38_full_gatk_HPV_HBV_HCV_spike-ins_v2.fa.gz.alt`.
- Extracted all 261 `_alt` references and `chr6:30000000-31000000` using the BAM index.
- Prepared input: **1,042,930 alignment records**, 3,380 header lines, 909,053,220 bytes of SAM.
- Postalt output: **3,628,954 alignment records**, with the same headers.

The BAM header reports bwa-mem3 alignment followed by sorting/duplicate marking. This tests compatibility on the available ALT-rich data, not an untouched bwa-mem2 stream. Name sorting and read-end regrouping cannot recover the original aligner order. The extraction includes regional records, not all distant records sharing a selected query name.

## Results

- Installed-k8 mode: Rust SAM output matched unmodified k8/reference output **byte for byte** over the full prepared input.
- Corrected default: Rust SAM output matched a temporary reference with only its sequence/quality reversal repaired **byte for byte** over the full prepared input.
- BAM input → BAM output: passed `samtools quickcheck`; decoding and normalization with samtools 1.21 matched the normalized default SAM output **byte for byte**. No filtering or sorting was applied during comparison.
- **61 Rust tests passed**: three focused algorithm tests and 58 integration tests, including 49 golden cases exercised in both modes.
- Live differential comparison passed all 49 targeted cases and **10,000 randomized templates**, seed `20261027`, in both modes. This includes HLA FASTQ comparisons.
- `cargo fmt --all --check`, `cargo clippy --all-targets -- -D warnings`, and release build passed.

The k8 runtime has two relevant defects. Indexed `Bytes` access was removed in k8 1.x, so the reference does not actually reverse generated opposite-strand sequence/quality. The default intentionally corrects this, as requested; `--k8-1.2-compat` retains the installed runtime behavior. Also, k8 1.2's `File.write(string)` allocates the string length while V8 writes an additional NUL, causing heap corruption on some HLA test inputs. A test-only wrapper copies ASCII bytes into a `Uint8Array` and uses the ArrayBuffer overload for those cases. The reference's HLA selection and output formatting are unchanged. See [k8's source](https://github.com/attractivechaos/k8/blob/v1.2/k8.cc) and [breaking changes](https://github.com/attractivechaos/k8/releases/tag/v1.0).

## Performance

Host: AMD EPYC 9555P, Linux x86-64, glibc 2.34. Rust 1.98.1, release build with thin LTO and one codegen unit, mimalloc allocator. `-j` sets processing/decompression workers (default 1) and `-@` BAM compression workers. No CPU affinity or exclusive host allocation was used.

Each row is the median of three sequential runs against the prepared input, discarding stdout to `/dev/null`. Input was already cached following correctness checks; no cold-cache claim is made. Startup, ALT-index construction and I/O are included. `/usr/bin/time` records wall/CPU time and peak RSS. Full individual measurements are in [benchmark.json](validation/benchmark.json) (re-measured 2026-10-28).

| Path | Wall time | CPU time | Peak RSS |
| --- | ---: | ---: | ---: |
| k8 reference, SAM → SAM | 62.24 s | 62.62 s | 109.7 MiB |
| Rust installed-k8 mode, SAM → SAM (`-j 1`) | 1.67 s | 1.86 s | 46.3 MiB |
| Rust corrected default, SAM → SAM (`-j 1`) | 1.66 s | 1.85 s | 48.3 MiB |
| Rust corrected default, SAM → SAM, `-j 0` | 1.83 s | 1.82 s | 26.3 MiB |
| Rust corrected default, SAM → SAM, `-j 4` | 0.50 s | 2.12 s | 88.4 MiB |
| Rust BAM → BAM, `-j 1 -@ 2` | 3.72 s | 11.35 s | 72.4 MiB |
| Rust BAM → BAM, `-j 4 -@ 4` | 1.70 s | 11.86 s | 130.5 MiB |
| `samtools view -b -@ 4` (copy only, for scale) | 2.12 s | 8.52 s | 8.5 MiB |

On this workload, SAM processing is **37× faster** than k8 with the default single processing worker. The BAM rows include decompression and compression and should not be compared directly with the SAM rows.

### Optimization of 2026-10-28

The ALT-rich dataset above is about 85% XA records, unlike whole-genome data, where most reads pass through unchanged. A second, pass-through-heavy dataset was prepared from the same BAM with `prepare_real_data.py --regions chr1:100000000-115000000`: 2,893,003 records (1.98 GB SAM), 5% with XA tags.

Profiling the previous build showed about 40% of time in glibc malloc/free (one heap string per SAM field), and BAM records went through a full BAM → text → BAM round trip. The optimized build keeps each record in one reused buffer, keeps BAM records in wire form unless they change, skips hit parsing for records that cannot involve an ALT contig, uses mimalloc, and processes batches of whole read-end groups in parallel, writing them in input order.

Equivalence: on both datasets, output was byte-identical to the previous build for SAM → SAM (both reversal modes), BAM → SAM, and samtools-normalized BAM → BAM and SAM → BAM, at `-j 0`, `1`, `4` and `8`. Truncated BAM input, an invalid FLAG and a malformed XA tag mid-stream produced identical partial output and messages. `validate_real.py` passed with hashes identical to [real-validation.json](validation/real-validation.json). The k8 differential suite (49 cases and 10,000 randomized templates, both modes) passed at `-j 0`, `1` and `4`. Rust tests: 4 unit, 62 integration (golden cases at several `-j` values, multi-batch streams, error semantics) and 1 doctest.

Old versus new build, median of three warm runs (full data in [optimization-benchmark.json](validation/optimization-benchmark.json)):

| Dataset | Path | Old | New | Peak RSS old → new |
| --- | --- | ---: | ---: | ---: |
| pass-through | SAM → SAM | 2.91 s | 0.83 s (`-j 1`), 0.45 s (`-j 8`) | 6 → 35 / 89 MiB |
| pass-through | BAM → SAM | 5.22 s | 2.73 s (`-j 0`), 0.44 s (`-j 4`) | 7 → 55 / 115 MiB |
| pass-through | BAM → BAM, `-@ 4` | 8.52 s | 2.36 s (`-j 4`); samtools copy 4.29 s | 9 → 108 MiB |
| pass-through | BAM → BAM, `-@ 8` | — | 1.20 s (`-j 4`); samtools copy 2.18 s | 133 MiB |
| ALT-rich | SAM → SAM | 4.97 s | 1.68 s (`-j 1`), 0.30 s (`-j 8`) | 6 → 49 / 139 MiB |
| ALT-rich | BAM → BAM, `-@ 4` | 9.93 s | 1.70 s (`-j 4`); samtools copy 2.12 s | 10 → 125 MiB |

Higher memory use comes from the batches in flight (at most `2 × jobs + 2`, about 512 KiB of input each) and per-thread allocator caches. It is bounded and does not grow with input size.

## Reproduce

```bash
cargo build --release --locked
module load samtools bwa-mem2

python3 scripts/prepare_real_data.py "$BAM" target/validation/real
python3 scripts/validate_real.py "$ALT" target/validation/real/input.sam \
  --bam target/validation/real/input.bam \
  --json target/validation/real-validation.json
python3 scripts/benchmark.py "$ALT" target/validation/real/input.sam \
  --bam target/validation/real/input.bam --repeats 3 \
  --json target/validation/benchmark.json
```

Set `BAM` and `ALT` to the paths above. `validate_real.py` streams outputs into SHA-256/byte/line comparisons instead of writing multi-gigabyte SAM outputs. The HLA checks use the synthetic differential suite because the supplied ALT file's contig names do not exercise the script's HLA allele naming convention.
