# Validation record — 2026-10-27

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

Host: AMD EPYC 9555P, Linux x86-64, glibc 2.34. Rust 1.98.1, release build with thin LTO and one codegen unit. Processing uses one thread; the BAM case uses two compression workers. No CPU affinity or exclusive host allocation was used.

Each row is the median of three sequential runs against the prepared input, discarding stdout to `/dev/null`. Input was already cached following correctness checks; no cold-cache claim is made. Startup, ALT-index construction and I/O are included. `/usr/bin/time` records wall/CPU time and peak RSS. Full individual measurements are in [benchmark.json](validation/benchmark.json).

| Path | Wall time | Peak RSS |
| --- | ---: | ---: |
| k8 reference, SAM → SAM | 62.13 s | 109.12 MiB |
| Rust installed-k8 mode, SAM → SAM | 4.85 s | 5.50 MiB |
| Rust corrected default, SAM → SAM | 4.90 s | 5.50 MiB |
| Rust corrected default, BAM → BAM, 2 workers | 9.79 s | 8.50 MiB |

On this workload, SAM processing is **12.8× faster** in compatibility mode, or **12.7× faster** with corrected reversal. These are dataset- and hardware-specific measurements. The BAM row includes decompression and compression and should not be interpreted as an identical-I/O comparison to the SAM rows.

The initial correctness runs writing full outputs to disk took 64.52 s for k8, 5.58 s for Rust compatibility mode, 5.75 s for corrected default, and 9.84 s for BAM → BAM. They are separate from the repeated `/dev/null` measurements.

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
