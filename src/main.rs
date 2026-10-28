// SPDX-License-Identifier: GPL-3.0-only
// Part of the BWA/bwakit postalt Rust port; see NOTICE.md for provenance.

//! CLI entry point and streaming read-end orchestration.
//!
//! Retains r985's `-p`, `-r` and `-v` meanings while adding explicit output/BAM
//! options and an installed-k8 reversal mode. Argument validation intentionally
//! reports errors for unsupported options and invalid numeric values instead of
//! reproducing JavaScript getopt/coercion behavior.
use anyhow::{Context, Result, bail};
use postalt_rs::{
    io::{Input, Output},
    pipeline::{self, Settings, Worker},
    postalt::AltIndex,
};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::Arc,
};

// Records are short-lived and numerous; mimalloc handles that pattern far
// better than glibc malloc (about 20-30% faster end to end when measured).
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
const HELP: &str = "Usage: postalt_rs [options] <alt.sam> [aln.sam|aln.bam|-]\n\n\
Postprocess BWA ALT alignments (bwa-postalt.js r985 compatibility).\n\
Input must be grouped by read name AND read end, in aligner order.\n\
Input defaults to stdin; output defaults to SAM on stdout.\n\n\
  -p STR                 Prefix for HLA gene FASTQ files\n\
  -r FLOAT               Minimum pa ratio (default 1, capped at 1)\n\
  -v                     Print reference compatibility version (r985)\n\
  -o, --output FILE      Output path (- for stdout)\n\
  -O, --output-format F  sam or bam (default sam)\n\
  -@, --threads INT      BAM compression workers (default 1)\n\
  -j, --jobs INT         Processing and BGZF decompression workers (default 1;\n\
                         0 runs everything on one thread)\n\
  --k8-1.2-compat        Reproduce installed k8's unchanged reverse sequences\n\
  --compression INT     BAM compression level 0..12 (default 6)\n\
  -h, --help             Show this help\n\
  --version              Show replacement version\n";
/// Validate CLI arguments, open streams and drive the processor to completion.
///
/// Errors retain their path/read context for `main` to report. Output may already
/// contain a header or earlier groups when a later record fails; processing is
/// streaming, not an atomic file replacement.
fn run() -> Result<()> {
    let mut args = std::env::args().skip(1).peekable();
    let mut prefix = None;
    let mut ratio: f64 = 1.0;
    let mut output: Option<PathBuf> = None;
    let mut k8_compat = false;
    let mut bam = false;
    let mut threads = 1;
    let mut jobs = 1;
    let mut compression = 6;
    let mut positional = Vec::new();
    while let Some(arg) = args.next() {
        let mut value = || {
            args.next()
                .with_context(|| format!("{arg} requires a value"))
        };
        match arg.as_str() {
            "-h" | "--help" => {
                print!("{HELP}");
                return Ok(());
            }
            // Preserve the reference's compatibility version; the replacement
            // package version is exposed separately as --version.
            "-v" => {
                println!("r985");
                return Ok(());
            }
            "--version" => {
                println!(
                    "postalt_rs {} (bwa-postalt.js r985)",
                    env!("CARGO_PKG_VERSION")
                );
                return Ok(());
            }
            "--k8-1.2-compat" => k8_compat = true,
            "-p" => prefix = Some(value()?),
            "-r" => ratio = parse_option(&arg, &value()?)?,
            "-o" | "--output" => output = Some(value()?.into()),
            "-O" | "--output-format" => {
                bam = match value()?.as_str() {
                    "sam" => false,
                    "bam" => true,
                    f => bail!("unsupported output format {f}"),
                };
            }
            "-@" | "--threads" => threads = parse_option(&arg, &value()?)?,
            "-j" | "--jobs" => jobs = parse_option(&arg, &value()?)?,
            "--compression" => compression = parse_option(&arg, &value()?)?,
            "--" => {
                positional.extend(args);
                break;
            }
            _ if arg.starts_with("-r") && arg.len() > 2 => ratio = parse_option("-r", &arg[2..])?,
            _ if arg.starts_with("-p") && arg.len() > 2 => prefix = Some(arg[2..].into()),
            _ if arg.starts_with('-') && arg != "-" => bail!("unknown option {arg}"),
            _ => positional.push(arg),
        }
    }
    if positional.is_empty() {
        eprint!("{HELP}");
        bail!("ALT file is required");
    }
    if positional.len() > 2 {
        bail!("expected an ALT file and optional alignment file");
    }
    // Deviation: reject nonfinite ratios instead of inheriting parseFloat/NaN
    // comparisons. Finite negative values remain allowed, matching r985.
    if !ratio.is_finite() {
        bail!("-r must be a finite number");
    }
    if !(1..=MAX_THREADS).contains(&threads) {
        bail!("threads must be between 1 and {MAX_THREADS}");
    }
    if jobs > MAX_THREADS {
        bail!("jobs must be at most {MAX_THREADS}");
    }
    if !(0..=12).contains(&compression) {
        bail!("compression must be between 0 and 12");
    }
    let index = Arc::new(AltIndex::from_path(Path::new(&positional[0]))?);
    // The FASTQ files are created (truncated) too: none may alias the output
    // or an input.
    let fastq_paths: Vec<String> = match prefix.as_deref().filter(|x| !x.is_empty()) {
        Some(prefix) => index
            .hla_genes()
            .iter()
            .map(|gene| format!("{prefix}.{gene}.fq"))
            .collect(),
        None => Vec::new(),
    };
    let created: Vec<(&str, &Path)> = output
        .as_deref()
        .filter(|p| *p != Path::new("-"))
        .map(|p| ("output", p))
        .into_iter()
        .chain(
            fastq_paths
                .iter()
                .map(|p| ("HLA FASTQ output", Path::new(p))),
        )
        .collect();
    let inputs = positional
        .iter()
        .filter(|p| *p != "-")
        .map(|p| ("input", Path::new(p)));
    for (i, &(role, path)) in created.iter().enumerate() {
        for (other, existing) in inputs.clone().chain(created[..i].iter().copied()) {
            if same_file(path, existing) {
                bail!(
                    "{role} {} is also the {other} {}",
                    path.display(),
                    existing.display()
                );
            }
        }
    }
    let input = Input::open(positional.get(1).map(Path::new), jobs)?;
    let mut out = Output::create(
        output.as_deref(),
        bam,
        input.header_text(),
        threads,
        compression,
    )?;
    // A nonempty prefix creates/truncates `<prefix>.<gene>.fq` even if no reads
    // eventually match that gene. Files already created are not rolled back.
    let mut fastq = Vec::new();
    for path in &fastq_paths {
        let file = File::create(path).with_context(|| format!("create {path}"))?;
        fastq.push(BufWriter::new(file));
    }
    let settings = Settings {
        ratio,
        k8_compat,
        fastq: !fastq.is_empty(),
    };
    let (in_header, out_header) = (input.bam_header().cloned(), out.bam_header().cloned());
    let make_worker = || {
        let (i, o) = (in_header.clone(), out_header.clone());
        Worker::new(Arc::clone(&index), &settings, i, o)
    };
    pipeline::run(input, &mut out, &mut fastq, jobs, make_worker)?;
    // Finalize both the HLA writers and alignment output explicitly so
    // buffered/compression failures surface.
    for file in &mut fastq {
        file.flush()?;
    }
    out.finish()
}
/// Most worker or compression threads accepted; more only exhausts the system.
const MAX_THREADS: usize = 1024;

/// Parse an option's value, naming the option on failure.
fn parse_option<T>(option: &str, value: &str) -> Result<T>
where
    T: std::str::FromStr,
    T::Err: std::error::Error + Send + Sync + 'static,
{
    value
        .parse()
        .with_context(|| format!("invalid value {value:?} for {option}"))
}

/// Whether two paths name the same file: identical, or resolving to the same
/// existing file (symlinks). A path that does not exist yet aliases only
/// itself.
fn same_file(a: &Path, b: &Path) -> bool {
    a == b
        || matches!(
            (std::fs::canonicalize(a), std::fs::canonicalize(b)),
            (Ok(x), Ok(y)) if x == y
        )
}

/// Report processing failures, treating a downstream closed pipe as completion.
fn main() -> std::process::ExitCode {
    match run() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(error) => {
            // Deviation from runtime-dependent k8 diagnostics: a consumer such
            // as `head` may close stdout early without making this CLI fail.
            if error.chain().any(|e| {
                e.downcast_ref::<std::io::Error>()
                    .is_some_and(|e| e.kind() == std::io::ErrorKind::BrokenPipe)
            }) {
                return std::process::ExitCode::SUCCESS;
            }
            eprintln!("postalt_rs: {error:#}");
            std::process::ExitCode::FAILURE
        }
    }
}
