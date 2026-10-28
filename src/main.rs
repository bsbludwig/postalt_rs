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
  -j, --jobs INT         Processing and BAM decompression workers (default 1;\n\
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
    let mut ratio = 1.0;
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
            "-r" => ratio = value()?.parse::<f64>()?,
            "-o" | "--output" => output = Some(value()?.into()),
            "-O" | "--output-format" => {
                bam = match value()?.as_str() {
                    "sam" => false,
                    "bam" => true,
                    f => bail!("unsupported output format {f}"),
                };
            }
            "-@" | "--threads" => threads = value()?.parse()?,
            "-j" | "--jobs" => jobs = value()?.parse()?,
            "--compression" => compression = value()?.parse()?,
            "--" => {
                positional.extend(args);
                break;
            }
            _ if arg.starts_with("-r") && arg.len() > 2 => ratio = arg[2..].parse()?,
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
    if threads == 0 {
        bail!("threads must be at least 1");
    }
    if !(0..=12).contains(&compression) {
        bail!("compression must be between 0 and 12");
    }
    // Explicit output paths are an extension to the stdout-only script. Guard
    // identical paths and resolvable symlink aliases before File::create truncates.
    // Stdin ("-") cannot alias a file, and inputs that cannot be resolved are
    // left for the open below to report.
    if let Some(out) = output.as_deref().filter(|p| *p != Path::new("-")) {
        for input in positional.iter().filter(|p| *p != "-") {
            let same = out == Path::new(input)
                || (out.exists()
                    && matches!(
                        (std::fs::canonicalize(out), std::fs::canonicalize(input)),
                        (Ok(a), Ok(b)) if a == b
                    ));
            if same {
                bail!("output must not overwrite an input file");
            }
        }
    }
    let index = Arc::new(AltIndex::from_path(Path::new(&positional[0]))?);
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
    if let Some(prefix) = prefix.filter(|x| !x.is_empty()) {
        for gene in index.hla_genes() {
            let path = format!("{prefix}.{gene}.fq");
            let file = File::create(&path).with_context(|| format!("create {path}"))?;
            fastq.push(BufWriter::new(file));
        }
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
