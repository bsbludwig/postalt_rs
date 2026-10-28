//! Malformed, unusual and adversarial input: each case fails with a clear
//! error or is handled as r985 handles it, never silently wrong output.
//! These tests need only the tracked `tests/fixtures/alt.sam`.
use seqair::bam::bgzf::BgzfStream;
use seqair::io::BgzfWriter;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn alt() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/alt.sam")
}

const HEADER: &str = "@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:altA\tLN:2000\n@SQ\tSN:altB\tLN:2000\n";

fn read(name: &str) -> String {
    let seq = "A".repeat(100);
    let qual = "I".repeat(100);
    format!("{name}\t0\tchr1\t1100\t60\t100M\t*\t0\t0\t{seq}\t{qual}\tNM:i:0\n")
}

/// Run the executable with `args`, the ALT file, and `input` on stdin.
fn run(args: &[&str], input: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
        .args(args)
        .arg(alt())
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    child.wait_with_output().unwrap()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// Decompress a BGZF file.
fn inflate(bgzf: &[u8]) -> Vec<u8> {
    let mut raw = Vec::new();
    BgzfStream::new(bgzf).read_to_end(&mut raw).unwrap();
    raw
}

/// Compress `raw` as BGZF.
fn deflate(raw: &[u8]) -> Vec<u8> {
    let mut writer = BgzfWriter::new(Vec::new());
    writer.write_all(raw).unwrap();
    writer.finish().unwrap()
}

/// A BAM of `body` (SAM lines) under `HEADER`, via the executable's own writer.
fn bam(body: &str) -> Vec<u8> {
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("in.bam");
    let input = format!("{HEADER}{body}");
    let out = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
        .args(["-O", "bam", "-o"])
        .arg(&path)
        .arg(alt())
        .arg("-")
        .stdin(Stdio::piped())
        .spawn()
        .and_then(|mut child| {
            child
                .stdin
                .take()
                .unwrap()
                .write_all(input.as_bytes())
                .unwrap();
            child.wait()
        })
        .unwrap();
    assert!(out.success());
    std::fs::read(path).unwrap()
}

/// Replace the header text of an uncompressed BAM, keeping its dictionary.
fn with_header_text(raw: &[u8], text: &[u8]) -> Vec<u8> {
    let l_text = i32::from_le_bytes(raw[4..8].try_into().unwrap()) as usize;
    let mut out = b"BAM\x01".to_vec();
    out.extend_from_slice(&i32::try_from(text.len()).unwrap().to_le_bytes());
    out.extend_from_slice(text);
    out.extend_from_slice(&raw[8 + l_text..]);
    out
}

#[test]
fn header_lines_among_alignments_pass_through_like_r985() {
    let input = format!("{HEADER}{}@CO\tmid-stream\n{}", read("r1"), read("r2"));
    let output = run(&["-j", "2"], input.as_bytes());
    assert!(output.status.success(), "{}", stderr(&output));
    // r985 prints a header line at once, while the group before it is
    // still buffered.
    let expected = format!("{HEADER}@CO\tmid-stream\n{}{}", read("r1"), read("r2"));
    assert_eq!(String::from_utf8(output.stdout).unwrap(), expected);

    let output = run(&["-O", "bam"], input.as_bytes());
    assert!(!output.status.success());
    assert!(
        stderr(&output).contains("cannot be written to BAM"),
        "{}",
        stderr(&output)
    );
}

#[test]
fn extreme_xa_values_are_errors_not_wrapped_mapqs() {
    let seq = "A".repeat(100);
    let qual = "I".repeat(100);
    for xa in [
        "altB,+101,100M,6823864795785271204;chr1,+5000,100M,0;",
        "altA,+9223372036854775807,100M,0;",
        "altA,+-50,100M,0;",
        "altA,101,100M,0;",
    ] {
        let input = format!(
            "{HEADER}r1\t0\taltA\t101\t60\t100M\t*\t0\t0\t{seq}\t{qual}\tNM:i:0\tXA:Z:{xa}\n"
        );
        let output = run(&[], input.as_bytes());
        assert!(!output.status.success(), "{xa}");
        assert!(
            stderr(&output).contains("read r1: invalid XA"),
            "{xa}: {}",
            stderr(&output)
        );
    }
}

#[test]
fn reserved_cigar_ops_in_xa_do_not_panic() {
    let seq = "A".repeat(100);
    let qual = "I".repeat(100);
    let input = format!(
        "{HEADER}r1\t0\taltA\t101\t60\t100M\t*\t0\t0\t{seq}\t{qual}\tNM:i:0\tXA:Z:chr1,+500,50M5B50M,1;\n"
    );
    let output = run(&[], input.as_bytes());
    assert!(
        output.status.code().is_some_and(|c| c <= 1),
        "{}",
        stderr(&output)
    );
    assert!(!stderr(&output).contains("panicked"), "{}", stderr(&output));
}

#[test]
fn bam_header_text_is_printed_as_samtools_prints_it() {
    let raw = inflate(&bam(&read("r1")));
    for text in [
        &b"@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:altA\tLN:2000\n@SQ\tSN:altB\tLN:2000"[..],
        &b"@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:altA\tLN:2000\n@SQ\tSN:altB\tLN:2000\n\0\0\0\0"[..],
    ] {
        let output = run(&[], &deflate(&with_header_text(&raw, text)));
        assert!(output.status.success(), "{}", stderr(&output));
        assert_eq!(
            String::from_utf8(output.stdout.clone()).unwrap(),
            format!("{HEADER}{}", read("r1"))
        );
    }
}

#[test]
fn truncated_and_padded_bgzf_input() {
    let body: String = (0..3000).map(|i| read(&format!("r{i}"))).collect();
    let file = bam(&body);
    // Cut off the EOF marker block: the input was truncated.
    let output = run(&[], &file[..file.len() - 28]);
    assert!(!output.status.success());
    assert!(
        stderr(&output).contains("EOF marker"),
        "{}",
        stderr(&output)
    );
    // An empty block between two data blocks is not the end.
    let raw = inflate(&file);
    let half = raw.len() / 2;
    // The first half's EOF marker block becomes an empty block mid-stream.
    let mut padded = deflate(&raw[..half]);
    padded.extend_from_slice(&deflate(&raw[half..]));
    let output = run(&["-j", "2"], &padded);
    assert!(output.status.success(), "{}", stderr(&output));
    let text = String::from_utf8(output.stdout).unwrap();
    assert_eq!(text, format!("{HEADER}{body}"));
}

#[test]
fn unsupported_inputs_are_named() {
    let output = run(&[], b"CRAM\x03\x00\x00\x00");
    assert!(
        stderr(&output).contains("CRAM input is not supported"),
        "{}",
        stderr(&output)
    );
    let raw = inflate(&bam(&read("r1")));
    let output = run(&[], &raw);
    assert!(
        stderr(&output).contains("uncompressed BAM"),
        "{}",
        stderr(&output)
    );
}

#[test]
fn early_stops_are_reported() {
    let input = format!("{HEADER}{}\n{}", read("r1"), read("r2"));
    let output = run(&[], input.as_bytes());
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout.clone()).unwrap(),
        format!("{HEADER}{}", read("r1"))
    );
    assert!(
        stderr(&output).contains("because it is blank"),
        "{}",
        stderr(&output)
    );
}

#[test]
fn outputs_must_not_alias_inputs_or_each_other() {
    let temp = tempfile::tempdir().unwrap();
    let prefix = temp.path().join("pre");
    let fastq = format!("{}.HLA-A.fq", prefix.display());
    let output = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
        .arg("-p")
        .arg(&prefix)
        .args(["-o", &fastq])
        .arg(alt())
        .arg("-")
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(
        stderr(&output).contains("is also the output"),
        "{}",
        stderr(&output)
    );
}

#[test]
fn thread_counts_are_bounded_and_options_named() {
    for args in [
        &["-j", "18446744073709551615"][..],
        &["-@", "100000"][..],
        &["-j", "x"][..],
    ] {
        let output = run(args, b"");
        assert!(!output.status.success(), "{args:?}");
        let err = stderr(&output);
        assert!(
            err.contains(args[0]) || err.contains("at most") || err.contains("between"),
            "{err}"
        );
    }
}
