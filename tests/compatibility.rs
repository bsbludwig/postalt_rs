// SPDX-License-Identifier: GPL-3.0-only
// Part of the BWA/bwakit postalt Rust port; see NOTICE.md for provenance.

//! Regression tests using independent k8-generated fixtures and I/O round trips.
//!
//! `expected.sam` uses repaired reversal; `k8.sam` retains installed-k8 behavior.
//! Golden FASTQs are shared because r985 extracts the first buffered input record,
//! not a generated reverse hit. Fixture regeneration is an explicit operation in
//! scripts/compatibility.py; these tests never update their own expected output.
use postalt_rs::io::{Input, Output};
use std::{
    fs,
    io::Write,
    path::PathBuf,
    process::{Command, Stdio},
};

/// Resolve fixtures independently of the test runner's working directory.
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}
/// Build a case-specific CLI invocation with isolated HLA destinations.
fn command(name: &str, temp: &std::path::Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_postalt_rs"));
    command.arg("-p").arg(temp.join("hla"));
    match name {
        "ratio_0_5" => {
            command.args(["-r", "0.5"]);
        }
        "ratio_capped" => {
            command.args(["-r", "2"]);
        }
        "ratio_negative" => {
            command.args(["-r", "-1"]);
        }
        _ => (),
    }
    command.arg(fixture("alt.sam"));
    command
}
/// Check both reversal modes and every independently generated HLA fixture,
/// single-threaded and with worker threads.
fn golden(name: &str) {
    for (compat, jobs) in [
        (false, "0"),
        (false, "1"),
        (false, "3"),
        (true, "0"),
        (true, "3"),
    ] {
        let temp = tempfile::tempdir().unwrap();
        let mut cmd = command(name, temp.path());
        cmd.args(["-j", jobs]);
        if compat {
            cmd.arg("--k8-1.2-compat");
        }
        let output = cmd.arg(fixture(name).join("input.sam")).output().unwrap();
        assert!(
            output.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let expected =
            fs::read(fixture(name).join(if compat { "k8.sam" } else { "expected.sam" })).unwrap();
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&expected),
            "{name}, compat={compat}, jobs={jobs}"
        );
        for entry in fs::read_dir(fixture(name)).unwrap() {
            let path = entry.unwrap().path();
            if path.extension().is_some_and(|x| x == "fq") {
                assert_eq!(
                    fs::read(temp.path().join(format!(
                        "hla.{}",
                        path.file_name().unwrap().to_str().unwrap()
                    )))
                    .unwrap(),
                    fs::read(&path).unwrap(),
                    "{name}: FASTQ {}",
                    path.display()
                );
            }
        }
    }
}
macro_rules! cases { ($($name:ident),* $(,)?) => { $(#[test] fn $name() { golden(stringify!($name)); })* }; }
cases!(
    passthrough,
    unmapped,
    hard_clip_reported,
    single_alt,
    single_unplaced_alt,
    same_locus,
    primary_capped,
    alt_floor,
    group_score_difference,
    other_group_best,
    tied_groups,
    reverse_alt,
    reverse_generated,
    reverse_reported,
    iupac_reverse,
    missing_sequence,
    missing_quality,
    multi_lift,
    insertion_projection,
    deletion_projection,
    skip_projection,
    clipped_start,
    clipped_end_quirk,
    reverse_clip,
    hard_clip_xa,
    nm_clamp_and_gaps,
    equal_x_ignored,
    no_nm,
    empty_xa,
    duplicate_tags,
    nonoverlap_pa,
    ratio_0_5,
    ratio_capped,
    ratio_negative,
    overlap_adjustment,
    opposite_strand,
    restore_om,
    zero_om_not_restored,
    current_cigar_overlap_quirk,
    read_end_boundary,
    name_boundary,
    hla_primary_hard_clip,
    hla_union_first_record,
    multiple_lifts_no_pa_reduction,
    transitive_overlap,
    touching_not_overlap,
    header_only,
    headerless,
    final_line_no_newline
);

#[test]
fn sam_stdin_and_explicit_stdout() {
    let mut child = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
        .arg(fixture("alt.sam"))
        .args(["-o", "-", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(fixture("same_locus/input.sam")).unwrap())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    assert_eq!(
        output.stdout,
        fs::read(fixture("same_locus/expected.sam")).unwrap()
    );
}

#[test]
fn bam_output_and_stream_input_preserve_records_and_aux_types() {
    use seqair::sam::stream::SamRecord;
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("input.bam");
    let source = fs::read_to_string(fixture("same_locus/input.sam")).unwrap();
    let header: String = source
        .lines()
        .filter(|l| l.starts_with('@'))
        .map(|s| format!("{s}\n"))
        .collect();
    let tags = "NM:i:0\tZA:A:x\tZZ:Z:a b\tZH:H:00FF\tZI:i:-2147483648\tZU:i:4294967295\tZF:f:0.5\tZD:d:0.125\tZE:B:i\tZQ:B:f\tBc:B:c,-128,0,127\tBC:B:C,0,255\tBs:B:s,-32768,32767\tBS:B:S,0,65535\tBi:B:i,-2147483648,2147483647\tBI:B:I,0,4294967295\tBf:B:f,0.5,-1.5";
    let records: Vec<SamRecord> = [
        format!(
            "types\t0\tchr1\t1\t255\t16M\tchr2\t2\t-100\t=ACMGRSVTWYHKDBN\tABCDEFGHIJKLMNOP\t{tags}"
        ),
        format!("unmapped\t4\t*\t0\t255\t*\t*\t0\t0\t=ACMGRSVTWYHKDBN\tABCDEFGHIJKLMNOP\t{tags}"),
        format!("missing\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\t{tags}"),
    ]
    .iter()
    .map(|line| SamRecord::parse(line.as_bytes()).unwrap())
    .collect();
    for threads in [1, 2] {
        let mut out = Output::create(Some(&path), true, &header, threads, 6).unwrap();
        for r in &records {
            out.write(r).unwrap();
        }
        out.finish().unwrap();
        let mut input = Input::open(Some(&path), threads - 1).unwrap();
        assert_eq!(input.header_text(), header);
        for original in &records {
            let actual = input.next_record().unwrap().unwrap();
            let mut a = Vec::new();
            let mut b = Vec::new();
            original.write_sam(&mut a).unwrap();
            actual.write_sam(&mut b).unwrap();
            assert_eq!(a, b);
        }
        assert!(input.next_record().unwrap().is_none());
    }
    // Also exercise the executable's BAM sniffing on a non-seekable pipe.
    let mut child = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
        .arg(fixture("alt.sam"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(path).unwrap())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    let mut expected = header.into_bytes();
    for r in &records {
        r.write_sam(&mut expected).unwrap();
    }
    assert_eq!(output.stdout, expected);
}

#[test]
fn bam_matches_sam_on_alt_processing() {
    for name in [
        "same_locus",
        "reverse_generated",
        "unmapped",
        "missing_sequence",
        "missing_quality",
        "duplicate_tags",
    ] {
        let temp = tempfile::tempdir().unwrap();
        let bam = temp.path().join("input.bam");
        let result = temp.path().join("output.bam");
        let mut input = Input::open(Some(&fixture(name).join("input.sam")), 0).unwrap();
        let mut out = Output::create(Some(&bam), true, input.header_text(), 1, 0).unwrap();
        while let Some(r) = input.next_record().unwrap() {
            out.write(&r).unwrap();
        }
        out.finish().unwrap();
        for jobs in ["0", "2"] {
            let output = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
                .args(["-O", "bam", "-@", "2", "-j", jobs, "-o"])
                .arg(&result)
                .arg(fixture("alt.sam"))
                .arg(&bam)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{name}: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            let mut input = Input::open(Some(&result), 0).unwrap();
            let mut actual = input.header_text().as_bytes().to_vec();
            while let Some(r) = input.next_record().unwrap() {
                r.write_sam(&mut actual).unwrap();
            }
            assert_eq!(
                String::from_utf8_lossy(&actual),
                fs::read_to_string(fixture(name).join("expected.sam")).unwrap(),
                "{name}, jobs={jobs}"
            );
            // BAM input decoded straight to SAM output must agree as well.
            let output = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
                .args(["-j", jobs])
                .arg(fixture("alt.sam"))
                .arg(&bam)
                .output()
                .unwrap();
            assert!(output.status.success());
            assert_eq!(
                String::from_utf8_lossy(&output.stdout),
                fs::read_to_string(fixture(name).join("expected.sam")).unwrap(),
                "{name}: BAM to SAM, jobs={jobs}"
            );
        }
    }
}

#[test]
fn invalid_inputs_fail_clearly() {
    for args in [
        vec![],
        vec!["--unknown"],
        vec!["-r"],
        vec!["--threads", "0", "missing.alt"],
        vec!["--compression", "13", "missing.alt"],
        vec!["-r", "NaN", "missing.alt"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
            .args(args)
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(!output.stderr.is_empty());
    }
    let temp = tempfile::tempdir().unwrap();
    let bad = temp.path().join("bad.sam");
    for text in [
        "broken\n",
        "q\tNaN\tchr1\t1\t10\t100M\t*\t0\t0\tA\tI\n",
        "q\t0\taltA\t1\t10\t100M\t*\t0\t0\tA\tI\tXA:Z:bad;\n",
    ] {
        fs::write(&bad, text).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
            .arg(fixture("alt.sam"))
            .arg(&bad)
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(!output.stderr.is_empty());
    }
}

#[test]
fn does_not_overwrite_input() {
    let temp = tempfile::tempdir().unwrap();
    let input = temp.path().join("reads.sam");
    fs::copy(fixture("same_locus/input.sam"), &input).unwrap();
    let before = fs::read(&input).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
        .arg("-o")
        .arg(&input)
        .arg(fixture("alt.sam"))
        .arg(&input)
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert_eq!(fs::read(&input).unwrap(), before);
}

#[test]
fn version_and_help() {
    let out = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
        .arg("-v")
        .output()
        .unwrap();
    assert!(out.status.success());
    assert_eq!(out.stdout, b"r985\n");
    let out = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
        .arg("--help")
        .output()
        .unwrap();
    assert!(out.status.success());
    assert!(String::from_utf8_lossy(&out.stdout).contains("--k8-1.2-compat"));
}

#[test]
fn gzip_sam_and_gzip_alt() {
    use flate2::{Compression, write::GzEncoder};
    let temp = tempfile::tempdir().unwrap();
    for (name, source) in [
        ("alt.gz", fixture("alt.sam")),
        ("input.gz", fixture("same_locus/input.sam")),
    ] {
        let mut encoder = GzEncoder::new(
            fs::File::create(temp.path().join(name)).unwrap(),
            Compression::default(),
        );
        encoder.write_all(&fs::read(source).unwrap()).unwrap();
        encoder.finish().unwrap();
    }
    let output = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
        .arg(temp.path().join("alt.gz"))
        .arg(temp.path().join("input.gz"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        fs::read(fixture("same_locus/expected.sam")).unwrap()
    );
}

#[test]
fn truncated_bam_is_an_error() {
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("bad.bam");
    let mut input = Input::open(Some(&fixture("same_locus/input.sam")), 0).unwrap();
    let mut output = Output::create(Some(&path), true, input.header_text(), 1, 0).unwrap();
    output
        .write(&input.next_record().unwrap().unwrap())
        .unwrap();
    output.finish().unwrap();
    let mut data = fs::read(&path).unwrap();
    data.truncate(data.len() / 2);
    fs::write(&path, data).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
        .arg(fixture("alt.sam"))
        .arg(path)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(!output.stderr.is_empty());
}

#[test]
fn format_detection_handles_single_byte_reads() {
    use std::io::{Cursor, Read};
    struct Fragmented(Cursor<Vec<u8>>);
    impl Read for Fragmented {
        fn read(&mut self, dst: &mut [u8]) -> std::io::Result<usize> {
            let len = dst.len().min(1);
            self.0.read(&mut dst[..len])
        }
    }
    let temp = tempfile::tempdir().unwrap();
    let bam = temp.path().join("input.bam");
    let sam = fs::read(fixture("same_locus/input.sam")).unwrap();
    let mut input = Input::from_reader(Box::new(Cursor::new(sam.clone())), 0).unwrap();
    let mut writer = Output::create(Some(&bam), true, input.header_text(), 1, 6).unwrap();
    writer
        .write(&input.next_record().unwrap().unwrap())
        .unwrap();
    writer.finish().unwrap();
    for bytes in [sam, fs::read(bam).unwrap()] {
        for threads in [0, 2] {
            let reader = Box::new(Fragmented(Cursor::new(bytes.clone())));
            let mut input = Input::from_reader(reader, threads).unwrap();
            assert_eq!(input.next_record().unwrap().unwrap().qname(), b"q");
            assert!(input.next_record().unwrap().is_none());
        }
    }
}

/// Run the executable on `input` with extra arguments, returning (success, stdout).
fn run_with(args: &[&str], input: &std::path::Path) -> (bool, Vec<u8>) {
    let output = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
        .args(args)
        .arg(fixture("alt.sam"))
        .arg(input)
        .output()
        .unwrap();
    (output.status.success(), output.stdout)
}

/// Prefix the QNAME of every alignment line (all lines in these bodies).
fn rename(body: &str, suffix: usize) -> String {
    body.lines()
        .map(|l| {
            let (name, rest) = l.split_once('\t').unwrap();
            format!("{name}_{suffix}\t{rest}\n")
        })
        .collect()
}

#[test]
fn many_batches_match_golden_output_for_every_job_count() {
    // Enough renamed copies of independent golden cases to span several
    // batches; each copy's expected output is the renamed golden output.
    let cases = [
        "same_locus",
        "reverse_generated",
        "overlap_adjustment",
        "restore_om",
        "read_end_boundary",
        "name_boundary",
        "multi_lift",
        "passthrough",
        "unmapped",
        "transitive_overlap",
    ];
    let split = |path: PathBuf| -> (String, String) {
        let text = fs::read_to_string(path).unwrap();
        let (header, body): (Vec<&str>, Vec<&str>) = text.lines().partition(|l| l.starts_with('@'));
        let join = |v: Vec<&str>| v.iter().map(|l| format!("{l}\n")).collect::<String>();
        (join(header), join(body))
    };
    let (header, _) = split(fixture("same_locus/input.sam"));
    let mut input = header.clone();
    let mut expected = header;
    let bodies: Vec<_> = cases
        .iter()
        .map(|case| {
            let input = split(fixture(case).join("input.sam")).1;
            (input, split(fixture(case).join("expected.sam")).1)
        })
        .collect();
    // Unique names per copy and case keep neighbouring groups separate.
    for copy in 0..1500 {
        for (i, (body, golden)) in bodies.iter().enumerate() {
            input += &rename(body, copy * cases.len() + i);
            expected += &rename(golden, copy * cases.len() + i);
        }
    }
    assert!(input.len() > 3 << 20, "input should span several batches");
    let temp = tempfile::tempdir().unwrap();
    let sam = temp.path().join("input.sam");
    fs::write(&sam, &input).unwrap();
    for jobs in ["0", "1", "4"] {
        let (ok, stdout) = run_with(&["-j", jobs], &sam);
        assert!(ok);
        assert!(stdout == expected.as_bytes(), "SAM, jobs={jobs}");
    }
    // The same through BAM input and output, decoded for comparison. BAM
    // normalizes some spellings (e.g. pa:f:0.500000), so expect the golden
    // output after the same round trip.
    let to_bam = |from: &std::path::Path, to: &std::path::Path| {
        let mut reader = Input::open(Some(from), 0).unwrap();
        let mut writer = Output::create(Some(to), true, reader.header_text(), 2, 1).unwrap();
        while let Some(r) = reader.next_record().unwrap() {
            writer.write(&r).unwrap();
        }
        writer.finish().unwrap();
    };
    let decode = |from: &std::path::Path| {
        let mut reader = Input::open(Some(from), 2).unwrap();
        let mut text = reader.header_text().as_bytes().to_vec();
        while let Some(r) = reader.next_record().unwrap() {
            r.write_sam(&mut text).unwrap();
        }
        text
    };
    let bam = temp.path().join("input.bam");
    to_bam(&sam, &bam);
    let golden = temp.path().join("expected.sam");
    let golden_bam = temp.path().join("expected.bam");
    fs::write(&golden, &expected).unwrap();
    to_bam(&golden, &golden_bam);
    let expected = decode(&golden_bam);
    for jobs in ["0", "4"] {
        let (ok, stdout) = run_with(&["-j", jobs], &bam);
        assert!(ok);
        assert!(stdout == expected, "BAM to SAM, jobs={jobs}");
        let out = temp.path().join(format!("out{jobs}.bam"));
        let (ok, _) = run_with(
            &[
                "-j",
                jobs,
                "-O",
                "bam",
                "-@",
                "2",
                "-o",
                out.to_str().unwrap(),
            ],
            &bam,
        );
        assert!(ok);
        assert!(decode(&out) == expected, "BAM to BAM, jobs={jobs}");
    }
}

#[test]
fn errors_keep_sequential_output_for_every_job_count() {
    // Groups flushed before a malformed record are written; the group before
    // it is not (its flush would have followed the record's successful parse).
    let body = fs::read_to_string(fixture("overlap_adjustment/input.sam")).unwrap();
    let (header, body): (Vec<&str>, Vec<&str>) = body.lines().partition(|l| l.starts_with('@'));
    let body: String = body.iter().map(|l| format!("{l}\n")).collect();
    let mut input: String = header.iter().map(|l| format!("{l}\n")).collect();
    for copy in 0..6000 {
        input += &rename(&body, copy);
    }
    let temp = tempfile::tempdir().unwrap();
    for bad in [
        "q\tNaN\tchr1\t1\t10\t100M\t*\t0\t0\tA\tI\n",
        "q\t0\taltA\t1\t10\t100M\t*\t0\t0\tA\tI\tXA:Z:bad;\n",
    ] {
        let path = temp.path().join("bad.sam");
        fs::write(&path, format!("{input}{bad}{}", rename(&body, 99999))).unwrap();
        let (ok, sequential) = run_with(&["-j", "0"], &path);
        assert!(!ok);
        assert!(sequential.len() > 1 << 20);
        for jobs in ["1", "4"] {
            let (ok, stdout) = run_with(&["-j", jobs], &path);
            assert!(!ok);
            assert!(stdout == sequential, "jobs={jobs}");
        }
    }
}

#[test]
fn stdin_input_with_existing_output_file() {
    let temp = tempfile::tempdir().unwrap();
    let out = temp.path().join("out.sam");
    fs::write(&out, "stale").unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_postalt_rs"))
        .arg("-o")
        .arg(&out)
        .arg(fixture("alt.sam"))
        .arg("-")
        .stdin(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(fixture("same_locus/input.sam")).unwrap())
        .unwrap();
    assert!(child.wait().unwrap().success());
    assert_eq!(
        fs::read(&out).unwrap(),
        fs::read(fixture("same_locus/expected.sam")).unwrap()
    );
}

#[test]
fn numeric_fields_are_written_in_canonical_spelling() {
    // r985 prints parseInt() values, so leading zeros and '+' disappear even
    // on records that pass through unchanged.
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("input.sam");
    fs::write(
        &path,
        "q\t+016\tchr1\t0100\t060\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n",
    )
    .unwrap();
    let (ok, stdout) = run_with(&[], &path);
    assert!(ok);
    assert_eq!(
        String::from_utf8_lossy(&stdout),
        "q\t16\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
    );
}
