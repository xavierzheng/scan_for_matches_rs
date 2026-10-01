//! One test per fix of version 0.1.0 (TODO groups A and B).  The expected
//! results follow from the pattern language, not from a recording.

use std::io::Write;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

struct Out {
    code: Option<i32>,
    stdout: String,
    stderr: String,
}

fn run_with(args: &[&str], pattern: &[u8], input: &[u8], files: &[(&str, &[u8])]) -> Out {
    let dir = std::env::temp_dir().join(format!(
        "sfm_fix_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    for (name, body) in files {
        std::fs::write(dir.join(name), body).unwrap();
    }
    let pat = dir.join("pattern");
    std::fs::write(&pat, pattern).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_scan_for_matches"))
        .current_dir(&dir)
        .args(args)
        .arg(&pat)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut stdin = child.stdin.take().unwrap();
    let input = input.to_vec();
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&input);
    });
    // the original hangs on some inputs: never wait for ever
    let start = Instant::now();
    while child.try_wait().unwrap().is_none() {
        if start.elapsed() > Duration::from_secs(60) {
            child.kill().unwrap();
            panic!("scan_for_matches did not finish");
        }
        std::thread::sleep(Duration::from_millis(5));
    }
    let out = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    Out {
        code: out.status.code(),
        stdout: String::from_utf8_lossy(&out.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&out.stderr).into_owned(),
    }
}

fn run(args: &[&str], pattern: &str, input: &[u8]) -> Out {
    run_with(args, format!("{pattern}\n").as_bytes(), input, &[])
}

fn random_dna(n: usize, seed: u64) -> String {
    let mut x = seed.max(1);
    (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"acgt"[(x % 4) as usize] as char
        })
        .collect()
}

fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'a' => 't',
            'c' => 'g',
            'g' => 'c',
            _ => 'a',
        })
        .collect()
}

/// Hit lines: (id, from, to, pieces).
fn hits(out: &Out) -> Vec<(String, i64, i64, String)> {
    let lines: Vec<&str> = out.stdout.lines().collect();
    lines
        .chunks(2)
        .map(|h| {
            let head = h[0].strip_prefix('>').unwrap();
            let (id, pos) = head.rsplit_once(":[").unwrap();
            let (a, b) = pos.trim_end_matches(']').split_once(',').unwrap();
            (
                id.to_string(),
                a.parse().unwrap(),
                b.parse().unwrap(),
                h[1].to_string(),
            )
        })
        .collect()
}

// ---- A: crashes, hangs, wrong results -------------------------------------

#[test]
fn option_n_stops_after_failed_sequences() {
    let out = run(&["-n", "2"], "GGGGGGGG", b">a\nACGT\n>b\nACGT\n>c\nACGT\n");
    assert_eq!(out.code, Some(1));
    assert_eq!(out.stderr, "exceeded limit of lines failing to match\n");
}

#[test]
fn option_m_limits_hits() {
    let out = run(&["-m", "2"], "A", b">x\nAAAAA\n>y\nAAAAA\n");
    assert_eq!(out.code, Some(0));
    assert_eq!(out.stdout, ">x:[1,1]\nA \n>x:[2,2]\nA \n");
}

#[test]
fn option_values_are_checked() {
    let out = run(&["-n", "x"], "A", b">x\nA\n");
    assert_eq!(out.code, Some(2));
    assert!(out.stderr.contains("invalid value on -n option"));
}

#[test]
fn header_at_end_of_file_without_newline() {
    let out = run(&[], "ACGT", b">x\nACGT\n>y");
    assert_eq!(out.code, Some(0));
    assert_eq!(out.stdout, ">x:[1,4]\nACGT \n");
}

#[test]
fn exactly_100_units() {
    let seq = random_dna(400, 7);
    let pattern = vec!["N"; 100].join(" ");
    let out = run(&[], &pattern, format!(">x\n{seq}\n").as_bytes());
    let h = hits(&out);
    assert_eq!(h.len(), 4);
    for (k, (_, a, b, pieces)) in h.iter().enumerate() {
        assert_eq!((*a, *b), (100 * k as i64 + 1, 100 * k as i64 + 100));
        let want: String = seq[100 * k..100 * k + 100]
            .chars()
            .map(|c| format!("{c} "))
            .collect();
        assert_eq!(pieces, &want);
    }
}

#[test]
fn more_than_100_units() {
    let seq = random_dna(400, 8);
    let pattern = vec!["1...1"; 150].join(" ");
    let out = run(&[], &pattern, format!(">x\n{seq}\n").as_bytes());
    let h = hits(&out);
    assert_eq!(h.len(), 2);
    assert_eq!((h[0].1, h[0].2, h[1].1, h[1].2), (1, 150, 151, 300));
}

#[test]
fn long_pattern_words_and_weights() {
    // more than 600 pattern code bytes, then a weight matrix
    let seq = random_dna(700, 9);
    let word = &seq[10..660];
    let out = run(
        &[],
        &format!("{word} {{(10,10,10,10)}} > 5"),
        format!(">x\n{seq}\n").as_bytes(),
    );
    let h = hits(&out);
    assert_eq!(h.len(), 1);
    assert_eq!((h[0].1, h[0].2), (11, 661));
}

#[test]
fn long_stem_with_mismatch() {
    let stem = random_dna(120, 10);
    let mut rc: Vec<char> = revcomp(&stem).chars().collect();
    rc[5] = if rc[5] == 'a' { 'c' } else { 'a' };
    let rc: String = rc.into_iter().collect();
    // flanks of "a" cannot pair with each other
    let seq = format!("{}{}tttttt{}{}", "a".repeat(20), stem, rc, "a".repeat(20));
    let out = run(
        &[],
        "p1=120...120 3...8 ~p1[1,0,0]",
        format!(">x\n{seq}\n").as_bytes(),
    );
    let h = hits(&out);
    assert_eq!(h.len(), 1);
    assert_eq!((h[0].1, h[0].2), (21, 266));
}

#[test]
fn many_open_choices() {
    // 130 mismatches, each one leaves a "delete" choice open
    let pattern = format!("{}[130,0,1]", "a".repeat(130));
    let seq = "c".repeat(130);
    let out = run(&[], &pattern, format!(">x\n{seq}\n").as_bytes());
    assert_eq!(out.code, Some(0));
    assert_eq!(out.stdout, format!(">x:[1,130]\n{seq} \n"));
}

#[test]
fn long_sequence_id() {
    let id = "i".repeat(1500);
    let out = run(&[], "cgt", format!(">{id}\nacgtacgt\n").as_bytes());
    assert_eq!(
        out.stdout,
        format!(">{id}:[2,4]\ncgt \n>{id}:[6,8]\ncgt \n")
    );
}

#[test]
fn large_ignore_list() {
    let mut ids: String = (0..25_000).map(|i| format!("id{i}\n")).collect();
    ids.push_str("skip\n");
    let out = run_with(
        &["-i", "ignore"],
        b"cgt\n",
        b">skip\nacgt\n>keep\nacgt\n",
        &[("ignore", ids.as_bytes())],
    );
    assert_eq!(out.code, Some(0));
    assert_eq!(out.stderr, "ignoring 25001 id(s)\n");
    assert_eq!(out.stdout, ">keep:[2,4]\ncgt \n");
}

#[test]
fn undefined_or_self_referring_names_are_pattern_errors() {
    for p in [
        "~p3",
        "p1=3...3 ~p2",
        "length(p4) < 3 p1=1...1",
        "p1=~p1",
        "p1=3...3 r7~p1",
        "p1=~p2 p2=~p1",
    ] {
        let out = run(&[], p, b">x\nacgt\n");
        assert_eq!(out.code, Some(1), "{p}");
        assert_eq!(
            out.stderr,
            format!("failed to parse pattern: {p} \n"),
            "{p}"
        );
    }
}

#[test]
fn name_and_rule_50_work_like_others() {
    let seq = format!(">x\n{}\n>y\n{}\n", random_dna(300, 13), random_dna(300, 14));
    let a = run(&["-c"], "p50=3...3 0...3 ~p50", seq.as_bytes());
    let b = run(&["-c"], "p1=3...3 0...3 ~p1", seq.as_bytes());
    assert!(!a.stdout.is_empty());
    assert_eq!(a.stdout, b.stdout);
    let a = run(
        &[],
        "r50={au,ua,gc,cg,gu,ug} p50=4...4 2...4 r50~p50",
        seq.as_bytes(),
    );
    let b = run(
        &[],
        "r1={au,ua,gc,cg,gu,ug} p1=4...4 2...4 r1~p1",
        seq.as_bytes(),
    );
    assert!(!a.stdout.is_empty());
    assert_eq!(a.stdout, b.stdout);
}

#[test]
fn negative_lengths_and_counts_are_pattern_errors() {
    for p in [
        "-3...5 AA",
        "2...-1",
        "AC[-1,0,0] GG",
        "p1=3...3 ~p1[0,-1,0]",
    ] {
        let out = run(&[], p, b">x\nacgt\n");
        assert_eq!(out.code, Some(1), "{p}");
    }
}

#[test]
fn more_inserts_than_letters_stay_inside_the_sequence() {
    let out = run(&[], "AC[0,5,0]", b">x\nacgtacgt\n");
    for (_, a, b, _) in hits(&out) {
        assert!((1..=8).contains(&a) && b <= 8, "hit [{a},{b}]");
    }
}

#[test]
fn mismatch_insert_and_delete_together() {
    // pattern acgtac against acgaac with t->a ... one mismatch, then a
    // deleted character: the search must try "delete" after "insert"
    let out = run(&[], "acgtacgt[1,1,1]", b">x\nggacgaacggtgg\n");
    let h = hits(&out);
    assert!(!h.is_empty(), "{}", out.stdout);
}

// ---- B: input problems ----------------------------------------------------

#[test]
fn windows_line_ends() {
    let unix = run(
        &["-c"],
        "p1=3...4 2...5 ~p1",
        b">x desc\nacgtacgtaacc\nggttaaccgg\n>y\nccccgggg\n",
    );
    let dos = run(
        &["-c"],
        "p1=3...4 2...5 ~p1",
        b">x desc\r\nacgtacgtaacc\r\nggttaaccgg\r\n>y\r\nccccgggg\r\n",
    );
    assert!(!unix.stdout.is_empty());
    assert_eq!(unix.stdout, dos.stdout);
    let pat_dos = run_with(&[], b"p1=3...4\r\n2...5 ~p1\r\n", b">y\nccccgggg\n", &[]);
    assert_eq!(pat_dos.code, Some(0));
    assert_eq!(pat_dos.stdout, ">y:[1,8]\nccc cg ggg \n");
}

#[test]
fn text_before_first_record() {
    let plain = run(&[], "cgt", b">x\nacgt\n>y\ncgtt\n");
    let lead = run(&[], "cgt", b"\n\n; a comment\n  \n>x\nacgt\n>y\ncgtt\n");
    assert_eq!(plain.stdout, ">x:[2,4]\ncgt \n>y:[1,3]\ncgt \n");
    assert_eq!(plain.stdout, lead.stdout);
}

#[test]
fn lower_case_protein() {
    let out = run(&["-p"], "ACD any(KR)", b">x\nmacdkw\n>y\nMACDKW\n");
    assert_eq!(out.stdout, ">x:[2,5]\nacd k \n>y:[2,5]\nACD K \n");
}

#[test]
fn complement_with_protein_is_rejected() {
    let out = run(&["-p", "-c"], "ACD", b">x\nACD\n");
    assert_eq!(out.code, Some(2));
    assert!(
        out.stderr
            .contains("-c (complementary strand) cannot be used with -p")
    );
}

#[test]
fn complement_keeps_lower_case_s() {
    // the reverse complement of "ast" is "ast"; C printed "aSt"
    let out = run(&["-c"], "a 1...1 t", b">x\nast\n");
    assert_eq!(out.stdout, ">x:[1,3]\na s t \n>x:[3,1]\na s t \n");
}

#[test]
fn bytes_above_0x7f_are_unknown() {
    // 0xA8 read a global in C (and matched "c"); here it is unknown, like "x"
    let a = run(&[], "c", b">x\naa\xa8tt\n");
    assert_eq!(a.stdout, "");
    let mut pattern = b"ac".to_vec();
    pattern.push(0xA8);
    pattern.push(b'\n');
    let b = run_with(&[], &pattern, b">x\nacc\n", &[]);
    assert_eq!(b.code, Some(1));
}
