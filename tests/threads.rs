//! `-t N`: the output with several threads must equal the output with one.

use std::io::Write;
use std::process::{Command, Stdio};

type Out = (Option<i32>, Vec<u8>, Vec<u8>);

fn scan(args: &[&str], pattern: &str, input: &[u8], files: &[(&str, &[u8])]) -> Out {
    scan_env(args, pattern, input, files, None)
}

/// As `scan`; `piece`: SFM_PIECE, the size of the pieces long records are
/// searched in (normally 1 Mb).
fn scan_env(
    args: &[&str],
    pattern: &str,
    input: &[u8],
    files: &[(&str, &[u8])],
    piece: Option<usize>,
) -> Out {
    let dir = std::env::temp_dir().join(format!(
        "sfm_thr_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    for (name, body) in files {
        std::fs::write(dir.join(name), body).unwrap();
    }
    let pat = dir.join("pattern");
    std::fs::write(&pat, format!("{pattern}\n")).unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_scan_for_matches"));
    match piece {
        Some(n) => cmd.env("SFM_PIECE", n.to_string()),
        None => cmd.env_remove("SFM_PIECE"),
    };
    let mut child = cmd
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
    let out = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    (out.status.code(), out.stdout, out.stderr)
}

/// Records of very different lengths (as in a genome: some long, many
/// short), some without hits.
fn fasta(seed: u64, records: usize) -> Vec<u8> {
    let mut x = seed;
    let mut rnd = move |n: u64| {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x % n
    };
    let mut out = String::new();
    for r in 0..records {
        let len = match rnd(10) {
            0 => 20_000 + rnd(30_000),
            1..=3 => rnd(30),
            _ => 200 + rnd(3000),
        } as usize;
        out.push_str(&format!(">r{r} record {r}\n"));
        let s: String = (0..len)
            .map(|_| {
                let k = if rnd(50) == 0 { 9 } else { 8 };
                b"acgtACGTn"[rnd(k) as usize] as char
            })
            .collect();
        for line in s.as_bytes().chunks(60) {
            out.push_str(std::str::from_utf8(line).unwrap());
            out.push('\n');
        }
    }
    out.into_bytes()
}

fn same_for_all_thread_counts(
    args: &[&str],
    pattern: &str,
    input: &[u8],
    files: &[(&str, &[u8])],
) -> Out {
    let one = scan(args, pattern, input, files);
    for t in ["2", "4", "8"] {
        let mut a = vec!["-t", t];
        a.extend_from_slice(args);
        let many = scan(&a, pattern, input, files);
        assert!(one == many, "-t {t} {args:?} {pattern}: output differs");
    }
    one
}

#[test]
fn patterns_and_options() {
    let input = fasta(11, 80);
    let patterns = [
        "p1=4...6 2...6 ~p1",
        "p1=3...3 0...4 p1[1,0,0]",
        "ACGT[1,0,1] 0...3 TTG",
        "{(10,0,0,10),(0,10,10,0),(10,10,0,0)} > 25",
        "(GAATTC | p1=4...4 1...3 ~p1) 2...5 RRY",
        "r1={au,ua,gc,cg,gu,ug} p1=5...5 3...6 r1~p1[1,0,0]",
    ];
    for p in patterns {
        for args in [
            &[][..],
            &["-c"][..],
            &["-o", "1"][..],
            &["-c", "-o", "1"][..],
        ] {
            let out = same_for_all_thread_counts(args, p, &input, &[]);
            assert_eq!(out.0, Some(0));
        }
    }
}

#[test]
fn hit_limit_and_miss_limit() {
    let input = fasta(12, 60);
    let p = "p1=4...6 2...6 ~p1";
    for m in ["0", "1", "7", "50", "100000"] {
        same_for_all_thread_counts(&["-m", m, "-c"], p, &input, &[]);
    }
    for n in ["1", "2", "5", "1000"] {
        let out = same_for_all_thread_counts(&["-n", n], "GGGGGGGGGG", &input, &[]);
        if n != "1000" {
            assert_eq!(out.0, Some(1));
        }
    }
    same_for_all_thread_counts(&["-n", "3", "-m", "20"], "p1=5...5 2...4 ~p1", &input, &[]);
}

#[test]
fn ignore_list() {
    let input = fasta(13, 40);
    let ids: String = (0..40).step_by(3).map(|r| format!("r{r}\n")).collect();
    same_for_all_thread_counts(
        &["-i", "ids"],
        "p1=4...6 2...6 ~p1",
        &input,
        &[("ids", ids.as_bytes())],
    );
}

#[test]
fn long_records_in_pieces() {
    // pieces of a few bases: hits across piece ends, hits over several
    // pieces, hits of length 0 (also just past the end), both strands
    let input = fasta(16, 30);
    let patterns = [
        "p1=4...6 2...6 ~p1",
        "p1=3...3 0...40 p1",
        "p1=0...1",
        "0...3 ACG",
        "ACGT[1,0,1] 0...3 TTG",
        "{(10,0,0,10),(0,10,10,0),(10,10,0,0)} > 25",
        "p1=4...4 0...30 ~p1 2...5 $",
        "p1=5...5 3...600 ~p1",
    ];
    for p in patterns {
        for args in [
            &[][..],
            &["-c"][..],
            &["-o", "1"][..],
            &["-c", "-m", "40"][..],
        ] {
            let one = scan(args, p, &input, &[]);
            for (t, piece) in [("2", 1), ("4", 3), ("8", 17), ("3", 400)] {
                let mut a = vec!["-t", t];
                a.extend_from_slice(args);
                let many = scan_env(&a, p, &input, &[], Some(piece));
                assert!(
                    one == many,
                    "-t {t}, pieces of {piece}, {args:?} {p}: output differs"
                );
            }
        }
    }
}

#[test]
fn damaged_gzip_in_the_middle() {
    let input = fasta(14, 60);
    let mut child = Command::new("gzip")
        .arg("-c")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut stdin = child.stdin.take().unwrap();
    let data = input.clone();
    let feeder = std::thread::spawn(move || stdin.write_all(&data).unwrap());
    let mut gz = child.wait_with_output().unwrap().stdout;
    feeder.join().unwrap();
    let n = gz.len();
    for b in &mut gz[n / 2..n / 2 + 50] {
        *b = 0x55;
    }
    let out = same_for_all_thread_counts(&["-c"], "p1=4...6 2...6 ~p1", &gz, &[]);
    assert_eq!(out.0, Some(1));
    assert!(!out.1.is_empty(), "hits before the damage are printed");
}

#[test]
fn pattern_using_earlier_sequences() {
    // ~p1 reads p1 before p1 is matched: the last match of p1 (from earlier
    // sequences) is used, so this pattern must run on one thread
    let input = fasta(15, 50);
    same_for_all_thread_counts(&[], "~p1 2...4 p1=3...3", &input, &[]);
    same_for_all_thread_counts(&["-c"], "(p1=4...4 | AC) 1...3 ~p1", &input, &[]);
}

#[test]
fn bad_thread_count() {
    for t in ["0", "x", "-2"] {
        let out = scan(&["-t", t], "ACGT", b">x\nACGT\n", &[]);
        assert_eq!(out.0, Some(2), "-t {t}");
        assert!(String::from_utf8_lossy(&out.2).contains("invalid value on -t option"));
    }
}

/// With warm-up, the writer searched a piece again from the end of the
/// last hit, which can be before the piece: a hit of length 0 at the last
/// position of the earlier piece was printed twice.
#[test]
fn warm_up_and_empty_hits() {
    let input = b">s0\nGGACCGaGcauaGuaAaGgcuAGcgVAGcgVuTgA\n>s1\ntGGgcgCtGGGacgtTTAAcg\n";
    let pat = "p1=0...1 ~p1\n";
    let dir = std::env::temp_dir().join(format!("sfm_thr_warm_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("pattern"), pat).unwrap();
    let run = |env: &[(&str, &str)], t: &str| {
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_scan_for_matches"));
        cmd.env_remove("SFM_PIECE").env_remove("SFM_WARM");
        for (k, v) in env {
            cmd.env(k, v);
        }
        let mut child = cmd
            .current_dir(&dir)
            .args(["-t", t, "-c", "pattern"])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        child.stdin.take().unwrap().write_all(input).unwrap();
        child.wait_with_output().unwrap().stdout
    };
    let one = run(&[], "1");
    for piece in ["1", "2", "3", "7"] {
        for warm in ["0", "1", "2", "5", "1000"] {
            let many = run(&[("SFM_PIECE", piece), ("SFM_WARM", warm)], "4");
            assert_eq!(one, many, "piece {piece} warm {warm}");
        }
    }
}
