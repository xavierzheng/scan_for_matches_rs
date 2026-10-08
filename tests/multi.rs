//! Several pattern files in one run (`--format` only): the hits of each
//! pattern equal the hits of a run with that pattern alone, for any `-t`.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

type Out = (Option<i32>, Vec<u8>, Vec<u8>);

fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("sfm_multi_{}_{name}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(dir: &Path, args: &[&str], input: &[u8], env: &[(&str, &str)]) -> Out {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_scan_for_matches"));
    cmd.env_remove("SFM_PIECE").env_remove("SFM_WARM");
    for (k, v) in env {
        cmd.env(k, v);
    }
    let mut child = cmd
        .current_dir(dir)
        .args(args)
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

/// Records of very different lengths, some long (searched in pieces).
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
        let len = match rnd(8) {
            0 => 30_000 + rnd(40_000),
            1 | 2 => rnd(40),
            _ => 200 + rnd(3000),
        } as usize;
        out.push_str(&format!(">r{r}\n"));
        let s: String = (0..len)
            .map(|_| {
                let k = if rnd(60) == 0 { 9 } else { 8 };
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

const PATS: [(&str, &str); 4] = [
    (
        "a.pat",
        "%@element te Name=A\np1=3...3 %@ target_site_duplication\np2=4...5 %@ terminal_inverted_repeat\n10...60\n~p2 %@ terminal_inverted_repeat\np1 %@ target_site_duplication\n",
    ),
    ("b.pat", "p1=4...6 2...8 ~p1\n"),
    // does not split (alternative first)
    ("c.pat", "( ACG | TTA ) 0...3 GG\n"),
    ("d.pat", "p1=3...3 CA 5...30 TG p1\n"),
];

fn setup(name: &str) -> PathBuf {
    let d = dir(name);
    for (f, t) in PATS {
        std::fs::write(d.join(f), t).unwrap();
    }
    d
}

/// JSON lines without the "name" and "pattern" fields, per pattern file.
fn per_pattern(out: &[u8]) -> Vec<(String, String)> {
    let text = String::from_utf8(out.to_vec()).unwrap();
    text.lines()
        .map(|l| {
            let pat = l
                .split("\"pattern\":\"")
                .nth(1)
                .map(|r| r[..r.find('"').unwrap()].to_string())
                .unwrap_or_default();
            let mut rest = String::new();
            for part in l.split(',') {
                if part.starts_with("\"name\":") || part.starts_with("\"pattern\":") {
                    continue;
                }
                rest.push_str(part);
                rest.push(',');
            }
            (pat, rest)
        })
        .collect()
}

#[test]
fn same_hits_as_single_runs_for_any_threads() {
    let d = setup("same");
    let input = fasta(11, 60);
    let files: Vec<&str> = PATS.iter().map(|(f, _)| *f).collect();
    for strand in [
        &[][..],
        &["-c"][..],
        &["-c", "--dedup"][..],
        &["-o", "1"][..],
    ] {
        let mut args: Vec<&str> = vec!["--format", "jsonl"];
        args.extend_from_slice(strand);
        // single runs, in order
        let mut want: Vec<(String, String)> = Vec::new();
        for f in &files {
            let mut a = args.clone();
            a.push(f);
            let (st, out, err) = run(&d, &a, &input, &[]);
            assert_eq!(st, Some(0), "{}", String::from_utf8_lossy(&err));
            want.extend(
                per_pattern(&out)
                    .into_iter()
                    .map(|(_, h)| (f.to_string(), h)),
            );
        }
        let mut a = args.clone();
        a.extend_from_slice(&files);
        let (st, one, _) = run(&d, &a, &input, &[]);
        assert_eq!(st, Some(0));
        let mut got = per_pattern(&one);
        let mut w = want.clone();
        // order: record by record, pattern by pattern; compare as sets
        // per pattern and check the full order across threads below
        got.sort();
        w.sort();
        assert_eq!(got, w, "{strand:?}");
        let mut t = vec!["-t", "4"];
        t.extend_from_slice(&a);
        for piece in ["700", "5000", "1000000"] {
            for warm in [None, Some("300")] {
                let mut env = vec![("SFM_PIECE", piece)];
                if let Some(w) = warm {
                    env.push(("SFM_WARM", w));
                }
                let (st, many, _) = run(&d, &t, &input, &env);
                assert_eq!(st, Some(0));
                assert_eq!(one, many, "{strand:?} piece {piece} warm {warm:?}");
            }
        }
    }
}

#[test]
fn names_unique_and_header() {
    let d = setup("names");
    let input = fasta(5, 40);
    let (st, out, _) = run(
        &d,
        &[
            "-c", "--dedup", "--format", "gff3", "a.pat", "b.pat", "d.pat",
        ],
        &input,
        &[],
    );
    assert_eq!(st, Some(0));
    let text = String::from_utf8(out).unwrap();
    assert!(text.starts_with(&format!(
        "##gff-version 3\n# scan_for_matches {} pattern=a.pat pattern=b.pat pattern=d.pat\n",
        env!("CARGO_PKG_VERSION")
    )));
    // b and d share the prefix sfm: one counter
    let mut tops = std::collections::HashSet::new();
    let mut regions = std::collections::HashSet::new();
    for l in text.lines() {
        if let Some(r) = l.strip_prefix("##sequence-region ") {
            assert!(regions.insert(r.to_string()), "region twice: {r}");
        }
        if l.starts_with('#') || l.contains("Parent=") {
            continue;
        }
        let id = l.split("ID=").nth(1).unwrap().split(';').next().unwrap();
        assert!(tops.insert(id.to_string()), "Name twice: {id}");
    }
    assert!(tops.iter().any(|n| n.starts_with('A')));
    assert!(tops.iter().any(|n| n.starts_with("sfm")));
}

#[test]
fn names_unique_when_a_prefix_ends_in_a_digit() {
    let d = setup("digit");
    std::fs::write(d.join("x.pat"), "%@element te Name=X\np1=4...6 2...8 ~p1\n").unwrap();
    std::fs::write(
        d.join("x1.pat"),
        "%@element te Name=X1\np1=4...6 2...8 ~p1\n",
    )
    .unwrap();
    let input = fasta(9, 30);
    let (st, out, _) = run(&d, &["--format", "bed6", "x.pat", "x1.pat"], &input, &[]);
    assert_eq!(st, Some(0));
    let text = String::from_utf8(out).unwrap();
    let mut names = std::collections::HashSet::new();
    for l in text.lines() {
        let n = l.split('\t').nth(3).unwrap();
        assert!(names.insert(n.to_string()), "Name twice: {n}");
    }
    assert!(names.len() > 20);
}

#[test]
fn limits_and_errors() {
    let d = setup("limits");
    let input = fasta(3, 30);
    // -m counts the hits of all patterns
    for t in ["1", "3"] {
        let (st, out, _) = run(
            &d,
            &["-t", t, "-m", "7", "--format", "bed6", "b.pat", "d.pat"],
            &input,
            &[("SFM_PIECE", "900")],
        );
        assert_eq!(st, Some(0));
        assert_eq!(out.iter().filter(|&&c| c == b'\n').count(), 7);
    }
    // plain output cannot tell the patterns apart
    let (st, _, err) = run(&d, &["b.pat", "d.pat"], &input, &[]);
    assert_eq!(st, Some(2));
    assert!(String::from_utf8_lossy(&err).contains("need --format"));
    // a missing second pattern file
    let (st, _, err) = run(&d, &["--format", "bed6", "b.pat", "nope.pat"], &input, &[]);
    assert_eq!(st, Some(2));
    assert!(String::from_utf8_lossy(&err).contains("nope.pat"));
    // a pattern that does not parse names its file
    std::fs::write(d.join("bad.pat"), "p1=3...3 ~p9\n").unwrap();
    let (st, _, err) = run(&d, &["--format", "bed6", "b.pat", "bad.pat"], &input, &[]);
    assert_eq!(st, Some(1));
    assert!(String::from_utf8_lossy(&err).contains("bad.pat"));
}
