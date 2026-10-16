//! Gap skipping for units after a wide range that were tried
//! length by length before: `pN[m,d,i]`, an inexact word, a weight unit
//! and a second range.  The hits must be the same as without skipping
//! (SFM_NO_SKIP=1).  `tests/fuzz_skip.py --slow` checks many more.

use std::io::Write;
use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_scan_for_matches");

/// stdout of `pat` on `fasta`.
fn run(pat: &str, args: &[&str], fasta: &str, env: &[(&str, &str)]) -> String {
    let d = std::env::temp_dir().join(format!("sfm_skipu_{}_{:x}", std::process::id(), {
        use std::hash::{BuildHasher, RandomState};
        RandomState::new().hash_one((pat, args, env))
    }));
    std::fs::create_dir_all(&d).unwrap();
    let pp = d.join("p.pat");
    std::fs::write(&pp, format!("{pat}\n")).unwrap();
    let mut child = Command::new(BIN)
        .args(args)
        .arg(&pp)
        .envs(env.iter().copied())
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(fasta.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    std::fs::remove_file(&pp).unwrap();
    std::fs::remove_dir(&d).unwrap();
    assert!(out.status.success());
    String::from_utf8(out.stdout).unwrap()
}

/// A random sequence with planted elements: a TSD, a gap, a unit to find
/// (with one change), the TSD again; some N and IUPAC codes.
fn fasta() -> String {
    let mut x: u64 = 11;
    let mut rnd = |n: u64| {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x % n
    };
    let mut seq = String::new();
    for k in 0..16 {
        let word = |n: u64, rnd: &mut dyn FnMut(u64) -> u64| -> String {
            (0..n).map(|_| b"ACGT"[rnd(4) as usize] as char).collect()
        };
        seq += &word(200 + rnd(400), &mut rnd);
        let tsd = word(4, &mut rnd);
        let rep = word(12, &mut rnd);
        let mut unit: Vec<char> = match k % 4 {
            0 => rep.chars().collect(),
            1 => "CACTACAAGAAA".chars().collect(),
            2 => "GACACCGC".chars().collect(),
            _ => "ACGTACGT".chars().collect(),
        };
        // one change: a mismatch, an N, a deletion or an insertion
        let at = rnd(unit.len() as u64) as usize;
        match rnd(5) {
            0 => unit[at] = 'T',
            1 => unit[at] = 'N',
            2 => {
                unit.remove(at);
            }
            3 => unit.insert(at, 'G'),
            _ => {}
        }
        let mid = word(50 + rnd(600), &mut rnd);
        let gap = if k % 7 == 0 { "NNNNN" } else { "" };
        seq += &format!(
            "{tsd}{rep}{mid}{gap}{}{tsd}",
            unit.iter().collect::<String>()
        );
    }
    format!(
        ">s1\n{seq}\n>s2\n{}\n",
        seq[seq.len() / 3..].replace("GATC", "GRYC")
    )
}

#[test]
fn new_skips_same_as_no_skip() {
    let fa = fasta();
    let pats = [
        // A: a direct repeat with errors (4-letter seeds, indel seeds)
        "p1=4...4 p2=12...12 50...700 p2[1,0,0] p1",
        "p1=4...4 p2=8...8 50...700 p2[1,0,0]",
        "p2=12...12 50...700 p2[1,1,0] p1=4...4",
        "p2=12...12 50...700 p2[2,1,1]",
        // B: an inexact word
        "p1=4...4 p2=12...12 50...700 CACTACAAGAAA[1,0,0] p1",
        "p1=4...4 50...700 CACTACAAGAAA[1,1,1] p1",
        "p1=4...4 50...700 ACGTACGT[0,0,2] p1",
        // C: a weight unit (and with an upper bound)
        "p1=4...4 50...700 {(16,0,84,0),(57,10,29,4),(0,80,0,20),(95,0,0,5),\
         (0,100,0,0),(18,60,20,2),(0,0,100,0),(0,50,50,0)} > 450 p1",
        "p1=4...4 50...700 600 > {(16,0,84,0),(57,10,29,4),(0,80,0,20),(95,0,0,5),\
         (0,100,0,0),(18,60,20,2),(0,0,100,0),(0,50,50,0)} > 450",
        // D: a second range, then a unit found by skipping
        "p1=4...4 50...700 0...30 ACGTACGT p1",
        "p1=4...4 p2=12...12 40...700 3...9 p2[1,0,0] p1",
        "p1=4...4 50...700 0...20 CACTACAAGAAA[1,0,0]",
    ];
    let mut total = 0;
    for pat in pats {
        // (-o with indels lists every alignment: slow without skipping)
        let o: &[&str] = if pat.contains(",1,") || pat.contains(",0,2]") {
            &["-c", "--strict-n"]
        } else {
            &["-o", "1"]
        };
        for args in [&[][..], &["-c"], o] {
            let a = run(pat, args, &fa, &[]);
            total += a.lines().filter(|l| l.starts_with('>')).count();
            assert_eq!(
                a,
                run(pat, args, &fa, &[("SFM_NO_SKIP", "1")]),
                "{pat} {args:?}"
            );
        }
    }
    assert!(total > 40, "{total}");
}

/// Pieces of a long record (`-t`): the starts of a word are found per
/// piece; the hits are those of one thread.
#[test]
fn new_skips_with_threads() {
    let fa = fasta();
    for pat in [
        "p1=4...4 50...700 CACTACAAGAAA[1,0,0] p1",
        "p1=4...4 p2=12...12 50...700 p2[1,0,0] p1",
        "p1=4...4 50...700 0...30 ACGTACGT[1,0,0] p1",
    ] {
        let one = run(pat, &["-c"], &fa, &[]);
        assert!(!one.is_empty(), "{pat}");
        assert_eq!(
            one,
            run(pat, &["-c", "-t", "3"], &fa, &[("SFM_PIECE", "700")]),
            "{pat}"
        );
    }
}
