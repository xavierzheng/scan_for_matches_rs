//! --align: `~pN[m,d,i]` with deletions or insertions tries every length
//! it can align with; and the index seeds of such units (default search,
//! same hits as without skipping).

use std::io::Write;
use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_scan_for_matches");

/// Hits (`id:[a,b]`) of `pat` on `fasta`.
fn hits(pat: &str, args: &[&str], fasta: &str, env: &[(&str, &str)]) -> Vec<String> {
    let d = std::env::temp_dir().join(format!("sfm_align_{}_{:x}", std::process::id(), {
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
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter_map(|l| l.strip_prefix('>').map(str::to_string))
        .collect()
}

fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|c| match c {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect()
}

/// TSD GTA, TIR CCATTAGTCAGA, its reverse complement without the last G:
/// the greedy alignment takes the G of the TSD and the TSD then fails.
fn element() -> String {
    let tir = "CCATTAGTCAGA";
    let r = rc(tir);
    format!(
        ">s\nTTTTTGTA{tir}{}{}GTACCCCC\n",
        "A".repeat(30),
        &r[..r.len() - 1]
    )
}

#[test]
fn align_finds_what_greedy_misses() {
    let pat = "p1=3...3 p2=12...12 20...40 ~p2[0,1,1] p1";
    assert!(hits(pat, &[], &element(), &[]).is_empty());
    assert_eq!(hits(pat, &["--align"], &element(), &[]), ["s:[6,64]"]);
}

#[test]
fn align_same_without_indels() {
    let fa = element();
    for pat in [
        "p1=3...3 p2=12...12 20...40 ~p2[1,0,0] p1",
        "p1=3...3 p2=11...11 20...40 ~p2 p1",
    ] {
        assert_eq!(hits(pat, &[], &fa, &[]), hits(pat, &["--align"], &fa, &[]));
    }
}

/// Seeds for `~p2[2,1,1]` after a wide gap (blocks of 4 letters, also at
/// the start of a run after N), and Myers' filter when p2 is too short for
/// seeds (19 letters): the same hits as without skipping, with and without
/// --align.
#[test]
fn indel_seeds_same_as_no_skip() {
    for len in [20, 19] {
        indel_skip_case(len);
    }
}

fn indel_skip_case(len: usize) {
    let mut x: u64 = 7;
    let mut rnd = |n: u64| {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x % n
    };
    let mut seq = String::new();
    for k in 0..6 {
        let tir: String = (0..len).map(|_| b"ACGT"[rnd(4) as usize] as char).collect();
        let mut end: Vec<char> = rc(&tir).chars().collect();
        end.remove(3 + rnd(10) as usize);
        end.insert(5 + rnd(10) as usize, 'A');
        end[rnd(19) as usize] = 'C';
        let mid: String = (0..100 + rnd(500))
            .map(|_| b"ACGT"[rnd(4) as usize] as char)
            .collect();
        // a run of N right before a block
        let gap = if k % 2 == 0 { "NNNN" } else { "" };
        seq += &format!(
            "TTAGCATTG{tir}{mid}{gap}{}TTAGCATTG",
            end.iter().collect::<String>()
        );
    }
    let fa = format!(">s\n{seq}\n");
    let pat = format!("p1=9...9 p2={len}...{len} 50...700 ~p2[2,1,1] p1");
    for args in [&[][..], &["--align"], &["-c"], &["-c", "--align"]] {
        let a = hits(&pat, args, &fa, &[]);
        assert!(!a.is_empty(), "{len} {args:?}");
        assert_eq!(
            a,
            hits(&pat, args, &fa, &[("SFM_NO_SKIP", "1")]),
            "{len} {args:?}"
        );
    }
}
