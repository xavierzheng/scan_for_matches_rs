# scan_for_matches (Rust port)

[![CI](https://github.com/xavierzheng/scan_for_matches_rs/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/xavierzheng/scan_for_matches_rs/actions/workflows/ci.yml)

A Rust port of **scan_for_matches** by Ross Overbeek (Argonne National
Laboratory): scan nucleotide or protein sequences in FASTA format for
patterns (ranges, hairpins/reverse complements with pairing rules, repeats,
inexact matches, weight matrices, alternatives, length limits, `^`/`$`,
`any()`/`notany()` for proteins).

The pattern language is described in [`README.original`](README.original)
(the original documentation).

## Versions

* **0.0.0** – exact reproduction of the original C program, including its
  bugs, crashes and fixed limits.
* **0.1.0** (this version) – the bugs and input problems of the original
  are fixed (see [`CHANGELOG.md`](CHANGELOG.md)).  For all other input the
  results are byte-identical to the original C program: same hits, same
  output format, same exit status, same messages.

## Build

```sh
cargo build --release
# binary: target/release/scan_for_matches
```

No third-party crates are used.

## Usage

Same command line as the original:

```sh
scan_for_matches [-c] [-p] [-n N] [-m N] [-o N] [-i ids_to_ignore] [-t N] pattern_file < fasta_input > hits
```

| option | meaning |
|---|---|
| `-c` | also search the opposite strand (reverse complement); not with `-p` |
| `-p` | protein sequences |
| `-n N` | stop (exit status 1) after N sequences without a hit |
| `-m N` | report at most N hits |
| `-o N` | show overlapping hits (the value is not used) |
| `-i file` | file of sequence ids to skip |
| `-t N` | use N threads: several records at the same time, and long records in pieces; default 1.  The output is the same as with one thread |

The FASTA input on stdin can be plain text or compressed with gzip or
bgzip (for example NCBI `*.fna.gz` files); compression is detected
automatically:

```sh
scan_for_matches -c pat_file < genome.fna.gz
```

Example (from the original README):

```sh
$ echo 'p1=4...7 3...8 ~p1' > pat_file
$ scan_for_matches pat_file < tests/data/test_dna_input
>tst1:[6,27]
cguaacc ggttaacc gguuacg 
>tst2:[6,27]
CGUAACC GGTTAACC GGUUACG 
```

## Differences from the original C program

Only these (details in `CHANGELOG.md`):

* no fixed limits: any number of pattern units, pattern codes and weights
  that fit on a pattern line, long hairpin stems and many mismatches,
  sequence ids of any length, any number of ids to ignore, sequences up to
  2 147 483 645 characters;
* `-n` and `-m` work (the original always crashed);
* undefined names (`~p3` without `p3`), names that refer to themselves
  (`p1=~p1`), undefined rule sets, negative lengths and negative
  mismatch/insert/delete counts are pattern errors ("failed to parse
  pattern", exit status 1) instead of crashes or reads outside memory;
* `p50` and `r50` work like the other names and rule sets;
* inexact matches with mismatches, inserts and deletes together try every
  choice (the original lost the "delete" choice after a mismatch, so it
  missed some hits);
* inexact words with more inserts than letters do not read past the end of
  the sequence;
* a header line at the end of the input without a newline ends the input
  (the original looped for ever);
* Windows line ends (`\r\n`) in sequences and patterns, and text or blank
  lines before the first `>`, are accepted;
* protein data is matched in upper case (lowercase data did not match);
  the output still shows the data as it is;
* `-c` with `-p` is an error (exit status 2);
* with `-c`, a lowercase `s` is shown as `s` (the original showed `S`);
* bytes 0x80 and above are unknown characters, like `x`;
* (after 0.1.0) gzip / bgzip compressed FASTA input is accepted;
* (after 0.1.0) `-t N` searches records in parallel.  A pattern that reads
  a name before it is certainly matched (defined later, or only in one
  branch of `( | )`) uses matches of earlier sequences, so it is always
  searched on one thread;
* (after 0.1.0) gap lengths that cannot lead to a hit are skipped (see
  Speed); the output does not change.

Unchanged on purpose (see `TODO.md`, group C): a pattern whose longest
match is 0 characters is a pattern error, the value of `-o` is not used, a
`length()` unit can see the last match of a name that is not matched yet.

## Tests

```sh
cargo test --release
```

* `tests/original_suite.rs` – the test suite shipped with the C program
  (`run_tests`, `test_output`).
* `tests/golden.rs` – cases recorded from the C reference program
  (`tests/make_golden.py`, split with `tests/split_golden.py`):
  `golden.tsv` (550 cases where this version equals C) and `fixed.tsv`
  (234 cases changed by the fixes).
* `tests/fixes.rs` – one test per fix, with expected results worked out
  from the pattern language.
* `tests/gzip_input.rs` – gzip and multi-member (bgzip) input, damaged
  input.
* `tests/threads.rs` – `-t 2/4/8` give the same output as `-t 1`
  (options, hit and miss limits, ignore list, damaged input, patterns
  that use earlier sequences); `tests/fuzz_threads.py` does the same with
  random patterns and input.
* `tests/fuzz_skip.py OLD NEW N SEED` – compares a build without gap
  skipping (for example 0.1.0) with a new one on random patterns with
  wide ranges followed by words, reverse complements and repeats.
* unit test in `src/engine.rs` – the inexact matcher against a separate
  recursive implementation of its search order (200 000 random cases).

Differential testing against the C program:

```sh
tests/build_reference.sh ../scan_for_matches_original target/reference/scan_for_matches
python3 tests/fuzz_compare.py --compat target/reference/scan_for_matches target/release/scan_for_matches 10000 1
```

With `--compat`, `fuzz_compare.py` generates only input that no fix
touches, so stdout, stderr and exit status must equal the C program's
(40 000 cases: no difference).  Without `--compat` it also generates
over-limit patterns, odd bytes and so on, as used for version 0.0.0.

The reference program is the original built on macOS arm64 with Apple
clang, `cc -std=gnu89 -O2` (the original does not compile with the default
flags of current compilers).

## Speed

Typical patterns run at about the speed of the C program (from 0.8× to
1.3× of its time, depending on the pattern).  Inexact matches with
mismatches, inserts and deletes together can take longer than in the
original, because the choices it skipped are now tried.

(After 0.1.0) A range followed by an exact word, an exact reverse
complement or an exact repeat, for example `p1=8...12 50...30000 ~p1`,
is much faster: only the gap lengths where the next unit can match are
tried (same order, same output).  For long gaps a 5-mer position index
of the record is used; it takes about 4 bytes per base for each thread
(a 74 Mb chromosome: about 300 MB per thread).  On 1 Mb of B. napus,
22 TIR patterns ran 100 to 500 times faster than before.

(After 0.1.0) With `-t N`, a long record (2 Mb or more) is cut into
pieces of 1 Mb of start positions, and the threads search the pieces at
the same time; the hits are joined in order, so the output is the same.
B. napus chromosome A1 (31 Mb) with `-t 8`: 3 to 4 times faster than
one thread.  More is not possible on that machine (Apple M4: 4 fast and
6 slow cores; the search is limited by memory speed).  Each thread keeps
a copy of the record (about 2 bytes per base); the 5-mer index of a
record is shared by all threads.  Whole B. napus genome, `-t 8 -c`:
at most 4.3 GB of memory.

## Layout

* `src/main.rs` – port of `scan_for_matches.c` (options, FASTA, output)
* `src/engine.rs` – port of `ggpunit.c` (pattern parser and matcher)
* `src/gz.rs` – gzip / bgzip input (system zlib)
* `src/sys.rs` – the few C library calls used (`getopt`, `sscanf`, stdio
  output, signals, `mmap`)
