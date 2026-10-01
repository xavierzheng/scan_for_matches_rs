# scan_for_matches (Rust port)

A Rust port of **scan_for_matches** by Ross Overbeek (Argonne National
Laboratory): scan nucleotide or protein sequences in FASTA format for
patterns (ranges, hairpins/reverse complements with pairing rules, repeats,
inexact matches, weight matrices, alternatives, length limits, `^`/`$`,
`any()`/`notany()` for proteins).

The goal of the port is **byte-identical results** with the original C
program: same hits, same output format, same exit status, same error
messages, same handling of odd input. The pattern language is described in
[`README.original`](README.original) (the original documentation).

## Build

```sh
cargo build --release
# binary: target/release/scan_for_matches
```

No third-party crates are used.

## Usage

Same command line as the original:

```sh
scan_for_matches [-c] [-p] [-o N] [-i ids_to_ignore] pattern_file < fasta_input > hits
```

| option | meaning |
|---|---|
| `-c` | also search the opposite strand (reverse complement) |
| `-p` | protein sequences |
| `-o N` | show overlapping hits (the value is not used) |
| `-i file` | file of sequence ids to skip |
| `-n`, `-m` | listed in the usage text, but in the original `getopt` string they take no argument, so the C program crashes (SIGSEGV) when they are given; the port does the same |

Example (from the original README):

```sh
$ echo 'p1=4...7 3...8 ~p1' > pat_file
$ scan_for_matches pat_file < tests/data/test_dna_input
>tst1:[6,27]
cguaacc ggttaacc gguuacg 
>tst2:[6,27]
CGUAACC GGTTAACC GGUUACG 
```

## How exact is it?

The port follows `scan_for_matches.c` and `ggpunit.c` statement by
statement, including their quirks, for example:

* option parsing uses the C library `getopt` (so messages and argument
  order rules are the platform's), `-n`/`-m` crash as described above;
* FASTA reading: the first record must start with `>`; only spaces and
  newlines are removed from sequences (`\r`, tabs, digits are kept); a
  header line at end of file without a newline makes the C program loop
  forever, and the port too;
* a pattern whose total maximum length is 0 (for example `^` alone) is
  reported as "failed to parse pattern";
* state that the C code keeps between sequences (last hit of each named
  unit, backtracking register, `past_last`) is kept the same way;
* `compl()` maps `s` to `S`, `-c` with `-p` runs the DNA complement code on
  protein data, bytes ≥ 0x80 index the code table with a signed `char`.

The original has fixed-size arrays (100 pattern units, 600 pattern code
bytes, 10500 weights, 100-entry stacks in the matcher, a 250 MB sequence
buffer). When a pattern or input is too big for them, the C program writes
past the end of the arrays: sometimes it silently changes the pattern,
sometimes it aborts (stack protector) or crashes. The port reproduces this
too: it keeps the C program's static data in one memory block laid out
exactly like the reference build, uses C-sized stack limits, and the same
buffer bounds. So results also match for over-limit patterns, `p50`/`r50`
(which alias other globals in C), negative ranges, and similar cases.

**Reference build.** Where C behaviour depends on memory layout, the port
follows the original built on macOS arm64 with Apple clang:
`cc -std=gnu89 -O2` (the README's `-O` and the Makefile's `-g -O2` give the
same layout). `tests/build_reference.sh` builds it.

**What cannot be identical.** In a few undefined-behaviour cases the C
program itself gives different results from run to run, because they depend
on address-space randomisation. No port can match those; they are:

* sequence ids of 1000 or more characters (the C `id[1000]` buffer
  overflows into a pointer, and the printed id then ends in random bytes);
* data bytes 0xA3/0xA4 in DNA mode, and bytes 0x80–0x9F when pattern names
  `p46`–`p49` are used (they read bytes of heap/static pointers);
* patterns with more than 100 units whose overflowed units overlay the
  pattern codes, in some cases (pointer bytes become pattern codes).

## Tests

```sh
cargo test --release
```

* `tests/original_suite.rs` – the test suite shipped with the C program
  (`run_tests`, `test_output`).
* `tests/golden.rs` – 788 cases recorded from the C reference
  (`tests/make_golden.py`): edge cases, over-limit patterns, options, odd
  bytes and random patterns; stdout and exit status must match.

Differential testing against the C program:

```sh
tests/build_reference.sh ../scan_for_matches_original target/reference/scan_for_matches
python3 tests/fuzz_compare.py target/reference/scan_for_matches target/release/scan_for_matches 10000 1
```

`fuzz_compare.py` generates random patterns (including stress patterns near
the C limits), FASTA input and options, and compares stdout, stderr and
exit status. A difference is reported as `C-unstable` when the C program
itself does not repeat its result. During development more than 100 000
cases were run with no difference outside that class.

## Speed

Typical patterns run at about the speed of the C program (from 0.8× to
1.3× of its time, depending on the pattern; exact words, inexact words,
weight matrices and alternatives are faster, heavy range backtracking is
slightly slower).

## Layout

* `src/main.rs` – port of `scan_for_matches.c` (options, FASTA, output)
* `src/engine.rs` – port of `ggpunit.c` (pattern parser and matcher)
* `src/sys.rs` – the few C library calls used for exact behaviour
  (`getopt`, `sscanf`, stdio output, signals, `mmap`)
