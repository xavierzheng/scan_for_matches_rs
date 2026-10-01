# Changelog

## Unreleased (0.2.0)

- FASTA input compressed with gzip or bgzip is read directly (detected by
  its first bytes; uses the system zlib library).  Damaged or truncated
  input stops with "gzip input: ..." and exit status 1.

## 0.1.0 — 2026-10-01

Fixes the bugs and input problems of the original (`TODO.md`, groups A
and B).  There is no option to get the old behaviour; use version 0.0.0 for
that.  For all input that these fixes do not touch, results are still
byte-identical to the original C program (40 000 random cases compared).

### Bugs of the original (group A)

- `-n N` and `-m N` work.  The original's option string lacked `:`, so
  both options always crashed.
- A header line at the end of the input without a newline ends the input
  (the original looped for ever).
- No fixed limits: pattern units (C: 100; with exactly 100 the original
  silently changed the pattern), pattern code bytes (600), weights
  (10 500), units in one hit (100), length of a complemented stem with
  mismatches (100), open choices in an inexact match (100), sequence id
  length (1000), ids to ignore (20 000), sequence length (250 000 000; now
  2 147 483 645).  The original wrote past these arrays: wrong pattern,
  random output, abort or crash.
- Undefined names (`~p3` without `p3`), names that refer to themselves
  directly or through other names (`p1=~p1`), and undefined rule sets
  (`r7~p1` without `r7={...}`) are pattern errors (the original crashed).
- `p50` and `r50` work like the other names and rule sets (the original
  accepted them but had only 50 entries and overwrote other variables).
- Negative lengths (`-3...5`) and negative mismatch/insert/delete counts
  are pattern errors (the original read outside the sequence, or treated
  the count as unlimited).  An inexact word with more inserts than
  letters does not start past the end of the sequence.
- Inexact matches with mismatches, inserts and deletes together try the
  "delete" choice after a mismatch too.  The original lost that choice, so
  it missed some hits.  These patterns can find other or more hits, and
  can take longer.

### Input problems of the original (group B)

- Windows line ends: `\r` is removed from sequences and treated like a
  newline in pattern files.
- Text and blank lines before the first `>` are skipped (the original
  gave no output).
- Protein data is matched in upper case, so lowercase data matches; the
  output shows the data as it is.
- `-c` together with `-p` is an error (exit status 2).  The original ran
  the DNA reverse complement on protein data.
- With `-c`, a lowercase `s` stays `s` (the original printed `S`).
- Bytes 0x80 and above are unknown characters, in data and in patterns
  (the original read other variables through a negative table index).

### Tests

- `tests/fixes.rs`: one test per fix.
- Unit test of the inexact matcher against an independent recursive
  implementation of its search order.
- Golden cases split into `golden.tsv` (equal to C) and `fixed.tsv`
  (changed by the fixes); `tests/split_golden.py`.
- `fuzz_compare.py --compat`: random input that no fix touches, compared
  with the C program.

## 0.0.0 — 2026-10-01

First version. A complete Rust port of `scan_for_matches` (Ross Overbeek,
Argonne National Laboratory) that **reproduces the original C program
exactly**.

- Same pattern language, command line options, output format, exit
  status and error messages as the original.
- Byte-identical results to the reference build of the C program
  (`cc -std=gnu89 -O2`, Apple clang, macOS arm64), including its bugs and
  its behaviour past its fixed limits (crashes, aborts, changed patterns).
  See `README.md`.
- Only cases where the C program itself gives different results from run
  to run (address-space randomisation) can differ; they are listed in
  `README.md`.
- Tests: the original test suite and 788 cases recorded from the C program
  (`cargo test --release`); differential fuzzer against the C program
  (`tests/fuzz_compare.py`), more than 100 000 cases run.
- The known problems of the original are kept on purpose and listed in
  `TODO.md`.
