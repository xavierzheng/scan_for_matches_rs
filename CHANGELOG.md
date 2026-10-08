# Changelog

## Unreleased

Less memory with `-t N`; the output does not change.

### Faster
- `-t N`: the two strands of a long record (letters and coded sequence)
  are made once and shared, read only, by all threads.  Before, each
  thread that searched a piece kept its own copy (about 2 bytes per base
  per thread: 300 Mb x 2 x 20 threads = 12 GB).
- `--format`: a worker sends the writer only the text of the units the
  format uses (BED: none; GFF3: the TSD/TIR parts; JSON lines: all).
  Before, every hit carried the text of every unit, also a 30 kb gap.
  At most 4 x N jobs are between the reader and the writer (finished
  pieces waited in memory behind a slow one).

### Fixed
- `--explain` / `--lint`: the 2nd and 3rd numbers of `[m,d,i]` were
  called "inserts" and "deletes"; they are deletions (a pattern letter
  missing in the data) and insertions (an extra letter in the data), as
  in README.original.  The hit length range of `--explain` used them the
  wrong way round.  The search itself did not change.

Maize (2.3 Gb, bgzip), `-t 20`, BED6, max RSS, same output as 0.4.0:
00_DTC 13.3 → 2.5 GB; 08_DTH_seed12 14.8 → 3.4 GB; 02_DTA_short5to7
16.0 → 5.1 GB; 02_DTA_short5to7 with `-o` (40 million hits) 97 GB, 222 s
→ 3.9 GB, 107 s.

## 0.4.0 — 2026-10-08

Several patterns in one run, pattern checks, a merge tool.  Without the
new options the output is byte-identical to 0.3.0 (only the version in
the GFF3 header line changes).

### New
- Several pattern files in one run (with `--format`): the input is read
  once, each record is searched by each pattern in turn.  The hits of
  each pattern equal a run with that pattern alone, for every `-t N`.
  Names: one counter for each Name prefix; a Name is never given twice
  (also when a prefix ends in a digit: `X` + 11 and `X1` + 1).  GFF3 header: `pattern=` for
  each file; JSON lines: a `"pattern"` field.  `-m` counts all hits.
  Before, extra file arguments were ignored; now they are patterns, and
  without `--format` they are an error (exit status 2).  With `-t N`,
  a thread holds a long record in one engine at a time (the others give
  the memory back).
- `--lint`: errors, warnings and notes for the traps of the pattern
  language (reversed range, name used before it is defined, symmetric
  pattern with `-c` and no `--dedup`, `-o -c`, slow wide gaps, text cut
  at 31 999 bytes, mismatches >= word length, bad labels, ...).  Exit
  status 1 when there is an error.
- `--explain`: each unit in plain words, the hit length range, labels.
- `--merge FILE.gff3 ...`: join the GFF3 of several runs (in priority
  order), remove the same element found twice (same span; or
  `--overlap F`, reciprocal), sort (seqid in natural order, then start),
  and give new genome-wide unique Names (`--name-prefix`,
  `--name-start`); `ID`/`Parent`/`Name` renamed, the kept element gets
  `Merged=file:old_Name`.  gzip input is accepted.  An ID used by two
  elements of one file (files joined with `cat`) is an error.
- `tools/jaspar2sfm.py`: a JASPAR matrix to a weight unit, with the
  cutoff for a p-value (exact distribution of the integer scores).
  Python standard library only.

### Fixed
- `-t N` with warm-up (0.3.0): a hit of length 0 at the last position of
  a piece could be printed twice (the writer searched the next piece
  again from the end of the last hit, before the piece).  Found by
  `fuzz_threads.py` (`p1=0...1 ~p1`, small pieces, `SFM_WARM`); TIR
  patterns have no hits of length 0.

B. napus genome (1.0 Gb, bgzip), 22 TIR patterns, `-c --dedup --format
gff3 -t 20`, Intel Xeon Gold 6238R: one run with the 22 files 648 s (max
RSS 6.4 GB); 22 runs (0.3.0) 762 s (max 4.9 GB).  The hits of each
pattern are the same.
`--merge` of that output (806 397 elements, 1 thread): 12.7 s, 2.5 GB;
53 890 elements with the same span removed (202 311 with `--overlap
0.9`).

Tests (HPC, Linux x86_64): `cargo test --release`; `fuzz_threads.py`
2 x 1000 and 500 with `--format`: no difference; `fuzz_compare.py
--compat` against the C program built on Linux: 20 000 cases, no
difference; `fuzz_skip.py` (0.1.0 vs 0.4.0) 800 + 1000 `--chain`: no
difference.

## 0.3.0 — 2026-10-08

Output formats, new options, faster threads.  Without the new options
the output (stdout, stderr, exit status) is byte-identical to 0.2.0.

### New
- `--format gff3|bed6|bed12|jsonl`: hits as GFF3, BED or JSON lines.
  The output is the same for every `-t N` (Names are given by the
  writer, in input order).
- Labels in `%` comments of the pattern file name the element and its
  parts: `%@element TYPE key=value ...` and `%@ TYPE key=value ...` at the
  end of a line (all units of that line = one feature).  The C program
  and older versions ignore them.  Bad labels: message and exit status 2
  (only with `--format`).
- GFF3: one `repeat_region` (when there are TSD labels), the element and
  one line for each labelled part, linked by `ID`/`Parent`; every line
  has the hit's unique `Name` (prefix + number, no fixed width),
  `Classification`, `Method` (default `structural`); `Sequence_ontology`
  when given or when the type is a known SO name; `TSD=`/`TIR=` on the
  element as in EDTA.
- `--name-prefix P`, `--name-start N`, `--type T`.
- `--dedup`: with `-c`, the reverse-strand copy of an element already
  found on the forward strand (same span) is dropped and the element gets
  strand `.`; `-m` counts the elements written.
- `-h` / `--help`: all options, with an example.
- `--input FILE` and `--output FILE` (`-`: stdin / stdout).  Without them
  the program reads stdin and writes stdout as before.  `-i` and `-o`
  keep their original meaning.
- `tir_scan_patterns/`: the 22 TIR candidate patterns (after Wicker et
  al. 2007), their labelled copies for `--format` (`labelled/`), and
  their documentation in English (`README.md`) and Traditional Chinese
  (`README_zh-TW.txt`).

### Faster
- bgzip input with `-t N` (N > 1): the bgzip blocks are decompressed by
  N threads (system zlib, no new library).  Before, one thread
  decompressed all the input, and the other threads waited for it.  A
  member that is not a bgzip block (plain gzip) and all the input after
  it are read by one thread, as before.  Output is the same for good
  input.  For damaged bgzip input, the hits printed before the error can
  differ from `-t 1` (a damaged block gives no data with N > 1).
- `-t N`: the k-mer index of a long record (one per strand, shared by
  the threads that search its pieces) is built by N threads.  Before,
  one thread built it and the others waited.  The index is the same.
- `-t N` without `-o`: patterns with long, dense hits (a wide gap and a
  hit almost everywhere) are fast with threads now.  Before, the writer
  thread searched most pieces again by itself: the chain of
  non-overlapping hits of a worker met the one-thread chain only after
  about 300 kb (measured), and pieces are 1 Mb.  Now, when the writer
  has to search more than 1/20 of the positions again, each worker
  starts its chain up to 2 Mb before its piece (warm-up; those hits are
  not printed), and new long records get pieces 4 times the warm-up.
  The writer still checks every piece, so the output is the same.
  Other patterns do not change (no warm-up).  `SFM_WARM=N` (tests)
  fixes the warm-up.

B. napus genome (1.0 Gb, bgzip), 22 TIR patterns, `-c --dedup --format
gff3 -t 20`, Intel Xeon Gold 6238R, outputs byte-identical (md5) to the
first 0.3.0 code (output formats only):

| | total of 22 patterns | 00_DTC | 02_DTA_short5to7 | 08_DTH_Tourist_seed8 |
|---|---|---|---|---|
| output formats only | 1546 s | 25.2 s | 337 s | 155 s |
| + bgzip threads, index threads | 1087 s | 9.6 s | 308 s | 160 s |
| + warm-up (0.3.0) | 762 s | 9.6 s | 106 s | 35 s |

Peak memory of the two dense patterns: 10 GB → 5 GB.

## 0.2.0 — 2026-10-03

Faster.  For all input that 0.1.0 reads, the output is byte-identical
to 0.1.0 with any number of threads; only the time changes.  22 TIR patterns of B. napus,
one thread: on 1 Mb of chromosome A1, 3963 s with 0.1.0 (C: 3162 s) and
3.2 s with 0.2.0, the same hits.  The whole genome (1.0 Gb, 22 patterns,
`-t 8`): 9.6 min, at most 5.6 GB of memory (0.1.0: estimated 44 CPU-days).

### Speed
- Gap lengths that cannot lead to a hit are skipped.  When a range
  (`a...b`) is followed by an exact word, an exact reverse complement
  (`~pN`), a reverse complement with mismatches only (`~pN[m,0,0]`) or an
  exact repeat (`pN`), only the lengths where that unit can match are
  tried, in the same order as before.  The fixed units after it (words,
  exact `~pN`, `pN` of names matched before the range) are checked at
  the same time, 8 bases at a time.  Long gaps use a 5-mer position index
  of each record (about 4 bytes per base, built once per record and
  strand).  Used for DNA patterns that do not use matches of earlier
  sequences.
- Reverse complement and the coding of the sequence use lookup tables.

### Threads (`-t N`, new option)
- Records are searched by N threads; the output is written in input
  order and is identical to one thread, including `-m`, `-n`, `-i` and
  read errors.  B. napus CDS (120 351 records): 42.8 s with 1 thread,
  9.2 s with 8 (before the gap skipping).
- A long record (2 Mb or more) is cut into pieces of 1 Mb of start
  positions that several threads search; the pieces are joined in order
  with the one-thread rule (the next hit starts at or after the end of
  the last one).  Not used for patterns that use matches of earlier
  sequences, that start with an alternative `( | )`, or that contain `^`.
  Each thread keeps a copy of the record it searches (about 2 bytes per
  base); the 5-mer index is shared.  Whole genome, 02_DTA_short5to7:
  349 s with 1 thread, 184 s with 2, 125 s with 4, 110 s with 8 (Apple
  M4, 4 fast and 6 slow cores; the search is limited by memory speed).
- Patterns that read a name before it is certainly matched (it is
  defined later, or only in one branch of an alternative) depend on
  earlier sequences; they are detected and searched on one thread.  The
  engine has no global state any more.

### Input
- FASTA input compressed with gzip or bgzip is read directly (detected by
  its first bytes; uses the system zlib library).  Damaged or truncated
  input stops with "gzip input: ..." and exit status 1; the output before
  the error does not depend on how the input arrives through a pipe.

### Checks
- `cargo test`; random patterns and sequences compared with the C program
  (`fuzz_compare.py --compat`, 28 000 cases in this version), with a build
  without gap skipping (`fuzz_skip.py`, also `--chain`, 8 100 cases), and
  between thread counts and piece sizes (`fuzz_threads.py`, 8 000 cases):
  no difference.

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
