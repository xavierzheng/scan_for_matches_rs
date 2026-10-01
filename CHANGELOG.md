# Changelog

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
