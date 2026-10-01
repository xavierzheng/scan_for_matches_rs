#!/bin/sh
# Build the original C scan_for_matches as the reference for the
# differential tests.  The 1990s K&R code needs -std=gnu89 with current
# compilers; the other flags are the ones of its Makefile.
#
# usage: tests/build_reference.sh [path/to/scan_for_matches_original] [output]
set -e
SRC=${1:-../scan_for_matches_original}
OUT=${2:-target/reference/scan_for_matches}
mkdir -p "$(dirname "$OUT")"
${CC:-cc} -std=gnu89 -w -g -O2 -o "$OUT" "$SRC/scan_for_matches.c" "$SRC/ggpunit.c"
echo "$OUT"
