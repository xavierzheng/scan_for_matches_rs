#!/usr/bin/env python3
"""Split golden cases recorded from the C program (make_golden.py) by the
behaviour of this version.

usage: split_golden.py RUST_BINARY < recorded.tsv

Cases where RUST_BINARY gives the same stdout and status as the C program
are written to tests/data/golden.tsv.  The others are cases changed by a fix
(see CHANGELOG.md); they are written to tests/data/fixed.tsv with this
version's result, after you have checked that every difference comes from a
fix.  Cases that take longer than 10 s (for example more than 100 optional
units, where C aborted and the full search is exponential) are left out.
"""
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

BIN = sys.argv[1]
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
lines = [l for l in sys.stdin.read().split("\n") if l]


def run(line):
    f = line.split("\t")
    pat, a, inp = bytes.fromhex(f[0]), bytes.fromhex(f[1]), bytes.fromhex(f[2])
    args = [x.decode() for x in a.split(b"\0")] if a else []
    with tempfile.TemporaryDirectory() as d:
        pp = os.path.join(d, "pat")
        with open(pp, "wb") as fh:
            fh.write(pat)
        try:
            p = subprocess.run([BIN] + args + [pp], input=inp, capture_output=True, timeout=10)
        except subprocess.TimeoutExpired:
            return None
    st = str(p.returncode) if p.returncode >= 0 else "sig%d" % -p.returncode
    return f, st, p.stdout


with ThreadPoolExecutor(min(8, os.cpu_count())) as ex:
    results = list(ex.map(run, lines))

same, fixed, slow = [], [], 0
for r in results:
    if r is None:
        slow += 1
        continue
    f, st, out = r
    if st == f[3] and out.hex() == f[4]:
        same.append("\t".join(f))
    else:
        fixed.append("\t".join(f[:3] + [st, out.hex()]))
with open(os.path.join(DATA, "golden.tsv"), "w") as fh:
    fh.write("\n".join(same) + "\n")
with open(os.path.join(DATA, "fixed.tsv"), "w") as fh:
    fh.write("\n".join(fixed) + "\n")
print("same as C: %d, changed by fixes: %d, left out (slow): %d" % (len(same), len(fixed), slow), file=sys.stderr)
