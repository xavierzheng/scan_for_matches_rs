#!/usr/bin/env python3
"""Random check of -t N: output with N threads must equal output with one.

usage: fuzz_threads.py BINARY [N_CASES] [SEED] [THREADS]

FUZZ_STRICT_N=1: DNA cases also get --strict-n.

Random patterns (from fuzz_compare.py, normal and stress), FASTA with many
records of very different lengths, random options (-c, -o, -m, -n, -i,
-p).  The run with threads also gets small pieces (SFM_PIECE), so long
records are searched in pieces, and a random warm-up before each piece
(SFM_WARM).  Runs two cases at a time (so 2 x THREADS
search threads).
"""
import os
import random
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fuzz_compare as fz  # noqa: E402

BIN = os.path.abspath(sys.argv[1])
N = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 1
T = sys.argv[4] if len(sys.argv) > 4 else "4"
FMT = os.environ.get("FUZZ_FORMAT") == "1"
STRICT = os.environ.get("FUZZ_STRICT_N") == "1"


def run(args, pat, inp, d, piece=None, warm=None):
    pp = os.path.join(d, "pat")
    with open(pp, "wb") as f:
        f.write(pat)
    try:
        env = dict(os.environ)
        if piece:
            env["SFM_PIECE"] = str(piece)
        if warm is not None:
            env["SFM_WARM"] = str(warm)
        p = subprocess.run([BIN] + args + [pp], input=inp, capture_output=True, timeout=60, cwd=d, env=env)
    except subprocess.TimeoutExpired:
        return None
    return (p.returncode, p.stdout, p.stderr)


def one(case):
    r = random.Random(SEED * 1000003 + case)
    protein = r.random() < 0.2
    pat = (fz.Pat(r, protein).build() if r.random() > 0.1 else fz.stress_pattern(r, protein)).encode("latin1") + b"\n"
    fasta = "".join(fz.rnd_fasta(r, protein) for _ in range(r.randint(1, 8)))
    ids = [l.split()[0][1:] for l in fasta.split("\n") if l.startswith(">")]
    args = []
    if protein:
        args.append("-p")
    elif r.random() < 0.4:
        args.append("-c")
    if r.random() < 0.3:
        args += ["-o", "1"]
    if r.random() < 0.2:
        args += ["-m", str(r.randint(0, 10))]
    if r.random() < 0.2:
        args += ["-n", str(r.randint(1, 5))]
    if STRICT and not protein:
        args.append("--strict-n")
    if FMT:
        # FUZZ_FORMAT=1: the output formats (numbering, --dedup) must not
        # depend on -t either
        args += ["--format", r.choice(["gff3", "bed6", "bed12", "jsonl"])]
        if r.random() < 0.5:
            args.append("--dedup")
    with tempfile.TemporaryDirectory() as d:
        if r.random() < 0.2 and ids:
            with open(os.path.join(d, "ign"), "w") as f:
                f.write("\n".join(r.sample(ids, r.randint(1, len(ids)))) + "\n")
            args += ["-i", "ign"]
        inp = fasta.encode("latin1")
        a = run(args, pat, inp, d)
        piece = r.choice([None, 1, 2, 3, 7, 30, 200, 1000])
        warm = r.choice([None, 0, 1, 2, 5, 17, 100, 1000])
        b = run(["-t", T] + args, pat, inp, d, piece, warm)
    if a is None or b is None:
        return ("timeout", None)
    return ("ok" if a == b else "DIFF", (pat, args, fasta[:300], a[0], b[0]))


with ThreadPoolExecutor(2) as ex:
    res = list(ex.map(one, range(N)))
stats = {}
for k, _ in res:
    stats[k] = stats.get(k, 0) + 1
print("results:", stats)
for k, info in res:
    if k == "DIFF":
        print(info)
        break
sys.exit(1 if stats.get("DIFF") else 0)
