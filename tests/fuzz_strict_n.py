#!/usr/bin/env python3
"""Random check of --strict-n against a small separate reference.

usage: fuzz_strict_n.py BINARY [N_CASES] [SEED]

Patterns of the TSD/TIR kind: names caught by ranges (`p1=3...6`), exact
words, gaps (some wide, so gap skipping is used), and `p1`, `<p1`, `~p1`
with 0-2 mismatches.  The data has runs of N and some IUPAC codes.  Each
case is run without and with --strict-n (forward strand, no -o); the hits
must equal the reference.  The run without the option checks the
reference itself.
"""
import os
import random
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

BIN = os.path.abspath(sys.argv[1])
N = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 1

BASES = "ACGT"
IUPAC = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def same(named, d, strict):
    """Does data letter d match the letter a name caught?"""
    if named in BASES:
        return named == d
    return not strict and d in IUPAC.get(named, "")


def unit_ends(u, seq, i, names, strict):
    """Ends of unit u starting at i, in the order the engine tries them."""
    kind = u[0]
    if kind == "range":
        _, name, a, b = u
        for e in range(i + a, min(i + b, len(seq)) + 1):
            if name:
                names[name] = seq[i:e]
            yield e
        return
    if kind == "word":
        w = u[1]
        if seq[i:i + len(w)] == w:
            yield i + len(w)
        return
    _, op, name, mis = u
    t = names[name]
    if i + len(t) > len(seq):
        return
    if op == "~":
        if any(c not in BASES for c in t):
            return  # as in C: a caught N never matches here
        t = "".join(COMP[c] for c in reversed(t))
    elif op == "<":
        t = t[::-1]
    for k, c in enumerate(t):
        d = seq[i + k]
        if d not in BASES:
            return
        if not (c == d if op == "~" else same(c, d, strict)):
            mis -= 1
            if mis < 0:
                return
    yield i + len(t)


def match_at(units, seq, i, names, strict):
    if not units:
        return i
    for e in unit_ends(units[0], seq, i, names, strict):
        r = match_at(units[1:], seq, e, names, strict)
        if r is not None:
            return r
    return None


def reference(units, recs, strict):
    out = []
    for rid, seq in recs:
        pos = 0
        while pos < len(seq):
            for s in range(pos, len(seq)):
                e = match_at(units, seq, s, {}, strict)
                if e is not None:
                    out.append(f"{rid}:[{s + 1},{e}]")
                    pos = e
                    break
            else:
                break
    return out


def rnd_pattern(r):
    units = []
    defined = []
    for k in range(r.randint(1, 2)):
        name = f"p{k + 1}"
        a = r.randint(1, 5)
        units.append(("range", name, a, a + r.choice([0, 0, 1, 3])))
        defined.append(name)
        if r.random() < 0.3:
            units.append(("word", "".join(r.choice(BASES) for _ in range(r.randint(1, 2)))))
    a = r.randint(0, 15)
    units.append(("range", None, a, a + r.choice([0, 5, 30, 200])))
    for name in reversed(defined):
        if r.random() < 0.3:
            units.append(("word", "".join(r.choice(BASES) for _ in range(r.randint(1, 2)))))
        units.append(("use", r.choice(["", "<", "~"]), name, r.choice([0, 0, 1, 2])))
        if name != defined[0]:
            a = r.randint(0, 8)
            units.append(("range", None, a, a + r.choice([0, 3, 40])))
    return units


def pattern_text(units):
    t = []
    for u in units:
        if u[0] == "range":
            t.append((f"{u[1]}=" if u[1] else "") + f"{u[2]}...{u[3]}")
        elif u[0] == "word":
            t.append(u[1])
        else:
            t.append(f"{u[1]}{u[2]}" + (f"[{u[3]},0,0]" if u[3] else ""))
    return " ".join(t) + "\n"


def rnd_seq(r):
    s = []
    n = r.randint(30, 600)
    while len(s) < n:
        x = r.random()
        if x < 0.04:
            s += "N" * r.randint(1, 12)
        elif x < 0.06:
            s.append(r.choice("RYSWKMBDHV"))
        elif x < 0.10 and len(s) > 12:
            # copy or reverse complement a few letters (TSD/TIR-like)
            k = r.randint(2, 8)
            part = s[-r.randint(k, 12):][:k]
            if r.random() < 0.5:
                part = [COMP.get(c, c) for c in reversed(part)]
            s += part
        else:
            s.append(r.choice(BASES))
    return "".join(s[:n])


def run(args, pat, fasta, d):
    pp = os.path.join(d, "pat")
    with open(pp, "w") as f:
        f.write(pat)
    try:
        p = subprocess.run([BIN] + args + [pp], input=fasta.encode(), capture_output=True, timeout=60)
    except subprocess.TimeoutExpired:
        return None
    return [l[1:] for l in p.stdout.decode().split("\n") if l.startswith(">")]


def one(case):
    r = random.Random(SEED * 1000003 + case)
    units = rnd_pattern(r)
    pat = pattern_text(units)
    recs = [(f"s{k}", rnd_seq(r)) for k in range(r.randint(1, 3))]
    fasta = "".join(f">{i}\n{s}\n" for i, s in recs)
    res = []
    with tempfile.TemporaryDirectory() as d:
        for strict in (False, True):
            got = run(["--strict-n"] if strict else [], pat, fasta, d)
            if got is None:
                return ("timeout", None)
            want = reference(units, recs, strict)
            if got != want:
                return ("DIFF", (strict, pat, fasta, got, want))
            res.append(len(want))
    return ("ok" if res[0] == res[1] else "ok-fewer", None)


with ThreadPoolExecutor(int(os.environ.get("FUZZ_WORKERS", "4"))) as ex:
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
