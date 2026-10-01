#!/usr/bin/env python3
"""Differential test: run the C original and the Rust port on random
patterns / FASTA inputs / options and compare stdout, stderr and exit code.

usage: fuzz_compare.py C_BINARY RUST_BINARY [N_CASES] [SEED]
"""
import os
import random
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

C_BIN, R_BIN = sys.argv[1], sys.argv[2]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 1
TIMEOUT = 10

DNA_SEQ = "acgtACGTuU"
DNA_AMB = "nNrRyYmMkKsSwWbBdDhHvV"
DNA_PAT = "acgtuACGTUnNrRyYmMkKsSwWbBdDhHvV"
AA = "ACDEFGHIKLMNPQRSTVWY"


def rnd_seq(r, protein):
    n = r.choice([0, 1, 3, 10, 30, 60, 120, 300])
    out = []
    for _ in range(n):
        x = r.random()
        if protein:
            if x < 0.85:
                out.append(r.choice(AA))
            elif x < 0.92:
                out.append(r.choice("XBZaclk*"))
            else:
                out.append(r.choice(" \n"))
        else:
            if x < 0.85:
                out.append(r.choice(DNA_SEQ))
            elif x < 0.93:
                out.append(r.choice(DNA_AMB))
            elif x < 0.97:
                out.append(r.choice(" \n"))
            else:
                out.append(r.choice("x-*.\r\t0"))
    # occasionally build an internal repeat / hairpin so complex patterns hit
    s = "".join(out)
    if s and r.random() < 0.5 and not protein:
        k = r.randint(2, 8)
        st = r.randint(0, max(0, len(s) - k))
        frag = s[st:st + k]
        comp = {"a": "t", "c": "g", "g": "c", "t": "a", "u": "a", "A": "T", "C": "G", "G": "C", "T": "A", "U": "A"}
        rc = "".join(comp.get(c, c) for c in reversed(frag))
        ins = r.choice([frag, rc, frag[::-1]])
        pos = r.randint(0, len(s))
        s = s[:pos] + "".join(r.choice(DNA_SEQ) for _ in range(r.randint(0, 6))) + ins + s[pos:]
    return s


def rnd_fasta(r, protein):
    recs = []
    for k in range(r.randint(1, 5)):
        hdr = ">" + r.choice(["s", "seq", "tst", "id_"]) + str(k)
        if r.random() < 0.3:
            hdr += " some description here"
        s = rnd_seq(r, protein)
        # wrap lines
        w = r.choice([10, 30, 60, 1000])
        lines = [s[i:i + w] for i in range(0, len(s), w)] or [""]
        recs.append(hdr + "\n" + "\n".join(lines) + "\n")
    txt = "".join(recs)
    if r.random() < 0.05:
        txt = txt.rstrip("\n")  # missing final newline (body)
    if r.random() < 0.03:
        txt = "\n" + txt  # leading blank line: C stops immediately
    return txt


class Pat:
    def __init__(self, r, protein):
        self.r = r
        self.protein = protein
        self.defined = []
        self.next_name = 1
        self.rules = []

    def mid(self):
        r = self.r
        if r.random() < 0.6:
            return ""
        return "[%d,%d,%d]" % (r.choice([0, 0, 1, 1, 2]), r.choice([0, 0, 0, 1, 2]), r.choice([0, 0, 0, 1, 2]))

    def word(self, n):
        if self.protein:
            return "".join(self.r.choice(AA + "X" + "acd") for _ in range(n))
        return "".join(self.r.choice(DNA_PAT if self.r.random() < 0.3 else "acgtACGT") for _ in range(n))

    def unit(self, depth):
        r = self.r
        kinds = ["range", "range", "exact", "sim", "repeat", "inv", "weight", "llim", "or", "start", "end"]
        if not self.protein:
            kinds += ["compl", "compl", "rcompl"]
        else:
            kinds += ["any", "notany", "any"]
        k = r.choice(kinds)
        name = ""
        if k in ("range", "exact", "sim", "weight", "any", "notany") and r.random() < 0.5 and self.next_name < 8:
            name = "p%d=" % self.next_name
            self.defined.append(self.next_name)
            self.next_name += 1
        if k == "range":
            a = r.randint(0, 6)
            b = a + r.randint(0, 6)
            return name + "%d...%d" % (a, b)
        if k == "exact":
            return name + self.word(r.randint(1, 4))
        if k == "sim":
            return name + self.word(r.randint(1, 6)) + "[%d,%d,%d]" % (r.randint(0, 2), r.randint(0, 2), r.randint(0, 2))
        if k in ("repeat", "inv", "compl", "rcompl", "llim"):
            if not self.defined:
                a = r.randint(1, 5)
                name = "p%d=" % self.next_name
                self.defined.append(self.next_name)
                self.next_name += 1
                return name + "%d...%d" % (a, a + r.randint(0, 3))
            n = r.choice(self.defined)
            if k == "repeat":
                return "p%d%s" % (n, self.mid())
            if k == "inv":
                return "<p%d%s" % (n, self.mid())
            if k == "compl":
                return "~p%d%s" % (n, self.mid())
            if k == "rcompl":
                rn = r.randint(1, 3)
                if rn not in self.rules:
                    self.rules.append(rn)
                return "r%d~p%d%s" % (rn, n, self.mid())
            if k == "llim":
                names = r.sample(self.defined, min(len(self.defined), r.randint(1, 3)))
                return "length(%s) < %d" % ("+".join("p%d" % x for x in names), r.randint(1, 15))
        if k == "weight":
            tup = 4 if not self.protein or r.random() < 0.3 else 20
            L = r.randint(1, 4)
            vecs = ["(" + ",".join(str(r.randint(0, 100)) for _ in range(tup)) + ")" for _ in range(L)]
            body = "{" + ",".join(vecs) + "}"
            cut = r.randint(0, 60 * L)
            if r.random() < 0.3:
                return name + "%d > %s > %d" % (cut + r.randint(10, 200), body, cut)
            return name + "%s > %d" % (body, cut)
        if k == "or":
            if depth > 1:
                return "%d...%d" % (1, 2)
            return "(" + self.lst(depth + 1, 1, 2) + " | " + self.lst(depth + 1, 1, 2) + ")"
        if k == "start":
            return "^"
        if k == "end":
            return "$"
        if k == "any":
            return name + "any(" + "".join(r.sample(AA, r.randint(1, 4))) + ")"
        if k == "notany":
            return name + "notany(" + "".join(r.sample(AA, r.randint(1, 4))) + ")"
        return "1...2"

    def lst(self, depth, lo, hi):
        return " ".join(self.unit(depth) for _ in range(self.r.randint(lo, hi)))

    def build(self):
        body = self.lst(0, 1, 5)
        rules = []
        for rn in self.rules:
            pairs = self.r.sample(["au", "ua", "gc", "cg", "gu", "ug", "ga", "ag", "AT", "TA", "GC", "CG"], self.r.randint(1, 6))
            rules.append("r%d={%s}" % (rn, ",".join(pairs)))
        txt = " ".join(rules + [body])
        r = self.r
        if r.random() < 0.1:
            # comments and newlines
            txt = "% comment line\n" + txt.replace(" ", "\n", 1) + "\n% trailing"
        if r.random() < 0.05:
            # random corruption
            i = r.randint(0, len(txt))
            txt = txt[:i] + r.choice(["(", ")", "[", "]", "|", ",", "~", "x", "p9", "..", "{"]) + txt[i:]
        return txt


def run(binary, args, pat_path, inp):
    try:
        p = subprocess.run([binary] + args + [pat_path], input=inp, capture_output=True, timeout=TIMEOUT)
        return (p.returncode, p.stdout, p.stderr)
    except subprocess.TimeoutExpired:
        return ("TIMEOUT", b"", b"")


def one(case):
    r = random.Random(SEED * 1000003 + case)
    protein = r.random() < 0.25
    pat = Pat(r, protein).build()
    fasta = rnd_fasta(r, protein)
    args = []
    if protein:
        args.append("-p")
    if r.random() < 0.4:
        args.append("-c")
    if r.random() < 0.3:
        args += ["-o", "1"]
    with tempfile.TemporaryDirectory() as d:
        pp = os.path.join(d, "pat")
        with open(pp, "w") as f:
            f.write(pat + "\n")
        a = run(C_BIN, args, pp, fasta.encode())
        b = run(R_BIN, args, pp, fasta.encode())
    if a[0] == "TIMEOUT" and b[0] == "TIMEOUT":
        return (case, "both-timeout", pat, args, fasta, a, b)
    if a != b:
        return (case, "DIFF", pat, args, fasta, a, b)
    return (case, "ok", pat, args, fasta, a, b)


def main():
    stats = {}
    bad = []
    hits = 0
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as ex:
        for res in ex.map(one, range(N)):
            stats[res[1]] = stats.get(res[1], 0) + 1
            if res[1] == "ok" and res[5][1]:
                hits += 1
            if res[1] == "DIFF":
                bad.append(res)
    print("results:", stats, "cases with hits:", hits)
    for case, _, pat, args, fasta, a, b in bad[:5]:
        print("=" * 60)
        print("case", case, "args", args)
        print("pattern:", repr(pat))
        print("fasta:", repr(fasta[:500]))
        print("C   :", a[0], a[1][:600], a[2][:300])
        print("Rust:", b[0], b[1][:600], b[2][:300])
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
