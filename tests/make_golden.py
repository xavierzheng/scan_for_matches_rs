#!/usr/bin/env python3
"""Record golden cases from the C reference program.

usage: make_golden.py C_BINARY > tests/data/golden.tsv

Each line: pattern<TAB>args<TAB>stdin<TAB>status<TAB>stdout, all byte
strings hex encoded (args joined by a NUL byte).  status is the exit code,
or "sigN" when the program was killed by signal N.  Only cases where the C
program gives the same result on repeated runs are kept (a few kinds of
undefined behaviour in C depend on address-space randomisation).
"""
import os
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fuzz_compare as fz  # noqa: E402

C_BIN = sys.argv[1]

IN3 = b""">s1 desc
acgtacguaaccggttaaccgguuacgtacgu
>s2
ACGTACGUAACCGGTTAACCGGUUACGTACGUnnnnaaaaaaaaaaaaRYMKacgt
>s3
GGGGAAAACCCCTTTTGGGGAAAACCCCTTTTACGTACGTAGCTAGCTAGCTTTTTT
"""
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
DNA_IN = open(os.path.join(DATA, "test_dna_input"), "rb").read()
PROT_IN = open(os.path.join(DATA, "test_prot_input"), "rb").read()

EDGE = [
    "~p3", "p1=3...3 ~p2", "length(p4) < 3 p1=1...1", "p1=~p1", "p1=3...3 r7~p1", "^", "0...0",
    "length(p1) < 3", "5...3 AC", "p1=3...3 p1=4...4", "(AC | GG", "(AC | GG) (TT | (CC | AA))",
    "p1=2...4 (p1 | ~p1) 0...2 $", "^ ACG", "AA $", "p1=2...2 <p1[1,0,0]", "AC[0,5,0]",
    "AC[-1,0,0] GG", "p1=2...3 1...2 ~p1[0,2,2]", "p1=3...3 p1[1,1,1]", "p1=3...3 length(p1) < -1",
    "r50={au,ua} p1=3...3 r50~p1", "r0={gc,cg,au,ua} p0=3...3 r0~p0", "p49=3...3 p49", "p51=3...3 p51",
    "{(1,2,3,4),(1,2,3)} > 1", "{(1,2,3)} > 1", "1000 > {(10,20,30,40),(10,20,30,40)} > 60",
    "p1=AAAA[1,0,0] 0...5 p1[0,1,1]", "p1=1...3 r1~p1 r1={aa}", "  ", "%only comment",
    "p1=3...3 p1 % comment\n ~p1", "p1=2...2 p2=1...1 length(p1+p2+p7) < 5",
    "p1=1...1 (~p1 | p1) p2=1...1 (p2 | AC)", "p50=3...3 p50", "p50=3...3 ~p50",
    "p50=2...4 0...3 p50[1,0,0]", "p1=2...2 p50=1...1 length(p1+p50) < 4", "-3...5 AA",
    "p1=4...7 3...8 ~p1", "p1=6...6 3...8 p1", "TATAA[1,0,0]", "RRRRYYYY",
    "p1=4...8 0...3 p2=6...8 p1 0...3 p2",
    "ARRYYTT p1=0...5 GCA[1,0,0] p2=1...6 ~p1 4...8 ~p2 p3=4...10 CCT length(p1+p2+p3) < 9",
    "600 > {(16,0,84,0),(57,10,29,4),(0,80,0,20),(95,0,0,5),(0,100,0,0),(18,60,20,2),(0,0,100,0),(0,50,50,0)} > 450",
    "(GAGA | (GCGCA | TTCGA))", "p1=6...6 <p1",
]
PROT_EDGE = ["p1=2...2 0...3 p1", "any(AC) notany(A) X", "p1=3...3 ~p1", "ACGT", "A",
             "{(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)} > 10", "{(1,2,3,4)} > 10",
             "p1=0...4 any(HQD) 1...3 notany(HK) p1", "CGXXYWG[1,0,0]", "<p1[1,0,0]"]
LIMITS = [
    " ".join(["1...1"] * 100), " ".join(["1...1"] * 101), " ".join(["1...1"] * 400), " ".join(["N"] * 100),
    " ".join(["N"] * 101), " ".join(["1...1"] * 100 + ["A"]), "N" * 650, "N" * 650 + " {(10,10,10,10)} > 5",
    "N" * 590 + " {(10,10,10,10)} > 5", "N" * 400 + " " + "N" * 400, "(A | " * 30 + "C" + ")" * 30,
    "p1=101...101 3...8 ~p1[1,0,0]", "p1=100...100 3...8 ~p1[1,0,0]", "p1=101...101 3...8 ~p1",
    "a" * 150 + "[150,150,1]", "a" * 101 + "[101,101,1]",
]
OPTS = [[], ["-c"], ["-o", "1"], ["-o1", "-c"], ["-co", "1"], ["-n"], ["-m"], ["-o", "1", "-m", "2"],
        ["-x"], ["--", "-c"], ["-p"], ["-pc"], ["-c", "--"]]


def run(args, pat, inp):
    with tempfile.TemporaryDirectory() as d:
        pp = os.path.join(d, "pat")
        with open(pp, "wb") as f:
            f.write(pat)
        res = set()
        for _ in range(3):
            try:
                p = subprocess.run([C_BIN] + args + [pp], input=inp, capture_output=True, timeout=10)
            except subprocess.TimeoutExpired:
                return None
            st = str(p.returncode) if p.returncode >= 0 else "sig%d" % -p.returncode
            res.add((st, p.stdout))
        if len(res) != 1:
            return None
        return res.pop()


def emit(args, pat, inp):
    r = run(args, pat, inp)
    if r is None:
        return 0
    st, out = r
    if len(out) > 2_000_000:
        return 0
    print("\t".join([pat.hex(), b"\0".join(a.encode() for a in args).hex(), inp.hex(), st, out.hex()]))
    return 1


def main():
    n = 0
    for pat in EDGE:
        for args in ([], ["-c"], ["-o", "1"]):
            n += emit(args, pat.encode() + b"\n", IN3)
    for pat in PROT_EDGE:
        n += emit(["-p"], pat.encode() + b"\n", PROT_IN)
    for pat in LIMITS:
        seq = b">x\n" + bytes(random.Random(1).choice(b"acgt") for _ in range(400)) + b"\n>y\n" + b"A" * 300 + b"\n"
        n += emit([], pat.encode() + b"\n", seq)
    for args in OPTS:
        n += emit(args, b"p1=4...7 3...8 ~p1\n", DNA_IN)
    weird = b">a\nacgt\x00acgt\n>b\r\nac gt\racgt\xa8\xffacgt\n>c\n" + bytes(range(0x80, 0xa0)) + b"acgt\n"
    for pat in ["n", "c", "p1=1...1 ~p1", "{(10,20,30,40)} > 5"]:
        for args in ([], ["-c"]):
            n += emit(args, pat.encode() + b"\n", weird)
    # random cases from the differential fuzzer (normal and stress mode)
    r0 = random.Random(20261001)
    for case in range(600):
        r = random.Random(r0.random())
        protein = r.random() < 0.25
        pat = fz.Pat(r, protein).build() if r.random() > 0.2 else fz.stress_pattern(r, protein)
        fasta = fz.rnd_fasta(r, protein)
        args = (["-p"] if protein else []) + (["-c"] if r.random() < 0.4 else []) + (["-o", "1"] if r.random() < 0.3 else [])
        n += emit(args, pat.encode("latin1") + b"\n", fasta.encode("latin1"))
    print("golden cases:", n, file=sys.stderr)


if __name__ == "__main__":
    main()
