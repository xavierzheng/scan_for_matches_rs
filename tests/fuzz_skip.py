#!/usr/bin/env python3
"""Random check of gap-length skipping: compare two builds (one without
skipping, for example 0.1.0, and one with it) on patterns with wide ranges
followed by exact words, exact reverse complements and exact repeats.

usage: fuzz_skip.py [--chain | --tir] OLD_BINARY NEW_BINARY [N_CASES] [SEED]

--chain: every pattern has a range followed by several units that check
fixed strings (words, named words, ~pN and pN, also of names in the chain
or of the range), which `next_start` checks together.

--tir: TIR patterns whose reverse complement allows mismatches
(`p2=10...24 60...8000 ~p2[1,0,0]`, 1 to 3 mismatches, few exact letters
after it): the gap skipping with tolerant masks.
"""
import os
import random
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

CHAIN = "--chain" in sys.argv
TIR = "--tir" in sys.argv
ARGS = [a for a in sys.argv[1:] if a not in ("--chain", "--tir")]
OLD, NEW = os.path.abspath(ARGS[0]), os.path.abspath(ARGS[1])
N = int(ARGS[2]) if len(ARGS) > 2 else 1000
SEED = int(ARGS[3]) if len(ARGS) > 3 else 1
COMP = {"a": "t", "c": "g", "g": "c", "t": "a"}


def rc(s):
    return "".join(COMP.get(c, c) for c in reversed(s))


def pattern(r):
    names = []
    units = []
    nxt = 1

    def name():
        nonlocal nxt
        n = nxt
        nxt += 1
        names.append(n)
        return "p%d=" % n

    for _ in range(r.randint(2, 6)):
        k = r.random()
        if k < 0.35 or not names:
            a = r.randint(0, 8)
            b = a + r.choice([0, 1, 3, 10, 60, 300, 2000])
            if r.random() < 0.05:
                a, b = b, a  # reversed range
            units.append((name() if r.random() < 0.6 else "") + "%d...%d" % (a, b))
        elif k < 0.55:
            n = r.choice(names)
            units.append("~p%d" % n + ("[1,0,0]" if r.random() < 0.1 else ""))
        elif k < 0.7:
            n = r.choice(names)
            units.append("p%d" % n + ("[0,1,0]" if r.random() < 0.1 else ""))
        elif k < 0.85:
            w = "".join(r.choice("acgt" if r.random() < 0.8 else "acgtnry") for _ in range(r.randint(1, 8)))
            units.append((name() if r.random() < 0.2 else "") + w)
        elif k < 0.9:
            units.append("<p%d" % r.choice(names))
        elif k < 0.95:
            units.append("(%s | %d...%d)" % (r.choice(["acg", "tt", "p1=2...3"]), r.randint(0, 3), r.randint(3, 50)))
        else:
            units.append("length(p%d) < %d" % (r.choice(names), r.randint(2, 30)))
    if r.random() < 0.05:
        units.insert(0, "^")
    if r.random() < 0.05:
        units.append("$")
    return " ".join(units)


def chain_pattern(r):
    units = []
    names = []
    nxt = 1
    for _ in range(r.randint(1, 2)):
        a = r.randint(0, 6)
        units.append("p%d=%d...%d" % (nxt, a, a + r.choice([0, 0, 1, 2])))
        names.append(nxt)
        nxt += 1
    if r.random() < 0.15:
        units.insert(0, r.choice(["acg", "t", "1...3"]))
    a = r.randint(0, 5)
    b = a + r.choice([10, 60, 300, 2000])
    rname = None
    if r.random() < 0.3:
        rname = nxt
        nxt += 1
        units.append("p%d=%d...%d" % (rname, a, b))
    else:
        units.append("%d...%d" % (a, b))
    chain_names = []
    for _ in range(r.randint(2, 4)):
        k = r.random()
        pool = names + chain_names + ([rname] if rname and r.random() < 0.3 else [])
        if k < 0.35:
            errs = r.choice(["", "", "[1,0,0]", "[2,0,0]", "[3,0,0]", "[0,1,0]", "[1,0,1]"])
            units.append("~p%d%s" % (r.choice(pool), errs))
        elif k < 0.6:
            units.append("p%d" % r.choice(pool))
        elif k < 0.85:
            w = "".join(r.choice("acgt" if r.random() < 0.9 else "acgtnry") for _ in range(r.randint(1, 6)))
            if r.random() < 0.3:
                units.append("p%d=%s" % (nxt, w))
                chain_names.append(nxt)
                nxt += 1
            else:
                units.append(w)
        elif k < 0.92:
            units.append("~p%d[1,0,0]" % r.choice(pool))
        else:
            units.append("(acg | tt)")
    if r.random() < 0.3:
        units.append("0...%d" % r.randint(0, 20))
    if r.random() < 0.2:
        units.append(r.choice(["acgt", "~p1", "p1"]))
    if r.random() < 0.15:
        # the range and the chain inside an alternative
        k = len(names) + (1 if units[0][0] not in "p" else 0)
        units = units[:k] + ["(" + " ".join(units[k:]) + " | tt)", r.choice(["acg", "p1", "~p1"])]
    return " ".join(units)


def tir_pattern(r):
    units = []
    tsd = r.random() < 0.6
    if tsd:
        a = r.randint(2, 9)
        units.append("p1=%d...%d" % (a, a + r.choice([0, 0, 1])))
    if r.random() < 0.3:
        units.append(r.choice(["ta", "cac", "g"]))
    k = r.randint(10, 24)
    units.append("p2=%d...%d" % (k, k + r.choice([0, 0, 2])))
    a = r.randint(0, 60)
    units.append("%d...%d" % (a, a + r.choice([60, 300, 2000, 8000])))
    units.append("~p2[%d,0,0]" % r.choice([1, 1, 2, 3]))
    if r.random() < 0.3:
        units.append(r.choice(["ta", "gtg", "c", "acgt"]))
    if tsd:
        units.append(r.choice(["p1", "p1", "p1[1,0,0]"]))
    return " ".join(units)


def sequence(r, n):
    s = [r.choice("acgt") for _ in range(n)]
    # planted hairpins and repeats
    for _ in range(r.randint(0, 6)):
        k = r.randint(3, 12)
        i = r.randint(0, max(0, n - 1))
        stem = "".join(r.choice("acgt") for _ in range(k))
        gap = r.randint(0, 400)
        ins = stem + "".join(r.choice("acgt") for _ in range(gap)) + (rc(stem) if r.random() < 0.6 else stem)
        s[i:i] = list(ins)
    # planted TIR-like elements: TSD, TIR, gap, reverse complement, TSD
    for _ in range(r.randint(0, 4) if CHAIN or TIR else 0):
        tsd = "".join(r.choice("acgt") for _ in range(r.randint(0, 9 if TIR else 6)))
        tir = "".join(r.choice("acgt") for _ in range(r.randint(10, 26) if TIR else r.randint(1, 10)))
        mid = "".join(r.choice("acgt") for _ in range(r.randint(0, 3000 if TIR else 400)))
        i = r.randint(0, max(0, len(s) - 1))
        end = list(rc(tir))
        for _ in range(r.choice([0, 0, 1, 2])):
            end[r.randrange(len(end))] = r.choice("acgtn")
        s[i:i] = list(tsd + tir + mid + "".join(end) + tsd)
    if r.random() < 0.2:
        for _ in range(r.randint(1, 10)):
            s[r.randrange(len(s))] = r.choice("nNryRY")
    return "".join(s)


# keep memory low: outputs are compared by hash, and a case whose output
# grows past OUT_LIMIT bytes is stopped and counted as "big"
OUT_LIMIT = 20_000_000
WORKERS = 4


def run(b, args, pat, inp):
    import hashlib
    import threading
    with tempfile.TemporaryDirectory() as d:
        pp = os.path.join(d, "pat")
        with open(pp, "w") as f:
            f.write(pat + "\n")
        p = subprocess.Popen([b] + args + [pp], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        def feed():
            try:
                p.stdin.write(inp)
                p.stdin.close()
            except BrokenPipeError:
                pass  # the program was stopped (time or output limit)

        threading.Thread(target=feed, daemon=True).start()
        timer = threading.Timer(60, p.kill)
        timer.start()
        h, n = hashlib.sha1(), 0
        while True:
            chunk = p.stdout.read(1 << 16)
            if not chunk:
                break
            h.update(chunk)
            n += len(chunk)
            if n > OUT_LIMIT:
                p.kill()
                p.wait()
                timer.cancel()
                return "big"
        err = p.stderr.read(1 << 16)
        rc = p.wait()
        timer.cancel()
        if rc == -9:
            return None
    return (rc, h.hexdigest(), n, err)


def one(case):
    r = random.Random(SEED * 1000003 + case)
    pat = chain_pattern(r) if CHAIN else tir_pattern(r) if TIR else pattern(r)
    sizes = [300, 2000, 6000, 20000] if TIR else [50, 300, 2000, 6000]
    recs = "".join(">r%d\n%s\n" % (k, sequence(r, r.choice(sizes))) for k in range(r.randint(1, 3)))
    args = (["-c"] if r.random() < 0.4 else []) + (["-o", "1"] if r.random() < 0.3 else [])
    a = run(OLD, args, pat, recs.encode())
    b = run(NEW, args, pat, recs.encode())
    if b is None and a is not None and a != "big":
        return ("NEW-timeout", (pat, args, len(recs)))
    if a is None or b is None:
        return ("timeout", None)
    if a == "big" or b == "big":
        return ("big", None)
    if a != b:
        return ("DIFF", (pat, args, len(recs)))
    return ("ok-hits" if a[2] else "ok", None)


with ThreadPoolExecutor(WORKERS) as ex:
    res = list(ex.map(one, range(N)))
stats = {}
for k, _ in res:
    stats[k] = stats.get(k, 0) + 1
print("results:", stats)
for k, info in res:
    if k in ("DIFF", "NEW-timeout"):
        print(info)
sys.exit(1 if stats.get("DIFF") or stats.get("NEW-timeout") else 0)
