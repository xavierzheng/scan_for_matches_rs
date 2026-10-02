#!/usr/bin/env python3
"""Random check of gap-length skipping: compare two builds (one without
skipping, for example 0.1.0, and one with it) on patterns with wide ranges
followed by exact words, exact reverse complements and exact repeats.

usage: fuzz_skip.py OLD_BINARY NEW_BINARY [N_CASES] [SEED]
"""
import os
import random
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

OLD, NEW = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
N = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 1
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
        threading.Thread(target=lambda: (p.stdin.write(inp), p.stdin.close()), daemon=True).start()
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
    pat = pattern(r)
    recs = "".join(">r%d\n%s\n" % (k, sequence(r, r.choice([50, 300, 2000, 6000]))) for k in range(r.randint(1, 3)))
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
