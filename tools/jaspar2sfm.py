#!/usr/bin/env python3
"""Convert a JASPAR position frequency matrix to a scan_for_matches weight unit.

usage: jaspar2sfm.py MATRIX [--pvalue 1e-4] [--scale 100] [--pseudocount 0.25]
                     [--background A,C,G,T] [--revcomp] [--name NAME]

Weights: w = round_half_up(scale * log2(((n + pc*bg) / (N + pc)) / bg)).
Cutoff T: smallest integer with P(score > T) <= pvalue under the background
(exact DP over the integer score distribution).  The engine tests score > T.
Several matrices in one file: each is converted (one block each); --name
keeps only the matrix whose id or name equals NAME.
Python standard library only.
"""
import argparse
import math
import re
import sys

BASES = "ACGT"
UNIT_WARN_LEN = 30000


class JasparError(Exception):
    pass


def parse_jaspar(text):
    """Return a list of (id, name, columns); columns[i] = [a, c, g, t] counts."""
    mats = []
    cur = None

    def finish():
        if cur is None:
            return
        rows = cur["rows"]
        if not rows:
            return
        missing = [b for b in BASES if b not in rows]
        if missing:
            raise JasparError("matrix %s: missing row(s) %s" % (cur["id"], ",".join(missing)))
        lens = {len(v) for v in rows.values()}
        if len(lens) != 1 or lens == {0}:
            raise JasparError("matrix %s: rows have different or zero length" % cur["id"])
        n = lens.pop()
        cols = [[rows[b][i] for b in BASES] for i in range(n)]
        mats.append((cur["id"], cur["name"], cols))

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            parts = line[1:].split(None, 1)
            if not parts:
                raise JasparError("empty header line")
            cur = {"id": parts[0], "name": parts[1].strip() if len(parts) > 1 else "",
                   "rows": {}}
            continue
        m = re.match(r"^([ACGTacgt])\s*\[?([^\]]*)\]?\s*$", line)
        if not m:
            raise JasparError("cannot parse line: %r" % line)
        if cur is None:
            cur = {"id": "matrix", "name": "", "rows": {}}
        b = m.group(1).upper()
        if b in cur["rows"]:
            raise JasparError("matrix %s: duplicate row %s" % (cur["id"], b))
        try:
            vals = [float(x) for x in m.group(2).split()]
        except ValueError:
            raise JasparError("non-numeric count in line: %r" % line)
        if any(v < 0 or math.isnan(v) or math.isinf(v) for v in vals):
            raise JasparError("bad (negative or non-finite) count in line: %r" % line)
        cur["rows"][b] = vals
    finish()
    if not mats:
        raise JasparError("no matrix found in input")
    return mats


def revcomp_columns(cols):
    """Reverse the columns; swap A<->T and C<->G."""
    return [[c[3], c[2], c[1], c[0]] for c in reversed(cols)]


def parse_background(s):
    try:
        v = [float(x) for x in s.split(",")]
    except ValueError:
        raise JasparError("bad --background %r" % s)
    if len(v) != 4 or any(x <= 0 for x in v):
        raise JasparError("--background needs 4 positive numbers A,C,G,T")
    t = sum(v)
    return [x / t for x in v]


def int_weights(cols, bg, scale, pseudocount):
    out = []
    for c in cols:
        n = sum(c)
        row = []
        for k in range(4):
            f = (c[k] + pseudocount * bg[k]) / (n + pseudocount)
            if f <= 0:
                raise JasparError("zero frequency: use --pseudocount > 0")
            row.append(int(math.floor(scale * math.log2(f / bg[k]) + 0.5)))
        out.append(row)
    return out


def score_distribution(weights, bg):
    """Exact {integer score: probability} over all words, by DP."""
    dist = {0: 1.0}
    for col in weights:
        nd = {}
        for s, p in dist.items():
            for k in range(4):
                key = s + col[k]
                nd[key] = nd.get(key, 0.0) + p * bg[k]
        dist = nd
    return dist


def find_cutoff(dist, pvalue):
    """Smallest integer T with P(score > T) <= pvalue; returns (T, P(score > T))."""
    keys = sorted(dist)
    tail = {}
    acc = 0.0
    for s in reversed(keys):
        tail[s] = acc  # P(score > s)
        acc += dist[s]
    eps = 1e-12 * pvalue
    for s in keys:
        if tail[s] <= pvalue + eps:
            # an integer between two score values has the tail of the lower
            # one (> pvalue), so the smallest qualifying integer is s itself
            return s, tail[s]
    return keys[-1], 0.0


def convert(mid, name, cols, args, bg):
    if args.revcomp:
        cols = revcomp_columns(cols)
    w = int_weights(cols, bg, args.scale, args.pseudocount)
    dist = score_distribution(w, bg)
    t, p = find_cutoff(dist, args.pvalue)
    maxs = sum(max(c) for c in w)
    unit = "{" + ",".join("(%d,%d,%d,%d)" % tuple(c) for c in w) + "} > %d" % t
    lines = [
        "%% matrix: %s %s%s" % (mid, name, " (reverse complement)" if args.revcomp else ""),
        "%% length: %d" % len(w),
        "%% scale: %g  pseudocount: %g" % (args.scale, args.pseudocount),
        "%% background A,C,G,T: %s" % ",".join("%.4g" % x for x in bg),
        "%% p-value asked: %g" % args.pvalue,
        "%% cutoff: %d (score > %d)" % (t, t),
        "%% exact p-value of cutoff: %.4e" % p,
        "%% max score: %d" % maxs,
        unit,
    ]
    if len(unit) > UNIT_WARN_LEN:
        sys.stderr.write("warning: unit is %d characters; the pattern file is cut at "
                         "31 999 bytes\n" % len(unit))
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("matrix", help="JASPAR file, or - for stdin")
    ap.add_argument("--pvalue", type=float, default=1e-4)
    ap.add_argument("--scale", type=float, default=100.0)
    ap.add_argument("--pseudocount", type=float, default=0.25)
    ap.add_argument("--background", default="0.25,0.25,0.25,0.25")
    ap.add_argument("--revcomp", action="store_true")
    ap.add_argument("--name", default=None,
                    help="keep only the matrix with this id or name")
    args = ap.parse_args(argv)
    try:
        if not (0 < args.pvalue < 1):
            raise JasparError("--pvalue must be between 0 and 1")
        if args.scale <= 0:
            raise JasparError("--scale must be > 0")
        if args.pseudocount < 0:
            raise JasparError("--pseudocount must be >= 0")
        bg = parse_background(args.background)
        try:
            text = sys.stdin.read() if args.matrix == "-" else open(args.matrix).read()
        except OSError as e:
            raise JasparError("cannot read %s: %s" % (args.matrix, e))
        mats = parse_jaspar(text)
        if args.name is not None:
            mats = [m for m in mats if args.name in (m[0], m[1])]
            if not mats:
                raise JasparError("no matrix with id or name %r" % args.name)
        blocks = [convert(i, n, c, args, bg) for i, n, c in mats]
    except JasparError as e:
        sys.stderr.write("jaspar2sfm: error: %s\n" % e)
        return 1
    sys.stdout.write("\n".join(blocks) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
