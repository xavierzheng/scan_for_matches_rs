import itertools
import math
import os
import random
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import jaspar2sfm as J  # noqa: E402

BIN = "/nfs/project_ssd/project3/pxzhe/TE_software/001_scan_for_match_rs/scan_for_matches_original_rs/target/release/scan_for_matches"
BG = [0.25] * 4

JASPAR = """>MA0000.1 TEST
A [ 10  0  0  5 ]
C [  0 10  0  5 ]
G [  0  0 10  0 ]
T [  0  0  0  0 ]
>MA0001.1 SECOND
A  [ 1 2 ]
C  [ 3 4 ]
G  [ 5 6 ]
T  [ 7 8 ]
"""


class Args:
    pvalue = 1e-4
    scale = 100.0
    pseudocount = 0.25
    revcomp = False


def brute(weights, bg):
    d = {}
    for w in itertools.product(range(4), repeat=len(weights)):
        s = sum(weights[i][k] for i, k in enumerate(w))
        p = math.prod(bg[k] for k in w)
        d[s] = d.get(s, 0.0) + p
    return d


class ParseTest(unittest.TestCase):
    def test_parse(self):
        m = J.parse_jaspar(JASPAR)
        self.assertEqual(len(m), 2)
        self.assertEqual(m[0][:2], ("MA0000.1", "TEST"))
        self.assertEqual(m[0][2], [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 10, 0], [5, 5, 0, 0]])
        self.assertEqual(m[1][2][1], [2, 4, 6, 8])

    def test_no_brackets(self):
        m = J.parse_jaspar(">x\nA 1 2\nC 3 4\nG 5 6\nT 7 8\n")
        self.assertEqual(m[0][2], [[1, 3, 5, 7], [2, 4, 6, 8]])

    def test_bad(self):
        for t in ["", ">x\nA [1 2]\nC [1 2]\nG [1 2]\n", ">x\nA [1 2]\nC [1]\nG [1 2]\nT [1 2]\n",
                  ">x\nA [1 z]\nC [1 2]\nG [1 2]\nT [1 2]\n", "hello\n"]:
            with self.assertRaises(J.JasparError):
                J.parse_jaspar(t)


class WeightTest(unittest.TestCase):
    def test_hand(self):
        # one column, counts 3,1,0,0; N=4; pc=0 -> not allowed for zero; use pc=1
        # f = (c + 0.25)/(4+1) ; w = round(100*log2(f/0.25))
        w = J.int_weights([[3, 1, 0, 0]], BG, 100, 1.0)
        exp = [int(math.floor(100 * math.log2((c + 0.25) / 5 / 0.25) + 0.5)) for c in (3, 1, 0, 0)]
        self.assertEqual(w[0], exp)
        self.assertEqual(w[0], [138, 0, -232, -232])

    def test_uniform_is_zero(self):
        self.assertEqual(J.int_weights([[5, 5, 5, 5]], BG, 100, 0.25), [[0, 0, 0, 0]])


class DPTest(unittest.TestCase):
    def test_dp_vs_brute(self):
        rnd = random.Random(1)
        for k in (4, 5, 6):
            cols = [[rnd.randint(0, 20) for _ in range(4)] for _ in range(k)]
            bg = [0.3, 0.2, 0.2, 0.3]
            w = J.int_weights(cols, bg, 100, 0.25)
            dp = J.score_distribution(w, bg)
            br = brute(w, bg)
            self.assertEqual(set(dp), set(br))
            for s in br:
                self.assertAlmostEqual(dp[s], br[s], places=12)

    def test_cutoff_rule(self):
        rnd = random.Random(2)
        cols = [[rnd.randint(0, 20) for _ in range(4)] for _ in range(5)]
        w = J.int_weights(cols, BG, 100, 0.25)
        d = brute(w, BG)
        for pv in (1e-1, 1e-2, 1e-3, 5e-3):
            t, p = J.find_cutoff(J.score_distribution(w, BG), pv)
            tail = lambda x: sum(v for s, v in d.items() if s > x)
            self.assertLessEqual(tail(t), pv * (1 + 1e-9))
            self.assertAlmostEqual(tail(t), p, places=12)
            self.assertGreater(tail(t - 1), pv)  # smallest integer

    def test_strict(self):
        # two words equally likely at top: P(score > 5)=0.5 only if scores 5/10, check strict >
        dist = {0: 0.5, 5: 0.25, 10: 0.25}
        self.assertEqual(J.find_cutoff(dist, 0.25), (5, 0.25))
        self.assertEqual(J.find_cutoff(dist, 0.24), (10, 0.0))


class RevcompTest(unittest.TestCase):
    def test_rc(self):
        cols = [[1, 2, 3, 4], [5, 6, 7, 8]]
        self.assertEqual(J.revcomp_columns(cols), [[8, 7, 6, 5], [4, 3, 2, 1]])
        self.assertEqual(J.revcomp_columns(J.revcomp_columns(cols)), cols)


class CliTest(unittest.TestCase):
    def run_cli(self, text, *extra):
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, "m.jaspar")
            with open(f, "w") as fh:
                fh.write(text)
            return subprocess.run([sys.executable, "-I", os.path.join(HERE, "jaspar2sfm.py"), f] + list(extra),
                                  capture_output=True, text=True)

    def test_output(self):
        r = self.run_cli(JASPAR, "--name", "MA0000.1", "--pvalue", "0.01")
        self.assertEqual(r.returncode, 0, r.stderr)
        lines = r.stdout.splitlines()
        self.assertTrue(all(l.startswith("%") for l in lines[:-1]))
        self.assertRegex(lines[-1], r"^\{\(.*\)\} > -?\d+$")

    def test_multi(self):
        r = self.run_cli(JASPAR)
        self.assertEqual(r.stdout.count("} > "), 2)

    def test_bad_exit(self):
        self.assertEqual(self.run_cli("junk\n").returncode, 1)
        self.assertEqual(self.run_cli(JASPAR, "--name", "nope").returncode, 1)


@unittest.skipUnless(os.access(BIN, os.X_OK), "scan_for_matches binary not found")
class EndToEnd(unittest.TestCase):
    def test_planted(self):
        counts = {"A": [0, 30, 0, 0, 25, 2, 0, 1, 20, 0],
                  "C": [0, 0, 30, 0, 3, 28, 0, 0, 5, 0],
                  "G": [30, 0, 0, 30, 2, 0, 30, 29, 3, 30],
                  "T": [0, 0, 0, 0, 0, 0, 0, 0, 2, 0]}
        text = ">T.1 T\n" + "".join("%s [ %s ]\n" % (b, " ".join(map(str, counts[b]))) for b in "ACGT")
        cols = J.parse_jaspar(text)[0][2]
        bg = BG
        w = J.int_weights(cols, bg, 100, 0.25)
        t, p = J.find_cutoff(J.score_distribution(w, bg), 1e-4)
        unit = "{" + ",".join("(%d,%d,%d,%d)" % tuple(c) for c in w) + "} > %d" % t
        cons = "".join("ACGT"[max(range(4), key=lambda k: c[k])] for c in cols)
        rnd = random.Random(7)
        n = 200000
        seq = [rnd.choice("ACGT") for _ in range(n)]
        pos = sorted(rnd.sample(range(0, n - 20, 10000), 20))
        for q in pos:
            seq[q:q + len(cons)] = cons
        with tempfile.TemporaryDirectory() as d:
            pat = os.path.join(d, "p.pat")
            with open(pat, "w") as fh:
                fh.write("% test\n" + unit + "\n")
            fa = os.path.join(d, "s.fa")
            with open(fa, "w") as fh:
                fh.write(">s\n" + "".join(seq) + "\n")
            with open(fa) as fin:
                r = subprocess.run([BIN, pat], stdin=fin, capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        starts = set()
        nhits = 0
        for line in r.stdout.splitlines():
            if line.startswith(">"):
                nhits += 1
                starts.add(int(line.split(":[")[1].split(",")[0]))
        # the engine reports non-overlapping hits: a random hit overlapping a
        # planted site may hide it, so "found" = some hit overlaps the site
        L = len(cons)
        missing = [q for q in pos if not any(q + 1 - L < s <= q + L for s in starts)]
        self.assertEqual(missing, [])
        exact = sum(1 for q in pos if q + 1 in starts)
        self.assertGreaterEqual(exact, 18)
        random_hits = nhits - 20
        expected = n * p
        self.assertLess(random_hits, 10 * expected + 10)


if __name__ == "__main__":
    unittest.main()
