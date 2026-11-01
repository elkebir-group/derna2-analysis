"""Shared helpers for the E2 figure pipeline.

Every script in this pipeline used to obtain these by exec-ing a prefix of pick_example3.py,
which meant an edit anywhere above that file's first statement silently changed what the other
scripts computed, and nothing recorded which version had run. They live here instead, imported
normally, so a change is visible in one place and the module's own checksum goes into each
figure's provenance.

Stdlib only: no ViennaRNA, no matplotlib. The interpreter that computes a figure's data and the
one that draws it are different, and this module has to import under both.
"""
import glob
import json
import os
import re

from e3paths import OUT, E  # noqa: F401

CUTOFFS = [0.80, 0.85, 0.90]
SLACK = 10.0
NESTS = ["F0", "F1", "F2", "F3"]
IUPAC = {"R": "[AG]", "Y": "[CU]", "N": "[ACGU]", "S": "[CG]", "W": "[AU]", "K": "[GU]",
         "M": "[AC]"}
SETS = {"I": "set1-manufacturing", "II": "set2-processing", "III": "set3-repeats"}
PATS = {k: [re.compile("".join(IUPAC.get(c, c) for c in
                                   line.split("#")[0].strip().upper().replace("T", "U")))
            for line in open(f"{E}/motifs/{f}.txt") if line.split("#")[0].strip()]
        for k, f in SETS.items()}

NAMED = {}
for line in open(f"{E}/FINISHED.tsv"):
    if line.startswith("#") or line.startswith("cell\t"):
        continue
    f = line.rstrip("\n").split("\t")
    NAMED[f[0]] = f[3]


def census(seq):
    return {k: sum(1 for r in pats for _ in r.finditer(seq)) for k, pats in PATS.items()}


def pairs(s):
    st, out = [], set()
    for i, c in enumerate(s):
        if c == "(":
            st.append(i)
        elif c == ")":
            out.add((st.pop(), i))
    return out


def dist(a, b):
    return len(a ^ b) / max(1, len(a) + len(b))


def path(cell, nest):
    if nest == "F0":
        for d in ("f0", "f0-l", "f0-five", "f0-xl"):
            p = f"{E}/screen/{d}/{cell}.jsonl"
            if os.path.exists(p):
                return p
        return None
    p = NAMED.get(f"{cell}_{nest}")
    if p and os.path.exists(p):
        return p
    c = glob.glob(f"{E}/*/{cell}_{nest}.r*.jsonl")
    if c:
        return max(c, key=lambda q: int(q.rsplit(".r", 1)[1].split(".")[0]))
    hit = glob.glob(f"{E}/*/{cell}_{nest}.jsonl")
    return hit[0] if hit else None


def load(p):
    rows = [json.loads(line) for line in open(p) if line.strip()]
    return rows[0], [r for r in rows[1:] if r.get("record") == "point"]


