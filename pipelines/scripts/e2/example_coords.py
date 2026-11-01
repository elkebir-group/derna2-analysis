"""Stage one of the Q92734 example plate: pick the five designs and lay out their folds.

Run under the ViennaRNA interpreter. P0's row drawn on a different instance and a different
bound, with the fifth panel chosen for a different reason.

The five designs, all at CAI >= 0.90:

  F0..F3   the min-MFE design at that bound on each rung's front, as P0's (c)-(f).
  F3'      another design on the F3 front that satisfies the same bound -- so an F3-feasible
           design, not an optimum -- chosen as the one whose fold is closest to F2's.

P0's F3' is an off-front design found by searching slack above the F3 optimum. This one is
on the front and costs nothing to find, which makes it a weaker search and a stronger claim:
the alternative was already among the designs the walk had enumerated.
"""
import glob
import json
import os
import sys

import RNA
from e3paths import OUT, E  # noqa: F401

CELL = sys.argv[1] if len(sys.argv) > 1 else "Q92734"
BOUND = float(sys.argv[2]) if len(sys.argv) > 2 else 0.90
NESTS = ["F0", "F1", "F2", "F3"]
# The rung the second design is drawn at, and the rung it is measured back toward: the
# break is not always at the deepest rung, and the alternative is only interesting at the
# rung that broke.
ALT = sys.argv[3] if len(sys.argv) > 3 else "F3"
# "none" draws the ladder alone -- four rungs, no second design -- and measures every fold
# against F0, the design the constraints are removed from.
# Where the second design comes from. "front" takes the closest point the walk already
# enumerated; "constructed" and "refined" take the design built by construct_alt.py and
# improved by refine_alt.py, which are off the front and so are a wider search at a price.
SRC = sys.argv[4] if len(sys.argv) > 4 else "front"
# The fold every panel is coloured and measured against. It defaults to the rung just left of
# ALT, which is the right anchor when the question is "did THIS rung's rearrangement have to
# happen". It is the wrong one when the second design was reached through an earlier rung's
# second design: an F3 design constructed from the F2 alternative keeps F1's shape, and
# scoring it against F2 -- the rung that destroyed that shape -- reports 0.666 for a fold
# that is 0.032 from the one it was built to keep. A pinned anchor can only answer the
# question the default asks, so it is an argument.
REFI = 0 if ALT == "none" else (NESTS.index(sys.argv[5]) if len(sys.argv) > 5
                                else NESTS.index(ALT) - 1)
REF = NESTS[REFI]
IUPAC = {"R": "[AG]", "Y": "[CU]", "N": "[ACGU]", "S": "[CG]", "W": "[AU]", "K": "[GU]",
         "M": "[AC]"}
SETS = {"I": "set1-manufacturing", "II": "set2-processing", "III": "set3-repeats"}
PATS = {k: [__import__("re").compile("".join(IUPAC.get(c, c) for c in
                                             line.split("#")[0].strip().upper().replace("T", "U")))
            for line in open(f"{E}/motifs/{f}.txt") if line.split("#")[0].strip()]
        for k, f in SETS.items()}

NAMED = {}
for line in open(f"{E}/FINISHED.tsv"):
    if line.startswith("#") or line.startswith("cell\t"):
        continue
    f = line.rstrip("\n").split("\t")
    NAMED[f[0]] = f[3]


def sites(seq):
    return sorted((m.start(), m.end(), roman)
                  for roman, pats in PATS.items() for r in pats for m in r.finditer(seq))


def pairs(s):
    st, out = [], []
    for i, c in enumerate(s):
        if c == "(":
            st.append(i)
        elif c == ")":
            out.append((st.pop(), i))
    return out


def dist(a, b):
    a, b = set(map(tuple, a)), set(map(tuple, b))
    return len(a ^ b) / max(1, len(a) + len(b))


def path(nest):
    if nest == "F0":
        for d in ("f0", "f0-l", "f0-five", "f0-xl"):
            p = f"{E}/screen/{d}/{CELL}.jsonl"
            if os.path.exists(p):
                return p
    p = NAMED.get(f"{CELL}_{nest}")
    if p and os.path.exists(p):
        return p
    c = glob.glob(f"{E}/*/{CELL}_{nest}.r*.jsonl")
    if c:
        return max(c, key=lambda q: int(q.rsplit(".r", 1)[1].split(".")[0]))
    return glob.glob(f"{E}/*/{CELL}_{nest}.jsonl")[0]


def load(p):
    rows = [json.loads(line) for line in open(p) if line.strip()]
    return rows[0], [r for r in rows[1:] if r.get("record") == "point"]


# P0 draws its structures in ViennaRNA's naive radial layout, not naview: RNApuzzler lays the
# constrained nests out several times flatter than the unconstrained one off a few base pairs,
# so the panels would differ for a reason unrelated to the claim. Both layouts are stored; the
# figure chooses.
PLOT_TYPE = 0


def layout0(dot):
    """The same coordinates P0 uses, normalised the same way as the naview ones below."""
    RNA.cvar.rna_plot_type = PLOT_TYPE
    c = RNA.get_xy_coordinates(dot)
    xs = [c.get(i).X for i in range(len(dot))]
    ys = [c.get(i).Y for i in range(len(dot))]
    cx, cy = sum(xs) / len(xs), sum(ys) / len(ys)
    step = sum(((xs[i + 1] - xs[i]) ** 2 + (ys[i + 1] - ys[i]) ** 2) ** .5
               for i in range(len(xs) - 1)) / max(1, len(xs) - 1)
    return [[(x - cx) / step, (y - cy) / step] for x, y in zip(xs, ys, strict=True)]


def layout(dot):
    c = RNA.naview_xy_coordinates(dot)
    xs = [c[i].X for i in range(len(dot))]
    ys = [c[i].Y for i in range(len(dot))]
    cx, cy = sum(xs) / len(xs), sum(ys) / len(ys)
    step = sum(((xs[i + 1] - xs[i]) ** 2 + (ys[i + 1] - ys[i]) ** 2) ** .5
               for i in range(len(xs) - 1)) / max(1, len(xs) - 1)
    return [[(x - cx) / step, (y - cy) / step] for x, y in zip(xs, ys, strict=True)]


def design(r, label, src):
    return {"label": label, "cai": r["cai"], "mfe": r["mfe"], "dot": r["structure"],
            "xy": layout(r["structure"]), "xy0": layout0(r["structure"]),
            "pairs": pairs(r["structure"]),
            "sites": sites(r["sequence"]), "nt": len(r["sequence"]), "seq": r["sequence"],
            "source": src}


fronts = {n: load(path(n)) for n in NESTS}
picks = {}
for n in NESTS:
    at = [r for r in fronts[n][1] if r["cai"] >= BOUND]
    picks[n] = min(at, key=lambda r: r["mfe"])

panels = [design(picks[n], n, os.path.basename(path(n))) for n in NESTS]
refp = panels[REFI]["pairs"]
if ALT == "none":
    pass
elif SRC == "front":
    feas = [r for r in fronts[ALT][1]
            if r["cai"] >= BOUND and r["mfe"] > picks[ALT]["mfe"]]
    alt = min(feas, key=lambda r: dist(refp, pairs(r["structure"])))
    src = os.path.basename(path(ALT))
else:
    aslug = CELL + ("" if ALT == "F3" else f"-{ALT.lower()}")
    f = f"{OUT}/alt-{aslug}.json" if SRC == "constructed" else f"{OUT}/refined-{CELL}.json"
    rec = json.load(open(f))
    r0 = rec["best"][0] if SRC == "constructed" else rec
    alt = {"cai": r0["cai"], "mfe": r0["mfe"], "sequence": r0["sequence"],
           "structure": r0.get("dot") or r0["structure"]}
    src = os.path.basename(f)
if ALT != "none":
    panels.append(design(alt, ALT + "′", src))

# The codons each design changes against the unconstrained one, by position, for the figure's
# codon strip. The designs encode one protein, so the codons correspond one to one.
seqs = [picks[n]["sequence"] for n in NESTS] + ([alt["sequence"]] if ALT != "none" else [])
codons = [[q[i:i + 3] for i in range(0, len(q), 3)] for q in seqs]
for p, c in zip(panels, codons, strict=True):
    p["n_codons"] = len(c)
    p["codons_differ"] = [i for i, (x, y) in enumerate(zip(codons[0], c, strict=True)) if x != y]

for a, b in zip(panels, panels[1:], strict=False):   # consecutive pairs
    b["dist_from_left"] = dist(a["pairs"], b["pairs"])
for p in panels:
    p["dist_to_F2"] = dist(refp, p["pairs"])
    p["rise"] = p["mfe"] - panels[0]["mfe"]

h0 = fronts["F0"][0]
name = h0["protein"]["name"]
rec = {"cell": CELL, "bound": BOUND, "ref": REFI, "alt_rung": ALT,
       "gene": name.split("GN=")[1].split()[0] if "GN=" in name else CELL,
       "L": len(h0["protein"]["sequence"]), "panels": panels,
       "fronts": {n: sorted((r["cai"], r["mfe"]) for r in fronts[n][1]) for n in NESTS},
       "vienna": RNA.__version__, "plot_type": PLOT_TYPE,
       "npts": {n: len(fronts[n][1]) for n in NESTS}}
SLUG = (f"{CELL}-c{round(BOUND * 100):02d}" + ("" if ALT == "F3" else f"-{ALT.lower()}")
        + ("" if SRC == "front" or ALT == "none" else f"-{SRC[:4]}"))
json.dump(rec, open(f"{OUT}/example-{SLUG}.json", "w"))
json.dump(rec, open(f"{OUT}/example-{CELL}.json", "w"))

print(f"{CELL} {rec['gene']} {rec['L']} aa / {panels[0]['nt']} nt, CAI >= {BOUND}")
for p in panels:
    print(f"  {p['label']:4} CAI {p['cai']:.4f}  MFE {p['mfe']:8.1f}  rise {p['rise']:+5.1f}  "
          f"{len(p['pairs']):4} pairs  {len(p['sites']):2} sites  "
          f"d(from left) {p.get('dist_from_left', float('nan')):.3f}  "
          f"d(to {REF}) {p['dist_to_F2']:.3f}")
if ALT != "none":
    print(f"  {ALT}′ premium over the {ALT} optimum: "
          f"{panels[4]['mfe'] - picks[ALT]['mfe']:+.2f} kcal/mol")
print(f"  (distances above are to {REF})")
