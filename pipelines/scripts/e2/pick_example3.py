"""Screen for an example whose ladder removes all three motif sets visibly.

An example plate is a staircase only if each rung has something to take away. The
requirement here is on the drawn designs, not on the price: the F0 pick must spell a site of
set I, of set II and of set III; the F1 pick must still carry set II and set III, since set I
is now forbidden; the F2 pick must still carry set III; and the F3 pick must carry none. Then
each panel differs from the last by a colour that leaves, and the final panel is empty.

That is strictly stronger than "set II binds". A set can be priced above zero on a protein and
still be absent from the one design a bound selects, in which case the figure shows nothing
happening at that rung however real the cost.

The structural signature from the earlier screen is kept as a secondary ranking: F0-F2 folding
alike, F3 not, and the closest F3-feasible design at the same cutoff back toward F2.
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


cells = sorted({os.path.basename(p)[:-6] for d in ("f0", "f0-l", "f0-five", "f0-xl")
                for p in glob.glob(f"{E}/screen/{d}/*.jsonl")})
hits, near = [], []
for cell in cells:
    paths = {n: path(cell, n) for n in NESTS}
    if not all(paths.values()):
        continue
    try:
        fronts = {n: load(paths[n]) for n in NESTS}
    except Exception:
        continue
    L = len(fronts["F0"][0]["protein"]["sequence"])

    for cut in CUTOFFS:
        picks = {}
        for n in NESTS:
            at = [r for r in fronts[n][1] if r["cai"] >= cut]
            if not at:
                picks = None
                break
            picks[n] = min(at, key=lambda r: r["mfe"])
        if not picks:
            continue
        cen = {n: census(picks[n]["sequence"]) for n in NESTS}
        staircase = (all(cen["F0"][k] for k in ("I", "II", "III"))
                     and cen["F1"]["II"] and cen["F1"]["III"] and cen["F2"]["III"]
                     and sum(cen["F3"].values()) == 0)
        pr = {n: pairs(picks[n]["structure"]) for n in NESTS}
        d01, d12, d23 = (dist(pr["F0"], pr["F1"]), dist(pr["F1"], pr["F2"]),
                         dist(pr["F2"], pr["F3"]))
        best = None
        for r in fronts["F3"][1]:
            if r["cai"] < cut:
                continue
            prem = r["mfe"] - picks["F3"]["mfe"]
            if prem <= 0 or prem > SLACK:
                continue
            d = dist(pr["F2"], pairs(r["structure"]))
            if best is None or d < best[0]:
                best = (d, prem)
        rec = {"cell": cell, "L": L, "cut": cut, "census": cen, "d01": d01, "d12": d12,
               "d23": d23, "d_alt": None if best is None else best[0],
               "premium": None if best is None else best[1],
               "mfe": {n: picks[n]["mfe"] for n in NESTS}}
        (hits if staircase else near).append(rec)

json.dump({"staircase": hits, "other": near}, open(f"{OUT}/pick_example3.json", "w"), indent=1)


def row(x):
    c = x["census"]
    alt = "   —" if x["d_alt"] is None else f"{x['d_alt']:.2f}"
    prem = "     —" if x["premium"] is None else f"{x['premium']:+5.1f}"
    return (f"  {x['cell']:8} {x['L']:4} {x['cut']:4.2f} | "
            + " ".join(f"{c[n]['I']}/{c[n]['II']}/{c[n]['III']}" for n in NESTS)
            + f" | {x['d01']:4.2f} {x['d12']:4.2f} {x['d23']:4.2f} | {alt} {prem} | "
            + f"{x['mfe']['F3'] - x['mfe']['F0']:+5.1f}")


print("sites in the drawn design, as set I / II / III, at F0  F1  F2  F3\n")
print(f"  {'cell':8} {'aa':>4} {'cut':>4} | F0    F1    F2    F3    | "
      f"{'d01':>4} {'d12':>4} {'d23':>4} | d_alt  prem | F0→F3")
print(f"\nFULL STAIRCASE — all three sets drawn at F0, each leaving in turn ({len(hits)}):")
for x in sorted(hits, key=lambda x: -(x["d23"] - (x["d_alt"] or x["d23"]))):
    print(row(x))
print("\nAll three sets present at F0 but the staircase breaks somewhere:")
for x in sorted(near, key=lambda x: -sum(1 for k in ("I", "II", "III") if x["census"]["F0"][k])):
    if all(x["census"]["F0"][k] for k in ("I", "II", "III")):
        print(row(x))
