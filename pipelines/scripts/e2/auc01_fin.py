"""The unit-square AUC on E2's finished set, read from the band owner's FINISHED.tsv.

A cell is finished when it is `complete` (the front is exact) or `retired` (it passed the
gate: share <= 0.80, hull_gap <= 0.20 and unvisited span <= 0.10 — bounded, not exact). The
TSV names the exact round file to read, which is the newest round with a clean exit rather
than whatever file is newest on disk: an in-flight checkpoint omits the boxes the parking
wrapper set aside, so its header reads falsely settled.

A retired front can only gain points, so its AUC here is a LOWER bound on the true one and
the drop from F0 an UPPER bound.

Finished is decided **per rung**: a protein enters this set once every rung the band owner
reports it over has finished, and it is absent from the rungs CELL-SET.md excludes. A
protein-level filter instead dropped a protein from the two rungs the owner retained
because a third was excluded, so the finished count read one lower than the owner's. Every
figure over these rows therefore carries its own rung denominator, and `rungs(row)` is the
only way to ask which nests a row has.
"""
import glob
import json
import math
import os

from cellset import EXCLUDE
from cellset import rungs as _rungs
from e3paths import OUT, E  # noqa: F401

NESTS = ["F0", "F1", "F2", "F3"]
N = 2001
U = [i / (N - 1) for i in range(N)]

rows = [line.rstrip("\n").split("\t") for line in open(f"{E}/FINISHED.tsv")
        if not line.startswith("#") and line.strip()]
hdr, tsv = rows[0], {}
for r in rows[1:]:
    d = dict(zip(hdr, r, strict=True))
    cell, nest = d["cell"].rsplit("_", 1)
    tsv.setdefault(cell, {})[nest] = d
# The rungs each protein is reported over, and the finished set: every reported rung done.
REPORTED = {c: _rungs(c, NESTS[1:]) for c in tsv}
FINISHED = sorted(c for c, v in tsv.items()
                  if REPORTED[c] and all(v.get(n, {}).get("status") in ("complete", "retired")
                                         for n in REPORTED[c]))
DROPPED = {(c, n): why for (c, n), why in EXCLUDE.items() if c in FINISHED and n in NESTS}


def rungs(row):
    """The nests a finished row carries, F0 first. Excluded rungs are simply absent."""
    return ["F0"] + REPORTED[row["cell"] if isinstance(row, dict) else row]


def load(path):
    pts = []
    for line in open(path):
        if not line.strip():
            continue
        r = json.loads(line)
        if r.get("record") == "point":
            pts.append((r["cai"], r["mfe"]))
    return pts


F0PATH = {os.path.basename(p)[:-6]: p for d in ("f0", "f0-l", "f0-five", "f0-xl")
          for p in glob.glob(f"{E}/screen/{d}/*.jsonl")}
FR = {c: {"F0": {"pts": load(F0PATH[c]), "status": "complete"}} for c in FINISHED}
for c in FINISHED:
    for n in REPORTED[c]:
        d = tsv[c][n]
        FR[c][n] = {"pts": load(d["round_file"]), "status": d["status"],
                    "share": float(d["share"]), "span": float(d["unvisited_span"])}


def gcurve(fr, xs):
    """min MFE among points with CAI >= x, for x ascending; inf where infeasible."""
    pts = sorted(fr["pts"], key=lambda p: -p[0])
    out, i, best = [], 0, math.inf
    for x in sorted(xs, reverse=True):
        while i < len(pts) and pts[i][0] >= x:
            best = min(best, pts[i][1])
            i += 1
        out.append(best)
    return out[::-1]


def ycurve(fr, xs, box):
    _, _, e_ideal, e_nadir = box
    return [0.0 if math.isinf(v) else max(0.0, min(1.0, (e_nadir - v) / (e_nadir - e_ideal)))
            for v in gcurve(fr, xs)]


def main():
    out = []
    for c in FINISHED:
        f0 = FR[c]["F0"]["pts"]
        c_lo, c_hi = min(p[0] for p in f0), max(p[0] for p in f0)
        e_ideal = min(p[1] for p in f0)
        e_nadir = max(m for a, m in f0 if a >= c_hi - 1e-12)
        box = (c_lo, c_hi, e_ideal, e_nadir)
        xs = [c_lo + u * (c_hi - c_lo) for u in U]
        r = {"cell": c, "cai_lo": c_lo, "cai_hi": c_hi, "e_ideal": e_ideal, "e_nadir": e_nadir,
             "band": tsv[c]["F1"]["band"],
             "rungs": ["F0"] + REPORTED[c],
             "exact": all(FR[c][n]["status"] == "complete" for n in REPORTED[c])}
        for n in r["rungs"]:
            y = ycurve(FR[c][n], xs, box)
            r[n] = sum(y) / len(y)
            r[n + "_status"] = FR[c][n]["status"]
            r[n + "_pts"] = len(FR[c][n]["pts"])
        out.append(r)
    json.dump(out, open(f"{OUT}/auc01_fin.json", "w"), indent=1)
    return out


if __name__ == "__main__":
    import statistics as st
    res = main()
    ex = [r for r in res if r["exact"]]
    print(f"{len(res)} finished proteins ({len(ex)} exact on every rung they are reported "
          f"over), {sum(1 for r in res if r['band'] == 'short')} short band")
    for (c, n), why in sorted(DROPPED.items()):
        print(f"  {c} is reported over {', '.join(REPORTED[c])} only: {n} {why}")
    for n in NESTS:
        v = [r[n] for r in res if n in r["rungs"]]
        print(f"  {n}: median AUC {st.median(v):.4f}  range {min(v):.4f}-{max(v):.4f}"
              f"  over {len(v)}")
    for a, b in zip(NESTS, NESTS[1:], strict=False):   # consecutive pairs
        at = [r for r in res if a in r["rungs"] and b in r["rungs"]]
        d = [r[a] - r[b] for r in at]
        de = [r[a] - r[b] for r in at if r["exact"]]
        print(f"  {a}->{b}: median drop {st.median(d):.4f} over {len(d)}  "
              f"| {st.median(de):.4f} over the exact {len(de)}  max {max(d):.4f}")
