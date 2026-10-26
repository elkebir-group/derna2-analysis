"""Figure 2's data, rebuilt from a sweep row file.

    python fig2_data.py SWEEP_JSONL PIN OUT_JSON

Only rows stamped with PIN are read: a sweep file accumulates across re-pins, and rows of
two builds are not one experiment.

Panel (a)'s λ-unreachable share is computed from the walk's own front points, not read
from the run's `vertices` field: that field is what the dichotomic search accepted at its
own tolerance, and the energies store to 0.1 kcal/mol.
"""

import json
import math
import statistics as st
import sys
from collections import defaultdict

from derna2_analysis.fronts import hull_vertices, pareto

ARMS = [
    ("derna2-walk", "DERNA2 walk"),
    ("derna2-hull", "DERNA2 hull"),
    ("mrnafold-hull", "mRNAfold"),
    ("lineardesign-hull", "LinearDesign"),
    ("lcds-exact-hull", "LinearCDSfold"),
    ("derna1-own", "DERNA1 own λ"),
]
BANDS = [(0, 150), (150, 300), (300, 600), (600, 10**6)]


def quartiles(xs):
    xs = sorted(xs)
    q = st.quantiles(xs, n=4) if len(xs) > 1 else [xs[0], xs[0], xs[0]]
    return {"lo": xs[0], "q1": q[0], "med": st.median(xs), "q3": q[2], "hi": xs[-1], "n": len(xs)}


def read(path, pin):
    """Cells, the walk's front points, and per-cell design sets and peak RSS.

    Only the walk's cell row carries `front_points` and `maxrss_kb`; for the weight arms
    both have to come from the per-λ rows, which is also where a design is actually
    recorded. Designs are counted as distinct rescored (CAI, MFE) pairs, so two λs that
    return the same design count once.
    """
    cells, points = {}, defaultdict(list)
    designs, rss = defaultdict(set), defaultdict(int)
    for line in open(path):
        if not line.strip():
            continue
        r = json.loads(line)
        if r.get("pin") != pin:
            continue
        kind, key = r.get("kind"), (r.get("arm"), r.get("protein"))
        if kind == "cell":
            cells[key] = r
            continue
        cai = r.get("cai_rescored") or r.get("cai") or r.get("cai_reported")
        mfe = r.get("mfe_rescored")
        if mfe is None:
            mfe = r.get("mfe") if r.get("mfe") is not None else r.get("mfe_reported")
        if kind == "front_point":
            if r["arm"] == "derna2-walk" and cai and cai > 0 and mfe is not None:
                points[r["protein"]].append((math.log(cai), mfe))
        else:
            if r.get("maxrss_kb"):
                rss[key] = max(rss[key], r["maxrss_kb"])
            if cai is not None and mfe is not None and r.get("exit", "ok") == "ok":
                designs[key].add((round(cai, 9), round(mfe, 4)))
    for key, cell in cells.items():
        cell["_designs"] = cell.get("front_points") or len(designs.get(key, ()))
        cell["_rss_kb"] = cell.get("maxrss_kb") or rss.get(key, 0)
    return cells, points


def unreachable(points, cells):
    """Share of each exact front no weight reaches, and the disagreement with the field."""
    rows, mismatch = [], []
    for protein, pts in points.items():
        cell = cells.get(("derna2-walk", protein))
        if not cell or cell.get("exit") != "ok" or len(pts) < 3:
            continue
        # Rescoring can dominate a point that was on the front under reported energies,
        # so the denominator is the front of the rescored points, not the rows returned.
        front = len(pareto(pts))
        computed = len(hull_vertices(pts))
        stored = cell.get("vertices")
        if stored is not None and stored != computed:
            mismatch.append({"protein": protein, "stored": stored, "computed": computed})
        rows.append(
            {
                "protein": protein,
                "n": cell["n_aa"],
                "returned": len(pts),
                "points": front,
                "vertices": computed,
                "share": 100.0 * (front - computed) / front,
            }
        )
    return rows, mismatch


def main():
    sweep, pin, out = sys.argv[1:]
    cells, points = read(sweep, pin)
    if not cells:
        raise SystemExit(f"{sweep}: no cell rows stamped with pin {pin}")
    reach, mismatch = unreachable(points, cells)

    ok = {
        a: {p for (x, p), c in cells.items() if x == a and c.get("exit") == "ok"} for a, _ in ARMS
    }
    paired = sorted(set.intersection(*ok.values()))

    rows = []
    for arm, label in ARMS:
        got = [cells[(arm, p)] for p in paired]
        pts = [c["_designs"] for c in got]
        rows.append(
            {
                "arm": arm,
                "label": label,
                "wall": quartiles([c["wall_s"] for c in got]),
                "gb": quartiles([c["_rss_kb"] / 1048576 for c in got]),
                "points": quartiles(pts),
                "per": quartiles([c["wall_s"] / n for c, n in zip(got, pts, strict=True) if n]),
                "completed": len(ok[arm]),
                "timeouts": sum(
                    1 for (x, _), c in cells.items() if x == arm and c.get("exit") == "timeout"
                ),
                "oom": sum(1 for (x, _), c in cells.items() if x == arm and c.get("exit") == "oom"),
                "attempted": sum(1 for (x, _) in cells if x == arm),
                "shortest": min((c["n_aa"] for (x, _), c in cells.items() if x == arm), default=0),
                "longest": max(
                    (
                        c["n_aa"]
                        for (x, _), c in cells.items()
                        if x == arm and c.get("exit") == "ok"
                    ),
                    default=0,
                ),
            }
        )

    bands = []
    for lo, hi in BANDS:
        got = [r["share"] for r in reach if lo < r["n"] <= hi]
        if got:
            bands.append({"band_lo": lo, "band_hi": hi, **quartiles(got)})

    n_paired = [cells[("derna2-walk", p)]["n_aa"] for p in paired]
    json.dump(
        {
            "rows": rows,
            "bands": bands,
            "reach": reach,
            "mismatch": mismatch,
            "meta": {
                "paired": len(paired),
                "n_lo": min(n_paired),
                "n_hi": max(n_paired),
                "fronts": len(reach),
                "cells": len(cells),
                "pin": pin,
            },
        },
        open(out, "w"),
        indent=1,
    )
    print(
        f"paired {len(paired)} ({min(n_paired)}-{max(n_paired)} aa), fronts {len(reach)}, "
        f"cells {len(cells)}, vertex-field mismatches {len(mismatch)}"
    )
    for r in rows:
        print(
            f"  {r['label']:<14} wall {r['wall']['med']:>8.1f}  GB {r['gb']['med']:>6.3f}  "
            f"pts {r['points']['med']:>7.1f}  s/pt {r['per']['med']:>8.3f}  "
            f"ok {r['completed']:>3d}  t/o {r['timeouts']:>3d}  max {r['longest']}"
        )


if __name__ == "__main__":
    main()
