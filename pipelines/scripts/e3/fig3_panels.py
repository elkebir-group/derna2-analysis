"""Figure 3's panels, from the extracted cells.

    python fig3_panels.py CELLS_JSON OUT_JSON

Three rules decide what is comparable here. A cell stopped on its budget still measured a
real peak RSS and a real wall, so cost panels keep it and mark it. A scaling exponent is
fitted per build and only where that build spans enough length to carry one, because the
builds sit in different length bands and a pooled fit reads their offsets as slope. And a
nest comparison is made inside one protein on one build, never across either.
"""

import json
import math
import statistics as st
import sys
from collections import defaultdict

NESTS = ["F1", "F2", "F3", "F4"]
SPAN = 1.5  # a build must cover this length ratio before an exponent is fitted
MIN_N = 8


def fit(pts):
    xs = [math.log10(a) for a, _ in pts]
    ys = [math.log10(b) for _, b in pts]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=True)) / sxx
    icept = my - slope * mx
    ss = sum((y - my) ** 2 for y in ys)
    res = sum((y - (icept + slope * x)) ** 2 for x, y in zip(xs, ys, strict=True))
    return {
        "exp": slope,
        "icept": icept,
        "r2": 1 - res / ss if ss else 1.0,
        "n": n,
        "lo": min(x for x, _ in pts),
        "hi": max(x for x, _ in pts),
    }


def main():
    cells = json.load(open(sys.argv[1]))["cells"]

    grid = defaultdict(dict)
    meta = {}
    for c in cells:
        key = (c["band"], c["protein"])
        meta[key] = c["aa"]
        rank = {"done": 3, "budget": 2, "stopped": 1, "failed": 0}
        prev = grid[key].get(c["nest"])
        now = "complete" if c["complete"] else c["state"]
        if prev is None or rank.get(c["state"], 0) >= rank.get(prev["raw"], 0):
            grid[key][c["nest"]] = {"state": now, "raw": c["state"], "points": c["points"]}
    coverage = [
        {
            "band": b,
            "protein": p,
            "aa": meta[(b, p)],
            "nests": {n: grid[(b, p)].get(n, {"state": "none"})["state"] for n in NESTS},
        }
        for (b, p) in sorted(grid, key=lambda k: meta[k])
    ]

    mem, wall = [], []
    by_build = defaultdict(list)
    for c in cells:
        if c["rss_kb"]:
            mem.append(
                {
                    "aa": c["aa"],
                    "gb": c["rss_kb"] / 1048576,
                    "build": c["build"],
                    "state": c["state"],
                }
            )
            by_build[c["build"]].append((c["aa"], c["rss_kb"] / 1048576))
        if c["wall_s"]:
            wall.append({"aa": c["aa"], "s": c["wall_s"], "state": c["state"]})
    fits = {}
    for build, pts in by_build.items():
        lo, hi = min(a for a, _ in pts), max(a for a, _ in pts)
        if len(pts) >= MIN_N and hi / lo >= SPAN:
            fits[build] = fit(pts)
    pooled = fit([(m["aa"], m["gb"]) for m in mem])

    series = []
    paired = defaultdict(dict)
    for c in cells:
        if c["points"] is None:
            continue
        key = (c["protein"], c["build"])
        prev = paired[key].get(c["nest"])
        if prev is None or c["points"] > prev["points"]:
            paired[key][c["nest"]] = c
    for (protein, build), got in sorted(paired.items()):
        if "F1" not in got or len(got) < 3:
            continue
        base = got["F1"]["points"]
        series.append(
            {
                "protein": protein,
                "build": build,
                "aa": got["F1"]["aa"],
                "ratio": {n: got[n]["points"] / base for n in NESTS if n in got},
                "points": {n: got[n]["points"] for n in NESTS if n in got},
                "complete": sum(1 for c in got.values() if c["complete"]),
            }
        )
    medians = {}
    for n in NESTS:
        got = [s["ratio"][n] for s in series if n in s["ratio"]]
        if got:
            medians[n] = {"med": st.median(got), "n": len(got), "lo": min(got), "hi": max(got)}

    json.dump(
        {
            "coverage": coverage,
            "mem": mem,
            "wall": wall,
            "fits": fits,
            "pooled": pooled,
            "series": series,
            "medians": medians,
            "meta": {
                "cells": len(cells),
                "complete": sum(1 for c in cells if c["complete"]),
                "with_front": sum(1 for c in cells if c["points"]),
                "censored": sum(1 for c in cells if c["state"] in ("budget", "stopped")),
                "builds": len(by_build),
                "aa_lo": min(c["aa"] for c in cells),
                "aa_hi": max(c["aa"] for c in cells),
            },
        },
        open(sys.argv[2], "w"),
        indent=1,
    )

    print(f"coverage rows {len(coverage)}  series {len(series)}  fits {len(fits)}")
    print(
        f"pooled exponent {pooled['exp']:.2f} (R2 {pooled['r2']:.2f}, n={pooled['n']}) "
        f"— across {len(by_build)} builds, not a law"
    )
    for b, f in sorted(fits.items(), key=lambda kv: -kv[1]["n"]):
        print(f"  {b}: exp {f['exp']:.2f}  R2 {f['r2']:.2f}  n={f['n']}  {f['lo']}-{f['hi']} aa")
    print(
        "median points found, relative to F1:", {n: round(v["med"], 3) for n, v in medians.items()}
    )


if __name__ == "__main__":
    main()
