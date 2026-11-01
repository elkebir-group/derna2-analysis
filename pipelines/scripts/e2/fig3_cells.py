"""Figure 3's data: what E2 has measured so far, read from the three bands' summaries.

    python fig3_cells.py E2_DIR OUT_JSON

Read-only. The driver writes a four-core mask unquoted, which splits one field into four
and shifts every later column left by three; such a row is repaired here the way
drivers/e2/fix_summary.py describes, and only where it is unambiguous. A row whose result
columns are empty is a cell that was stopped before it wrote anything: it keeps its wall
and its peak RSS, which are real, and contributes no front.
"""

import csv
import json
import math
import sys
from pathlib import Path

HEAD = [
    "protein",
    "aa",
    "nest",
    "sets",
    "build",
    "threads",
    "core",
    "exit",
    "wall_s",
    "peak_rss_kb",
    "points",
    "hull",
    "complete",
    "open_boxes",
    "load",
]
# What a non-zero exit means for the cell, and whether its wall is a runtime or a bound.
STOPPED = {"124": "budget", "15": "stopped", "143": "stopped"}


def repair(row):
    """A core mask split across k+1 fields, rejoined — or None if the row is not that."""
    extra = len(row) - len(HEAD)
    if extra <= 0:
        return row if len(row) == len(HEAD) else None
    mask = row[6 : 7 + extra]
    if not all(f.isdigit() and 0 <= int(f) < 256 for f in mask):
        return None
    return row[:6] + [",".join(mask)] + row[7 + extra :]


def band(path, header):
    out = []
    with open(path) as f:
        for i, raw in enumerate(csv.reader(f)):
            if (header and i == 0) or not raw:
                continue
            row = repair(raw)
            if row is None:
                print(f"  skipped unrepairable row in {path}: {len(raw)} fields", file=sys.stderr)
                continue
            r = dict(zip(HEAD, row, strict=True))
            yield_row = {
                "protein": r["protein"],
                "aa": int(r["aa"]),
                "nest": r["nest"],
                "build": r["build"],
                "exit": r["exit"],
                "wall_s": float(r["wall_s"]) if r["wall_s"] else None,
                "rss_kb": int(r["peak_rss_kb"]) if r["peak_rss_kb"] else None,
                "points": int(r["points"]) if r["points"] else None,
                "hull": int(r["hull"]) if r["hull"] else None,
                "complete": r["complete"] == "True",
                "state": "done" if r["exit"] == "0" else STOPPED.get(r["exit"], "failed"),
            }
            out.append(yield_row)
    return out


def long_band(path):
    """The long band writes its result fields as a `;`-joined blob in the last column."""
    out = []
    for raw in csv.reader(open(path)):
        if not raw:
            continue
        protein, aa, nest, build, _threads, _core, exit_, wall, rss, blob = raw[:10]
        kv = dict(p.split("=", 1) for p in blob.split(";") if "=" in p)
        out.append(
            {
                "protein": protein,
                "aa": int(aa),
                "nest": nest,
                "build": build,
                "exit": exit_,
                "wall_s": float(wall) if wall else None,
                "rss_kb": int(rss) if rss else None,
                "points": int(kv["points"]) if "points" in kv else None,
                "hull": int(kv["hull"]) if "hull" in kv else None,
                "complete": kv.get("complete") == "True",
                "state": "done" if exit_ == "0" else STOPPED.get(exit_, "failed"),
                "share": float(kv["share"]) if "share" in kv else None,
                "gaps": int(kv["gaps"]) if "gaps" in kv else None,
            }
        )
    return out


def fit(points):
    """Least squares on log10 x, log10 y: the exponent and its R²."""
    xs = [math.log10(x) for x, _ in points]
    ys = [math.log10(y) for _, y in points]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=True)) / sxx
    a = my - b * mx
    ss = sum((y - my) ** 2 for y in ys)
    res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys, strict=True))
    return {"exponent": b, "intercept": a, "r2": 1 - res / ss if ss else 1.0, "n": n}


def main():
    root, out = Path(sys.argv[1]), sys.argv[2]
    cells = []
    for name, path, header in [
        ("short", root / "grid" / "summary.csv", True),
        ("mid", root / "mid" / "summary.csv", False),
    ]:
        for r in band(path, header):
            cells.append({**r, "band": name})
    for r in long_band(root / "long" / "summary.csv"):
        cells.append({**r, "band": "long"})

    # Peak RSS is measured whether or not a cell finished, so the memory law uses them all.
    mem = [(c["aa"], c["rss_kb"] / 1048576) for c in cells if c["rss_kb"]]
    json.dump({"cells": cells, "mem_fit": fit(mem), "mem": sorted(mem)}, open(out, "w"), indent=1)

    done = sum(1 for c in cells if c["complete"])
    print(
        f"cells {len(cells)}  complete {done}  with a front {sum(1 for c in cells if c['points'])}"
    )
    law = fit(mem)
    print(f"memory law: peak GB ~ aa^{law['exponent']:.2f}  R2 {law['r2']:.3f}  n={len(mem)}")
    for b in ("short", "mid", "long"):
        got = [c for c in cells if c["band"] == b]
        print(
            f"  {b:6s} {len(got):3d} cells, {len({c['protein'] for c in got})} proteins, "
            f"{min(c['aa'] for c in got)}-{max(c['aa'] for c in got)} aa, "
            f"complete {sum(1 for c in got if c['complete'])}"
        )


if __name__ == "__main__":
    main()
