"""The planned Figure 3: nested fronts, the price at matched CAI, and the cost ratio.

    python fig3_planned.py E2_DIR OUT_JSON

The price of a blacklist is read against the unconstrained front of the same protein, at
equal CAI: for a constrained point, the baseline's best MFE among its points of at least
that CAI. A constrained front that stopped early is missing points, not holding wrong
ones, so a price at a CAI it did reach is exact while its coverage of CAI is partial.

A cell whose blacklist removes nothing from the baseline front is dropped from the price
and kept in the cost, which is the plan's rule and two different filters over one table.
"""

import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

SETS = {
    "F1": ["manufacturing"],
    "F2": ["manufacturing", "processing"],
    "F3": ["manufacturing", "processing", "repeats"],
    "F4": ["manufacturing", "processing", "repeats", "initiation"],
}
NESTS = list(SETS)


def front(path):
    """(header, points sorted by increasing CAI)."""
    lines = Path(path).read_text().splitlines()
    head = json.loads(lines[0])
    pts = []
    for line in lines[1:]:
        if not line.strip():
            continue
        r = json.loads(line)
        if r.get("record") == "point":
            pts.append({"cai": float(r["cai"]), "mfe": float(r["mfe"])})
    pts.sort(key=lambda p: p["cai"])
    return head, pts


def best_round(root, protein, nest):
    """A cell's rounds accumulate; the one with the most points is the current front."""
    best = None
    for path in Path(root).glob(f"{protein}_{nest}*.jsonl"):
        head, pts = front(path)
        if best is None or len(pts) > len(best[1]):
            best = (head, pts, path.name)
    return best


def price(base, pts):
    """ΔMFE at matched CAI: how much worse the constrained design is at the same adaptation."""
    out = []
    for p in pts:
        # the baseline's best energy among designs at least as adapted
        cands = [b["mfe"] for b in base if b["cai"] >= p["cai"] - 1e-12]
        if cands:
            out.append(p["mfe"] - min(cands))
    return out


def main():
    root = Path(sys.argv[1])
    binding = [
        json.loads(x)
        for x in (root / "screen" / "binding.jsonl").read_text().splitlines()
        if x.strip()
    ]
    occurs = defaultdict(int)
    for r in binding:
        occurs[(r["protein"], r["set"])] += r["designs_carrying"]

    proteins = sorted(
        {p.stem for p in (root / "screen" / "f0").glob("*.jsonl")}
        & {f.name.split("_")[0] for f in (root / "grid").glob("*_F*.jsonl")}
    )

    cells, fronts = [], {}
    for protein in proteins:
        base_head, base = front(root / "screen" / "f0" / f"{protein}.jsonl")
        fronts[protein] = {
            "F0": {
                "points": [(round(p["cai"], 4), round(p["mfe"], 1)) for p in base],
                "n": len(base),
                "wall": base_head["wall_s"],
                "complete": base_head.get("stop") is None,
            }
        }
        for nest in NESTS:
            got = best_round(root / "grid", protein, nest)
            if not got:
                continue
            head, pts, name = got
            removed = sum(occurs[(protein, s)] for s in SETS[nest])
            deltas = price(base, pts)
            fronts[protein][nest] = {
                "points": [(round(p["cai"], 4), round(p["mfe"], 1)) for p in pts],
                "n": len(pts),
                "wall": head["wall_s"],
                "complete": head.get("complete") is True,
                "stopped_by": head.get("stopped_by"),
            }
            cells.append(
                {
                    "protein": protein,
                    "nest": nest,
                    "aa": head.get("protein_aa"),
                    "file": name,
                    "points": len(pts),
                    "base_points": len(base),
                    "wall": head["wall_s"],
                    "base_wall": base_head["wall_s"],
                    "ratio": head["wall_s"] / base_head["wall_s"],
                    "complete": head.get("complete") is True,
                    "stopped_by": head.get("stopped_by"),
                    "binds": removed > 0,
                    "carrying": removed,
                    "dmfe_med": st.median(deltas) if deltas else None,
                    "dmfe_lo": min(deltas) if deltas else None,
                    "dmfe_hi": max(deltas) if deltas else None,
                    "dmfe_n": len(deltas),
                    "cai_lo": min(p["cai"] for p in pts) if pts else None,
                    "cai_hi": max(p["cai"] for p in pts) if pts else None,
                }
            )

    # The instance for panel a is the one whose cells are furthest along: most certified
    # complete, then most blacklists run. Naming the rule keeps the choice out of taste.
    shown = max(
        proteins,
        key=lambda p: (
            sum(1 for c in cells if c["protein"] == p and c["complete"]),
            sum(1 for c in cells if c["protein"] == p),
        ),
    )
    json.dump(
        {
            "cells": cells,
            "fronts": {shown: fronts[shown]},
            "instance": shown,
            "proteins": proteins,
            "meta": {
                "cells": len(cells),
                "proteins": len(proteins),
                "binding": sum(1 for c in cells if c["binds"]),
                "complete": sum(1 for c in cells if c["complete"]),
            },
        },
        open(sys.argv[2], "w"),
        indent=1,
    )

    print(f"{len(proteins)} proteins, {len(cells)} cells, instance for panel a: {shown}")
    for c in cells:
        tag = "binds" if c["binds"] else "VACUOUS"
        d = f"{c['dmfe_med']:.2f}" if c["dmfe_med"] is not None else "  -  "
        print(
            f"  {c['protein']:8s} {c['nest']} {tag:8s} pts {c['points']:5d}/{c['base_points']:5d} "
            f"dMFE med {d}  ratio {c['ratio']:6.2f}  complete {c['complete']}  "
            f"stop {c['stopped_by']}"
        )


if __name__ == "__main__":
    main()
