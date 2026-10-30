"""Figure 4's run table: one cell row per (arm, protein) of one pin, chosen rather than taken last.

    python fig4_cells.py SWEEP_JSONL PIN OUT_CSV

A pin can carry two attempts at the same cell (a re-run, or two binaries of the same build
series), and keeping whichever line came last would make the choice an accident of file order.
Keep instead the best the method achieved on this pin: a run that delivered over one that did
not, and among equals the faster. `attempts` records how many rows the cell had.
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


def rank(r):
    """Lower is better: a delivered run first, then the shorter wall time."""
    return (r.get("exit") != "ok", r.get("wall_s") or 0.0)


def main():
    sweep, pin, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
    out.parent.mkdir(parents=True, exist_ok=True)
    seen = defaultdict(list)
    for line in open(sweep):
        r = json.loads(line)
        if r.get("pin") == pin and r.get("kind") == "cell":
            seen[(r["arm"], r["protein"])].append(r)
    rows = []
    for (arm, protein), rs in seen.items():
        best = min(rs, key=rank)
        rows.append(
            {
                "protein": protein,
                "arm": arm,
                "n_aa": best.get("n_aa"),
                "exit": best.get("exit"),
                "wall_s": best.get("wall_s"),
                "attempts": len(rs),
            }
        )
    t = pd.DataFrame(rows).sort_values(["arm", "protein"])
    t.to_csv(out, index=False)
    dup = t[t.attempts > 1]
    print(len(t), "cells,", len(dup), "with more than one attempt")
    print(dup.groupby("arm").size().to_dict())


if __name__ == "__main__":
    main()
