"""Figure 4's per-protein table: each arm's exit and wall time, and the dichotomic hull search's
cost with and without Carrillo-Lipman pruning.

    python fig4_per_protein.py SWEEP_JSONL HULL_DIR OUT_CSV

The columns Figure 4 reads from E2's analysis table (`e2_per_protein.csv`), computed as that
table computes them, so the output equals those columns of it:

- one cell row per (arm, protein): the later row in the file, whatever its pin. The arm
  `mrnafold-hull-c8` is left out, and so are the four proteins whose band was stopped partway
  (P05160, Q06187, Q8TAA9, Q9UNN5);
- `dichotomic_pruned_min`: the walk's own hull sweep, `sweep_s` in HULL_DIR/<P>.walk.json;
- `dichotomic_unpruned_projected_min`: the unpruned search runs the exact DP 2V-1 times for V
  hull vertices, at derna2-hull's wall per λ (over the λs it finished, where it timed out).
"""

import csv
import json
import sys
from pathlib import Path

ARMS = [
    "derna2-walk",
    "derna2-hull",
    "mrnafold-hull",
    "lineardesign-hull",
    "lcds-exact-hull",
    "derna1-own",
]
EXCLUDE = {"mrnafold-hull-c8"}
STOPPED = {"Q8TAA9", "Q9UNN5", "Q06187", "P05160"}


def cells(sweep):
    """(arm, protein) -> the later cell row, and protein -> length."""
    cell, n_aa = {}, {}
    for line in open(sweep):
        try:
            r = json.loads(line)
        except ValueError:
            continue
        if r.get("kind") == "cell" and r["arm"] not in EXCLUDE:
            cell[(r["arm"], r["protein"])] = r
            n_aa[r["protein"]] = r["n_aa"]
    return cell, n_aa


def dichotomic(cell, hull_dir, p):
    """Minutes with pruning, DP runs without, projected minutes without; None if unmeasured."""
    r = cell.get(("derna2-hull", p))
    if not r:
        return None
    try:
        walk = json.load(open(hull_dir / f"{p}.walk.json"))
    except (OSError, ValueError):
        return None
    done = int(str(r.get("lambdas_done") or 0).split("/")[0])
    if not done:
        return None
    runs = 2 * len(walk["vertices"]) - 1
    return round(walk["sweep_s"] / 60, 2), runs, round(r["wall_s"] / done * runs / 60, 2)


def main():
    sweep, hull_dir, out = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    cell, n_aa = cells(sweep)
    prots = sorted({p for (a, p) in cell if a == "derna2-walk"}, key=lambda p: (n_aa[p], p))
    prots = [p for p in prots if p not in STOPPED]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            ["protein", "n_aa"]
            + [f"{a}_{x}" for a in ARMS for x in ("exit", "wall_min")]
            + [
                "dichotomic_pruned_min",
                "dichotomic_unpruned_dp_runs",
                "dichotomic_unpruned_projected_min",
            ]
        )
        for p in prots:
            row = [p, n_aa[p]]
            for a in ARMS:
                r = cell.get((a, p))
                row += [r["exit"], round(r["wall_s"] / 60, 2)] if r else ["", ""]
            row += list(dichotomic(cell, hull_dir, p) or ("", "", ""))
            w.writerow(row)
    print(len(prots), "proteins")


if __name__ == "__main__":
    main()
