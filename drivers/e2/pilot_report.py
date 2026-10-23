#!/usr/bin/env python3
"""Read a concurrency pilot's levels and say what each cost.

Two kinds of cell, and a wall ratio answers only one of them:

* A cell that **finished** is the same work at every level, so its wall is
  comparable directly: the ratio against the serial baseline is the slowdown.
* A cell the limit **stopped** has the same wall at every level by
  construction — the limit — so its wall says nothing. What varies is how far
  it got, so its slowdown is the work it did at the baseline over the work it
  did here. A capped arm is otherwise invisible to the pilot while still
  contending for the host, which is the half a wall-only comparison misses.

    python3 pilot_report.py --out pilot/ --levels 1,4,8 --slack 0.05

The verdict is the highest level whose median slowdown, over both kinds, is
within `--slack` of 1.
"""

import argparse
import json
import statistics
from pathlib import Path


def cells_of(path):
    """Each cell row of a sweep file, by (arm, protein)."""
    out = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("kind") == "cell":
            out[(row["arm"], row["protein"])] = row
    return out


def slowdown(ref, row):
    """How much slower `row` ran than `ref`, and on which measure.

    `None` where the pair cannot be compared: a cell that finished at one
    level and was stopped at the other is not the same work either way, and
    reading it as a slowdown would turn a changed workload into a cost.
    """
    done, ref_done = row.get("exit") == "ok", ref.get("exit") == "ok"
    if done and ref_done:
        return (row["wall_s"] / ref["wall_s"], "wall") if ref.get("wall_s") else None
    stopped = row.get("exit") == "timeout" and ref.get("exit") == "timeout"
    if stopped and row.get("lambdas_done") and ref.get("lambdas_done"):
        # Same wall by construction, so the work done is what moved.
        return ref["lambdas_done"] / row["lambdas_done"], "work"
    return None


def compare(baseline, got, level, slack, log):
    """One level against the serial baseline, over both kinds of cell."""
    rows = []
    for key, row in sorted(got.items()):
        ref = baseline.get(key)
        if ref is None:
            continue
        measured = slowdown(ref, row)
        if measured:
            rows.append((measured[0], measured[1], key))
    if not rows:
        log(f"c{level}: no comparable cell")
        return False
    median = statistics.median([r for r, _, _ in rows])
    worst, how, where = max(rows)
    by_work = sum(1 for _, how, _ in rows if how == "work")
    log(f"c{level}: {len(rows)} cells ({by_work} by work done), "
        f"median {median:.2f}x, worst {worst:.2f}x "
        f"({where[0]} {where[1]}, by {how})")
    for ratio, how, key in sorted(rows, reverse=True):
        log(f"    {key[0]:20} {key[1]:16} {ratio:5.2f}x  by {how}")
    return median <= 1 + slack


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="the pilot's output root")
    ap.add_argument("--levels", default="1,4,8")
    ap.add_argument("--slack", type=float, default=0.05)
    args = ap.parse_args()
    base = Path(args.out)

    def log(msg):
        print(msg, flush=True)

    levels = [int(x) for x in args.levels.split(",")]
    baseline, passed = None, []
    for level in levels:
        path = base / f"c{level}" / "sweep.jsonl"
        if not path.exists():
            log(f"c{level}: no rows yet ({path})")
            continue
        got = cells_of(path)
        if baseline is None:
            baseline = got
            passed.append(level)
            log(f"c{level}: baseline of {len(got)} cells")
            continue
        if compare(baseline, got, level, args.slack, log):
            passed.append(level)
    if passed:
        log(f"within {args.slack:.0%} of the baseline: {passed}")
        log(f"concurrency to run the grid at: {max(passed)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
