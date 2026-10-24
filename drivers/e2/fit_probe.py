#!/usr/bin/env python3
"""The hull probe's table, fits and frontier — every number the plan quotes.

Reads the rung records written by `e2_hull_probe*.sh` and prints the rung table,
the log-log fits for wall and peak RSS, the projections to the long arm's
lengths, and the frontier: the largest n whose hull completes inside a budget.

Peak RSS fits a clean power law; the wall does not, so its fit is printed with
its r² and with the consecutive-rung exponents, which are what the discontinuity
shows up in. Walls are E-core walls under load — load at launch is printed per
rung, and the n = 200 control pins what load costs (+3%).

Usage: `python3 -I fit_probe.py`. Stdlib only.
"""

import csv
import glob
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
BUDGETS = ((3600, "1 h"), (43200, "12 h"))
LONG_ARM = (984, 1629)
CAP_GB = 52.0


def rungs():
    """Every completed rung, as (motif set, n, wall s, rss kB, load, jsonl path)."""
    out = []
    for summary in sorted(glob.glob(os.path.join(HERE, "e3-probe*", "summary.csv"))):
        root = os.path.dirname(summary)
        for row in csv.DictReader(open(summary)):
            if row.get("exit") != "0":
                continue
            label = row.get("label") or f"n{row['rung']}"
            n = int(row.get("n") or row["rung"])
            path = os.path.join(root, f"{label}.jsonl")
            out.append({
                "set": row.get("motifs", "25-pattern"),
                "label": label,
                "n": n,
                "wall": float(row["wall_s"]),
                "rss_kb": float(row["peak_rss_kb"]),
                "load": float(row["load_at_launch"]),
                "path": path,
            })
    return out


def probe_stats(path):
    """Per-fold cost, split at λ = 0.01, plus probe and fold counts."""
    if not os.path.exists(path):
        return {}
    ps = [json.loads(l) for l in open(path) if '"lambda"' in l]
    if not ps:
        return {}
    lo = [p for p in ps if p["lambda"] < 0.01]
    hi = [p for p in ps if p["lambda"] >= 0.01]
    fold_s = lambda g: sum(p["s"] for p in g) / max(sum(p["folds"] for p in g), 1)
    return {
        "probes": len(ps),
        "folds": sum(p["folds"] for p in ps),
        "s_fold_ord": fold_s(hi),
        "s_fold_low": fold_s(lo),
        "pruned": sum(p["pruned"] for p in ps),
        "vertices": sum(1 for l in open(path) if "walk-vertex" in l),
    }


def fit(xs, ys):
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    k = len(lx)
    mx, my = sum(lx) / k, sum(ly) / k
    b = sum((x - mx) * (y - my) for x, y in zip(lx, ly)) / sum((x - mx) ** 2 for x in lx)
    a = my - b * mx
    res = sum((y - (a + b * x)) ** 2 for x, y in zip(lx, ly))
    tot = sum((y - my) ** 2 for y in ly)
    return b, a, (1 - res / tot if tot else float("nan"))


def report(label, rows):
    print(f"\n=== {label}  ({len(rows)} completed rungs)")
    head = ("n", "wall s", "RSS GB", "vtx", "probes", "f/probe", "s/fold ord", "s/fold 1e-4", "load")
    print("{:>5} {:>9} {:>7} {:>4} {:>6} {:>7} {:>10} {:>11} {:>5}".format(*head))
    for r in rows:
        st = probe_stats(r["path"])
        fp = st.get("folds", 0) / st["probes"] if st.get("probes") else float("nan")
        print(f"{r['n']:5d} {r['wall']:9.1f} {r['rss_kb'] / 1048576:7.3f} "
              f"{st.get('vertices', 0):4d} {st.get('probes', 0):6d} {fp:7.2f} "
              f"{st.get('s_fold_ord', float('nan')):10.2f} "
              f"{st.get('s_fold_low', float('nan')):11.2f} {r['load']:5.1f}")

    if len(rows) < 3:
        print("  (fewer than three rungs: no fit — a two-point line is exact by construction)")
        return
    ns = [r["n"] for r in rows]
    for name, ys, unit in (("wall", [r["wall"] for r in rows], "s"),
                           ("RSS", [r["rss_kb"] / 1048576 for r in rows], "GB")):
        b, a, r2 = fit(ns, ys)
        proj = "   ".join(f"n={n}: {math.exp(a + b * math.log(n)):.1f} {unit}" for n in LONG_ARM)
        print(f"  {name:4s} exponent {b:.2f}  r2 {r2:.4f}   {proj}")
        if name == "RSS":
            at = math.exp(a + b * math.log(LONG_ARM[0]))
            print(f"       against a {CAP_GB:.0f} GB cap: {'inside' if at < CAP_GB else 'over'}")
    b, a, _ = fit(ns, [r["wall"] for r in rows])
    for budget, lbl in BUDGETS:
        print(f"  frontier: hull completes within {lbl} up to n = "
              f"{math.exp((math.log(budget) - a) / b):.0f}")
    print("  consecutive-rung wall exponents: " + ", ".join(
        f"{math.log(rows[i + 1]['wall'] / rows[i]['wall']) / math.log(rows[i + 1]['n'] / rows[i]['n']):.2f}"
        for i in range(len(rows) - 1)))


def main():
    rows = rungs()
    by_set = {}
    for r in rows:
        by_set.setdefault(r["set"], []).append(r)
    for label, group in by_set.items():
        # One rung per n. Where a length was measured twice, the control re-run
        # wins: it ran under load comparable to the rungs above it, so the fit
        # is not mixing a quiet measurement with loaded ones.
        best = {}
        for r in sorted(group, key=lambda r: r["label"].startswith("control")):
            best[r["n"]] = r
        report(label, sorted(best.values(), key=lambda r: r["n"]))

    if len(by_set) > 1:
        print("\n=== blacklist size at equal n (the density-versus-length discriminator)")
        sets = list(by_set)
        common = set.intersection(*({r["n"] for r in by_set[s]} for s in sets))
        for n in sorted(common):
            cells = []
            for s in sets:
                r = next(x for x in by_set[s] if x["n"] == n)
                st = probe_stats(r["path"])
                cells.append(f"{s}: {r['wall']:.1f} s, {st.get('s_fold_ord', float('nan')):.2f} s/fold")
            print(f"  n={n}: " + "   |   ".join(cells))


if __name__ == "__main__":
    main()
