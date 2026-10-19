#!/usr/bin/env python3
"""The hull probe's fit, and the pre-registered verdict on the long arm.

Reads `e3-probe/summary.csv`, fits log(cost) against log(n) over the rungs that
completed, extrapolates to the long arm's lengths and applies the decision rule
fixed in E3-PLAN.md § 3a before the probe ran:

    if predicted peak RSS at 984 aa exceeds the cap at --workers 1, the Cas arm
    is not exact and runs budgeted with anytime partial fronts; if it is inside
    the cap, the arm stays open pending Phase B's gap-solve cost.

Peak RSS is the primary readout. Walls are load-qualified lower bounds, so their
fit is reported and carries that clause. Stdlib only.
"""

import csv
import math
import os
import sys

CAP_GB = 52.0          # tamazirt's memory cap for this account
LONG_ARM = (984, 1629)  # the shortest and longest Cas ortholog
HERE = os.path.dirname(os.path.abspath(__file__))


def fit(xs, ys):
    """Least squares on log-log: returns (exponent, intercept, r2)."""
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    n = len(lx)
    mx, my = sum(lx) / n, sum(ly) / n
    sxx = sum((x - mx) ** 2 for x in lx)
    sxy = sum((x - mx) * (y - my) for x, y in zip(lx, ly))
    b = sxy / sxx
    a = my - b * mx
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(lx, ly))
    ss_tot = sum((y - my) ** 2 for y in ly)
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")
    return b, a, r2


def predict(a, b, n):
    return math.exp(a + b * math.log(n))


def main():
    path = os.path.join(HERE, "e3-probe", "summary.csv")
    rows = [r for r in csv.DictReader(open(path)) if r["exit"] == "0" and r["rung"]]
    if len(rows) < 2:
        print(f"{len(rows)} completed rung(s): not enough to fit", file=sys.stderr)
        return 1

    ns = [int(r["rung"]) for r in rows]
    rss = [float(r["peak_rss_kb"]) / 1048576 for r in rows]   # GB
    wall = [float(r["wall_s"]) for r in rows]

    print(f"completed rungs: {ns}")
    print(f"peak RSS (GB):   {[round(v, 3) for v in rss]}")
    print(f"wall (s):        {[round(v, 1) for v in wall]}")
    print()

    for name, ys, unit in (("peak RSS", rss, "GB"), ("wall", wall, "s")):
        b, a, r2 = fit(ns, ys)
        print(f"{name}: exponent {b:.2f}, r2 {r2:.4f}")
        for n in LONG_ARM:
            print(f"    n = {n:4d}: {predict(a, b, n):.1f} {unit}")
        if name == "wall":
            print("    (walls are load-qualified lower bounds; see the plan)")
        print()

    b, a, _ = fit(ns, rss)
    at984 = predict(a, b, LONG_ARM[0])
    verdict = "NOT EXACT: budgeted, anytime partial fronts" if at984 > CAP_GB \
        else "OPEN: inside the cap, pending Phase B's gap-solve cost"
    print(f"predicted RSS at 984 aa = {at984:.1f} GB against a {CAP_GB:.0f} GB cap")
    print(f"verdict on the long arm: {verdict}")
    print("note: the hull lower-bounds the front, so a NOT EXACT verdict here is final")
    return 0


if __name__ == "__main__":
    sys.exit(main())
