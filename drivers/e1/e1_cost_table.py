#!/usr/bin/env python3
"""Per-lambda cost per arm, fitted on E1's timing rows and extrapolated to E1's
instances.

E1 measured one lambda per cell, which is exactly the unit E1 multiplies by a
vertex count. For each arm this fits log(wall) and log(peak RSS) against
log(n) over the single-core `ok` rows, reports the exponent, and extrapolates
to the lengths E1 adds.

Memory is extrapolated with more confidence than time: a peak-RSS fit held to
0.4% at 1.8x the fitted length on this substrate, while instruction counts
broke at 2.1x. Both are printed with the extrapolation factor so the reader
can see how far past the data each number sits.

Given E1's own `sweep.jsonl` as a second argument it also reports what E1 has
measured for itself, which is the only source for two arms the core grid never
ran: the hull-gap walk, whose unit of work is a whole front rather than a
lambda, and DERNA v1's own adaptive sweep, whose lambda count is chosen per
instance and so cannot be predicted from a cost model at all. That count is
deterministic and survives a loaded host; the walls beside it do not, so the
load average each cell ran under is printed with them.

Usage: cost_table.py [TIME_JSONL [E1_SWEEP_JSONL]]
"""
import json
import math
import sys
from collections import defaultdict

# E1's instance lengths beyond the ladder: the Cas orthologs.
CAS = [984, 1053, 1082, 1368, 1629]
# The 12 h per-cell limit, in seconds.
LIMIT_S = 12 * 3600


def fit_loglog(xs, ys):
    """Least-squares slope and intercept of log(y) = a*log(x) + b."""
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    n = len(lx)
    mx = sum(lx) / n
    my = sum(ly) / n
    sxx = sum((v - mx) ** 2 for v in lx)
    sxy = sum((lx[i] - mx) * (ly[i] - my) for i in range(n))
    a = sxy / sxx
    b = my - a * mx
    resid = max(abs(ly[i] - (a * lx[i] + b)) for i in range(n))
    return a, b, resid


def main(path):
    # arm variant -> list of (n, wall, rss_kb), single core, ok only,
    # one row per (n, corner) at its best attempt.
    by_arm = defaultdict(dict)
    for line in open(path):
        r = json.loads(line)
        if r["exit"] != "ok" or r["ncores"] != 1:
            continue
        key = f"{r['arm']}" + (f" {r['variant']}" if r["variant"] != "-" else "")
        # Keep the slowest corner per length: the cost table is for admission,
        # so the pessimistic corner is the one that decides.
        prev = by_arm[key].get(r["n_aa"])
        if prev is None or r["wall_s"] > prev[0]:
            by_arm[key][r["n_aa"]] = (r["wall_s"], r["maxrss_kb"])

    print(f"{'arm':22} {'pts':>4} {'n range':>12} "
          f"{'wall exp':>9} {'rss exp':>8}   per-lambda cost at n")
    print(f"{'':22} {'':>4} {'':>12} {'':>9} {'':>8}   "
          f"{'1273':>10} {'1368':>10} {'1629':>10}")
    print("-" * 110)

    rows = {}
    for arm in sorted(by_arm):
        pts = sorted(by_arm[arm].items())
        if len(pts) < 5:
            continue
        ns = [n for n, _ in pts]
        walls = [v[0] for _, v in pts]
        rsses = [v[1] for _, v in pts]
        wa, wb, wres = fit_loglog(ns, walls)
        ra, rb, _ = fit_loglog(ns, rsses)

        def wall_at(n):
            return math.exp(wb) * n ** wa

        def rss_gb_at(n):
            return math.exp(rb) * n ** ra / 1e6

        rows[arm] = (wa, wb, ra, rb, max(ns))
        cells = []
        for n in (1273, 1368, 1629):
            w = wall_at(n)
            cells.append(f"{w / 60:8.0f}m" if w < 3600 else f"{w / 3600:8.1f}h")
        print(f"{arm:22} {len(pts):4d} {min(ns):5d}-{max(ns):<6d} "
              f"{wa:9.2f} {ra:8.2f}   " + " ".join(f"{c:>10}" for c in cells)
              + f"   (fit resid {wres:.2f} in log)")

    print()
    print("Lambdas affordable inside the 12 h cell limit, by arm and length")
    print("(limit / per-lambda wall; the vertex count E1 needs is typically "
          "tens to low hundreds)")
    print()
    hdr = [1012, 1273] + CAS
    print(f"{'arm':22} " + " ".join(f"{n:>7}" for n in hdr))
    print("-" * 90)
    for arm, (wa, wb, ra, rb, nmax) in sorted(rows.items()):
        cells = []
        for n in hdr:
            w = math.exp(wb) * n ** wa
            k = LIMIT_S / w
            cells.append(f"{k:7.0f}" if k >= 1 else f"{k:7.2f}")
        print(f"{arm:22} " + " ".join(cells)
              + f"   [data to n={nmax}, so n=1629 is {1629 / nmax:.2f}x out]")

    print()
    print("Peak RSS extrapolated, GB per concurrent process")
    print(f"{'arm':22} " + " ".join(f"{n:>7}" for n in hdr))
    print("-" * 90)
    for arm, (wa, wb, ra, rb, nmax) in sorted(rows.items()):
        cells = [f"{math.exp(rb) * n ** ra / 1e6:7.1f}" for n in hdr]
        print(f"{arm:22} " + " ".join(cells))


def e1_measured(path):
    """What E1's own cells cost, per arm and length.

    Three arms have no per-lambda cost to report and are printed on their own
    terms: the walk returns a whole front in one run, and DERNA v1's adaptive
    sweep picks its own lambda count, which is the column the cost model cannot
    supply. For the replay arms the per-lambda mean over the cell's own lambda
    rows is the comparable number.
    """
    cells, per_lambda = [], defaultdict(list)
    for line in open(path):
        r = json.loads(line)
        if r.get("kind") == "cell":
            if r["exit"] == "ok":
                cells.append(r)
        elif r.get("wall_s") is not None and r.get("exit") == "ok":
            per_lambda[(r["arm"], r["protein"])].append(r["wall_s"])
    if not cells:
        return
    print()
    print("E1's own measured cells")
    print(f"{'arm':22} {'instance':14} {'n':>5} {'count':>6} {'wall':>9} "
          f"{'per unit':>9} {'peak':>8} {'load':>6}")
    print("-" * 92)
    for c in sorted(cells, key=lambda c: (c["n_aa"], c["arm"])):
        lams = per_lambda.get((c["arm"], c["protein"]), [])
        unit = (sum(lams) / len(lams) if lams
                else (c["wall_s"] / c["lambdas_done"]
                      if c["lambdas_done"] else None))
        rss = c.get("maxrss_kb")
        print(f"{c['arm']:22} {c['protein']:14} {c['n_aa']:5d} "
              f"{c['lambdas_done']:6d} {c['wall_s'] / 60:8.2f}m "
              f"{(f'{unit:8.1f}s' if unit else '        -'):>9} "
              f"{(f'{rss / 1e6:7.2f}G' if rss else '       -'):>8} "
              f"{c.get('load', float('nan')):6.1f}")
    own = sorted(((c["n_aa"], c["lambdas_done"], c["wall_s"])
                  for c in cells if c["arm"] == "derna1-own"))
    if len(own) > 1:
        print()
        print("DERNA v1's own lambda count against length "
              "(deterministic; the walls beside it are not)")
        for n, k, w in own:
            print(f"  n={n:5d}  {k:3d} lambdas  {w / 60:7.2f} min  "
                  f"{w / k:7.1f} s/lambda")
        a, b, res = fit_loglog([n for n, _, _ in own],
                               [w for _, _, w in own])
        print(f"  wall exponent {a:.2f} over n={own[0][0]}-{own[-1][0]}, "
              f"fit resid {res:.2f} in log")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1
         else "run/paper-e1/time.jsonl")
    e1 = (sys.argv[2] if len(sys.argv) > 2
          else "run/paper-e2/sweep.jsonl")
    try:
        e1_measured(e1)
    except OSError:
        pass
