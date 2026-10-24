#!/usr/bin/env python3
"""The per-λ price of motif-freedom, read off two fronts rather than re-solved.

The price at a trade-off weight λ is the difference between the best design the
constraint allows at λ and the best design at λ with no constraint. Both are
already in the front records: a weighted-sum optimum can only ever be a vertex
of the convex hull of the feasible set, so the optimum at λ is the point of the
front minimizing

    f(λ) = λ·mfe − (1−λ)·cai_log

and the price is that argmin over the constrained front against the argmin over
the baseline front. No design call is needed, and none is made here.

The λs reported are the baseline front's own hull λs. A point record carries no
λ — it carries `mfe`, `cai_log` and `on_hull` — so the λs are recovered from the
hull's geometry: consecutive hull vertices exchange optimality at a breakpoint
λ, and each vertex is reported at the midpoint of the interval it wins.

**A partial front is not admissible.** Where a walk stopped on its budget the
front may be missing the vertex that is optimal at some λ, and the argmin then
returns a point that is suboptimal with nothing to say so. A front whose header
reports `complete: false` is refused rather than priced.

The hull recovered here is checked against the `on_hull` flags the walk itself
wrote, which is what caught an ordering fault during development: slope descends
along ascending λ, and building the envelope the other way left 2 vertices of 61
while still pricing a front against itself at exactly zero. Over all 61 baseline
fronts of the screen the recovered hull, the flagged hull and the `hull` column of
`f0/summary.csv` now agree exactly, and every self-price is 0.

Usage:
    python3 price_per_lambda.py BASELINE.jsonl CONSTRAINED.jsonl [-o OUT.jsonl]
    python3 price_per_lambda.py --self-test BASELINE.jsonl

Stdlib only. Reads front records and nothing from the cloned corpus, so it
needs no `-I`.
"""

import argparse
import json
import sys

# A hull vertex's own λ interval can be narrower than this where two vertices
# nearly coincide; such an interval is reported at its midpoint like any other.
TIES = 1e-12


class PartialFront(Exception):
    """A front whose walk did not finish, and so cannot be priced."""


def read_front(path):
    """A front's points as (mfe, cai_log) pairs, with its hull flags.

    Raises PartialFront if the header says the walk stopped short. A header
    without a `complete` field predates the budget flags and is taken as
    complete: plain `front` has no budget to stop on, and its completion is
    attested by the exit code in the screen's summary.
    """
    points, flagged = [], []
    for line in open(path):
        rec = json.loads(line)
        if rec.get("record") == "header" and rec.get("complete") is False:
            raise PartialFront(path)
        if "mfe" in rec and "cai_log" in rec:
            points.append((rec["mfe"], rec["cai_log"]))
            flagged.append(bool(rec.get("on_hull")))
    if not points:
        raise ValueError(f"no points in {path}")
    return points, flagged


def objective(point, lam):
    """The weighted objective at λ, which the front minimizes."""
    mfe, cai_log = point
    return lam * mfe - (1 - lam) * cai_log


def best_at(points, lam):
    """The front's optimum at λ: the point minimizing the weighted objective."""
    return min(points, key=lambda p: objective(p, lam))


def as_line(point):
    """A point as the line its weighted objective traces in λ.

    f(λ) = λ·mfe − (1−λ)·cai_log is affine in λ, so each design is a line
    `slope·λ + intercept`. The optimum at λ is then the lowest line at λ, and
    the λ-attainable designs are the vertices of the lower envelope — which is
    what makes the hull recoverable without any geometry in (mfe, cai).
    """
    mfe, cai_log = point
    return mfe + cai_log, -cai_log


def crosses_at(left, right):
    """The λ where two lines meet, or None when they are parallel."""
    (slope_l, int_l), (slope_r, int_r) = left, right
    if abs(slope_l - slope_r) <= TIES:
        return None
    return (int_r - int_l) / (slope_l - slope_r)


def envelope(points):
    """The λ-attainable points and the λ interval each one wins, over (0, 1).

    The lower envelope of the lines, built by the usual monotone scan. Slope
    descends along ascending λ — a steeper line is lowest while λ is small and
    is overtaken as λ grows — so the scan runs over descending slope and the
    stack comes out in the order the λs do. A line is dropped when its two
    neighbours already meet below it, and one of two lines of equal slope is
    dropped outright.
    """
    lines = sorted(((as_line(p), p) for p in points), key=lambda it: -it[0][0])
    stack = []
    for line, point in lines:
        if stack and abs(stack[-1][0][0] - line[0]) <= TIES:
            if stack[-1][0][1] <= line[1]:
                continue
            stack.pop()
        while len(stack) >= 2 and _hidden(stack[-2][0], stack[-1][0], line):
            stack.pop()
        stack.append((line, point))

    bounds = [0.0]
    for (left, _), (right, _) in zip(stack, stack[1:]):
        at = crosses_at(left, right)
        bounds.append(at if at is not None else bounds[-1])
    bounds.append(1.0)

    out = []
    for point, lo, hi in zip((p for _, p in stack), bounds, bounds[1:]):
        lo, hi = max(lo, 0.0), min(hi, 1.0)
        if hi - lo > TIES:
            out.append((point, lo, hi))
    return out


def _hidden(left, middle, right):
    """Whether `middle` is above where its neighbours already cross."""
    at = crosses_at(left, right)
    if at is None:
        return True
    return middle[0] * at + middle[1] >= left[0] * at + left[1] - TIES


def hull(points):
    """The λ-attainable designs, in ascending λ."""
    return [p for p, _, _ in envelope(points)]


def hull_lambdas(points):
    """One λ per attainable design: the midpoint of the interval it wins."""
    return [(lo + hi) / 2 for _, lo, hi in envelope(points)]


def price(baseline, constrained, lambdas):
    """Δobjective, ΔMFE and ΔCAI at each λ, constrained against baseline."""
    rows = []
    for lam in lambdas:
        free = best_at(baseline, lam)
        held = best_at(constrained, lam)
        rows.append({
            "lambda": lam,
            "mfe_baseline": free[0],
            "mfe_constrained": held[0],
            "cai_log_baseline": free[1],
            "cai_log_constrained": held[1],
            "d_objective": objective(held, lam) - objective(free, lam),
            "d_mfe": held[0] - free[0],
            "d_cai_log": held[1] - free[1],
        })
    return rows


def check_hull_agrees(points, flagged):
    """How often the recovered hull matches the front's own `on_hull` flags.

    The records already say which points the walk put on the hull, so this is a
    cross-check of the geometry here against the walk's, not a definition.
    """
    recovered = set(hull(points))
    marked = {p for p, on in zip(points, flagged) if on}
    return len(recovered & marked), len(recovered), len(marked)


def self_test(path):
    """Pricing a front against itself must give exactly zero at every λ.

    This is the test that the λ arithmetic, the hull recovery and the argmin
    agree: if any of the three disagreed, a front would appear to have a price
    against itself.
    """
    points, flagged = read_front(path)
    lambdas = hull_lambdas(points)
    rows = price(points, points, lambdas)
    worst = max(abs(r["d_objective"]) for r in rows)
    agreed, recovered, marked = check_hull_agrees(points, flagged)
    on_hull = sum(1 for lam in lambdas
                  if any(best_at(points, lam) == p
                         for p, on in zip(points, flagged) if on))

    print(f"{path}: {len(points)} point(s), {marked} flagged on_hull")
    print(f"  hull recovered here: {recovered} vertex(es), "
          f"{agreed} also flagged by the walk")
    print(f"  hull lambdas: {len(lambdas)}")
    print(f"  argmin lands on a flagged hull vertex: {on_hull}/{len(lambdas)}")
    print(f"  worst |d_objective| against itself: {worst:.3e}")
    ok = worst == 0.0 and on_hull == len(lambdas)
    print("  PASS" if ok else "  FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(
        description="the per-lambda price of motif-freedom, read off two fronts")
    ap.add_argument("baseline")
    ap.add_argument("constrained", nargs="?")
    ap.add_argument("-o", "--out")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test(args.baseline)
    if not args.constrained:
        ap.error("a constrained front is required unless --self-test")

    try:
        baseline, flagged = read_front(args.baseline)
        constrained, _ = read_front(args.constrained)
    except PartialFront as exc:
        print(f"refused: {exc} reports complete: false, so its argmin at some "
              f"lambda may be suboptimal with nothing to say so", file=sys.stderr)
        return 2

    lambdas = hull_lambdas(baseline)
    rows = price(baseline, constrained, lambdas)

    if args.out:
        with open(args.out, "w") as fh:
            for row in rows:
                fh.write(json.dumps(row) + "\n")
        print(f"wrote {args.out}")

    agreed, recovered, marked = check_hull_agrees(baseline, flagged)
    print(f"baseline: {len(baseline)} point(s), hull {recovered} recovered / "
          f"{marked} flagged / {agreed} agreeing")
    print(f"priced at {len(lambdas)} hull lambda(s)")
    worst = max(rows, key=lambda r: r["d_objective"])
    print(f"worst d_objective {worst['d_objective']:+.4f} at "
          f"lambda={worst['lambda']:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
