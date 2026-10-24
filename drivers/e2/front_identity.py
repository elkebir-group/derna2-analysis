"""Compare two fronts of the same molecule without letting float noise decide.

Exact set equality on a key containing `cai_log` does not test what it looks like
it tests. `cai_log` is a sum of per-codon logarithms, so two routes to one point
differ in the last bits whenever the summation order differs -- measured here at
1.42e-14 between a resumed walk and an uninterrupted one, with `mfe` identical to
the bit. An exact gate reports that as a different point and a front of 479 as
"DIFFERS", which is a statement about IEEE addition, not about the walk.

The opposite error is a tolerance chosen to make a comparison pass. A slack wide
enough to absorb noise is also wide enough to admit a genuinely different point,
and then the gate cannot fail. So this does not take a tolerance as an input. It
builds a bijection between the two fronts, greedily by nearest point, and
*reports* the largest deviation the bijection had to accept. The caller reads that
number: 1e-14 is two routes to one point, 1e-2 is two different points, and
nothing in between is silently absorbed.

`mfe` is held exact. It is a sum over a lookup table with few enough terms that
the bit pattern is reproducible, which the measurement confirms, so a deviation
there is a finding rather than noise.

Designs are reported, never gated. One point can be spelled by more than one
co-optimal design -- a resumed front and an uninterrupted one chose different
sequences at a tied hull vertex here -- and which one a walk reports is a
tie-breaking detail, in the same way that one sequence can carry more than one
co-optimal structure.
"""

import json
import sys

# Far above the measured summation noise (1.4e-14) and far below the energy
# quantum (1e-2). A deviation in between is unexplained, so it fails.
NOISE_CEILING = 1e-9


def load(path):
    with open(path) as fh:
        records = [json.loads(line) for line in fh]
    header = next(r for r in records if "points" in r)
    return header, [r for r in records if "sequence" in r]


def bijection(left, right):
    """Pair each left point with its nearest unused right point.

    Returns the pairs and whatever could not be paired. Greedy on a front this
    size is exact in practice: the points are far apart relative to the noise
    being tolerated, so the nearest unused partner is the only candidate.
    """
    pool = list(range(len(right)))
    pairs, unpaired = [], []
    for i, a in enumerate(left):
        best, best_d = None, None
        for j in pool:
            b = right[j]
            if a["mfe"] != b["mfe"]:
                continue
            d = abs(a["cai_log"] - b["cai_log"])
            if best_d is None or d < best_d:
                best, best_d = j, d
        if best is None:
            unpaired.append(i)
        else:
            pool.remove(best)
            pairs.append((i, best, best_d))
    return pairs, unpaired, pool


def main():
    if len(sys.argv) != 4:
        print("usage: front_identity.py LABEL_A FRONT_A FRONT_B")
        return 2
    label, pa, pb = sys.argv[1], sys.argv[2], sys.argv[3]
    ha, a = load(pa)
    hb, b = load(pb)

    print("=== %s ===" % label)
    print("%-10s %7s %6s %6s %10s %6s" % ("", "points", "hull", "gaps", "complete", "open"))
    for name, h, pts in (("A", ha, a), ("B", hb, b)):
        print("%-10s %7d %6d %6s %10s %6d"
              % (name, len(pts), sum(1 for r in pts if r.get("on_hull")),
                 h.get("gaps"), h.get("complete"), len(h.get("open", []))))

    pairs, only_a, only_b = bijection(a, b)
    print()
    print("paired %d of %d; unpaired A %d, unpaired B %d"
          % (len(pairs), max(len(a), len(b)), len(only_a), len(only_b)))

    if only_a or only_b:
        print("  *** an unpaired point has no partner at the same exact mfe, which")
        print("  *** float noise cannot explain. These are different fronts.")
        for i in only_a[:5]:
            print("      only A: mfe %.17g cai_log %.17g" % (a[i]["mfe"], a[i]["cai_log"]))
        for j in only_b[:5]:
            print("      only B: mfe %.17g cai_log %.17g" % (b[j]["mfe"], b[j]["cai_log"]))
        return 1

    worst = max((d for _, _, d in pairs), default=0.0)
    exact = sum(1 for _, _, d in pairs if d == 0.0)
    print("cai_log deviation across the bijection: worst %.3e, exact on %d of %d"
          % (worst, exact, len(pairs)))
    if worst == 0.0:
        print("  bit-identical on every paired point")
    elif worst < NOISE_CEILING:
        print("  below NOISE_CEILING = %.0e, so consistent with two summation" % NOISE_CEILING)
        print("  orders over the same points and not with two different fronts")
        print("  (a real difference is bounded below by the 1e-2 energy quantum,")
        print("  and mfe is exact on every pair)")
    else:
        print("  *** at or above NOISE_CEILING = %.0e. Too large for summation" % NOISE_CEILING)
        print("  *** order, so treat it as a real difference.")
        return 1

    differing = [(i, j) for i, j, _ in pairs if a[i]["sequence"] != b[j]["sequence"]]
    print()
    print("designs differing at a paired point: %d of %d (reported, not gated)"
          % (len(differing), len(pairs)))
    for i, j in differing[:5]:
        print("   mfe %.17g cai_log %.17g  on_hull %s/%s"
              % (a[i]["mfe"], a[i]["cai_log"], a[i].get("on_hull"), b[j].get("on_hull")))
    if differing:
        print("  a point reachable by more than one co-optimal design; which one a")
        print("  walk reports is a tie-break, so it is not evidence of divergence")

    order = [(r["mfe"], r["sequence"]) for r in a] == [(r["mfe"], r["sequence"]) for r in b]
    print()
    print("emission order: %s (reported, not gated)"
          % ("identical" if order else "differs"))
    print("VERDICT: the two fronts are the same front")
    return 0


if __name__ == "__main__":
    sys.exit(main())
