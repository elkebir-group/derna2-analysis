"""The long band's selection, from the lattice census rather than from fronts.

The band is ten instances at 1000 aa and above, published as anytime curves: a
complete constrained front is out of reach there, so each cell runs under a
fixed budget and reports the front it reached.

Selection cannot use a baseline front. One costs 12 h at 1000 aa and 65 h at
1629 on the measured `n^3.53`, which would spend longer screening the band than
running it. The census is used instead: `spellable` is a lattice query, cheap at
any length and immune to host load, and `occurs` is read off the partial front
each cell produces anyway.

ELIGIBILITY is that every set of the nest has at least one spellable motif.
The plan pre-registered the stricter `every entry of the nest is spellable`;
that was amended before any front ran, with its reason recorded, because a nest
step's price is non-vacuous when one motif of its set binds and the strict form
therefore dropped instances for a property the claim does not rest on. The
amendment restores two: Q6UXY8, missing one of seventeen manufacturing motifs,
and Q96HA7, missing one of three processing motifs. For set IV the two forms
coincide, that set holding exactly one motif, so the instances that lose F4 lose
it outright rather than marginally.

The two forms are both reported, so a reader can see what the criterion cost.

THE RULE, fixed before any front in this band exists:

  1. Qualify: at least 1000 aa, a natural sequence, and at least one spellable
     motif in each of sets I-IV.
  2. The Cas orthologs first, by identifier. They are the band's reason for
     reaching this length at all, and a stratified draw could otherwise leave
     the band without one -- which is what the census found before the three
     missing fastas were sourced.
  3. Then stratify the remaining slots in LOG length. Cost grows as about
     n^3.53 and the band spans 2.81x, so equal linear strata would place most
     boundaries inside the cheap half; equal ratios spread them over the decade
     the exponent is claimed across.
  4. Final tie: alphabetical by identifier. Arbitrary, fixed, and independent of
     every measured quantity, so the data cannot steer it.

A reversed synthetic sequence in the mirror is excluded by rule 1 and named in
the report, because a band characterizing real proteins cannot reach its tenth
instance by admitting a sequence that is not one.
"""

import argparse
import collections
import json
import math

NEST = ("manufacturing", "processing", "repeats", "initiation")
SYNTHETIC = ("ReverseCovid",)
CAS_PREFIX = "PRT-"
LOW, WANTED = 1000, 10


def census(path):
    """protein -> {set -> {motif -> spellable}}, and protein -> length."""
    spellable = collections.defaultdict(dict)
    aa = {}
    with open(path) as fh:
        for line in fh:
            rec = json.loads(line)
            if "spellable" not in rec:
                continue
            per_set = spellable[rec["protein"]].setdefault(rec["set"], {})
            per_set[rec["motif"]] = rec["spellable"]
            aa[rec["protein"]] = rec["aa"]
    return spellable, aa


def deepest(per_set, strict):
    """The deepest nest this lattice admits, 0 if not even the first set."""
    depth = 0
    for step, name in enumerate(NEST, 1):
        motifs = per_set.get(name, {})
        live = sum(1 for v in motifs.values() if v > 0)
        admits = live == len(motifs) if strict else live > 0
        if not (admits and depth == step - 1):
            break
        depth = step
    return depth


def qualifiers(spellable, aa, strict=False):
    """Everything eligible for the full nest, as [(name, aa)]."""
    out = []
    for name in sorted(spellable):
        if name in SYNTHETIC or aa[name] < LOW:
            continue
        if deepest(spellable[name], strict) == len(NEST):
            out.append((name, aa[name]))
    return sorted(out, key=lambda r: r[1])


def log_strata(lo, hi, count):
    """`count` equal-ratio bands over [lo, hi], as (index, low, high)."""
    step = (math.log(hi) - math.log(lo)) / count
    edges = [math.exp(math.log(lo) + step * k) for k in range(count + 1)]
    return [(k, edges[k], edges[k + 1]) for k in range(count)]


def take_cas(pool, chosen, used):
    """Rule 2: the Cas orthologs, by identifier."""
    for name, length in sorted(pool, key=lambda r: r[0]):
        if name.startswith(CAS_PREFIX) and len(chosen) < WANTED:
            chosen.append((name, length, "Cas ortholog"))
            used.add(name)


def nearest(cands, mid):
    """Rule 4: closest to a stratum's geometric midpoint, then alphabetical."""
    return min(cands, key=lambda r: (abs(math.log(r[1]) - math.log(mid)), r[0]))


def fill_strata(pool, chosen, used, strata, inside_only):
    """Rule 3, one pass. `inside_only` fills strata that hold a candidate."""
    for _, lo, hi in strata:
        if len(chosen) >= WANTED:
            return
        mid = math.sqrt(lo * hi)
        if inside_only:
            cands = [r for r in pool if r[0] not in used and lo <= r[1] < hi]
            label = "stratum %d-%d aa" % (round(lo), round(hi) - 1)
        else:
            if any(lo <= n < hi for _, n, _ in chosen):
                continue
            cands = [r for r in pool if r[0] not in used]
            label = "stratum %d-%d aa, empty, nearest taken" % (round(lo), round(hi) - 1)
        if not cands:
            continue
        pick = nearest(cands, mid)
        chosen.append((pick[0], pick[1], label))
        used.add(pick[0])


def select(pool):
    """Rules 2-4 over [(name, aa)]."""
    chosen, used = [], set()
    take_cas(pool, chosen, used)
    remaining = WANTED - len(chosen)
    if remaining <= 0:
        return chosen
    rest = [r for r in pool if r[0] not in used]
    if not rest:
        return chosen
    strata = log_strata(min(n for _, n in rest), max(n for _, n in rest) + 1, remaining)
    fill_strata(pool, chosen, used, strata, inside_only=True)
    fill_strata(pool, chosen, used, strata, inside_only=False)
    return chosen


def report_tiers(spellable, aa):
    """What each criterion admits, so the amendment's cost stays visible."""
    for strict, label in ((True, "every entry of the nest"),
                          (False, "at least one entry per set")):
        tiers = collections.defaultdict(list)
        for name in spellable:
            if name in SYNTHETIC:
                continue
            tiers[deepest(spellable[name], strict)].append(name)
        counts = " ".join("F%d:%d" % (d, len(tiers[d]))
                          for d in sorted(tiers, reverse=True))
        print("  %-28s %s" % (label, counts))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--census",
                    default="run/paper-e3/screen/lattice-xl.jsonl")
    args = ap.parse_args()

    spellable, aa = census(args.census)
    print("census: %d instance(s), %d-%d aa"
          % (len(spellable), min(aa.values()), max(aa.values())))
    print("excluded as not a protein: %s" % ", ".join(SYNTHETIC))
    print("deepest nest admitted, by criterion:")
    report_tiers(spellable, aa)

    pool = qualifiers(spellable, aa)
    print()
    print("qualifiers for the full nest: %d, wanted %d" % (len(pool), WANTED))
    if len(pool) < WANTED:
        print("FEWER THAN WANTED. The band is those %d, and the supply is"
              % len(pool))
        print("re-examined rather than the criterion relaxed again.")

    chosen = select(pool)
    print()
    print("=== the long band, %d instance(s) ===" % len(chosen))
    print("%-12s %5s  %s" % ("instance", "aa", "selected as"))
    for name, length, why in chosen:
        print("%-12s %5d  %s" % (name, length, why))
    lens = [n for _, n, _ in chosen]
    print()
    print("length span: %d-%d aa (%.2fx)" % (min(lens), max(lens), max(lens) / min(lens)))
    print("cells: %d x 4 = %d" % (len(chosen), 4 * len(chosen)))
    return 0 if len(chosen) == WANTED else 1


if __name__ == "__main__":
    raise SystemExit(main())
