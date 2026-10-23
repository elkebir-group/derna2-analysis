"""The middle band's selection, applying a rule fixed before the data arrived.

The band is ten of the 32 proteins at 500-1000 aa. Its pool is 24: the 19 screened
at 501-700 aa, and 5 inherited E2 fronts at 719-962 aa. An earlier selection took
ten from the nineteen alone and was withdrawn -- it spanned 1.29x in length,
against within-host fits over a 1.08x range that returned R-squared 0.01, so it
could not have supported a claim in length.

The rule below was written down before the five were screened, because choosing it
afterwards is what pre-registration exists to prevent. Both of the original rule's
tie-breaks had already failed and neither is reused: length-tertile stratification,
because no protein at or below 211 aa binds all four sets, and a GC tie-break, for
want of any surplus to break ties over.

  1. Qualify: 500-1000 aa, a complete F0, and sets I, II and III all occurring.
  2. All-four binders first. A qualifier on which set IV also occurs is taken
     before anything else, because those are the only instances with an F4 price
     and a stratified draw could otherwise empty the band's four-nest tier --
     which is what the headline table is made of.
  3. Then stratify: ten equal 50 aa strata over 500-1000, each contributing the
     unused qualifier nearest its midpoint; an empty stratum passes to the nearest
     unused qualifier by absolute distance to that midpoint.
  4. Final tie: alphabetical by accession. Arbitrary, fixed, and independent of
     every measured quantity, so the data cannot steer it.

Step 1 keeps `occurs` on F0 as the criterion although it is known to understate --
a set occurred 7 times on one F0 and 121 times on the F1 beside it -- because the
criterion was pre-registered and changing it after seeing the data is the error
this module is shaped to avoid. The understatement is reported, not corrected.

If fewer than ten qualify, the band is however many do, and the condition for
screening the 13 unscreened proteins at 701-1000 aa is re-examined rather than
assumed.
"""

import argparse
import glob
import json
import os

SETS_I3 = ("manufacturing", "processing", "repeats")
SET_IV = "initiation"
LOW, HIGH, STRATA, WANTED = 500, 1000, 10, 10


def binding(path):
    """(protein, set) -> whether it occurs, and each protein's length."""
    occurs, aa = {}, {}
    with open(path) as fh:
        for line in fh:
            rec = json.loads(line)
            key = (rec["protein"], rec["set"])
            occurs[key] = occurs.get(key, False) or bool(rec["occurs"])
            if rec.get("aa"):
                aa[rec["protein"]] = rec["aa"]
    return occurs, aa


def front_length(path):
    """A front's protein length, for an inherited front with no screen row."""
    with open(path) as fh:
        for line in fh:
            rec = json.loads(line)
            if "points" in rec:
                return len(rec["protein"]["sequence"]), rec.get("complete")
    return None, None


def qualifiers(occurs, aa):
    """Step 1: in range, complete F0, sets I-III all occurring."""
    out = []
    for p in sorted({q for q, _ in occurs}):
        n = aa.get(p)
        if n is None or not (LOW <= n <= HIGH):
            continue
        if all(occurs.get((p, s)) for s in SETS_I3):
            out.append((p, n, bool(occurs.get((p, SET_IV)))))
    return out


def select(pool):
    """Steps 2-4. pool is [(accession, aa, binds_set_IV)]."""
    chosen, used = [], set()

    # 2. all-four binders first, alphabetically among themselves
    for p, n, iv in sorted(pool, key=lambda r: r[0]):
        if iv and len(chosen) < WANTED:
            chosen.append((p, n, iv, "all four sets"))
            used.add(p)

    # 3. ten equal strata, in TWO passes. A stratum with a candidate inside it is
    # filled first; only then do empty strata draw on what is left. A single pass
    # lets an empty stratum's fall-through consume a protein that a later stratum
    # would have matched exactly -- in a dry run, an empty 800-849 stratum took a
    # 962 aa protein that belonged to 950-999, costing the band its top stratum.
    width = (HIGH - LOW) / STRATA
    bounds = [(k, LOW + width * k, LOW + width * (k + 1)) for k in range(STRATA)]

    for k, lo, hi in bounds:
        if len(chosen) >= WANTED:
            break
        mid = (lo + hi) / 2
        inside = [r for r in pool if r[0] not in used and lo <= r[1] < hi]
        if not inside:
            continue
        # 4. alphabetical by accession breaks an exact distance tie
        pick = min(inside, key=lambda r: (abs(r[1] - mid), r[0]))
        chosen.append((pick[0], pick[1], pick[2],
                       "stratum %d (%d-%d aa)" % (k + 1, lo, hi - 1)))
        used.add(pick[0])

    for k, lo, hi in bounds:
        if len(chosen) >= WANTED:
            break
        mid = (lo + hi) / 2
        rest = [r for r in pool if r[0] not in used]
        if not rest:
            break
        if any(lo <= n < hi for _, n, _, _ in chosen):
            continue        # this stratum was already filled in pass one
        pick = min(rest, key=lambda r: (abs(r[1] - mid), r[0]))
        chosen.append((pick[0], pick[1], pick[2],
                       "stratum %d (%d-%d aa), empty, nearest taken"
                       % (k + 1, lo, hi - 1)))
        used.add(pick[0])
    return chosen


def options(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    S = "/scratch/xinyugu4/paper-e3/screen"
    ap.add_argument("--binding", default=os.path.join(S, "binding-l.jsonl"),
                    help="the screened half of the pool")
    ap.add_argument("--inherited", default=os.path.join(S, "binding-five.jsonl"),
                    help="the inherited half, once screen_occurs has run on it")
    ap.add_argument("--inherited-fronts", default="/scratch/xinyugu4/paper-e3/e2/out/hull",
                    help="where the five land, checked only to report what is missing")
    return ap.parse_args(argv)


def main():
    args = options()
    occurs, aa = binding(args.binding)
    screened = len({p for p, _ in occurs})

    have_inherited = os.path.exists(args.inherited)
    if have_inherited:
        o2, a2 = binding(args.inherited)
        occurs.update(o2)
        aa.update(a2)

    pool = qualifiers(occurs, aa)
    print("pool: %d screened at 501-700 aa%s" % (screened,
          ", plus the inherited five" if have_inherited else ""))
    print("qualifiers (I-III all occur, complete F0, 500-1000 aa): %d" % len(pool))
    for p, n, iv in sorted(pool, key=lambda r: r[1]):
        print("  %-8s %4d aa  %s" % (p, n, "binds set IV -> four nests" if iv else "three nests"))

    if not have_inherited:
        print()
        print("THE FIVE HAVE NOT BEEN SCREENED, so no selection is made.")
        print("Expected: Q96NI6 719, Q9BXB4 747, Q9HCI6 787, Q9H3R1 872, Q9Y5G9 962")
        missing = [f for f in ("Q96NI6", "Q9BXB4", "Q9HCI6", "Q9H3R1", "Q9Y5G9")
                   if not os.path.exists(os.path.join(args.inherited_fronts,
                                                      f + ".front.jsonl"))]
        print("Fronts not yet present: %s" % (", ".join(missing) if missing else "none"))
        print()
        print("Selecting from the nineteen alone spans 1.29x in length and was")
        print("withdrawn for that reason, so this is a deliberate refusal rather")
        print("than a missing input.")
        return 1

    if len(pool) < WANTED:
        print()
        print("Only %d qualify, fewer than the %d wanted. The band is those %d, and"
              % (len(pool), WANTED, len(pool)))
        print("the condition for screening the 13 at 701-1000 aa is now live again.")

    chosen = select(pool)
    print()
    print("=== the middle band, %d instance(s) ===" % len(chosen))
    print("%-8s %5s %-12s %s" % ("protein", "aa", "nests", "selected as"))
    for p, n, iv, why in chosen:
        print("%-8s %5d %-12s %s" % (p, n, "F1-F4" if iv else "F1-F3", why))
    four = sum(1 for _, _, iv, _ in chosen if iv)
    print()
    print("cells: %d x 4 + %d x 3 = %d"
          % (four, len(chosen) - four, four * 4 + (len(chosen) - four) * 3))
    lens = [n for _, n, _, _ in chosen]
    print("length span: %d-%d aa (%.2fx)" % (min(lens), max(lens), max(lens) / min(lens)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
