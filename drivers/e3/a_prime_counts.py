"""A': the counts table. How much of each motif set a free front actually spells.

Publishes counts only, never a price, which is why it was scheduled onto the
contended socket: `occurs` is a substring scan over emitted designs and `spellable`
is an exact lattice calculation, and neither is a timing, so contention cannot
change a number in this table.

Three quantities per (protein, set), and the distinctions between them carry the
whole of the table's meaning:

  spellable          start positions the codon lattice admits at all. Zero means
                     the set forbids nothing on this protein, so any price
                     measured against it is a price against no constraint.
  occurs             occurrences across every design on the free front.
  designs_carrying   how many of those designs hold at least one. This is the one
                     to read as a rate, because `occurs` counts an occurrence
                     twice when two designs share it and a front's design count
                     varies by an order of magnitude across proteins.

What the table cannot say, stated here because the number invites the stronger
reading: these are occurrences on **F0**. The step that prices set k is taken
against F(k-1), and the earlier sets have by then moved the designs. Measured: one
set occurred 7 times on an F0 and 121 times on the F1 beside it, because forbidding
the previous set pushed the walk into designs that spell this one. So a zero in
`occurs` is *absent from the free front*, not *forbids nothing*.

*Forbids nothing* is what `spellable = 0` means, and the two must not be read as
one. Over the 501-700 aa band, **23** of 76 cells have `occurs = 0` and only **7** of
those are unspellable -- every one of the seven being set IV, whose single motif
the lattice cannot place at all on those proteins. The other **16** are motifs the
lattice admits and the free optimum happens to avoid, which is exactly the class
that can reappear once an earlier set is forbidden.
"""

import argparse
import json
import os

SETS = ("manufacturing", "processing", "repeats", "initiation")
ROMAN = {"manufacturing": "I", "processing": "II", "repeats": "III",
         "initiation": "IV"}


def read(path):
    with open(path) as fh:
        return [json.loads(line) for line in fh]


def table(binding, lattice):
    """(protein, set) -> counts, plus each protein's aa length."""
    aa, cells = {}, {}
    for rec in read(lattice):
        if "pp_adjacent" in rec:
            aa[rec["protein"]] = rec["aa"]
    for rec in read(binding):
        key = (rec["protein"], rec["set"])
        cell = cells.setdefault(key, {"spellable": 0, "occurs": 0,
                                      "carrying": 0, "designs": 0,
                                      "motifs": 0, "unspellable": 0})
        cell["spellable"] += rec["spellable"] or 0
        cell["occurs"] += rec["occurs"]
        cell["carrying"] = max(cell["carrying"], rec["designs_carrying"])
        cell["designs"] = rec["front_designs"]
        cell["motifs"] += 1
        cell["unspellable"] += 0 if rec["spellable"] else 1
        aa.setdefault(rec["protein"], rec.get("aa"))
    return cells, aa


def options(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    S = "/scratch/xinyugu4/paper-e3/screen"
    ap.add_argument("--binding", default=os.path.join(S, "binding-l.jsonl"))
    ap.add_argument("--lattice", default=os.path.join(S, "lattice-l.jsonl"))
    ap.add_argument("--label", default="501-700 aa")
    ap.add_argument("--out", default=os.path.join(S, "a-prime-l.csv"))
    return ap.parse_args(argv)


def main():
    args = options()
    cells, aa = table(args.binding, args.lattice)
    proteins = sorted({p for p, _ in cells}, key=lambda p: (aa.get(p) or 0, p))

    with open(args.out, "w") as fh:
        fh.write("protein,aa,set,set_roman,motifs,motifs_unspellable,"
                 "spellable,occurs,designs_carrying,front_designs,"
                 "carrying_share\n")
        for p in proteins:
            for s in SETS:
                c = cells.get((p, s))
                if not c:
                    continue
                share = c["carrying"] / c["designs"] if c["designs"] else 0.0
                fh.write("%s,%s,%s,%s,%d,%d,%d,%d,%d,%d,%.4f\n"
                         % (p, aa.get(p), s, ROMAN[s], c["motifs"],
                            c["unspellable"], c["spellable"], c["occurs"],
                            c["carrying"], c["designs"], share))

    print("A' counts, %s, %d protein(s) x %d set(s)"
          % (args.label, len(proteins), len(SETS)))
    print()
    print("%-8s %5s  %s" % ("protein", "aa",
                            "  ".join("%-18s" % ("set " + ROMAN[s]) for s in SETS)))
    print("%-8s %5s  %s" % ("", "",
                            "  ".join("%-18s" % "carrying / designs" for _ in SETS)))
    for p in proteins:
        row = []
        for s in SETS:
            c = cells.get((p, s))
            if not c:
                row.append("%-18s" % "-")
            elif not c["spellable"]:
                row.append("%-18s" % "unspellable")
            else:
                row.append("%-18s" % ("%d / %d" % (c["carrying"], c["designs"])))
        print("%-8s %5s  %s" % (p, aa.get(p), "  ".join(row)))

    print()
    print("per set, over %d protein(s)" % len(proteins))
    print("%-14s %7s %9s %10s %12s %14s"
          % ("set", "motifs", "occurs>0", "spellable=0", "median share",
             "max share"))
    for s in SETS:
        cs = [cells[(p, s)] for p in proteins if (p, s) in cells]
        on = [c for c in cs if c["occurs"]]
        none = [c for c in cs if not c["spellable"]]
        shares = sorted(c["carrying"] / c["designs"] for c in cs if c["designs"])
        med = shares[len(shares) // 2] if shares else 0.0
        print("%-14s %7d %9s %10s %11.1f%% %13.1f%%"
              % ("%s (%s)" % (s, ROMAN[s]), cs[0]["motifs"] if cs else 0,
                 "%d/%d" % (len(on), len(cs)), "%d/%d" % (len(none), len(cs)),
                 100 * med, 100 * max(shares, default=0.0)))

    vacuous = [(p, s) for p in proteins for s in SETS
               if (p, s) in cells and not cells[(p, s)]["spellable"]]
    print()
    print("cells where the set forbids nothing at all (spellable = 0): %d of %d"
          % (len(vacuous), len(proteins) * len(SETS)))
    for p, s in vacuous:
        print("  %s / set %s -- a price measured here is a price against no"
              " constraint" % (p, ROMAN[s]))
    print()
    print("wrote %s" % args.out)


if __name__ == "__main__":
    main()
