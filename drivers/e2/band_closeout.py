"""Assemble the long band's 19 summary rows from the fronts, in length order.

The band accumulated its rows three ways: four went straight into `summary.csv`,
twelve were written as per-protein fragments under `summary.d/`, and three
(P07947, P40222, Q494X3) have a complete front on disk and no row at all,
their driver having been killed before it wrote one.

Every column except `exit` is recoverable from a front's own header, so the three
missing rows can be rebuilt. That rebuild is only trustworthy if the same
procedure reproduces the rows that were actually recorded, so it is run over all
nineteen fronts and checked against the sixteen recorded rows first. A
reconstruction that cannot reproduce a measured row is not evidence about the
rows it is standing in for.

Two columns are the driver's observations rather than the walk's, and neither can
be rebuilt:

`exit` is not carried by any front. For a front that reports `open: []` and a
settled gap for every hull interval, exit 0 is an inference from completeness
rather than a measurement, and the provenance column says so.

`wall_s` is the driver's `$(date +%s) - start` around the whole process, which
brackets the walk's own `wall_s` and exceeds it by 1-2 s on every one of the
sixteen recorded rows -- process startup, the front write and teardown. They are
different quantities, so they get different columns: `wall_s` stays the driver's
and is blank where no driver row survives, and `walk_wall_s` is the header's and
is present for all nineteen. Carrying the header's figure under the driver's name
would make a column that silently means one thing on sixteen rows and another on
three, and the band's length trend is read off that column.
"""

import json
import os
import sys

BAND = "run/paper-e3/screen/f0-l"
FRAG = os.path.join(BAND, "summary.d")
COLUMNS = ("protein", "aa", "exit", "wall_s", "points", "hull", "gaps")
# the driver's wall is not rebuildable; the walk's own is, and they differ
REBUILDABLE = ("protein", "aa", "exit", "points", "hull", "gaps")


def front_row(acc):
    """The summary row a front implies, and whether the front is complete."""
    with open(os.path.join(BAND, acc + ".jsonl")) as fh:
        records = [json.loads(line) for line in fh]
    header = next(r for r in records if "points" in r)
    points = [r for r in records if "sequence" in r]
    hull = sum(1 for r in points if r.get("on_hull"))
    complete = not header.get("open") and header.get("gaps") == hull - 1
    return {
        "protein": acc,
        "aa": len(header["protein"]["sequence"]),
        "exit": 0 if complete else None,
        "walk_wall_s": round(header["wall_s"]),
        "points": header["points"],
        "hull": hull,
        "gaps": header["gaps"],
    }, complete


def recorded_rows():
    """Every row already on disk, from summary.csv and from the fragments."""
    rows = {}
    path = os.path.join(BAND, "summary.csv")
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh.read().splitlines()[1:]:
                if line.strip():
                    rows[line.split(",")[0]] = (line, "summary.csv")
    for name in sorted(os.listdir(FRAG)):
        if name.endswith(".csv"):
            with open(os.path.join(FRAG, name)) as fh:
                line = fh.read().strip()
            if line:
                rows[line.split(",")[0]] = (line, "summary.d/" + name)
    return rows


def as_fields(row, columns):
    return ["" if row.get(c) is None else str(row.get(c)) for c in columns]


def main():
    accs = sorted(f[:-6] for f in os.listdir(BAND) if f.endswith(".jsonl"))
    recorded = recorded_rows()

    rebuilt, incomplete, mismatched = {}, [], []
    for acc in accs:
        row, complete = front_row(acc)
        rebuilt[acc] = row
        if not complete:
            incomplete.append(acc)

    print("=== control: does the rebuild reproduce the rows that were recorded? ===")
    print("    (on every column a front can carry; the driver's wall_s is compared")
    print("     separately below, being a different quantity)")
    idx = {c: COLUMNS.index(c) for c in REBUILDABLE}
    gaps_in_wall = []
    for acc in accs:
        if acc not in recorded:
            continue
        line, source = recorded[acc]
        got = line.split(",")
        mine = as_fields(rebuilt[acc], REBUILDABLE)
        theirs = [got[idx[c]] for c in REBUILDABLE]
        agree = mine == theirs
        if not agree:
            mismatched.append(acc)
        print("  %-8s %-22s %s" % (acc, source, "agrees" if agree else "DIFFERS"))
        if not agree:
            for c, a, b in zip(REBUILDABLE, theirs, mine):
                if a != b:
                    print("      %-8s recorded %s, rebuilt %s" % (c, a, b))
        driver = int(got[COLUMNS.index("wall_s")])
        rebuilt[acc]["wall_s"] = driver
        gaps_in_wall.append((acc, driver, rebuilt[acc]["walk_wall_s"]))
    print("  %d of %d recorded rows reproduced on every rebuildable column" % (
        len(recorded) - len(mismatched), len(recorded)))

    print()
    print("=== the two walls, on the rows that carry both ===")
    deltas = [d - w for _, d, w in gaps_in_wall]
    print("  driver wall minus walk wall: min %+d s, max %+d s, over %d rows"
          % (min(deltas), max(deltas), len(deltas)))
    if min(deltas) < 0:
        print("  *** a driver wall BELOW the walk's own is impossible: the driver")
        print("  *** brackets the walk. Investigate before using either column.")
        return 1
    print("  always non-negative, as the bracketing requires, and never more than")
    print("  %d s -- so the walk's wall is a sound stand-in for a length trend," % max(deltas))
    print("  but it is reported as its own column rather than as the driver's.")

    if mismatched:
        print("\nThe rebuild does not reproduce measured rows, so it cannot stand in")
        print("for the missing ones. Resolve the disagreement before writing.")
        return 1
    if incomplete:
        print("\nIncomplete fronts, which get no reconstructed exit code: %s"
              % ", ".join(incomplete))

    header = ("protein", "aa", "exit", "wall_s", "walk_wall_s",
              "points", "hull", "gaps")
    out = os.path.join(BAND, "summary.full.csv")
    with open(out, "w") as fh:
        fh.write(",".join(header) + ",provenance\n")
        for acc in sorted(accs, key=lambda a: rebuilt[a]["aa"]):
            if acc in recorded:
                source = recorded[acc][1]
            else:
                source = "rebuilt from front; exit inferred from completeness, no driver wall"
            fh.write("%s,%s\n" % (",".join(as_fields(rebuilt[acc], header)), source))

    print("\n=== %s, %d rows, shortest first ===" % (out, len(accs)))
    print("%-8s %5s %5s %9s %11s %8s %6s %6s  %s"
          % ("protein", "aa", "exit", "wall_s", "walk_wall_s", "points",
             "hull", "gaps", "provenance"))
    for acc in sorted(accs, key=lambda a: rebuilt[a]["aa"]):
        r = rebuilt[acc]
        source = recorded[acc][1] if acc in recorded else "REBUILT from front"
        print("%-8s %5d %5s %9s %11d %8d %6d %6d  %s"
              % (acc, r["aa"], r["exit"],
                 r.get("wall_s") if r.get("wall_s") is not None else "-",
                 r["walk_wall_s"], r["points"], r["hull"], r["gaps"], source))
    return 0


if __name__ == "__main__":
    sys.exit(main())
