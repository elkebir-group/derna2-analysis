"""Does the hull recomputation reproduce the flags the walk itself wrote?

`convert_e1_front` has to mark `on_hull` on inherited fronts that do not carry
it, so the routine that marks it is checked against fronts that DO -- every E2
front on disk, which the walk flagged itself. A converter validated only on its
own output would pass whatever it did.

Two numbers are asserted, not just the first:

  * agreement, which is the point of the check; and
  * how many collinear cases the orientation tolerance actually decided, which
    is the control. A tolerance that never fires and one that is swallowing
    real bends both produce agreement here, so the run reports the largest
    relative cross it absorbed and the smallest genuine bend it kept. Those two
    populations have to stay separated by orders of magnitude; if they ever
    approach each other, the cutoff is inside live data and the agreement above
    stops meaning anything.
"""

import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from convert_e1_front import COLLINEAR_RELATIVE, upper_hull  # noqa: E402

DEFAULT_GLOB = "run/paper-e3/grid/*.jsonl"


def points_of(path):
    """The point rows of an E2 front."""
    with open(path) as fh:
        records = [json.loads(line) for line in fh]
    return [r for r in records if r.get("record") == "point"]


def compare(points):
    """(recomputed vertices, the walk's own, the absorbed magnitudes)."""
    vertices, absorbed = upper_hull([(r["mfe"], r["cai_log"]) for r in points])
    mine = {r["sequence"] for r in points if (r["mfe"], r["cai_log"]) in vertices}
    theirs = {r["sequence"] for r in points if r.get("on_hull")}
    return mine, theirs, absorbed


def smallest_kept_bend(points):
    """The smallest relative cross the walk treated as a genuine bend."""
    chain, kept = [], []
    for point in sorted({(r["mfe"], r["cai_log"]) for r in points}):
        while len(chain) >= 2:
            cross, scale = orientation(chain[-2], chain[-1], point)
            if scale and cross < -COLLINEAR_RELATIVE * scale:
                kept.append(abs(cross) / scale)
                break
            chain.pop()
        chain.append(point)
    return min(kept) if kept else None


def orientation(a, b, c):
    """The cross product of a-b against a-c, and the magnitude it came from."""
    ux, uy = b[0] - a[0], b[1] - a[1]
    vx, vy = c[0] - a[0], c[1] - a[1]
    return ux * vy - uy * vx, abs(ux * vy) + abs(uy * vx)


def report(path, points):
    """One front's row, and whether it agreed."""
    mine, theirs, absorbed = compare(points)
    agreed = mine == theirs
    print("%-30s %5d pts  hull %3d/%3d  absorbed %2d (worst %s)  %s"
          % (os.path.basename(path), len(points), len(mine), len(theirs),
             len(absorbed), "-" if not absorbed else "%.1e" % max(absorbed),
             "MATCH" if agreed else "DIFFER +%d/-%d"
             % (len(mine - theirs), len(theirs - mine))))
    return agreed, absorbed


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fronts", default=DEFAULT_GLOB)
    args = ap.parse_args()

    disagreed, absorbed, kept, total = 0, [], [], 0
    for path in sorted(glob.glob(args.fronts)):
        points = points_of(path)
        if not points:
            continue
        total += len(points)
        agreed, cases = report(path, points)
        disagreed += not agreed
        absorbed.extend(cases)
        bend = smallest_kept_bend(points)
        if bend is not None:
            kept.append(bend)

    print()
    print("%d point(s) over %d front(s); disagreeing fronts: %d"
          % (total, len(kept), disagreed))
    if not absorbed:
        print("CONTROL FAILED: the orientation tolerance decided nothing, so this")
        print("run does not establish that collinear points are handled at all.")
        return 1
    worst, finest = max(absorbed), min(kept)
    print("collinear absorbed: %d, worst %.3e relative" % (len(absorbed), worst))
    print("genuine bends kept: smallest %.3e relative" % finest)
    print("cutoff %.0e sits %.0fx above the worst absorbed and %.3gx below the"
          % (COLLINEAR_RELATIVE, COLLINEAR_RELATIVE / worst, finest / COLLINEAR_RELATIVE))
    print("smallest kept, an empty gap of %.1e" % (finest / worst))
    return 1 if disagreed else 0


if __name__ == "__main__":
    raise SystemExit(main())
