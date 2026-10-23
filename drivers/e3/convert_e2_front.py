"""Convert an inherited E2 front into the record schema E3's readers expect.

The middle band's pool is 24: nineteen proteins screened here, and five fronts
inherited from E2 at 719-962 aa. The inherited five are exact whole fronts and
are not re-run, a complete front at 962 aa costing more than the band it would
join. They arrive in E2's format, which differs from E3's in four ways:

  * the energy key is `mfe_kcal_mol` rather than `mfe`;
  * there is no `record` discriminator, so a reader cannot tell a point from a
    header;
  * there is no header at all, so protein identity and length are absent;
  * there is no `on_hull`, so the lambda-reachable subset is not marked.

Only the first three are renames and additions. `on_hull` has to be recomputed,
and recomputing it from a rounded field is where this conversion can go wrong
quietly: E2 prints `cai_log` to six decimals, so each coordinate carries up to
5e-7 of rounding, which is the same order as the solver's own admission slacks.
Whether that flips a hull vertex is a property of the front, not of the format,
so it is MEASURED per front and reported rather than assumed. The quantity is
the smallest vertical distance in `cai_log` from an interior point up to the
hull chain; the conversion is refused when that clearance falls within
`CLEARANCE_FACTOR` times the rounding half-width, because below that margin the
printed precision cannot establish which points are vertices and the front has
to be re-emitted at full precision instead.

On the one E2 front available while this was written -- 32 aa, 31 points, 11
vertices, matching E2's own recorded count -- the clearance was 7.4e-3, four
orders clear of the floor, and two hundred perturbations at the rounding
half-width changed no vertex. One front is not the five, which is why the check
ships as a per-front column rather than as a note saying the rounding is fine.

THE RECOMPUTATION IS VALIDATED AGAINST E3's OWN FLAGS, not against itself: the
same routine run over 21 E3 fronts and 10,103 points reproduces every `on_hull`
exactly, which `check_hull_recompute.py` re-runs. It did not at first. An exact
orientation test marked 1 to 3 extra vertices on 10 of those 21 fronts and
missed none, and each extra sat between 0 and 1.2e-13 above the chain --
collinear to within float noise, which the walk drops and an exact test keeps.
The tolerance that separates the two is not tuned: over the same fronts the
collinear cases reach at most 1.2e-13 in relative terms while the smallest
genuine bend is 6.0e-6, an empty gap of 5.1e7, and any cutoff inside it gives
identical answers. Both gates carry a control arm that differs only in the
quantity under test -- an interior point 1e-5 under the hull is refused and the
same front at 1e-2 converts, and the check above fails if the tolerance decided
no case at all.

WHAT IS NOT RECONSTRUCTED. `wall_s`, `peak_rss_kb`, `threads` and the motif list
describe a run that did not happen here, and `gaps` and `open` describe the walk
that produced the front. All are emitted as null. In particular `gaps` is not
set to `vertices - 1`: E3's completeness test reads exactly that equality, so
filling it in from the vertex count would make the test pass by construction and
report completeness without having checked anything. Completeness is instead
declared by the caller, recorded as `complete` beside the declaration's grounds
in `converted_from`, and left for a consumer to trust or re-derive.
"""

import argparse
import json
import os

E2_KEYS = frozenset(("mfe_kcal_mol", "cai", "cai_log", "sequence", "structure"))
CAI_LOG_DECIMALS = 6
CLEARANCE_FACTOR = 100.0
# a cross product of differences of floats is zero only by luck; this is
# relative to the magnitudes multiplied; on 21 E3 fronts the collinear cases reach
# 1.2e-13 and the smallest genuine bend is 6.0e-6, so this sits in an empty gap
COLLINEAR_RELATIVE = 1e-12


def half_width(decimals):
    """The most a value printed to `decimals` places can have been moved."""
    return 0.5 * 10.0 ** -decimals


def read_e2(path):
    """The point rows of an E2 front, rejecting anything not that shape."""
    rows = []
    with open(path) as fh:
        for number, line in enumerate(fh, 1):
            if not line.strip():
                continue
            rec = json.loads(line)
            if frozenset(rec) != E2_KEYS:
                raise SystemExit(
                    "%s:%d is not an E2 point row: keys %s"
                    % (path, number, sorted(rec))
                )
            rows.append(rec)
    if not rows:
        raise SystemExit("%s holds no rows" % path)
    return rows


def check_one_protein(rows, aa):
    """Every row spells the same protein, at the length the fasta gives."""
    lengths = {len(r["sequence"]) for r in rows}
    if len(lengths) != 1:
        raise SystemExit("rows disagree on length: %s" % sorted(lengths))
    nt = lengths.pop()
    if nt != 3 * aa:
        raise SystemExit("front is %d nt, protein is %d aa (%d nt)" % (nt, aa, 3 * aa))
    sequences = {r["sequence"] for r in rows}
    if len(sequences) != len(rows):
        raise SystemExit(
            "%d rows spell only %d distinct sequences" % (len(rows), len(sequences))
        )


def coordinates(rows):
    """Each row as the (mfe, cai_log) pair the hull is taken over."""
    return [(r["mfe_kcal_mol"], r["cai_log"]) for r in rows]


def upper_hull(points):
    """The lambda-reachable vertices, and the collinear cases absorbed.

    A front trades energy against adaptation, and a scalarization reaches only
    the vertices of the hull in that plane, so membership is a geometric query
    on the pair and not a property any single coordinate carries.

    A point lying ON a segment is not a vertex: it has the same scalarized
    value as its neighbours and is reached only at a tie. In exact arithmetic
    the orientation test returns zero there and the point is dropped, but a
    cross product of two differences of floats returns zero only by luck, so
    the test admits a tolerance RELATIVE to the magnitudes it multiplies. The
    cases the tolerance decides are returned beside the vertices, because a
    tolerance that is never exercised and one that is swallowing real bends
    produce the same vertex count.
    """
    chain, absorbed = [], []
    for point in sorted(set(points)):
        while len(chain) >= 2 and not turns_right(chain[-2], chain[-1], point, absorbed):
            chain.pop()
        chain.append(point)
    return set(chain), absorbed


def turns_right(a, b, c, absorbed):
    """Whether a-b-c bends away from the reachable side, collinear counting as not."""
    ux, uy = b[0] - a[0], b[1] - a[1]
    vx, vy = c[0] - a[0], c[1] - a[1]
    cross = ux * vy - uy * vx
    scale = abs(ux * vy) + abs(uy * vx)
    floor = -COLLINEAR_RELATIVE * scale
    if floor <= cross < 0:
        absorbed.append(abs(cross) / scale if scale else 0.0)
    return cross < floor


def clearance(points, vertices):
    """The smallest gap in `cai_log` from an interior point up to the hull.

    None when no interior point lies under a segment, which is the degenerate
    case of a front that is all vertices -- nothing to misclassify.
    """
    chain = sorted(vertices)
    gaps = [
        gap
        for point in points
        if point not in vertices
        for gap in [distance_below(chain, point)]
        if gap is not None
    ]
    return min(gaps) if gaps else None


def distance_below(chain, point):
    """How far under the hull chain a point sits, in `cai_log`."""
    for left, right in zip(chain, chain[1:]):
        if left[0] <= point[0] <= right[0] and right[0] != left[0]:
            span = (point[0] - left[0]) / (right[0] - left[0])
            return left[1] + (right[1] - left[1]) * span - point[1]
    return None


def resolve(points, vertices):
    """Whether the printed precision can tell vertices from interior points."""
    floor = CLEARANCE_FACTOR * half_width(CAI_LOG_DECIMALS)
    gap = clearance(points, vertices)
    if gap is None:
        return None, floor, True
    return gap, floor, gap >= floor


def read_protein(path):
    """A fasta's name and residues."""
    with open(path) as fh:
        lines = [l.strip() for l in fh if l.strip()]
    if not lines or not lines[0].startswith(">"):
        raise SystemExit("%s is not a fasta" % path)
    return lines[0][1:], "".join(lines[1:])


def header_record(name, residues, rows, declared, source, build, gap, floor,
                  absorbed):
    """The header E3's readers look for, with what was not measured left null."""
    return {
        "record": "header",
        "format": "derna-result",
        "format_version": 1,
        "command": "front",
        "protein": {"name": name, "sequence": residues},
        "points": len(rows),
        "complete": declared == "exact",
        "gaps": None,
        "open": None,
        "wall_s": None,
        "peak_rss_kb": None,
        "threads": None,
        "motifs": [],
        "utr5": None,
        "utr3": None,
        "converted_from": {
            "format": "e2-front",
            "path": os.path.abspath(source),
            "build": build,
            "completeness": declared,
            "on_hull": "recomputed from (mfe, cai_log)",
            "cai_log_decimals": CAI_LOG_DECIMALS,
            "hull_clearance_cai_log": gap,
            "hull_clearance_floor": floor,
            "collinear_absorbed": len(absorbed),
            "collinear_worst_relative": max(absorbed) if absorbed else None,
            "fields_not_measured": [
                "gaps",
                "open",
                "wall_s",
                "peak_rss_kb",
                "threads",
                "motifs",
            ],
        },
    }


def point_record(row, vertices):
    """One E2 row in E3's point schema."""
    return {
        "record": "point",
        "mfe": row["mfe_kcal_mol"],
        "cai": row["cai"],
        "cai_log": row["cai_log"],
        "on_hull": (row["mfe_kcal_mol"], row["cai_log"]) in vertices,
        "sequence": row["sequence"],
        "structure": row["structure"],
    }


def write_front(path, header, rows, vertices):
    """The converted front, header first, points in the order they arrived."""
    with open(path, "w") as fh:
        fh.write(json.dumps(header) + "\n")
        for row in rows:
            fh.write(json.dumps(point_record(row, vertices)) + "\n")


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("source", help="the front in E2 format")
    ap.add_argument("out", help="where to write the converted front")
    ap.add_argument("--protein", required=True, help="the fasta the front designs")
    ap.add_argument("--build", required=True,
                    help="the build that produced the source front")
    ap.add_argument("--completeness", required=True, choices=("exact", "partial"),
                    help="what the source front is, on its own record")
    return ap.parse_args()


def main():
    args = parse_args()
    name, residues = read_protein(args.protein)
    rows = read_e2(args.source)
    check_one_protein(rows, len(residues))

    points = coordinates(rows)
    vertices, absorbed = upper_hull(points)
    gap, floor, resolved = resolve(points, vertices)

    print("%s: %d aa, %d points, %d hull vertices"
          % (name.split()[0], len(residues), len(rows), len(vertices)))
    print("hull clearance in cai_log: %s, floor %.3g"
          % ("none, every point is a vertex" if gap is None else "%.3e" % gap, floor))
    print("collinear cases the tolerance decided: %d, worst %s relative"
          % (len(absorbed), "none" if not absorbed else "%.2e" % max(absorbed)))
    if not resolved:
        print("REFUSED. At six printed decimals a coordinate carries up to %.1e of"
              % half_width(CAI_LOG_DECIMALS))
        print("rounding, and the nearest interior point sits %.3e under the hull," % gap)
        print("so which points are vertices is not established by this file.")
        print("Re-emit the source front with full-precision cai_log.")
        return 1

    header = header_record(name, residues, rows, args.completeness,
                           args.source, args.build, gap, floor, absorbed)
    write_front(args.out, header, rows, vertices)
    print("wrote %s" % os.path.abspath(args.out))
    print("declared %s; gaps, open, wall_s, peak_rss_kb, threads and motifs are null"
          % args.completeness)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
