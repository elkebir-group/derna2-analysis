"""Pareto fronts in (log-CAI, MFE): log-CAI is maximized, MFE (kcal/mol) minimized."""


def pareto(points):
    """The non-dominated (cai_log, energy) points, by decreasing cai_log.

    A point is dominated when another has cai_log at least as high and energy at least as
    low, one of them strictly. Of exact duplicates, one is kept.
    """
    out, best = [], float("inf")
    for a, e in sorted(set(points), key=lambda p: (-p[0], p[1])):
        if e < best:
            out.append((a, e))
            best = e
    return out


def hull_vertices(points):
    """The front points that are optimal for some weighted sum λ·E − (1 − λ)·cai_log with
    λ in [0, 1]: the lower-left convex hull of the front, by decreasing cai_log.

    Only these points are reachable by sweeping a scalarization weight; the rest of the
    front lies in the hull's concave pockets.
    """
    front = pareto(points)
    hull = []
    for p in front:
        # Moving to lower cai_log, energy must fall at an increasing rate: drop the last
        # vertex while it lies on or above the chord from its predecessor to p.
        while len(hull) >= 2:
            (a0, e0), (a1, e1) = hull[-2], hull[-1]
            if (e1 - e0) * (p[0] - a0) <= (p[1] - e0) * (a1 - a0):
                hull.pop()
            else:
                break
        hull.append(p)
    return hull


def hypervolume(points, ref):
    """Area dominated by the front and bounded by ref = (cai_log_ref, energy_ref), the worst
    corner. Points that do not dominate ref add nothing."""
    ra, re = ref
    # Sweep by decreasing cai_log: each point extends the dominated region down to its
    # energy, over the cai_log range from ref up to it.
    area, best = 0.0, re
    for a, e in pareto(points):
        if a <= ra or e >= best:
            continue
        area += (a - ra) * (best - e)
        best = e
    return area
