import itertools
import random

from derna2_analysis.fronts import hull_vertices, hypervolume, pareto


def brute_pareto(points):
    return sorted(
        {p for p in points if not any(q != p and q[0] >= p[0] and q[1] <= p[1] for q in points)},
        key=lambda p: -p[0],
    )


def brute_hull(points):
    """Points optimal for some λ on a fine grid of weights; ties broken toward both ends."""
    out = set()
    for k in range(2001):
        lam = k / 2000
        scores = [lam * e - (1 - lam) * a for a, e in points]
        best = min(scores)
        out |= {p for p, s in zip(points, scores, strict=True) if abs(s - best) < 1e-12}
    return out


def random_points(rng, n):
    return [(rng.randint(-40, 0), rng.randint(-60, 0)) for _ in range(n)]


def test_pareto_matches_brute_force():
    rng = random.Random(1)
    for _ in range(200):
        pts = random_points(rng, rng.randint(1, 25))
        assert pareto(pts) == brute_pareto(pts)


def test_hull_vertices_are_exactly_the_scalarization_optima():
    rng = random.Random(2)
    for _ in range(200):
        pts = random_points(rng, rng.randint(1, 25))
        hull = set(hull_vertices(pts))
        # Every hull vertex wins some weight, and every point that wins a weight on the
        # grid and is a front point is a vertex or lies on a hull edge (collinear).
        reached = brute_hull(pareto(pts))
        assert hull <= reached
        assert all(p in hull or on_hull_edge(p, sorted(hull)) for p in reached)


def on_hull_edge(p, hull):
    for (a0, e0), (a1, e1) in itertools.pairwise(hull):
        if min(a0, a1) <= p[0] <= max(a0, a1) and (e1 - e0) * (p[0] - a0) == (p[1] - e0) * (
            a1 - a0
        ):
            return True
    return False


def test_hull_drops_a_point_in_a_concave_pocket():
    # (−1, −5) lies above the chord from (0, 0) to (−2, −20): no weight picks it.
    assert hull_vertices([(0, 0), (-1, -5), (-2, -20)]) == [(0, 0), (-2, -20)]


def test_hypervolume_of_a_staircase():
    # Rectangles from ref (−10, 0): (0,−1) gives 10×1, (−5,−4) adds 5×3.
    assert hypervolume([(0, -1), (-5, -4)], ref=(-10, 0)) == 10 + 15


def test_hypervolume_ignores_dominated_and_out_of_range_points():
    base = hypervolume([(0, -1), (-5, -4)], ref=(-10, 0))
    assert hypervolume([(0, -1), (-5, -4), (-6, -2), (-11, -50), (1, 5)], ref=(-10, 0)) == base


def test_hypervolume_matches_grid_count():
    rng = random.Random(3)
    for _ in range(50):
        pts = random_points(rng, rng.randint(1, 10))
        ref = (-41, 1)
        cells = sum(
            1
            for x in range(-41, 0)
            for y in range(-60, 1)
            if any(x + 1 <= a and y >= e for a, e in pts)
        )
        assert hypervolume(pts, ref) == cells
