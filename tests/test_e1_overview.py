import numpy as np

from derna2_analysis import e1_overview as ov


def test_group_edges():
    assert [ov.group_of(n) for n in (78, 200, 201, 350, 351, 500, 501, 1000, 1001)] == [
        0, 0, 1, 1, 2, 2, 3, 3, 4,
    ]


def test_finished_slack_and_budget():
    ok = {"exit": "ok", "n_aa": 300, "wall_s": 7200 * 1.02}
    assert ov.finished(ok)
    assert not ov.finished({**ok, "wall_s": 7200 * 1.03})
    assert ov.finished({**ok, "n_aa": 1200, "wall_s": 20000})
    assert not ov.finished({**ok, "exit": "timeout"})
    assert not ov.finished(None)


def test_read_keeps_pin_and_last_cell():
    rows = [
        {"kind": "cell", "arm": "a", "protein": "P", "pin": "x", "wall_s": 1},
        {"kind": "cell", "arm": "a", "protein": "P", "pin": "x", "wall_s": 2},
        {"kind": "cell", "arm": "a", "protein": "P", "pin": "y", "wall_s": 9},
        {"kind": "front_point", "arm": "derna2-walk", "protein": "P", "pin": "x",
         "cai_rescored": 0.9, "mfe_rescored": -10.0},
        {"kind": "front_point", "arm": "mrnafold-hull", "protein": "P", "pin": "x",
         "cai_rescored": 0.8, "mfe_rescored": -20.0},
    ]
    cells, points = ov.read_e1(rows, "x")
    assert cells[("a", "P")]["wall_s"] == 2
    assert points["P"] == [(0.9, -10.0)]


def test_read_falls_back_to_reported_points():
    rows = [
        {"kind": "front_point", "arm": "derna2-walk", "protein": "P", "pin": "x", "cai": 0.9, "mfe": -10.0},
        {"kind": "front_point", "arm": "derna2-walk", "protein": "Q", "pin": "x", "cai": 0.8, "mfe": -9.0},
        {"kind": "front_point", "arm": "derna2-walk", "protein": "Q", "pin": "x",
         "cai": 0.8, "mfe": -9.0, "cai_rescored": 0.81, "mfe_rescored": -9.5},
    ]
    _, points = ov.read_e1(rows, "x")
    assert points["P"] == [(0.9, -10.0)] and points["Q"] == [(0.81, -9.5)]


def test_front_hull_and_steps():
    # (CAI, MFE): the second point lies above the chord of its neighbours, so off the hull;
    # the fifth is dominated.
    f = ov.Front([(0.5, -10.0), (0.7, -5.0), (0.9, -4.5), (0.95, -2.0), (0.6, 0.0)])
    assert list(f.cai.round(2)) == [0.5, 0.7, 0.9, 0.95]
    assert list(f.on_hull) == [True, False, True, True]
    d_cai, d_mfe = f.steps(hull_only=False)
    assert np.allclose(d_cai, [0.2, 0.2, 0.05]) and np.allclose(d_mfe, [5.0, 0.5, 2.5])
    assert len(f.steps(hull_only=True)[0]) == 2
    assert ov.mfe_at_cai(f, 0.8, hull_only=False) == -4.5
    assert ov.mfe_at_cai(f, 0.8, hull_only=True) == -4.5
    assert ov.mfe_at_cai(f, 0.99, hull_only=False) is None


def test_speedups():
    walk = {"wall_s": 50.0}
    assert ov.speedups(walk, {"exit": "ok", "wall_s": 38.0}, sweep_s=19.0) == (2.0, 19.0 / 50.0)
    # An unpruned search that timed out gives no pruning speed-up, but the walk ratio stays.
    assert ov.speedups(walk, {"exit": "timeout", "wall_s": 99.0}, 19.0) == (None, 19.0 / 50.0)
    assert ov.speedups(walk, None, 19.0) == (None, 19.0 / 50.0)
    assert ov.speedups(walk, None, None) == (None, None)


def test_hull_dp_is_measured_where_it_finished_and_projected_where_it_stopped():
    assert ov.hull_dp_s({"exit": "ok", "wall_s": 38.0}) == (38.0, False)
    stopped = {"exit": "timeout", "wall_s": 7200.0, "lambdas_done": 90, "lambdas_planned": 180}
    assert ov.hull_dp_s(stopped) == (14400.0, True)
    none_done = {"exit": "timeout", "wall_s": 7200.0, "lambdas_done": 0, "lambdas_planned": 9}
    assert ov.hull_dp_s(none_done) is None
    assert ov.hull_dp_s(None) is None
