import json

from derna2_analysis.fig1 import so_front, uncovered
from derna2_analysis.io import write_jsonl


def report(lo, hi, front, complete=True, open_=()):
    return {
        "lo": lo,
        "hi": hi,
        "complete": complete,
        "open": [list(x) for x in open_],
        "front": [{"cai_log": a, "energy_kcal_mol": e, "sequence": s} for a, e, s in front],
    }


def test_so_front_keeps_the_non_dominated_designs_of_seeds_and_windows(tmp_path):
    write_jsonl(
        tmp_path / "seed_stair.jsonl",
        [
            {"cai_log": -10.0, "mfe_kcal_mol": -50.0, "sequence": "S1"},
            {"cai_log": -1.0, "mfe_kcal_mol": -20.0, "sequence": "S2"},
        ],
    )
    (tmp_path / "r1").mkdir()
    rep = report(-10.0, -1.0, [(-5.0, -40.0, "W1"), (-4.0, -45.0, "W2")])
    (tmp_path / "r1" / "w0.jsonl").write_text(json.dumps(rep) + "\n")
    front, _ = so_front(tmp_path)
    # W1 is dominated by W2 (higher cai_log, lower energy).
    assert [p["sequence"] for p in front] == ["S1", "W2", "S2"]


def test_uncovered_counts_a_stopped_window_less_what_it_left_open():
    stopped = {"lo": -6.7, "hi": -3.6, "open": [[-4.8, -3.6]]}
    follow_up = {"lo": -4.8, "hi": -3.6, "open": []}
    assert uncovered([stopped, follow_up], -6.7, -3.6) == []
    assert uncovered([stopped], -6.7, -3.6) == [(-4.8, -3.6)]


def test_uncovered_finds_a_gap_between_windows():
    a = {"lo": -10.0, "hi": -6.0, "open": []}
    b = {"lo": -5.0, "hi": -1.0, "open": []}
    assert uncovered([a, b], -10.0, -1.0) == [(-6.0, -5.0)]
