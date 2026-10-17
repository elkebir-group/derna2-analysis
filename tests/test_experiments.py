import json

from derna2_analysis.e1 import best_time, cells
from derna2_analysis.fig1 import so_front, uncovered
from derna2_analysis.io import write_jsonl


def run(cell, exit="ok", contaminated=False, wall=1.0, rep=0):
    return {
        "cell": cell,
        "rep": rep,
        "exit": exit,
        "contaminated": contaminated,
        "wall_s": wall,
        "user_s": wall,
        "sys_s": 0.0,
    }


def test_best_time_takes_the_latest_clean_run():
    rows = [
        run("c", wall=1),
        run("c", contaminated=True, wall=2),
        run("c", wall=3),
        run("c", contaminated=True, wall=4),
    ]
    assert best_time(rows)["c"]["wall_s"] == 3


def test_best_time_falls_back_to_the_latest_run():
    rows = [run("c", contaminated=True, wall=1), run("c", contaminated=True, wall=2)]
    got = best_time(rows)["c"]
    assert got["wall_s"] == 2 and got["clean"] is False and got["attempts"] == 2


def work_row(arm, protein, mfe, exit="ok", rss=100):
    return {
        "cell": f"{arm}.-.{protein}.mfe",
        "arm": arm,
        "variant": "-",
        "protein": protein,
        "n_aa": 30,
        "corner": "mfe",
        "exit": exit,
        "maxrss_kb": rss,
        "instructions": 1,
        "translates": True,
        "mfe_rescored": mfe,
        "cai_rescored": 0.7,
    }


def test_cells_measure_the_gap_against_the_reference_on_the_same_instance(tmp_path):
    work = [
        work_row("derna-exact-30", "A", -10.0),
        work_row("lineardesign", "A", -9.5),
        work_row("derna-exact-30", "B", -20.0),
        work_row("lineardesign", "B", -20.0),
    ]
    time = [dict(run(r["cell"]), maxrss_kb=1) for r in work]
    write_jsonl(tmp_path / "w.jsonl", work)
    write_jsonl(tmp_path / "t.jsonl", time)
    gaps = {
        (c["arm"], c["protein"]): c["mfe_gap_kcal_mol"]
        for c in cells(tmp_path / "w.jsonl", tmp_path / "t.jsonl")
    }
    assert gaps[("lineardesign", "A")] == 0.5
    assert gaps[("lineardesign", "B")] == 0.0


def test_a_skipped_cell_reads_memory_and_energy_from_its_timed_run(tmp_path):
    ref = work_row("derna-exact-30", "A", -10.0)
    skipped = work_row("lcds-exact", "A", None, exit="skipped", rss=5)
    timed = dict(
        run(skipped["cell"]),
        maxrss_kb=900,
        mfe_rescored=-9.0,
        cai_rescored=0.7,
        instructions=7,
        translates=True,
    )
    write_jsonl(tmp_path / "w.jsonl", [ref, skipped])
    write_jsonl(tmp_path / "t.jsonl", [dict(run(ref["cell"]), maxrss_kb=1), timed])
    got = {c["arm"]: c for c in cells(tmp_path / "w.jsonl", tmp_path / "t.jsonl")}
    assert got["lcds-exact"]["maxrss_kb"] == 900
    assert got["lcds-exact"]["mfe_gap_kcal_mol"] == 1.0


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
