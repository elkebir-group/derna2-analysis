"""Figures 2 and 3's data builders: the rules a reader of those figures relies on without
seeing them — which rows are read, how a mangled summary row is repaired, and that hull
vertices come from the geometry, not from a stored flag."""

import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "pipelines" / "scripts"


def load(rel):
    spec = importlib.util.spec_from_file_location(Path(rel).stem, SCRIPTS / rel)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_rows(path, rows):
    path.write_text("".join(json.dumps(r) + "\n" for r in rows))
    return path


# A staircase whose middle point sits above the chord: on the front, off the hull.
STAIRCASE = [(0.9, 0.0), (0.8, -5.0), (0.6, -20.0)]


def walk_rows(pin, vertices):
    cell = {
        "kind": "cell",
        "arm": "derna2-walk",
        "protein": "P1",
        "n_aa": 50,
        "exit": "ok",
        "wall_s": 1.0,
        "vertices": vertices,
        "pin": pin,
    }
    points = [
        {
            "kind": "front_point",
            "arm": "derna2-walk",
            "protein": "P1",
            "cai": c,
            "mfe": e,
            "pin": pin,
        }
        for c, e in STAIRCASE
    ]
    return [cell, *points]


def test_fig2_reads_only_the_rows_of_its_pin(tmp_path):
    fig2 = load("e1/fig2_data.py")
    sweep = write_rows(tmp_path / "sweep.jsonl", walk_rows("aaaa", 3) + walk_rows("bbbb", 2))
    cells, points = fig2.read(sweep, "bbbb")
    assert cells[("derna2-walk", "P1")]["vertices"] == 2
    assert len(points["P1"]) == len(STAIRCASE)


def test_fig2_unreachable_share_is_computed_and_a_disagreeing_field_is_reported(tmp_path):
    fig2 = load("e1/fig2_data.py")
    # The stored field claims all three points are vertices; the geometry has two.
    sweep = write_rows(tmp_path / "sweep.jsonl", walk_rows("p", 3))
    cells, points = fig2.read(sweep, "p")
    rows, mismatch = fig2.unreachable(points, cells)
    assert [(r["points"], r["vertices"]) for r in rows] == [(3, 2)]
    assert rows[0]["share"] == pytest.approx(100 / 3)
    assert mismatch == [{"protein": "P1", "stored": 3, "computed": 2}]


def test_fig3_rejoins_a_core_mask_split_across_fields():
    cells = load("e2/fig3_cells.py")
    head = ["P1", "300", "F1", "s", "b", "4"]
    tail = ["0", "60.0", "1024", "10", "4", "True", "0", "1.0"]
    assert cells.repair(head + ["4", "5", "6", "7"] + tail) == head + ["4,5,6,7"] + tail
    assert cells.repair(head + ["4"] + tail) == head + ["4"] + tail


def test_fig3_refuses_a_long_row_that_is_not_a_split_mask():
    cells = load("e2/fig3_cells.py")
    row = ["P1", "300", "F1", "s", "b", "4", "x", "y"] + ["0"] * 8
    assert cells.repair(row) is None


def test_fig3_strip_hull_ignores_a_flag_that_disagrees_with_the_geometry(tmp_path):
    strip = load("e2/fig3_strip.py")
    points = [
        {"record": "point", "cai_log": math.log(c), "mfe": e, "on_hull": True} for c, e in STAIRCASE
    ]
    _, pts = strip.read_front(write_rows(tmp_path / "front.jsonl", [{"complete": True}, *points]))
    assert strip.hull_xy(pts) == [[round(math.log(0.6), 2), -20.0], [round(math.log(0.9), 2), 0.0]]


def test_fig3_price_reads_the_baseline_at_matched_or_higher_cai():
    planned = load("e2/fig3_planned.py")
    base = [{"cai": 0.9, "mfe": -10.0}, {"cai": 0.5, "mfe": -30.0}]
    # At CAI 0.7 the baseline's best among designs at least as adapted is -10, not -30.
    assert planned.price(base, [{"cai": 0.7, "mfe": -8.0}]) == [2.0]
    # A design more adapted than anything in the baseline has no price.
    assert planned.price(base, [{"cai": 0.95, "mfe": 0.0}]) == []
