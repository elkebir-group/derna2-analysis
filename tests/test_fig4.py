"""Figure 4's builders: which hull step panel (a) magnifies, which cell row stands for a run,
how the per-protein table picks its rows, and that the overlap check can fail."""

import importlib.util
import json
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "pipelines" / "scripts" / "e1"


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


overview = load("fig4_overview")
steps = load("fig4_steps")


def front(points):
    """(mfe, cai, on_hull) rows as fig4_steps writes them."""
    return np.array(points, float)


def test_the_marked_step_holds_the_most_designs_among_readable_slopes():
    # Two hull steps over ranges of 100 kcal/mol and 0.1 CAI: a steep one (slope 10) with three
    # designs inside, and one of slope 0.53 with two. The steep one cannot show both
    # differences, so the other is marked.
    pts = front(
        [
            (-100, 0.70, 1),
            (-80, 0.701, 0),
            (-70, 0.702, 0),
            (-60, 0.703, 0),
            (-50, 0.705, 1),
            (-40, 0.73, 0),
            (-30, 0.76, 0),
            (0, 0.80, 1),
        ]
    )
    lo, hi = overview.widest_step(pts)
    assert (lo[1], hi[1]) == (0.705, 0.80)


def test_with_no_readable_step_the_fullest_is_marked():
    pts = front(
        [
            (-100, 0.70, 1),
            (-80, 0.701, 0),
            (-70, 0.702, 0),
            (-50, 0.705, 1),
            (-49, 0.75, 0),
            (-48, 0.80, 1),
        ]
    )
    lo, hi = overview.widest_step(pts)
    assert (lo[1], hi[1]) == (0.70, 0.705)


def test_log_limits_always_show_two_decades():
    ticks, (lo, hi) = overview.decades(2e-4, 6e-4)
    shown = [t for t in ticks if lo <= t <= hi]
    assert len(shown) >= 2


def test_a_tie_in_the_rounded_cai_keeps_both_designs():
    # Equal 6-decimal CAI, different log CAI: the worse-energy design still adds CAI, so both
    # are Pareto-optimal on the raw values.
    rows = [(-1000, -16.069031, 0.959017), (-1003, -16.069100, 0.959017)]
    assert len(steps.pareto(rows)) == 2


def test_the_overlap_check_flags_touching_text_and_passes_apart_text():
    for x2, expect in ((0.31, True), (0.8, False)):
        fig, ax = plt.subplots(figsize=(3, 2))
        ax.set_axis_off()
        ax.text(0.3, 0.5, "label one")
        ax.text(x2, 0.5, "label two")
        found = [p for p in overview.check(fig) if "overlaps" in p]
        plt.close(fig)
        assert bool(found) is expect, found


def test_text_over_a_point_is_flagged():
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.scatter([0.5], [0.5], s=20)
    ax.text(0.5, 0.5, "here", ha="center", va="center")
    problems = overview.check(fig)
    plt.close(fig)
    assert any("covers a mark" in p for p in problems)


def write_rows(path, rows):
    path.write_text("".join(json.dumps(r) + "\n" for r in rows))
    return path


def test_the_run_table_keeps_the_delivered_then_faster_attempt(tmp_path):
    cells = load("fig4_cells")
    sweep = write_rows(
        tmp_path / "sweep.jsonl",
        [
            {
                "kind": "cell",
                "pin": "p",
                "arm": "a",
                "protein": "X",
                "n_aa": 10,
                "exit": "timeout",
                "wall_s": 5.0,
            },
            {
                "kind": "cell",
                "pin": "p",
                "arm": "a",
                "protein": "X",
                "n_aa": 10,
                "exit": "ok",
                "wall_s": 9.0,
            },
            {
                "kind": "cell",
                "pin": "p",
                "arm": "a",
                "protein": "X",
                "n_aa": 10,
                "exit": "ok",
                "wall_s": 7.0,
            },
            {
                "kind": "cell",
                "pin": "other",
                "arm": "a",
                "protein": "X",
                "n_aa": 10,
                "exit": "ok",
                "wall_s": 1.0,
            },
        ],
    )
    out = tmp_path / "cells.csv"
    import sys

    argv, sys.argv = sys.argv, ["fig4_cells.py", str(sweep), "p", str(out)]
    try:
        cells.main()
    finally:
        sys.argv = argv
    lines = out.read_text().splitlines()
    assert lines[1].split(",")[3:6] == ["ok", "7.0", "3"]


def test_the_per_protein_table_reads_the_later_cell_of_any_pin(tmp_path):
    per = load("fig4_per_protein")
    sweep = write_rows(
        tmp_path / "sweep.jsonl",
        [
            {
                "kind": "cell",
                "pin": "p",
                "arm": "derna2-walk",
                "protein": "X",
                "n_aa": 10,
                "exit": "ok",
                "wall_s": 60.0,
            },
            {
                "kind": "cell",
                "pin": "q",
                "arm": "derna2-walk",
                "protein": "X",
                "n_aa": 10,
                "exit": "ok",
                "wall_s": 120.0,
            },
            {
                "kind": "cell",
                "pin": "p",
                "arm": "mrnafold-hull-c8",
                "protein": "X",
                "n_aa": 10,
                "exit": "ok",
                "wall_s": 1.0,
            },
        ],
    )
    cell, n_aa = per.cells(sweep)
    assert cell[("derna2-walk", "X")]["wall_s"] == 120.0
    assert ("mrnafold-hull-c8", "X") not in cell
    assert per.dichotomic(cell, tmp_path, "X") is None


@pytest.mark.parametrize("name", ["fig4_cells", "fig4_per_protein", "fig4_steps", "fig4_overview"])
def test_the_scripts_name_no_host_path(name):
    text = (SCRIPTS / f"{name}.py").read_text()
    assert "/scratch/" not in text and "/home/" not in text
