"""The overview figure's hull marks: panel a claims a property of the front, so it is
computed from the points and not read from whatever flag the search that wrote the file
accepted at its own tolerance."""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DEPOSIT = ROOT / "data" / "results" / "fig1" / "a_cdsonly_front.jsonl"


def plot_module():
    spec = importlib.util.spec_from_file_location(
        "fig1_plot", ROOT / "pipelines" / "scripts" / "fig1" / "plot.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_hull_marks_ignore_a_flag_that_disagrees_with_the_geometry():
    # A staircase whose middle point sits above the chord: on the front, off the hull.
    rows = [
        {"cai_log": 0.0, "mfe_kcal_mol": 0.0, "hull_vertex": False},
        {"cai_log": -1.0, "mfe_kcal_mol": -5.0, "hull_vertex": True},
        {"cai_log": -2.0, "mfe_kcal_mol": -20.0, "hull_vertex": False},
    ]
    x, y = plot_module().hull_xy(rows, "mfe_kcal_mol")
    assert list(zip(x, y, strict=True)) == [(-2.0, -20.0), (0.0, 0.0)]


@pytest.mark.skipif(not DEPOSIT.exists(), reason="fig1 deposit not linked")
def test_the_computed_hull_reproduces_the_deposited_flag():
    # The flag and the geometry agree on this front; the panel computes it so that a front
    # where they disagree is drawn by the geometry.
    from derna2_analysis.io import read_jsonl

    rows = read_jsonl(DEPOSIT)
    x, _ = plot_module().hull_xy(rows, "mfe_kcal_mol")
    assert len(x) == sum(1 for r in rows if r.get("hull_vertex")) == 21
