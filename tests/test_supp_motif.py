"""The motif supplement's two load-bearing rules: the gate that refuses a colliding plate,
and what the pinned caches are allowed to depend on.

Both exist because they failed once. The gate is here because a panel letter printed on its
neighbour's y label and was found by cropping a PNG rather than by reading the figure, and
because an overlap check that has never been seen to fire is not a check. The digest tests
are here because the finished-set digest twice moved for a reason that was not a
measurement: once over cells that were still running, once over the directory the deposit
happened to be mounted in.
"""

import importlib.util
import itertools
import os
import random
import re
import sys
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "pipelines" / "scripts" / "e2"
# supp_common reads its paths from the environment at import, as the steps that run it as a
# prefix require. The motif sets and their names are committed, so pointing it at the
# repository's own inputs needs no deposit.
os.environ.setdefault("E2_DEPOSIT", str(ROOT / "data" / "results" / "e2"))
os.environ.setdefault("E2_SUPP", str(ROOT / "data" / "derived" / "e2" / "supp"))
os.environ.setdefault("E2_MOTIF_DIR", str(ROOT / "inputs" / "motifs"))
os.environ.setdefault("E2_MOTIF_SETS",
                      "set1-manufacturing,set2-processing,set3-repeats")


# supp_admits imports supp_common by name, the way the rules run it (PYTHONPATH).
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gate = load("supp_gate")
common = load("supp_common")
admits = load("supp_admits")


# ---- the gate ---------------------------------------------------------------------------

def two_labels(x2):
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.set_axis_off()
    ax.text(.3, .5, "label one")
    ax.text(x2, .5, "label two")
    return fig


def test_the_gate_flags_touching_text_and_passes_text_set_apart():
    # The whole point of the gate: it has to come back empty on one of these and not on the
    # other, or it is reporting the same answer whatever it is shown.
    for x2, expect in ((.31, True), (.80, False)):
        fig = two_labels(x2)
        found = gate.overlaps(fig)
        plt.close(fig)
        assert bool(found) is expect, found


def test_the_gate_refuses_to_pass_a_colliding_figure():
    fig = two_labels(.31)
    with pytest.raises(AssertionError):
        gate.check(fig)
    plt.close(fig)


def test_a_panel_letter_is_compared_with_its_own_panels_title():
    """Owning a panel excuses a letter from intruding on it, not from landing on its title.

    The letter is drawn on the figure and given its panel as owner, so the "strayed into a
    panel it is not part of" test skips that pair deliberately. Grouping the two tests under
    one owner also skipped the letter against the panel's own labels, and that is the pair
    that collided: `(a)` sat one pixel from the heading it labels.
    """
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.set_title("a heading")
    t = fig.text(*ax.title.get_window_extent().transformed(fig.transFigure.inverted()).p0,
                 "(a)", ha="left", va="bottom")
    found = gate.overlaps(fig, {id(t): f"{id(ax)}"})
    plt.close(fig)
    assert any("(a)" in pair[:2] and "a heading" in pair[:2] for pair in found), found


def test_two_texts_in_one_panel_are_compared_unless_they_declare_one_block():
    """Nothing but the author places two texts inside an axes, so they are compared.

    The exemption is opt-in and visible in the drawing code: two lines of one caption share
    a gid, which says their leading is deliberate. Sharing a gid by accident would hide a
    real collision, so both directions are pinned here.
    """
    for gid, expect in ((None, True), ("one-caption", False)):
        fig, ax = plt.subplots(figsize=(3, 2))
        ax.set_axis_off()
        ax.text(.5, .50, "first line", ha="center", gid=gid)
        ax.text(.5, .47, "second line", ha="center", gid=gid)
        found = gate.overlaps(fig)
        plt.close(fig)
        assert bool(found) is expect, found


def test_a_label_over_a_panel_it_does_not_belong_to_is_flagged():
    fig, (left, right) = plt.subplots(1, 2, figsize=(4, 2))
    for ax in (left, right):
        ax.set_axis_off()
    # Placed in the left axes, drawn over the right one.
    left.text(1.6, .5, "wandered", transform=left.transAxes)
    found = gate.overlaps(fig)
    plt.close(fig)
    assert any(p[1] == "<a panel it is not part of>" for p in found), found


def test_a_label_outside_the_page_margin_is_flagged():
    """Off the page, and also merely nearer the edge than the clearance.

    A label a pixel inside the canvas is not off it and still reads as falling off the
    page, so the margin is held to the same clearance as every other pair.
    """
    for x, expect in ((-.05, True), (.002, True), (.5, False)):
        fig, ax = plt.subplots(figsize=(3, 2))
        ax.set_axis_off()
        fig.text(x, .5, "edge")
        found = [p for p in gate.overlaps(fig) if p[1] == "<outside the page margin>"]
        plt.close(fig)
        assert bool(found) is expect, (x, found)


def test_the_clearance_is_a_distance_on_the_page_and_not_a_pixel_count():
    """Doubling the dpi must not halve what the gate asks for."""
    fig = two_labels(.36)
    fig.set_dpi(100)
    at_100 = gate.overlaps(fig)
    fig.set_dpi(200)
    at_200 = gate.overlaps(fig)
    plt.close(fig)
    assert bool(at_100) == bool(at_200), (at_100, at_200)


# ---- what the pinned caches depend on ----------------------------------------------------

ROW = "{cell}\tshort\t{status}\t{path}\t0.1\t0.0\t0.0\t12\n"


def finished_table(tmp_path, rows, stamp="2026-11-02T03:08:08-06:00"):
    text = f"# E2 finished set, written {stamp} by finished_tsv.py\n"
    text += "cell\tband\tstatus\tround_file\tshare\thull_gap\tunvisited_span\tpoints\n"
    text += "".join(ROW.format(**r) for r in rows)
    (tmp_path / "FINISHED.tsv").write_text(text)
    common.E = str(tmp_path)
    return common.finished()[1]


BASE = [{"cell": "P01860_F1", "status": "complete", "path": "/mnt/a/P01860_F1.r5.jsonl"},
        {"cell": "P01860_F2", "status": "open", "path": "/mnt/a/P01860_F2.r3.jsonl"}]


def test_the_finished_digest_ignores_where_the_deposit_is_mounted(tmp_path):
    here = finished_table(tmp_path, BASE)
    moved = [dict(r, path=r["path"].replace("/mnt/a/", "/scratch/someone/else/"))
             for r in BASE]
    assert finished_table(tmp_path, moved) == here


def test_the_finished_digest_ignores_a_round_of_a_cell_that_is_still_open(tmp_path):
    here = finished_table(tmp_path, BASE)
    advanced = [dict(r, path=r["path"].replace(".r3.", ".r4."))
                if r["status"] == "open" else r for r in BASE]
    assert finished_table(tmp_path, advanced) == here


def test_the_finished_digest_follows_a_finished_cells_round(tmp_path):
    here = finished_table(tmp_path, BASE)
    advanced = [dict(r, path=r["path"].replace(".r5.", ".r6."))
                if r["status"] == "complete" else r for r in BASE]
    assert finished_table(tmp_path, advanced) != here


def test_the_finished_digest_ignores_the_tables_write_time(tmp_path):
    assert finished_table(tmp_path, BASE, stamp="2026-01-01T00:00:00-06:00") \
        == finished_table(tmp_path, BASE)


# ---- the admissibility scan is exact ------------------------------------------------------

# Residues with at most four codons, so one peptide enumerates in at most 4**5 sequences.
ALPHABET = "GAPKWEFTV"
PATTERNS = ["GGGGG", "AAUAAA", "GCGGCCGC", "GAAGAC"]


def spellings(peptide):
    return ("".join(c) for c in itertools.product(*(common.CODONS[a] for a in peptide)))


@pytest.mark.parametrize("motif", PATTERNS)
def test_the_scan_agrees_with_every_synonymous_sequence(motif):
    """The window is as many residues as a motif can span at any phase.

    One residue too few turns an exact scan into a sample that still returns a plausible
    number, and nothing downstream could tell: a motif missed here is reported as
    inadmissible, which is the answer that says a constraint removed nothing.
    """
    rx = re.compile("".join(common.IUPAC.get(c, c) for c in motif))
    rng = random.Random(f"scan{motif}")
    peptides = ["".join(rng.choice(ALPHABET) for _ in range(5)) for _ in range(40)]
    # Every peptide that can spell the motif somewhere, found by brute force over the whole
    # sequence rather than over any window.
    truth = {p for p in peptides if any(rx.search(s) for s in spellings(p))}
    got = {p for p in peptides if admits.admits(p, motif)}
    assert got == truth, sorted(got ^ truth)


def test_the_scan_is_not_vacuous():
    """Both answers occur on the peptides the test above draws, for at least one motif."""
    seen = set()
    for motif in PATTERNS:
        rng = random.Random(f"scan{motif}")
        peptides = ["".join(rng.choice(ALPHABET) for _ in range(5)) for _ in range(40)]
        seen |= {admits.admits(p, motif) is not None for p in peptides}
    assert seen == {True, False}
