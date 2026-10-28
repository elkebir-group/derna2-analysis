"""The revised overview figure's drawing rules that a reader of the figure relies on: the protein
drawn is the protein designed, a cross marks exactly the point it forbids, and a zoom goes where
it covers the least of the molecule."""

import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("RNA")  # the structure layouts come from ViennaRNA

import matplotlib  # noqa: E402

matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
STRUCTURE = ROOT / "inputs" / "structures" / "AF-P15421-F1-model_v6.pdb"
FASTA = ROOT / "inputs" / "glpe_hbb" / "GLPE_HUMAN_78.fasta"
THREE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def overview2():
    spec = importlib.util.spec_from_file_location(
        "overview2", ROOT / "pipelines" / "scripts" / "fig1" / "overview2.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_drawn_structure_is_the_designed_protein():
    ca = [
        line
        for line in STRUCTURE.read_text().splitlines()
        if line.startswith("ATOM") and line[12:16].strip() == "CA"
    ]
    protein = "".join(
        line.strip() for line in FASTA.read_text().splitlines() if not line.startswith(">")
    )
    assert "".join(THREE[line[17:20]] for line in ca) == protein
    xyz, helix, plddt = overview2().ca_trace(STRUCTURE)
    assert len(xyz) == len(helix) == len(plddt) == len(protein)
    assert helix.any() and not helix.all()


def test_a_cross_is_centred_on_the_point_it_marks():
    ov2 = overview2()
    fig = plt.figure(figsize=(6.5, 2.95))
    before = len(fig.artists)
    ov2.forbid(fig, 0.4, 0.7, "red")
    strokes = fig.artists[before:]
    assert len(strokes) == 2
    for stroke in strokes:
        xs, ys = stroke.get_data()
        assert (sum(xs) / 2, sum(ys) / 2) == pytest.approx((0.4, 0.7))
    plt.close(fig)


def test_a_zoom_goes_to_the_emptier_corner():
    ov2 = overview2()
    fig, ax = plt.subplots()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    # every base in the left half: the right-hand candidate covers none of them
    xs, ys = [0.1, 0.2, 0.3, 0.1], [0.9, 0.8, 0.2, 0.5]
    left, right = [0.0, 0.6, 0.4, 0.4], [0.6, 0.6, 0.4, 0.4]
    assert ov2.emptiest(ax, xs, ys, [left, right]) == right
    plt.close(fig)


def test_every_front_states_every_property():
    ov2 = overview2()
    assert ov2.PROPERTIES == ("UTRs", "motifs", "complete", "modular")
    for *_, has in ov2.FRONTS:
        assert len(has) == len(ov2.PROPERTIES)
    # only DERNA2's front is modular, and it is the only front drawn with squares
    assert [m for _, _, m, *_ in ov2.FRONTS] == ["s", "o", "o"]
    # the circles show only their filled glyph; the legend under the table gives the fills
    assert [fills for *_, fills, _ in ov2.FRONTS] == [(True, False), (True,), (True,)]
    assert [has[-1] for *_, has in ov2.FRONTS] == [True, False, False]


def test_the_marker_legend_has_one_glyph_per_shape_and_fill():
    ov2 = overview2()
    ov2.ov.set_style()
    fig = plt.figure(figsize=(6.5, 2.95))
    before = len(fig.artists)
    ov2.marker_legend(fig, 0.03, 0.49, 0.5, 5.8)
    glyphs = {(a.get_marker(), a.get_markerfacecolor()) for a in fig.artists[before:]}
    assert len(fig.artists) - before == 4
    assert {m for m, _ in glyphs} == {"s", "o"}
    assert len(glyphs) == 4
    plt.close(fig)
