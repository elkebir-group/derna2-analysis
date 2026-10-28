"""Every pinned deposit says where it came from, and the fingerprint that ties rows to a
binary is the one the code computes."""

import json
import os
from pathlib import Path

import pytest

from derna2_analysis.provenance import fnv1a64

ROOT = Path(__file__).resolve().parents[1]
PINS = ROOT / "pins"
WINDOWS = ROOT / "data" / "results" / "fig1" / "so_windows"
# The window binary this test can check the stamp against, where the host has it. Named by
# the environment rather than written here: the path is one host's, and a test that spelled
# it out would skip everywhere else while looking as though it had run.
FIG1_BIN = Path(os.environ["FIG1_WINDOW_BIN"]) if "FIG1_WINDOW_BIN" in os.environ else None
STAMP = "1ff256ead61c92ba"


def names_files(manifest):
    """A manifest that is all comments is a freeze that has not happened yet; `verify`
    refuses it for the same reason, so it carries nothing to attribute."""
    return any(
        line.strip() and not line.startswith("#") for line in manifest.read_text().splitlines()
    )


def pinned_names():
    """The experiment each non-empty manifest under pins/ belongs to."""
    return sorted(
        {
            p.name.removeprefix("results-").removesuffix(".sha256")
            for p in PINS.glob("*.sha256")
            if names_files(p)
        }
    )


def test_there_are_pins_to_check():
    assert pinned_names()


@pytest.mark.parametrize("name", pinned_names())
def test_every_pinned_experiment_records_its_provenance(name):
    # A manifest says the files have not changed; it does not say what produced them, and
    # a deposit whose build is not written down anywhere cannot be rebuilt or compared.
    assert (PINS / f"{name}-PROVENANCE.md").is_file()


def test_fnv1a64_is_the_stamp_the_runs_wrote():
    # The window reports' binary_fnv64 identifies a pinned binary, which is what lets rows
    # from a build with no recorded commit still be attributed.
    assert fnv1a64(__file__) == fnv1a64(__file__)
    if FIG1_BIN is None or not FIG1_BIN.is_file():
        pytest.skip("set FIG1_WINDOW_BIN to the pinned scope37_window to check its stamp")
    assert fnv1a64(FIG1_BIN) == STAMP


@pytest.mark.skipif(not WINDOWS.is_dir(), reason="fig1 deposit not linked")
def test_the_window_reports_are_one_build_and_the_design_rows_are_unstamped():
    # Two kinds of row live here. A report stamps its binary; a design row is a sequence
    # and its two coordinates, and carries nothing about the run — so a design row cannot
    # be attributed by inspection, and reading the absence as a second build is wrong.
    stamped, bare = set(), 0
    for path in WINDOWS.rglob("*.jsonl"):
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("binary_fnv64"):
                stamped.add(row["binary_fnv64"])
            else:
                bare += 1
                assert set(row) == {"cai_log", "mfe_kcal_mol", "sequence"}, sorted(row)
    assert stamped == {STAMP}
    assert bare == 18
