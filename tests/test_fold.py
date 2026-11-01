"""`free_fold` against a design whose energy is recorded, and against a real binary.

The binary writes an integer and the rows everywhere else are kcal/mol, so what these
tests are about is the factor between them and the refusals that catch a binary whose
output does not line up with what it was asked to fold.
"""

import textwrap
from pathlib import Path

import pytest

from derna2_analysis.fold import free_fold, free_folds

# One deposited design and the energy data/results/fig1 records for it: a change in the
# factor, or in which column is read, shows up here as a number off by a power of ten.
DESIGN = (
    "AUGUAUGGCAAGAUCAUCUUCGUCCUGCUGCUCUCGGGGAUCGUGUCGAUCUCCGCGAGCAGCACGACGGGGGUGGCC"
    "AUGCAUACGAGUACCAGCAGUAGCGUGACUAAGAGUUAUAUAUCCUCACAGACCAACGGCAUCACCUUGAUAAAUUGG"
    "UGGGCGAUGGCCCGCGUAAUUUUCGAGGUGAUGCUGGUGGUCGUGGGGAUGAUAAUUCUUAUCAGCUACUGCAUUCGU"
)
DESIGN_MFE = -148.7


def fake_fold(tmp_path, body):
    """A stand-in for the binary: whatever `body` prints, per line of input."""
    path = tmp_path / "fake_fold"
    path.write_text("#!/usr/bin/env python3\n" + textwrap.dedent(body))
    path.chmod(0o755)
    return path


def test_the_energy_is_hundredths_of_a_kcal(tmp_path):
    fold = fake_fold(
        tmp_path,
        """
        import sys
        for line in sys.stdin:
            print(("." * len(line.strip())) + "\\t-14870")
        """,
    )
    assert free_fold(fold, DESIGN)["mfe_kcal_mol"] == DESIGN_MFE


def test_the_structure_comes_back_as_written(tmp_path):
    fold = fake_fold(tmp_path, """print("(((...)))\\t-120")""")
    assert free_fold(fold, "GGGAAACCC") == {"structure": "(((...)))", "mfe_kcal_mol": -1.2}


def test_a_structure_of_the_wrong_length_is_refused(tmp_path):
    fold = fake_fold(tmp_path, """print(".....\\t-100")""")
    with pytest.raises(ValueError, match="structure of 5"):
        free_fold(fold, "GGGAAACCC")


def test_fewer_folds_than_sequences_is_refused(tmp_path):
    fold = fake_fold(tmp_path, """print("...\\t-10")""")
    with pytest.raises(ValueError, match="1 fold"):
        free_folds(fold, ["GGG", "CCC"])


def test_a_sequence_that_is_not_one_is_refused(tmp_path):
    with pytest.raises(ValueError, match="not a sequence"):
        free_folds(tmp_path / "absent", ["GGG AAA"])


def test_nothing_to_fold_runs_nothing(tmp_path):
    assert free_folds(tmp_path / "absent", []) == []


def test_one_process_folds_every_sequence(tmp_path):
    """Thousands of designs are folded per figure, so the call is not per sequence."""
    counter = tmp_path / "calls"
    fold = fake_fold(
        tmp_path,
        f"""
        import pathlib, sys
        p = pathlib.Path({str(counter)!r})
        p.write_text(str(int(p.read_text()) + 1) if p.exists() else "1")
        for line in sys.stdin:
            print(("." * len(line.strip())) + "\\t-100")
        """,
    )
    assert len(free_folds(fold, ["GGG", "CCC", "AAA"])) == 3
    assert counter.read_text() == "1"


def built_mfe_fold():
    """A built `mfe_fold` under data/derived/bin, or None where none has been built."""
    return next(iter(sorted(Path("data/derived/bin").glob("derna2-*/mfe_fold"))), None)


@pytest.mark.skipif(built_mfe_fold() is None, reason="no built mfe_fold under data/derived/bin")
def test_the_real_binary_reproduces_a_deposited_energy():
    got = free_fold(built_mfe_fold(), DESIGN)
    assert got["mfe_kcal_mol"] == DESIGN_MFE
    assert len(got["structure"]) == len(DESIGN)
