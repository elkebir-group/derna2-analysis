"""E1 tables: one row per cell, and the correctness summary per method and corner.

    python tables.py WORK TIME CELLS_TSV CORRECTNESS_TSV EXPECTED_WORK_ROWS

Correctness is judged on the re-scored MFE (one folder re-folds every design), not on the
energy a tool reports. At the CAI corner every exact method must reach CAI = 1, since
every amino acid has a codon of weight 1.
"""

import sys

import pandas as pd

from derna2_analysis.e1 import cells
from derna2_analysis.io import expect_count, read_jsonl

TOL = 1e-6


def correctness(df):
    rows = []
    for (arm, variant, corner), g in df.groupby(["arm", "variant", "corner"]):
        ok = g[g.work_exit == "ok"]
        gap = ok.mfe_gap_kcal_mol.dropna()
        rows.append(
            {
                "arm": arm,
                "variant": variant,
                "corner": corner,
                "cells": len(g),
                "ok": len(ok),
                "translates": int(ok.translates.fillna(False).sum()),
                "mfe_equal_reference": int((gap.abs() <= TOL).sum()),
                "mfe_above_reference": int((gap > TOL).sum()),
                "mfe_below_reference": int((gap < -TOL).sum()),
                "max_gap_kcal_mol": gap.max() if len(gap) else None,
                "cai_one": int(((ok.cai_rescored - 1.0).abs() <= 1e-9).sum()),
                "time_unclean": int((g.time_clean == False).sum()),  # noqa: E712
            }
        )
    return pd.DataFrame(rows)


def main():
    work, time, cells_tsv, corr_tsv, expected = sys.argv[1:]
    expect_count(read_jsonl(work), int(expected), work)
    df = pd.DataFrame(cells(work, time))
    df.sort_values(["corner", "arm", "variant", "n_aa", "protein"]).to_csv(
        cells_tsv, sep="\t", index=False
    )
    correctness(df).to_csv(corr_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
