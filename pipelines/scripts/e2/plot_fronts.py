"""E2: what the whole front costs each method, from the sweep rows of one pin.

    python plot_fronts.py SWEEP_JSONL PIN OUT_PDF OUT_TSV

Left: a cell's wall time (its whole λ sweep, or derna2's walk) against protein length.
Right: Pareto points each cell recovered, against length. Cells that did not finish (no
hull, timeout, out of memory) are in the table and left out of the plot.

Only rows stamped with PIN are read: a sweep file accumulates across re-pins, and rows of two
builds are not one experiment. The table holds one row per (protein, arm) cell.
"""

import sys

import pandas as pd

from derna2_analysis.figures import SLOTS, save, style
from derna2_analysis.io import read_jsonl, rows_of_pin

ARMS = {
    "derna2-walk": ("derna2, whole front", SLOTS[0], "o"),
    "derna2-hull": ("derna2, hull sweep", SLOTS[6], "X"),
    "lineardesign-hull": ("LinearDesign", SLOTS[1], "s"),
    "lcds-exact-hull": ("LinearCDSfold", SLOTS[2], "^"),
    "mrnafold-hull": ("mRNAfold", SLOTS[3], "D"),
    "derna1-hull": ("DERNA v1, hull λs", SLOTS[4], "v"),
    "derna1-own": ("DERNA v1, own sweep", SLOTS[7], "h"),
}


def cell_table(rows):
    """One record per (protein, arm): the cell row's wall time and status, and the number
    of distinct non-dominated points among the designs it returned (re-scored, so every
    arm is counted on the same folder)."""
    out = []
    for c in (r for r in rows if r.get("kind") == "cell"):
        designs = [
            (r["cai_rescored"], r["mfe_rescored"])
            for r in rows
            if r["arm"] == c["arm"]
            and r["protein"] == c["protein"]
            and r.get("mfe_rescored") is not None
        ]
        front = {
            p for p in designs if not any(q != p and q[0] >= p[0] and q[1] <= p[1] for q in designs)
        }
        out.append(
            {
                "protein": c["protein"],
                "n_aa": c["n_aa"],
                "arm": c["arm"],
                "exit": c.get("exit"),
                "wall_s": c.get("wall_s"),
                "maxrss_kb": c.get("maxrss_kb"),
                "front_points": len(front),
            }
        )
    return pd.DataFrame(out)


def main():
    sweep, pin, out_pdf, out_tsv = sys.argv[1:]
    rows = rows_of_pin(read_jsonl(sweep), pin)
    df = cell_table(rows)
    df.to_csv(out_tsv, sep="\t", index=False)
    style()
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharex=True)
    for arm, (label, color, marker) in ARMS.items():
        g = df[(df.arm == arm) & (df.exit == "ok")].sort_values("n_aa") if len(df) else df
        if g.empty:
            continue
        for ax, field in zip(axes, ["wall_s", "front_points"], strict=True):
            ax.plot(g.n_aa, g[field], color=color, marker=marker, label=label)
    axes[0].set_ylabel("wall time per cell (s)")
    axes[0].set_yscale("log")
    axes[1].set_ylabel("Pareto points recovered")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("protein length (aa)")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7)
    fig.suptitle(
        f"derna2 pin {pin}: {df.protein.nunique() if len(df) else 0} instances", fontsize=8
    )
    save(fig, out_pdf)


if __name__ == "__main__":
    main()
