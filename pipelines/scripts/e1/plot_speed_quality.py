"""E1 approximate modes: what each beam width buys and costs against derna2's exact run.

    python plot_speed_quality.py CELLS_TSV OUT_PDF

One column per approximate method; the top row is the re-scored MFE gap to the exact
optimum (kcal/mol, 0 = optimal), the bottom row wall time as a multiple of derna2's exact
run on the same cell. MFE corner, P0DTC2 prefix ladder. Color steps from light to dark
with beam width.
"""

import sys

import pandas as pd

from derna2_analysis.figures import INK_2, WIDTH_RAMP, label_line_end, save, style

FAMILIES = [
    ("derna-open2", "derna2 beam, 2 open"),
    ("derna-open5", "derna2 beam, 5 open"),
    ("lcds-beam", "LinearCDSfold beam"),
]
WIDTHS = ["b100", "b500", "b2000"]


def main():
    cells_tsv, out = sys.argv[1:]
    df = pd.read_csv(cells_tsv, sep="\t")
    df = df[(df.corner == "mfe") & df.ladder]
    exact = df[(df.arm == "derna-exact-30") & (df.time_exit == "ok")].set_index("protein")
    style()
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.6), sharex=True, sharey="row")
    for col, (arm, title) in enumerate(FAMILIES):
        for width, color in zip(WIDTHS, WIDTH_RAMP, strict=True):
            g = df[(df.arm == arm) & (df.variant == width)].sort_values("n_aa")
            g = g[g.protein.isin(exact.index)]
            if g.empty:
                continue
            ratio = g.wall_s.values / exact.loc[g.protein, "wall_s"].values
            for row, y in enumerate([g.mfe_gap_kcal_mol.values, ratio]):
                ax = axes[row][col]
                ax.plot(g.n_aa, y, color=color, marker="o", markersize=2.5)
                if col == len(FAMILIES) - 1:
                    label_line_end(ax, g.n_aa.iloc[-1], y[-1], f"width {width[1:]}")
        axes[0][col].set_title(title)
        axes[1][col].axhline(1.0, color=INK_2, linewidth=0.8, linestyle="--")
        axes[1][col].set_xlabel("protein length (aa)")
        for row in range(2):
            axes[row][col].set_xscale("log")
    axes[1][0].set_yscale("log")
    axes[0][0].set_ylabel("MFE above exact (kcal/mol)")
    axes[1][0].set_ylabel("wall time / derna2 exact")
    save(fig, out)


if __name__ == "__main__":
    main()
