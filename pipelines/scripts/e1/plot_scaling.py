"""E1 scaling: wall-clock time and peak memory against protein length, exact serial methods.

    python plot_scaling.py CELLS_TSV OUT_PDF

Columns are the two corners (MFE and CAI), rows wall time and peak RSS, both log-log.
Lines follow the P0DTC2 prefix ladder; the UniProt proteins are the unconnected markers.
Five series, so the legend carries identity together with a distinct marker per method.
"""

import sys

import pandas as pd

from derna2_analysis.figures import METHODS, save, style

ARMS = ["derna-exact-30", "lineardesign", "lcds-exact", "mrnafold", "cdsfold"]
CORNERS = [("mfe", "MFE corner"), ("cai", "CAI corner")]
ROWS = [("wall_s", "wall-clock time (s)", 1.0), ("maxrss_kb", "peak RSS (GB)", 1 / 2**20)]


def main():
    cells_tsv, out = sys.argv[1:]
    df = pd.read_csv(cells_tsv, sep="\t")
    df = df[(df.variant == "-") & df.arm.isin(ARMS) & (df.time_exit == "ok")]
    style()
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4), sharex=True)
    for col, (corner, title) in enumerate(CORNERS):
        for row, (field, ylabel, scale) in enumerate(ROWS):
            ax = axes[row][col]
            for arm in ARMS:
                g = df[(df.arm == arm) & (df.corner == corner)].dropna(subset=[field])
                if g.empty:
                    continue
                label, color, marker = METHODS[arm]
                lad = g[g.ladder].sort_values("n_aa")
                ax.plot(
                    lad.n_aa,
                    lad[field] * scale,
                    color=color,
                    marker=marker,
                    markersize=3,
                    label=label,
                )
                other = g[~g.ladder]
                ax.plot(
                    other.n_aa,
                    other[field] * scale,
                    linestyle="none",
                    marker=marker,
                    markersize=3,
                    markerfacecolor="none",
                    markeredgecolor=color,
                    alpha=0.8,
                )
            ax.set_xscale("log")
            ax.set_yscale("log")
            if row == 0:
                ax.set_title(title)
            if row == 1:
                ax.set_xlabel("protein length (aa)")
            if col == 0:
                ax.set_ylabel(ylabel)
    axes[0][0].legend(loc="upper left", fontsize=7)
    save(fig, out)


if __name__ == "__main__":
    main()
