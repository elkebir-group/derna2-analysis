"""Shared figure style: one color, marker and label per method, fixed across every figure.

The categorical colors are assigned in a fixed order and follow the method, never its
rank in a panel, so a method reads the same in every figure. Three of them sit under
3:1 contrast on white, so every series also carries a distinct marker and a direct
label: identity is never carried by color alone.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"

# Categorical slots 1-8, light mode, in their validated order.
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
# One-hue ordinal ramp (blue 250 / 450 / 650) for beam widths, narrow to wide.
WIDTH_RAMP = ["#86b6ef", "#2a78d6", "#104281"]

# arm -> (label, color, marker). Exact serial methods take slots 1-4 in this order.
METHODS = {
    "derna-exact-30": ("derna2", SLOTS[0], "o"),
    "lineardesign": ("LinearDesign", SLOTS[1], "s"),
    "lcds-exact": ("LinearCDSfold", SLOTS[2], "^"),
    "mrnafold": ("mRNAfold", SLOTS[3], "D"),
    "derna-exact-20": ("derna2, loop cap 20", SLOTS[5], "P"),
    "derna-threaded": ("derna2, 8 threads", SLOTS[6], "X"),
    "mrnafold-threaded": ("mRNAfold, 8 threads", SLOTS[7], "h"),
}


def style():
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.edgecolor": INK_2,
            "axes.labelcolor": INK,
            "axes.titlesize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "xtick.color": INK_2,
            "ytick.color": INK_2,
            "lines.linewidth": 1.5,
            "lines.markersize": 4,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "savefig.bbox": "tight",
        }
    )


def save(fig, out_pdf):
    """Write the figure as the PDF the rule names and a PNG beside it."""
    out_pdf = Path(out_pdf)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.with_suffix(".png"), dpi=200)
    plt.close(fig)


def label_line_end(ax, x, y, text):
    """Direct label at a series' last point, in text ink rather than the series color."""
    ax.annotate(
        text, (x, y), xytext=(4, 0), textcoords="offset points", va="center", fontsize=7, color=INK
    )
