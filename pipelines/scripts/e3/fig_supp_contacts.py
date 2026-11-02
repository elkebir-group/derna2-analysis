"""Supplement figure: where the CDS-UTR contacts of the concatenation designs fall, per UTR pair.

    python fig_supp_contacts.py CONTACTS.json OUT.pdf

CONTACTS.json is contacts_data.py's output. One row per UTR pair, three panels per row: nt before
the start codon, percent along the CDS, nt after the stop codon. The two UTR panels keep one nt
scale in every row, so a short UTR is narrow. Bars are the share of the pair's non-modular
concatenation designs with a contact endpoint in the bin, coloured by the UTR of the contact.
"""

import importlib.util
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mp  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

PAIRS = ["hbb", "mrna1273", "hba", "curevac", "bnt"]
DISPLAY = {"hbb": "HBB", "mrna1273": "mRNA-1273", "hba": "HBA", "curevac": "CureVac", "bnt": "BNT"}
UTR5, UTR3 = "#56b4e9", "#9a9a9a"   # the two UTRs' colours in the E3 plate
INK = "#222222"


def style():
    """The paper's figure style, which the E1 overview defines: Times through usetex."""
    path = Path(__file__).resolve().parents[1] / "e1" / "fig4_overview.py"
    spec = importlib.util.spec_from_file_location("fig4_overview", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.set_style()
    return mod.FS, mod.FSS


def main():
    src, out = sys.argv[1:3]
    fs, fss = style()
    data = json.loads(Path(src).read_text())
    W = 7.3
    rh, rg, top, bot = 0.74, 0.33, 0.20, 0.46
    H = top + 5 * rh + 4 * rg + bot
    fig = plt.figure(figsize=(W, H))
    sc5, sc3 = 0.0105, 0.0063                   # inches per nt, the same in every row
    xc, wc = 1.52, 3.4                          # the CDS panel's left edge and width
    gap, lab_gap = 0.10, 0.30                   # between panels; y label to its axes
    x_lab = []
    for r, pair in enumerate(PAIRS):
        d = data[pair]
        l5, l3 = d["utr5_nt"], d["utr3_nt"]
        y0 = H - top - (r + 1) * rh - r * rg
        axs = [fig.add_axes([(xc - gap - l5 * sc5) / W, y0 / H, l5 * sc5 / W, rh / H]),
               fig.add_axes([xc / W, y0 / H, wc / W, rh / H]),
               fig.add_axes([(xc + wc + gap) / W, y0 / H, l3 * sc3 / W, rh / H])]
        n = max(d["designs"], 1)
        tops = []
        for ax, edges, parts in ((axs[0], d["utr5_edges"], [(d["utr5"], UTR5)]),
                                 (axs[1], d["cds_edges"], [(d["cds_with_5"], UTR5), (d["cds_with_3"], UTR3)]),
                                 (axs[2], d["utr3_edges"], [(d["utr3"], UTR3)])):
            edges = np.array(edges)
            bottom = np.zeros(len(edges) - 1)
            for hist, col in parts:
                share = 100 * np.array(hist) / n
                ax.bar(edges[:-1], share, bottom=bottom, width=np.diff(edges), align="edge", color=col, lw=0,
                       zorder=3)
                bottom = bottom + share
            tops.append(bottom.max())
        for ax in axs:
            ax.set_ylim(0, max(tops) * 1.08 if max(tops) > 0 else 1)
            ax.tick_params(labelbottom=True)
        axs[0].set_xlim(-l5, 0)
        axs[1].set_xlim(0, 100)
        axs[2].set_xlim(0, l3)
        for ax in axs[1:]:
            ax.spines["left"].set_visible(False)
            ax.tick_params(axis="y", labelleft=False, length=0)
        axs[0].set_xticks([-60, -30, 0] if l5 >= 60 else [-30, 0])
        axs[1].set_xticks([0, 25, 50, 75, 100])
        step = 50 if l3 <= 150 else 100
        axs[2].set_xticks(list(range(step, l3 + 1, step)))
        axs[0].set_ylabel("designs (\\%)")
        axs[0].yaxis.set_label_coords(-lab_gap / (l5 * sc5), 0.5)
        x_lab.append(xc - gap - l5 * sc5 - lab_gap - 0.07)
        if r == 0:
            axs[1].legend(handles=[mp.Patch(fc=UTR5, lw=0, label="5$'$ UTR contact"),
                                   mp.Patch(fc=UTR3, lw=0, label="3$'$ UTR contact")],
                          loc="upper center", ncol=2, frameon=False, fontsize=fss, handlelength=0.9,
                          columnspacing=1.0, handletextpad=0.4, borderaxespad=0.1)
        if r == len(PAIRS) - 1:
            axs[0].set_xlabel("nt before start codon")
            axs[1].set_xlabel("position in CDS (\\%)")
            axs[2].set_xlabel("nt after stop codon")
    # panel titles left-aligned in one column, at the leftmost y label
    for r, pair in enumerate(PAIRS):
        y0 = H - top - (r + 1) * rh - r * rg
        fig.text(min(x_lab) / W, (y0 + rh + 0.03) / H, f"\\textbf{{({'abcde'[r]})}} {DISPLAY[pair]}", ha="left",
                 va="bottom", fontsize=fs)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, metadata={"CreationDate": None})


if __name__ == "__main__":
    main()
