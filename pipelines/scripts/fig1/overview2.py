"""The revised overview figure: (a) the problem as input -> method -> output, (b) the fronts with a
table that is their legend, (c)-(e) the three designs' structures.

    python overview2.py RESULTS_DIR MRNAFOLD_JSONL DESIGNS_JSON FASTA UTR5_TXT UTR3_TXT PDB \
        OUT_PDF [ochre|nc|plddt]

Drawing helpers (fronts, zooms, RNApuzzler structure layouts) come from overview_draw.py beside
this file. PDB is the AlphaFold model of the protein; its C-alpha trace is drawn as a ribbon,
thicker where the backbone is helical.
"""

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

spec = importlib.util.spec_from_file_location("ov", Path(__file__).with_name("overview_draw.py"))
ov = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ov)
plt, COL, FS, FSS, tt = ov.plt, ov.COL, ov.FS, ov.FSS, ov.tt
GREEN, ORANGE, PURPLE = ov.GREEN, ov.ORANGE, ov.PURPLE

TOP = 0.735  # top of the (b)-(e) row, figure fraction; (a) sits above it


# ---- (a): input -> method -> output ------------------------------------------------------


def ca_trace(pdb):
    ca = [row for row in open(pdb) if row.startswith("ATOM") and row[12:16].strip() == "CA"]
    xyz = np.array([[float(row[30:38]), float(row[38:46]), float(row[46:54])] for row in ca])
    # an AlphaFold model stores its per-residue confidence (pLDDT) in the B-factor column
    plddt = np.array([float(row[60:66]) for row in ca])
    # helical where residue i and i+3 sit about one turn apart
    d3 = np.linalg.norm(xyz[3:] - xyz[:-3], axis=1)
    helix = np.zeros(len(xyz), bool)
    for i, d in enumerate(d3):
        if d < 6.0:
            helix[i : i + 4] = True
    return xyz, helix, plddt


def spline(p, n=8):
    # Catmull-Rom through the C-alpha positions, n samples per segment
    p = np.vstack([p[0], p, p[-1]])
    out = []
    for i in range(1, len(p) - 2):
        p0, p1, p2, p3 = p[i - 1], p[i], p[i + 1], p[i + 2]
        for t in np.linspace(0, 1, n, endpoint=False):
            out.append(
                0.5
                * (
                    (2 * p1)
                    + (-p0 + p2) * t
                    + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t**2
                    + (-p0 + 3 * p1 - 3 * p2 + p3) * t**3
                )
            )
    out.append(p[-2])
    return np.array(out)


PLDDT_BANDS = [(90, "#0053d6"), (70, "#65cbf3"), (50, "#ffdb13"), (0, "#ff7d45")]  # AlphaFold DB


def residue_colors(mode, plddt):
    # one RGB per residue: a single hue, N- to C-terminus, or AlphaFold's confidence bands
    n = len(plddt)
    if mode == "nc":
        cmap = plt.get_cmap("turbo")
        return [np.array(cmap(0.08 + 0.84 * i / (n - 1))[:3]) for i in range(n)]
    if mode == "plddt":
        hexes = [next(c for lo, c in PLDDT_BANDS if p > lo) for p in plddt]
        return [np.array([int(h[k : k + 2], 16) / 255 for k in (1, 3, 5)]) for h in hexes]
    return [np.array([0.722, 0.471, 0.098])] * n


def protein_cartoon(ax, pdb, mode="ochre"):
    xyz, helix, plddt = ca_trace(pdb)
    base = residue_colors(mode, plddt)
    xyz = xyz - xyz.mean(axis=0)
    # view down the axis of least spread, longest extent horizontal
    _, _, vt = np.linalg.svd(xyz, full_matrices=False)
    xyz = xyz @ vt.T
    n = 8
    s = spline(xyz, n)
    hel = np.repeat(helix, n)[: len(s) - 1]
    z = s[:, 2]
    zn = (z - z.min()) / (z.max() - z.min())
    segs = sorted(range(len(s) - 1), key=lambda k: z[k])
    for k in segs:
        shade = 0.55 + 0.45 * zn[k]  # nearer is darker
        c = base[min(k // n, len(base) - 1)] * shade + (1 - shade)
        lw = 3.2 if hel[k] else 1.1
        ax.plot(s[k : k + 2, 0], s[k : k + 2, 1], color=c, lw=lw, solid_capstyle="round", zorder=k)
    ax.set_aspect("equal")
    ax.axis("off")


def forbid(fig, x, y, color, half_in=0.035):
    # a cross centred on (x, y): this is not allowed. half_in is half its width in inches.
    w, h = fig.get_size_inches()
    dx, dy = half_in / w, half_in / h
    for sx in (-1, 1):
        fig.add_artist(
            Line2D(
                [x - dx, x + dx],
                [y - sx * dy, y + sx * dy],
                color=color,
                lw=1.0,
                solid_capstyle="round",
                transform=fig.transFigure,
                zorder=5,
            )
        )


def cartoon_key(fig, mode, x, y):
    # what the ribbon's colors mean, under the cartoon
    if mode == "nc":
        cmap = plt.get_cmap("turbo")
        for i in range(40):
            fig.add_artist(
                Rectangle(
                    (x + i * 0.0015, y),
                    0.0015,
                    0.012,
                    transform=fig.transFigure,
                    fc=cmap(0.08 + 0.84 * i / 39),
                    ec="none",
                )
            )
        fig.text(x - 0.004, y + 0.006, "N", ha="right", va="center", fontsize=5.6)
        fig.text(x + 0.064, y + 0.006, "C", ha="left", va="center", fontsize=5.6)
    elif mode == "plddt":
        labels = ["$>$90", "70--90", "50--70", "$<$50"]
        fig.text(x - 0.032, y + 0.003, "pLDDT", ha="left", va="center", fontsize=5.6)
        for k, ((_, c), lab) in enumerate(zip(PLDDT_BANDS, labels, strict=True)):
            xx, yy = x + 0.012 + (k % 2) * 0.038, y + 0.012 - (k // 2) * 0.02
            fig.add_artist(
                Rectangle(
                    (xx, yy - 0.006), 0.008, 0.012, transform=fig.transFigure, fc=c, ec="none"
                )
            )
            fig.text(xx + 0.011, yy, lab, va="center", fontsize=5.0)


def schematic(fig, w, u5, u3, motifs, pdb, mode="ochre"):
    # protein cartoon, then the fixed inputs as text, an arrow, and the designed mRNA
    cart = fig.add_axes([0.03, TOP + 0.075, 0.12, 0.15])
    protein_cartoon(cart, pdb, mode)
    yc = TOP + 0.048  # the baseline the four captions of (a) share
    fig.text(0.09, yc, r"protein $\mathbf{w}$ (78 aa)", ha="center", va="center", fontsize=FSS)
    cartoon_key(fig, mode, 0.06, TOP + 0.025)
    x, y, dy = 0.17, 0.92, 0.047
    rows = [
        (COL["u5"], r"5$'$ UTR $\mathbf{u}_5$", tt(u5[:20] + r"\ldots{}") + f" ({len(u5)} nt)"),
        (COL["u3"], r"3$'$ UTR $\mathbf{u}_3$", tt(u3[:20] + r"\ldots{}") + f" ({len(u3)} nt)"),
        (
            COL["m"],
            r"motifs $\mathcal{M}$",
            r"$\{$"
            + ", ".join(tt(m) + (f" ({ov.NAMES[m]})" if m in ov.NAMES else "") for m in motifs)
            + r"$\}$",
        ),
    ]
    for c, name, val in rows:
        fig.add_artist(
            Line2D(
                [x, x + 0.018],
                [y, y],
                color=c,
                lw=2,
                solid_capstyle="butt",
                transform=fig.transFigure,
            )
        )
        fig.text(x + 0.024, y, name, va="center", fontsize=FSS)
        fig.text(x + 0.115, y, val, va="center", fontsize=FSS)
        y -= dy
    fig.text(
        x,
        yc,
        r"fixed: $\mathbf{u}_5$, $\mathbf{u}_3$, $\mathcal{M}$;"
        r" designed: the CDS $\mathbf{v}$ encoding $\mathbf{w}$",
        va="center",
        fontsize=FSS,
        color="#555555",
    )
    # the method
    ya = 0.865
    fig.add_artist(
        FancyArrowPatch(
            (0.52, ya),
            (0.615, ya),
            transform=fig.transFigure,
            arrowstyle="-|>",
            mutation_scale=9,
            lw=0.9,
            color="k",
        )
    )
    fig.text(0.5675, ya + 0.02, r"\textbf{DERNA2}", ha="center", va="bottom", fontsize=FSS)
    fig.text(
        0.5675,
        yc + 0.015,
        "Pareto front of\nmodular designs",
        ha="center",
        va="center",
        fontsize=FSS,
        linespacing=1.0,
    )
    # the output: one mRNA as three boxes, with what it may not contain drawn on it
    x0, yb, h = 0.63, ya - 0.022, 0.044
    edges = []
    for c, lab, wd in (
        (COL["u5"], r"5$'$ UTR", 0.07),
        (COL["cds"], r"CDS $\mathbf{v}$", 0.19),
        (COL["u3"], r"3$'$ UTR", 0.10),
    ):
        fig.add_artist(
            FancyBboxPatch(
                (x0, yb),
                wd,
                h,
                boxstyle="round,pad=0,rounding_size=0.004",
                transform=fig.transFigure,
                fc=c,
                ec="white",
                lw=0.8,
            )
        )
        edges.append((x0, x0 + wd))
        # the CDS label sits left of the motif chip drawn below
        lx = x0 + (0.06 if c == COL["cds"] else wd / 2)
        fig.text(lx, yb + h / 2, lab, ha="center", va="center", fontsize=FSS, color="white")
        x0 += wd
    (_, _), (c0, c1), (_, _) = edges
    top = yb + h
    # no CDS-UTR base pair: an arc from each UTR into the CDS, both of one span and height,
    # each crossed out at its apex
    span, rise = 0.07, 0.042
    for a in (edges[0][0] + 0.035, c1 - 0.035):
        u = np.linspace(0, np.pi, 40)
        fig.add_artist(
            Line2D(
                a + span * (1 - np.cos(u)) / 2,
                top + 0.003 + rise * np.sin(u),
                color=COL["x"],
                lw=0.8,
                transform=fig.transFigure,
            )
        )
        forbid(fig, a + span / 2, top + 0.0025 + rise, COL["x"])
    # no motif of M inside the CDS: a motif block, crossed out where it sits
    mx = c0 + 0.13
    fig.add_artist(
        Rectangle(
            (mx, yb + 0.006),
            0.03,
            h - 0.012,
            transform=fig.transFigure,
            fc=COL["m"],
            ec="white",
            lw=0.5,
        )
    )
    forbid(fig, mx + 0.015, yb + h / 2, "white")
    fig.text(
        0.63,
        yc,
        r"\textbf{modular:} no motif of $\mathcal{M}$ in $\mathbf{v}$,"
        r" no CDS--UTR base pair",
        va="center",
        fontsize=FSS,
        linespacing=1.1,
    )


# ---- (b): the fronts, and the table that is their legend --------------------------------


def mark(c, m, filled=True):
    return Line2D([], [], marker=m, ls="", ms=3.0, mew=0.5, mec=c, mfc=c if filled else "white")


def table_legend(fig, left, right, bottom, top):
    # 2x2: rows modular / non-modular, columns complete / incomplete; each name with its marker.
    # A closed grid in light grey, so it frames the legend without outweighing the plot.
    ax = fig.add_axes([left, bottom, right - left, top - bottom])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    grid = "#9a9a9a"
    c1, c2 = 0.19, 0.405  # column edges: row label | complete | incomplete
    h1, h2 = 2 / 3, 1 / 3  # row edges: header | modular | non-modular
    ax.add_patch(Rectangle((c1, h1), 1 - c1, 1 - h1, fc="#efefef", ec="none"))
    ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, ec=grid, lw=0.5))
    for yy in (h1, h2):
        ax.plot([0, 1], [yy, yy], color=grid, lw=0.4)
    for xx in (c1, c2):
        ax.plot([xx, xx], [0, 1], color=grid, lw=0.4)
    pad = 0.012
    ax.text(pad, (1 + h1) / 2, r"\textbf{Pareto front}", va="center", fontsize=FSS)
    ax.text(c1 + pad, (1 + h1) / 2, "complete", va="center", fontsize=FSS)
    ax.text(c2 + pad, (1 + h1) / 2, "incomplete", va="center", fontsize=FSS)
    ax.text(pad, (h1 + h2) / 2, "modular", va="center", fontsize=FSS)
    ax.text(pad, h2 / 2, "non-modular", va="center", fontsize=FSS)

    def entry(x, y, color, m, text):
        # filled and open marker at fixed offsets, then the name at one left edge per column
        for k, filled in enumerate((True, False)):
            ax.add_line(
                Line2D(
                    [x + 0.02 + 0.03 * k],
                    [y],
                    marker=m,
                    ls="",
                    ms=3.0,
                    mew=0.5,
                    mec=color,
                    mfc=color if filled else "white",
                )
            )
        ax.text(x + 0.07, y, text, va="center", fontsize=FSS)

    entry(c1, (h1 + h2) / 2, PURPLE, "s", "DERNA2")
    entry(c2, (h1 + h2) / 2, ORANGE, "o", r"mRNAfold$^{*}$")
    ax.text((c1 + c2) / 2, h2 / 2, "--", ha="center", va="center", fontsize=FSS, color="#888888")
    entry(c2, h2 / 2, GREEN, "o", "LinearDesign, LinearCDSfold, DERNA")
    return ax


def panel_b(fig, rect, fronts, designs):
    ax = fig.add_axes(rect)
    ax.set_xlim(0.85, 1.003)
    shown = [r["free_mfe_kcal_mol"] for r in fronts.b + fronts.mf + fronts.cd if ov.cai(r) >= 0.85]
    shown += [r["mfe_kcal_mol"] for r in fronts.so if ov.cai(r) >= 0.85]
    ax.set_ylim(*designs.get("ylim", (min(shown) - 6.5, max(shown) + 3)))
    ov.plot_fronts(ax, fronts, 2.2, 1.8)
    xs = [ov.cai(designs[k]) for k in ("d1", "d2", "d3")]
    ys = [designs[k]["free_mfe_kcal_mol"] for k in ("d1", "d2", "d3")]
    ins = ov.zoom_inset(
        ax,
        fronts,
        [0.03, 0.56, 0.42, 0.42],
        min(xs) - 0.007,
        max(xs) + 0.0125,
        min(ys) - 4.0,
        max(ys) + 6.0,
        [((0, 1), (0, 0)), ((1, 1), (1, 0))],
    )
    pts = sorted((r["cai_log"], r["mfe_kcal_mol"]) for r in fronts.so)
    lo, hi = max(ov.pairwise(fronts.hull), key=lambda e: ov.gap(e, pts))
    ov.zoom_inset(
        ax,
        fronts,
        [0.62, 0.03, 0.37, 0.26],
        math.exp(lo[0] / 78) - 0.0025,
        math.exp(hi[0] / 78) + designs.get("inset2_pad", 0.0025),
        lo[1] - 6.5,
        hi[1] + 3.0,
        [((1, 0), (1, 1))],
    )
    ov.mark_designs(ins, designs, designs.get("label_at", {}))
    ax.set_xticks([0.85, 0.9, 0.95, 1.0])
    ax.set_xticklabels(["$0.85$", "$0.9$", "$0.95$", "$1$"])
    ax.set_xlabel("CAI", labelpad=1)
    ax.set_ylabel("MFE (kcal/mol)", labelpad=1)
    ax.tick_params(pad=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    return ax


# ---- (c)-(e): the structures ---------------------------------------------------------------


def emptiest(ax, xs, ys, candidates):
    # of the candidate inset bounds (axes fractions), the one covering the fewest bases
    ax.apply_aspect()
    pts = (ax.transData + ax.transAxes.inverted()).transform(np.column_stack([xs, ys]))

    def covered(b):
        x0, y0, w, h = b
        return sum(
            1
            for px, py in pts
            if x0 - 0.02 <= px <= x0 + w + 0.02 and y0 - 0.02 <= py <= y0 + h + 0.02
        )

    return min(candidates, key=covered)


def structures(fig, designs, left, right, bottom, top):
    gr = fig.add_gridspec(
        1,
        3,
        left=left,
        right=right,
        bottom=bottom,
        top=top,
        wspace=0.04,
        width_ratios=[1.0, 1.0, 1.0],
    )
    axs = [fig.add_subplot(gr[0, k]) for k in range(3)]
    d1, d2 = designs["d1"], designs["d2"]
    helix = ov.cds_utr_helix(d1["free_structure"])
    xs, ys = ov.draw_structure(
        fig, axs[0], d1["free_structure"], vroom=designs.get("vroom_c", 0.15)
    )
    axs[0].set_anchor("C")
    spot = designs.get("spot_c") or emptiest(
        axs[0], xs, ys, [[0.24, 0.66, 0.76, 0.34], [0.0, 0.66, 0.76, 0.34]]
    )
    ov.letter_zoom(
        axs[0],
        d1,
        xs,
        ys,
        helix,
        COL["x"],
        spot,
        margin=designs.get("margin_c", 12),
        cds_name=r"$\mathbf{v}''$",
    )
    m0, mlen, spelled, mot = ov.motif_sites(d2["sequence"], designs["motifs"])
    xs, ys = ov.draw_structure(
        fig, axs[1], d2["free_structure"], mot, vroom=0.6, align=(m0, m0 + mlen - 1)
    )
    axs[1].set_anchor("C")
    # a band under the molecule for the CDS-UTR zoom, so that zoom covers none of it
    y0, y1 = axs[1].get_ylim()
    axs[1].set_ylim(y0 - 0.3 * (y1 - y0), y1)
    ov.letter_zoom(
        axs[1],
        d2,
        xs,
        ys,
        set(range(m0, m0 + mlen)),
        COL["m"],
        [0.02, 0.76, 0.9, 0.24],
        margin=18,
        name=ov.NAMES.get(spelled, spelled),
    )
    spot = designs.get("spot_d") or emptiest(
        axs[1], xs, ys, [[0.0, 0.0, 0.36, 0.24], [0.64, 0.0, 0.36, 0.24]]
    )
    ov.letter_zoom(
        axs[1],
        d2,
        xs,
        ys,
        ov.cds_utr_helix(d2["free_structure"]),
        COL["x"],
        spot,
        margin=13,
        cds_name=r"$\mathbf{v}'$",
        corners=[((0, 0), (0, 1)), ((1, 0), (1, 1))],
    )
    ov.draw_structure(fig, axs[2], designs["d3"]["free_structure"])
    return gr, axs


def color_key(fig, left, right, y):
    items = [
        (COL["u5"], r"5$'$ UTR"),
        (COL["cds"], "CDS"),
        (COL["u3"], r"3$'$ UTR"),
        (COL["m"], r"motif of $\mathcal{M}$"),
        (COL["x"], "CDS--UTR base pair"),
    ]
    texts = [fig.text(0, y, lab, va="center", fontsize=FSS) for _, lab in items]
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    boxes = [inv.transform(tx.get_window_extent()) for tx in texts]
    widths = [b[1, 0] - b[0, 0] for b in boxes]
    sw, pad = 0.018, 0.005
    gap = (right - left - sum(w + sw + pad for w in widths)) / (len(items) - 1)
    x = left
    for (c, _), tx, w in zip(items, texts, widths, strict=True):
        fig.add_artist(
            Line2D(
                [x, x + sw], [y, y], color=c, lw=2, solid_capstyle="butt", transform=fig.transFigure
            )
        )
        tx.set_x(x + sw + pad)
        x += sw + pad + w + gap


def titles(fig, designs, xs_cde, fronts):
    d1, d2, d3 = designs["d1"], designs["d2"], designs["d3"]

    def stats(d):
        return f"CAI {ov.cai(d):.3f}, {ov.mfe_text(d['free_mfe_kcal_mol'])} kcal/mol"

    def found(d):
        return ", ".join(
            ov.NAMES.get(m, tt(m)) + (" site" if m in ov.NAMES else "") for m in d["motifs_found"]
        )

    fig.text(
        0.005,
        0.985,
        r"\textbf{(a)} glycophorin-E (GLPE\_HUMAN, UniProt"
        r" P15421) with " + designs.get("utr_label", r"the human $\alpha$-globin (HBA1) UTRs"),
        va="top",
        fontsize=FS,
    )
    fig.text(0.005, TOP - 0.005, r"\textbf{(b)}", va="top", fontsize=FS)
    heads = [
        (
            "(c)",
            GREEN,
            "o",
            r"$\mathbf{v}''$: UTR-unaware",
            d1,
            f"{ov.cds_utr_pairs(d1['free_structure'])} CDS--UTR pairs"
            f", {len(ov.cds_utr_helix(d1['free_structure'])) // 2} shown",
        ),
        (
            "(d)",
            ORANGE,
            "o",
            r"$\mathbf{v}'$: mRNAfold",
            d2,
            f"{ov.cds_utr_pairs(d2['free_structure'])} CDS--UTR pairs, {found(d2)}",
        ),
        ("(e)", PURPLE, "s", r"$\mathbf{v}$: DERNA2", d3, "no CDS--UTR pair, no XhoI"),
    ]
    # filled when the design is a hull vertex of its method, as in (b)
    filled = [
        fronts.on_b(d1),
        fronts.on_mf(d2),
        (d3["cai_log"], d3["mfe_kcal_mol"]) in set(fronts.hull),
    ]
    inv = fig.transFigure.inverted()
    for x, (letter, col, m, name, d, extra), full in zip(xs_cde, heads, filled, strict=True):
        y = TOP - 0.005
        lab = fig.text(x, y, r"\textbf{" + letter + "}", va="top", fontsize=FS)
        fig.canvas.draw()
        bb = inv.transform(lab.get_window_extent())
        xm, ym = bb[1, 0] + 0.009, (bb[0, 1] + bb[1, 1]) / 2
        fig.add_artist(
            Line2D(
                [xm],
                [ym],
                marker=m,
                ms=3.2,
                mew=0.6,
                mec=col,
                mfc=col if full else "white",
                ls="",
                transform=fig.transFigure,
            )
        )
        fig.text(xm + 0.008, y, name, va="top", fontsize=FS)
        fig.text(x, y - 0.048, stats(d) + "\n" + extra, va="top", fontsize=5.9, linespacing=1.15)


def main():
    results, mrnafold, designs_json, fasta, utr5, utr3, pdb, out = sys.argv[1:9]
    mode = sys.argv[9] if len(sys.argv) > 9 else "ochre"  # ochre, nc or plddt
    fronts = ov.load_fronts(results, mrnafold)
    designs = json.loads(Path(designs_json).read_text())
    ov.L5 = len(Path(utr5).read_text().strip())
    ov.L3 = len(Path(utr3).read_text().strip())
    w, u5, u3 = ov.load_inputs(fasta, utr5, utr3, designs)
    ov.set_style()
    fig = plt.figure(figsize=(6.5, 2.95))
    schematic(fig, w, u5, u3, designs["motifs"], pdb, mode)
    fig.add_artist(
        Line2D(
            [0.0, 1.0],
            [TOP + 0.015, TOP + 0.015],
            color="#cccccc",
            lw=0.5,
            transform=fig.transFigure,
        )
    )
    # (b): table on top, plot below
    table_legend(fig, 0.07, 0.49, TOP - 0.165, TOP - 0.035)
    fig.text(
        0.49,
        TOP - 0.172,
        r"filled: hull vertex; open: off the hull."
        r" $^{*}$modularity not guaranteed.",
        ha="right",
        va="top",
        fontsize=5.4,
        color="#444444",
    )
    panel_b(fig, [0.07, 0.09, 0.413, TOP - 0.305], fronts, designs)
    # (c)-(e)
    left, right = 0.505, 0.975
    gr, _ = structures(fig, designs, left, right, 0.07, TOP - 0.115)
    xs_cde = [gr[0, k].get_position(fig).x0 for k in range(3)]
    titles(fig, designs, xs_cde, fronts)
    color_key(fig, left, right, 0.03)
    fig.savefig(out)
    fig.savefig(out.rsplit(".", 1)[0] + ".png", dpi=300)


if __name__ == "__main__":
    main()
