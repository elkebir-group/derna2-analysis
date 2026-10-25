"""The overview figure: GLPE_HUMAN (78 aa) with the hbb UTRs and one restriction site forbidden.

    python overview.py RESULTS_DIR MRNAFOLD_JSONL DESIGNS_JSON FASTA UTR5_TXT UTR3_TXT OUT_PDF

(a) the modular front (DERNA2) against mRNAfold's designs with CDS-UTR pairs discouraged (one per
    hull weight of the UTR-masking front, from the fig1_fronts deposit) and the hull vertices of
    the concatenation front, the designs a weighted sum can return (a proxy for LinearDesign,
    LinearCDSfold and DERNA), each at the MFE of the refolded mRNA; a zoom on the three designs
    of (b)-(d) and one on the hull edge of the modular front with the most designs above it.
    Color: method; shape: modular (square) or not (circle); fill: on the convex hull of the
    modular front or above it.
(b)-(d) the structures of the three designs of overview_designs.py; (c) zooms on the motif.
Text is set by LaTeX in Times at the sizes of the paper's Fig. 2 (7 pt, small labels 6.3 pt), at
the paper's text width (6.5 in). Structure layouts come from ViennaRNA (NAView).
"""

import json
import math
import sys
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

import matplotlib
import RNA
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Rectangle

from derna2_analysis.io import read_jsonl

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

L5, LC = 50, 234   # hbb 5' UTR and the CDS of GLPE_HUMAN (78 aa)
# (b)-(d) use none of the method colors of (a)
COL = {"u5": "#56b4e9", "cds": "#1f3a93", "u3": "#9a9a9a", "x": "#d62728", "m": "#c51b8a"}
GREEN, ORANGE, PURPLE = "#2e9e5b", "#e8772e", "#6a4bc4"
FS, FSS = 7, 6.3
PREAMBLE = (r"\usepackage[T1]{fontenc}\usepackage{times}\usepackage{textcomp}"
            r"\usepackage{amsmath}")


@dataclass
class Fronts:
    """The three fronts of (a) and which of their designs lie on their methods' hulls."""

    b: list            # concatenation designs, refolded with the UTRs
    mf: list           # mRNAfold's designs, CDS-UTR pairs discouraged
    cd: list           # the UTR-masking front, which mRNAfold's weights sample
    so: list           # the modular front
    hull: list         # convex hull of the modular front in (log CAI, MFE)
    on_hull_b: set     # CDS sequences on the hull of the front without UTRs

    def on_b(self, r):
        return r["cds_sequence"] in self.on_hull_b

    def on_mf(self, r):
        return True


def cai(r):
    return math.exp(r["cai_log"] / 78)


def by_cai(rows):
    return sorted(rows, key=cai)


def region(i):
    if i < L5:
        return "u5"
    return "cds" if i < L5 + LC else "u3"


def pairs(s):
    stack, out = [], []
    for i, c in enumerate(s):
        if c == "(":
            stack.append(i)
        elif c == ")":
            out.append((stack.pop(), i))
    return out


def cds_utr_pairs(s):
    return sum((region(i) == "cds") != (region(j) == "cds") for i, j in pairs(s))


def mfe_text(v):
    return r"\textminus" + f"{-v:.1f}"


def tt(x):
    return r"\texttt{" + x + "}"


def lower_hull(pts):
    h = []
    for q in sorted(pts):
        while len(h) >= 2 and ((h[-1][0] - h[-2][0]) * (q[1] - h[-2][1])
                               - (h[-1][1] - h[-2][1]) * (q[0] - h[-2][0])) <= 0:
            h.pop()
        h.append(q)
    return h


def hull_curve(h):
    # the hull edges are straight in log CAI, so each is sampled before mapping to CAI
    xs, ys = [], []
    for (x0, y0), (x1, y1) in pairwise(h):
        for t in range(21):
            xs.append(math.exp((x0 + t / 20 * (x1 - x0)) / 78))
            ys.append(y0 + t / 20 * (y1 - y0))
    return xs, ys


def gap(edge, pts):
    # the largest height of a design above the hull edge
    (ex0, ey0), (ex1, ey1) = edge
    heights = [q[1] - (ey0 + (q[0] - ex0) / (ex1 - ex0) * (ey1 - ey0))
               for q in pts if ex0 < q[0] < ex1]
    return max(heights or [0])


def load_mrnafold(path):
    # the guarded arm, with its CAI under our codon weights
    return [dict(r, cai_log=r["theirs_cai_log"]) for r in read_jsonl(path) if r["arm"] == "guard"]


def load_fronts(results, mrnafold):
    res = Path(results)
    b = read_jsonl(res / "b_cdsonly_with_utrs.scored.jsonl")
    mf = load_mrnafold(mrnafold)
    so = read_jsonl(res / "d_so_front.jsonl")
    assert not any(r["symmetric_so"] for r in b + mf)   # no teal or orange point is modular
    on_hull_b = {r["sequence"] for r in read_jsonl(res / "a_cdsonly_front.jsonl")
                 if r.get("hull_vertex")}
    return Fronts(
        b=b, mf=mf, so=so, cd=read_jsonl(res / "cd_symmetric_guarded_front.scored.jsonl"),
        hull=lower_hull([(r["cai_log"], r["mfe_kcal_mol"]) for r in so]),
        on_hull_b=on_hull_b,
    )


def load_inputs(fasta, utr5, utr3, designs):
    w = "".join(line.strip() for line in Path(fasta).read_text().splitlines()
                if not line.startswith(">"))
    u5, u3 = Path(utr5).read_text().strip(), Path(utr3).read_text().strip()
    assert (len(u5), len(u3), len(w)) == (L5, 132, 78)
    assert all(designs[k]["sequence"] == u5 + designs[k]["sequence"][L5:L5 + LC] + u3
               for k in ("d1", "d2", "d3"))
    return w, u5, u3


def set_style():
    plt.rcParams.update({
        "text.usetex": True, "font.family": "serif", "text.latex.preamble": PREAMBLE,
        "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS,
        "xtick.labelsize": FSS, "ytick.labelsize": FSS, "legend.fontsize": FSS,
        "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    })


def zoom_box(parent, child, x0, x1, y0, y1, color, corners, rad=0.0):
    # the zoomed region on the parent and two dashed lines to the matching corners of the inset
    parent.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=color, lw=0.5,
                               zorder=6))
    for (px, py), (cx, cy) in corners:
        parent.add_artist(ConnectionPatch(xyA=((x0, x1)[px], (y0, y1)[py]),
                                          coordsA=parent.transData, xyB=(cx, cy),
                                          coordsB=child.transAxes, color=color, lw=0.5,
                                          ls=(0, (3, 2)), connectionstyle=f"arc3,rad={rad}",
                                          zorder=6))


def plot_others(a, fronts, ms):
    # not modular: mRNAfold's designs and the hull vertices of the concatenation front, the only
    # designs a weighted sum returns, as circles at the MFE of the refolded mRNA. Most of them
    # are the same design, so the concatenation circles are drawn larger, behind.
    # in the background: the designs of the two fronts that no weight reaches, as faint dots
    mf_seqs = {r["sequence"] for r in fronts.mf}
    hull_b = [r for r in fronts.b if not fronts.on_b(r)]
    hull_cd = [r for r in fronts.cd if r["sequence"] not in mf_seqs]
    for rows, col in ((hull_b, GREEN), (hull_cd, ORANGE)):
        a.plot([cai(r) for r in rows], [r["free_mfe_kcal_mol"] for r in rows], "o",
               ms=0.9 * ms, mec="none", mfc=col, alpha=0.3, zorder=1)
    for rows, on, col, size, z in ((fronts.b, fronts.on_b, GREEN, 1.7 * ms, 2.5),
                                   (fronts.mf, fronts.on_mf, ORANGE, ms, 2.6)):
        sel = [r for r in by_cai(rows) if on(r)]
        a.plot([cai(r) for r in sel], [r["free_mfe_kcal_mol"] for r in sel], "o", ms=size,
               mec=col, mfc=col, mew=0.5, zorder=z)


def plot_modular(a, fronts, ms):
    # the modular front: filled squares on its convex hull, open squares above it
    rows = by_cai(fronts.so)
    vertices = set(fronts.hull)
    xs = [cai(r) for r in rows]
    ys = [r["mfe_kcal_mol"] for r in rows]
    a.plot(xs, ys, color=PURPLE, lw=0.7, zorder=3)
    on = [(r["cai_log"], r["mfe_kcal_mol"]) in vertices for r in rows]
    a.plot([x for x, o in zip(xs, on, strict=True) if o],
           [y for y, o in zip(ys, on, strict=True) if o],
           "s", ms=ms, color=PURPLE, zorder=3.2)
    a.plot([x for x, o in zip(xs, on, strict=True) if not o],
           [y for y, o in zip(ys, on, strict=True) if not o],
           "s", ms=ms, mfc="white", mec=PURPLE, mew=0.5, zorder=3.2)


def plot_fronts(a, fronts, ms_others, ms_modular):
    plot_others(a, fronts, ms_others)
    plot_modular(a, fronts, ms_modular)
    hx, hy = hull_curve(fronts.hull)
    a.plot(hx, hy, "--", color="k", lw=0.45, zorder=2.8)


def zoom_inset(ax, fronts, bounds, x0, x1, y0, y1, corners, rad=0.0):
    ins = ax.inset_axes(bounds)
    plot_fronts(ins, fronts, 2.6, 2.6)
    ins.set_xlim(x0, x1)
    ins.set_ylim(y0, y1)
    ins.set_xticks([])
    ins.set_yticks([])
    for sp in ins.spines.values():
        sp.set_linewidth(0.4)
    zoom_box(ax, ins, x0, x1, y0, y1, "k", corners, rad)
    return ins


def mark_designs(ins, designs):
    # the designs of (b)-(d), circled and labeled in the zoom
    # each label at the upper right of its circle
    marks = (("d1", "free_mfe_kcal_mol", r"$\mathbf{v}''$", (0, 3), "left", "bottom"),
             ("d2", "free_mfe_kcal_mol", r"$\mathbf{v}'$", (2, 2), "left", "bottom"),
             ("d3", "mfe_kcal_mol", r"$\mathbf{v}$", (2, 2), "left", "bottom"))
    for key, field, label, at, ha, va in marks:
        d = designs[key]
        ins.plot(cai(d), d[field], "o", ms=6.0, mfc="none", mec="k", mew=0.7, zorder=5)
        ins.annotate(label, (cai(d), d[field]), xytext=at, textcoords="offset points",
                     ha=ha, va=va, zorder=6)


def add_legend(ax):
    style = {"ls": "", "ms": 2.8, "mew": 0.5}

    def method(c, m):
        return Line2D([], [], marker=m, mec=c, mfc=c, **style)

    def shape(m, filled):
        return Line2D([], [], marker=m, mec="k", mfc="k" if filled else "white", **style)

    faint = Line2D([], [], marker="o", ls="", ms=2.2, mec="none", mfc="k", alpha=0.3)
    derna2 = (method(PURPLE, "s"),
              Line2D([], [], marker="s", mec=PURPLE, mfc="white", **style))
    handles = [derna2, method(ORANGE, "o"), method(GREEN, "o"),
               shape("s", True), shape("s", False), shape("o", True), faint]
    labels = ["DERNA2", "mRNAfold", "LinearDesign,\nLinearCDSfold,\nDERNA, \\ldots",
              "modular, on hull", "modular, above hull", "not modular, on hull",
              "not modular, above hull"]
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0, pos.width, pos.height - 0.19])
    handles.insert(3, Line2D([], [], lw=0))
    labels.insert(3, "")
    ax.legend(handles, labels, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.4)},
              loc="lower left", bbox_to_anchor=(-0.2, 1.0), ncol=2,
              frameon=False, handletextpad=0.4, borderaxespad=0.1, labelspacing=0.15,
              handlelength=1.3, columnspacing=0.5)


def panel_a(fig, cell, fronts, designs):
    ax = fig.add_subplot(cell)
    ax.set_xlim(0.85, 1.003)
    ax.set_ylim(-187, -118)
    plot_fronts(ax, fronts, 2.2, 1.8)
    # the first zoom covers the three designs of (b)-(d)
    xs = [cai(designs[k]) for k in ("d1", "d2", "d3")]
    ys = [designs[k]["free_mfe_kcal_mol"] for k in ("d1", "d2", "d3")]
    ins = zoom_inset(ax, fronts, [0.03, 0.56, 0.47, 0.42], min(xs) - 0.006, max(xs) + 0.0085,
                     min(ys) - 4.0, max(ys) + 6.0, [((0, 1), (0, 0)), ((1, 1), (1, 0))])
    # the second zoom: the hull edge of the modular front with the most designs above it
    pts = sorted((r["cai_log"], r["mfe_kcal_mol"]) for r in fronts.so)
    lo, hi = max(pairwise(fronts.hull), key=lambda e: gap(e, pts))
    zoom_inset(ax, fronts, [0.62, 0.03, 0.37, 0.26],
               math.exp(lo[0] / 78) - 0.0015, math.exp(hi[0] / 78) + 0.0015,
               lo[1] - 1.0, hi[1] + 1.0, [((0, 0), (0, 1)), ((1, 0), (1, 1))])
    mark_designs(ins, designs)
    ax.set_xticks([0.85, 0.9, 0.95, 1.0])
    ax.set_xticklabels(["$0.85$", "$0.9$", "$0.95$", "$1$"])
    ax.set_xlabel("CAI", labelpad=1)
    ax.set_ylabel("MFE of the mRNA (kcal/mol)", labelpad=1)
    ax.tick_params(pad=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    add_legend(ax)
    return ax


def layout(s, aspect, room=0.0, vroom=0.0, align=None):
    # NAView coordinates, rotated to the angle at which the drawing fills a cell of the given
    # width/height best
    RNA.cvar.rna_plot_type = 1
    xy = RNA.get_xy_coordinates(s)
    x0 = [xy.get(i).X for i in range(len(s))]
    y0 = [xy.get(i).Y for i in range(len(s))]

    def rotate(a):
        c, sn = math.cos(a), math.sin(a)
        return ([c * x - sn * y for x, y in zip(x0, y0, strict=True)],
                [sn * x + c * y for x, y in zip(x0, y0, strict=True)])

    def scale(a):
        xs, ys = rotate(a)
        return min(aspect / ((max(xs) - min(xs)) * (1 + room)),
                   1 / ((max(ys) - min(ys)) * (1 + vroom)))

    best = max((math.radians(d) for d in range(0, 180, 3)), key=scale)
    if align:
        # rotate so that the segment between two positions is horizontal (either direction)
        i, j = align
        a0 = -math.atan2(y0[j] - y0[i], x0[j] - x0[i])
        best = max((a0, a0 + math.pi), key=scale)
    return rotate(best)


def draw_structure(fig, ax, s, motif=(), room=0.0, vroom=0.0, align=None):
    # backbone, bases colored by region, base pairs; CDS-UTR pairs in red, motif bases on top
    fig.canvas.draw()
    bb = ax.get_position()
    xs, ys = layout(s, bb.width * fig.get_figwidth() / (bb.height * fig.get_figheight()),
                    room, vroom, align)
    ax.plot(xs, ys, color="#d9d9d9", lw=0.4, zorder=1)
    red = {k for i, j in pairs(s) if (region(i) == "cds") != (region(j) == "cds")
           for k in (i, j)}
    ax.scatter(xs, ys, s=[3.0 if i in red else 1.3 for i in range(len(s))],
               c=[COL["x"] if i in red else COL[region(i)] for i in range(len(s))],
               linewidths=0, zorder=2)
    if motif:
        ax.scatter([xs[i] for i in motif], [ys[i] for i in motif], s=6, c=COL["m"],
                   linewidths=0, zorder=4)
    for i, j in pairs(s):
        x = i in red
        ax.plot([xs[i], xs[j]], [ys[i], ys[j]], color=COL["x"] if x else "#cfcfcf",
                lw=1.6 if x else 0.3, solid_capstyle="butt", zorder=3 if x else 1)
    w = max(xs) - min(xs)
    pad = 0.03 * w
    ax.set_xlim(min(xs) - pad - room * w, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad + vroom * (max(ys) - min(ys)))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_anchor("SW")
    return xs, ys


def motif_sites(seq, motif):
    # the positions of the motif and of its reverse complement in the CDS
    rc = motif[::-1].translate(str.maketrans("ACGU", "UGCA"))
    m0 = L5 + seq[L5:L5 + LC].index(motif)
    m1 = L5 + seq[L5:L5 + LC].index(rc)
    return m0, rc, set(range(m0, m0 + len(motif))) | set(range(m1, m1 + len(rc)))


def cds_utr_helix(s):
    # the longest stack of CDS-UTR base pairs: (i, j), (i - 1, j + 1), ...
    cross = sorted((i, j) for i, j in pairs(s) if (region(i) == "cds") != (region(j) == "cds"))
    best, run = [], []
    for p in cross:
        if run and (p[0] - run[-1][0], p[1] - run[-1][1]) in ((1, -1), (2, -1), (1, -2)):
            run.append(p)
        else:
            run = [p]
        best = max(best, run, key=len)
    return {k for p in best for k in p}


def strand_labels(xs, ys, show):
    # v at the CDS strand and u_3 (or u_5) at the UTR strand of a CDS-UTR helix: (label, color,
    # anchor base) per strand, and whether the strands lie above each other (rows)
    strands = {}
    for i in show:
        strands.setdefault(region(i), []).append(i)
    names = {"cds": r"$\mathbf{v}$", "u5": r"$\mathbf{u}_5$", "u3": r"$\mathbf{u}_3$"}
    centers = [(sum(xs[i] for i in v) / len(v), sum(ys[i] for i in v) / len(v))
               for v in strands.values()]
    (xa, ya), (xb, yb) = centers
    rows = abs(ya - yb) > abs(xa - xb)
    def pick(v):
        return min(v, key=lambda i: xs[i]) if rows else max(v, key=lambda i: ys[i])

    return [(names[k], COL[k], pick(v)) for k, v in strands.items()], rows


def letter_zoom(ax, d, xs, ys, focus, color, bounds, margin=9, corners=None, name=None):
    # the letters of the focus bases and their partners, magnified in an inset of ax
    s, seq = d["free_structure"], d["sequence"]
    show = set(focus)
    show |= {j for i, j in pairs(s) if i in show} | {i for i, j in pairs(s) if j in show}
    sx, sy = [xs[i] for i in show], [ys[i] for i in show]
    rx, ry = 0.5 * (max(sx) - min(sx)) + margin, 0.5 * (max(sy) - min(sy)) + margin
    cx, cy = (max(sx) + min(sx)) / 2, (max(sy) + min(sy)) / 2
    zi = ax.inset_axes(bounds)
    near = [i for i in range(len(xs)) if abs(xs[i] - cx) < rx and abs(ys[i] - cy) < ry]
    for k in range(len(xs) - 1):
        if k in near or k + 1 in near:
            zi.plot(xs[k:k + 2], ys[k:k + 2], color="#d9d9d9", lw=0.4, zorder=1)
    for i, j in pairs(s):
        if i in near and j in near:
            cross = (region(i) == "cds") != (region(j) == "cds")
            zi.plot([xs[i], xs[j]], [ys[i], ys[j]], color=COL["x"] if cross else "#b0b0b0",
                    lw=1.2 if cross else 0.4, zorder=4)
    for i in near:
        # motif letters in the motif color; around CDS-UTR pairs the letters keep their region's
        # color and the pairs themselves are red
        c = color if i in focus and color != COL["x"] else COL[region(i)]
        letter = r"\textbf{" + seq[i] + "}" if i in focus else seq[i]
        zi.text(xs[i], ys[i], letter, ha="center", va="center", fontsize=6.0, color=c,
                zorder=5, bbox={"boxstyle": "circle,pad=0.02", "fc": "white", "ec": "none"})
    # labels inside the inset: left of the rows, or above the columns, in room added there
    x0, x1, y0, y1 = cx - rx, cx + rx, cy - ry, cy + ry
    labels, rows = [], True
    if color == COL["x"]:
        labels, rows = strand_labels(xs, ys, show)
    if name:
        labels = [(name, color, min(focus, key=lambda i: xs[i]))]
    if labels and rows:
        x0 -= (0.35 if name else 0.12) * (x1 - x0)
        for text, c, i in labels:
            zi.text(x0 + 0.02 * (x1 - x0), ys[i], text, ha="left", va="center", fontsize=FSS,
                    color=c)
    elif labels:
        y1 += 0.3 * (y1 - y0)
        for text, c, i in labels:
            zi.text(xs[i], y1 - 0.03 * (y1 - y0), text, ha="center", va="top", fontsize=FSS,
                    color=c)
    zi.set_xlim(x0, x1)
    zi.set_ylim(y0, y1)
    # stretched to the inset's shape, so that the two strands of a helix sit apart
    zi.set_aspect("auto")
    zi.set_xticks([])
    zi.set_yticks([])
    for sp in zi.spines.values():
        sp.set_linewidth(0.4)
        sp.set_color(color)
    # by default the box's top corners join the inset's bottom ones; straight connectors
    corners = corners or [((0, 1), (0, 0)), ((1, 1), (1, 0))]
    zoom_box(ax, zi, x0, x1, y0, y1, color, corners, rad=0.0)


def structure_panels(fig, ax, cell, designs, key_height, key_last):
    # right side: the key's last row sits on the bottom of (a), the structures fill the rest
    fig.canvas.draw()
    abot = (ax.get_tightbbox(fig.canvas.get_renderer())
            .transformed(fig.transFigure.inverted()).y0)
    kb = abot - key_last * key_height + 0.012
    gr = fig.add_gridspec(2, 3, left=cell.get_position(fig).x0, right=0.995, top=0.80,
                          bottom=kb, height_ratios=[0.80 - kb - key_height, key_height],
                          wspace=0.02, hspace=0.0)
    axs = [fig.add_subplot(gr[0, k]) for k in range(3)]
    d1, d2 = designs["d1"], designs["d2"]
    helix = cds_utr_helix(d1["free_structure"])
    lo = min(i for i in helix if region(i) == "cds")
    hi = max(i for i in helix if region(i) == "cds")
    xs, ys = draw_structure(fig, axs[0], d1["free_structure"], vroom=1.4, align=(lo, hi))
    axs[0].set_anchor("C")
    letter_zoom(axs[0], d1, xs, ys, helix, COL["x"], [0.0, 0.62, 1.0, 0.38])
    m0, _, mot = motif_sites(d2["sequence"], designs["motif"])
    xs, ys = draw_structure(fig, axs[1], d2["free_structure"], mot, vroom=0.65,
                            align=(m0, m0 + len(designs["motif"]) - 1))
    axs[1].set_anchor("C")
    letter_zoom(axs[1], d2, xs, ys, mot, COL["m"], [0.0, 0.62, 1.0, 0.38],
                name=designs["motif_name"])
    letter_zoom(axs[1], d2, xs, ys, cds_utr_helix(d2["free_structure"]), COL["x"],
                [0.68, 0.0, 0.32, 0.30], margin=6,
                corners=[((1, 1), (0, 1)), ((1, 0), (0, 0))])
    draw_structure(fig, axs[2], designs["d3"]["free_structure"])
    return gr, axs


def add_titles(fig, ax, cell, gr, axs, designs):
    d1, d2, d3 = designs["d1"], designs["d2"], designs["d3"]

    def stats(d):
        return f"CAI {cai(d):.3f}, MFE {mfe_text(d['free_mfe_kcal_mol'])}"

    titles = [
        (ax, r"\textbf{(a)} GLPE with the hbb UTRs"),
        (axs[0], r"\textbf{(b)} $(\mathbf{v}'',P'')$" + "\n" + stats(d1) + "\n"
         + f"{cds_utr_pairs(d1['free_structure'])} CDS--UTR pairs"),
        (axs[1], r"\textbf{(c)} $(\mathbf{v}',P')$" + "\n" + stats(d2) + "\n"
         + f"{cds_utr_pairs(d2['free_structure'])} CDS--UTR pairs, "
         + designs["motif_name"] + " site"),
        (axs[2], r"\textbf{(d)} $(\mathbf{v},P)$" + "\n" + stats(d3) + "\n"
         + "no CDS--UTR pairs"),
    ]
    # (a)'s title starts at its y-axis label, the others at their cells
    cells = {id(a): gr[0, k].get_position(fig).x0 for k, a in enumerate(axs)}
    cells[id(ax)] = cell.get_position(fig).x0
    fig.canvas.draw()
    for a, t in titles:
        renderer = fig.canvas.get_renderer()
        if a is ax:
            x = ax.get_tightbbox(renderer).transformed(fig.transFigure.inverted()).x0
        else:
            x = cells[id(a)]
        fig.text(x, 0.975, t, ha="left", va="top", fontsize=FS)


def add_key(fig, cell, w, u5, u3, motif, rc, motif_name):
    lax = fig.add_subplot(cell)
    lax.axis("off")
    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    rows = ((COL["u5"], r"5$'$ UTR $\mathbf{u}_5$", [tt(u5) + r" (50 nt)"]),
            (COL["cds"], r"protein $\mathbf{w}$", [tt(w[:39]), tt(w[39:]) + r" (78 aa)"]),
            (COL["u3"], r"3$'$ UTR $\mathbf{u}_3$",
             [tt(u3[:22] + r"\ldots{}" + u3[-22:]) + r" (132 nt)"]),
            (COL["m"], r"motifs $\mathcal{M}$",
             [", ".join(tt(m) for m in dict.fromkeys((motif, rc)))
              + r" (forbidden, " + motif_name + ")"]))
    y, dy = 0.92, 0.165
    for c, name, seqs in rows:
        lax.plot([0.0, 0.035], [y, y], color=c, lw=2, solid_capstyle="butt")
        lax.text(0.045, y, name, va="center", fontsize=FSS)
        for k, sq in enumerate(seqs):
            lax.text(0.165, y - k * dy, sq, va="center", fontsize=FSS)
        y -= dy * max(1, len(seqs))
    lax.plot([0.0, 0.035], [y, y], color=COL["x"], lw=2, solid_capstyle="butt")
    lax.text(0.045, y, "CDS--UTR base pair", va="center", fontsize=FSS)


def main():
    results, mrnafold, designs_json, fasta, utr5, utr3, out = sys.argv[1:8]
    fronts = load_fronts(results, mrnafold)
    with open(designs_json) as f:
        designs = json.load(f)
    w, u5, u3 = load_inputs(fasta, utr5, utr3, designs)
    set_style()
    fig = plt.figure(figsize=(6.5, 2.25))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 3], left=0.075, right=0.995, top=0.80,
                          bottom=0.14, wspace=0.05)
    ax = panel_a(fig, gs[0], fronts, designs)
    # the key's height (figure fraction) and its last row (key coordinates)
    gr, axs = structure_panels(fig, ax, gs[1], designs, key_height=0.30, key_last=0.095)
    add_titles(fig, ax, gs[0], gr, axs, designs)
    _, rc, _ = motif_sites(designs["d2"]["sequence"], designs["motif"])
    add_key(fig, gr[1, :], w, u5, u3, designs["motif"], rc, designs["motif_name"])
    fig.savefig(out)
    fig.savefig(out.rsplit(".", 1)[0] + ".png", dpi=300)


if __name__ == "__main__":
    main()
