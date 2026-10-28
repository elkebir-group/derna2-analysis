"""Drawing helpers for the revised overview figure (overview2.py): the fronts of panel (b), its
zoom insets and the marks on the three designs, and the RNApuzzler layouts of the structures with
their letter-level zooms.

Positions are split into regions by the module constants L5, LC and L3 (5' UTR, CDS and 3' UTR
lengths in nt), which overview2.py sets from the UTRs of the figure it draws. Text is set by LaTeX
in Times at the sizes of the paper's figures (7 pt, small labels 6.3 pt).
"""

import math
import re
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

import matplotlib
import RNA
from matplotlib.patches import ConnectionPatch, Rectangle

from derna2_analysis.io import read_jsonl

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

NAMES = {"CUCGAG": "XhoI"}  # restriction sites named in the figure
# region lengths in nt; overview2.py sets L5 and L3 from the UTRs it is given
L5, LC, L3 = 37, 234, 111
# (b)-(d) use none of the method colors of (a)
COL = {"u5": "#56b4e9", "cds": "#1f3a93", "u3": "#9a9a9a", "x": "#d62728", "m": "#c51b8a"}
GREEN, ORANGE, PURPLE = "#2e9e5b", "#e8772e", "#6a4bc4"
FS, FSS = 7, 6.3
PREAMBLE = (
    r"\usepackage[T1]{fontenc}\usepackage{times}\usepackage{textcomp}"
    r"\usepackage{amsmath}\usepackage{pifont}"
)


@dataclass
class Fronts:
    """The three fronts of (a) and which of their designs lie on their methods' hulls."""

    b: list  # concatenation designs, refolded with the UTRs
    mf: list  # mRNAfold's designs, CDS-UTR pairs discouraged
    cd: list  # the UTR-masking front, which mRNAfold's weights sample
    so: list  # the modular front
    hull: list  # convex hull of the modular front in (log CAI, MFE)
    on_hull_b: set  # CDS sequences on the hull of the front without UTRs

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
        while (
            len(h) >= 2
            and (
                (h[-1][0] - h[-2][0]) * (q[1] - h[-2][1])
                - (h[-1][1] - h[-2][1]) * (q[0] - h[-2][0])
            )
            <= 0
        ):
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
    heights = [
        q[1] - (ey0 + (q[0] - ex0) / (ex1 - ex0) * (ey1 - ey0)) for q in pts if ex0 < q[0] < ex1
    ]
    return max(heights or [0])


def load_mrnafold(path):
    # the guarded arm, with its CAI under our codon weights
    return [dict(r, cai_log=r["theirs_cai_log"]) for r in read_jsonl(path) if r["arm"] == "guard"]


def load_fronts(results, mrnafold):
    res = Path(results)
    b = read_jsonl(res / "b_cdsonly_with_utrs.scored.jsonl")
    mf = load_mrnafold(mrnafold)
    so = read_jsonl(res / "d_so_front.jsonl")
    # no green or orange point drawn is modular (modular competitor designs lie below CAI 0.85)
    assert not any(r["symmetric_so"] for r in b + mf if cai(r) >= 0.85)
    on_hull_b = {
        r["sequence"] for r in read_jsonl(res / "a_cdsonly_front.jsonl") if r.get("hull_vertex")
    }
    return Fronts(
        b=b,
        mf=mf,
        so=so,
        cd=read_jsonl(res / "cd_symmetric_guarded_front.scored.jsonl"),
        hull=lower_hull([(r["cai_log"], r["mfe_kcal_mol"]) for r in so]),
        on_hull_b=on_hull_b,
    )


def load_inputs(fasta, utr5, utr3, designs):
    w = "".join(
        line.strip() for line in Path(fasta).read_text().splitlines() if not line.startswith(">")
    )
    u5, u3 = Path(utr5).read_text().strip(), Path(utr3).read_text().strip()
    assert (len(u5), len(u3), len(w)) == (L5, L3, 78)
    assert all(
        designs[k]["sequence"] == u5 + designs[k]["sequence"][L5 : L5 + LC] + u3
        for k in ("d1", "d2", "d3")
    )
    return w, u5, u3


def set_style():
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "text.latex.preamble": PREAMBLE,
            "font.size": FS,
            "axes.titlesize": FS,
            "axes.labelsize": FS,
            "xtick.labelsize": FSS,
            "ytick.labelsize": FSS,
            "legend.fontsize": FSS,
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
        }
    )


def zoom_box(parent, child, x0, x1, y0, y1, color, corners, rad=0.0, z=6):
    # the zoomed region on the parent and two dashed lines to the matching corners of the inset
    parent.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=color, lw=0.5, zorder=z))
    for (px, py), (cx, cy) in corners:
        parent.add_artist(
            ConnectionPatch(
                xyA=((x0, x1)[px], (y0, y1)[py]),
                coordsA=parent.transData,
                xyB=(cx, cy),
                coordsB=child.transAxes,
                color="#999999" if color == "k" else color,
                lw=0.5,
                ls="-",
                connectionstyle=f"arc3,rad={rad}",
                zorder=z,
            )
        )


def plot_others(a, fronts, ms):
    # not modular: mRNAfold's designs and the hull vertices of the concatenation front, the only
    # designs a weighted sum returns, as circles at the MFE of the refolded mRNA. Most of them
    # are the same design, so the concatenation circles are drawn larger, behind.
    # in the background: the designs of the two fronts that no weight reaches, as open circles
    mf_seqs = {r["sequence"] for r in fronts.mf}
    hull_b = [r for r in fronts.b if not fronts.on_b(r)]
    hull_cd = [r for r in fronts.cd if r["sequence"] not in mf_seqs]
    for rows, col in ((hull_b, GREEN), (hull_cd, ORANGE)):
        a.plot(
            [cai(r) for r in rows],
            [r["free_mfe_kcal_mol"] for r in rows],
            "o",
            ms=0.9 * ms,
            mec=col,
            mfc="none",
            mew=0.45,
            zorder=1,
        )
    for rows, on, col, size, z in (
        (fronts.b, fronts.on_b, GREEN, 1.7 * ms, 2.5),
        (fronts.mf, fronts.on_mf, ORANGE, ms, 2.6),
    ):
        sel = [r for r in by_cai(rows) if on(r)]
        a.plot(
            [cai(r) for r in sel],
            [r["free_mfe_kcal_mol"] for r in sel],
            "o",
            ms=size,
            mec=col,
            mfc=col,
            mew=0.5,
            zorder=z,
        )


def plot_modular(a, fronts, ms):
    # the modular front: filled squares on its convex hull, open squares above it
    rows = by_cai(fronts.so)
    vertices = set(fronts.hull)
    xs = [cai(r) for r in rows]
    ys = [r["mfe_kcal_mol"] for r in rows]
    a.plot(xs, ys, color=PURPLE, lw=0.7, zorder=3)
    on = [(r["cai_log"], r["mfe_kcal_mol"]) in vertices for r in rows]
    a.plot(
        [x for x, o in zip(xs, on, strict=True) if o],
        [y for y, o in zip(ys, on, strict=True) if o],
        "s",
        ms=ms,
        color=PURPLE,
        zorder=3.2,
    )
    a.plot(
        [x for x, o in zip(xs, on, strict=True) if not o],
        [y for y, o in zip(ys, on, strict=True) if not o],
        "s",
        ms=ms,
        mfc="white",
        mec=PURPLE,
        mew=0.5,
        zorder=3.2,
    )


def plot_fronts(a, fronts, ms_others, ms_modular):
    plot_others(a, fronts, ms_others)
    plot_modular(a, fronts, ms_modular)
    hx, hy = hull_curve(fronts.hull)
    a.plot(hx, hy, "--", color=PURPLE, lw=0.5, zorder=2.8)


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


def mark_designs(ins, designs, label_at=None):
    # the designs of (b)-(d), circled and labeled in the zoom
    # each label at the upper right of its circle
    marks = (
        ("d1", "free_mfe_kcal_mol", r"$\mathbf{v}''$", (0, 3), "left", "bottom"),
        ("d2", "free_mfe_kcal_mol", r"$\mathbf{v}'$", (-3, 2), "right", "bottom"),
        ("d3", "mfe_kcal_mol", r"$\mathbf{v}$", (2, 2), "left", "bottom"),
    )
    for key, field, label, at, ha, va in marks:
        # a design file may move a label: {"d2": [dx, dy, ha, va]}
        if label_at and key in label_at:
            *at, ha, va = label_at[key]
        d = designs[key]
        ins.plot(cai(d), d[field], "o", ms=5.0, mfc="none", mec="k", mew=0.6, zorder=5)
        ins.annotate(
            label, (cai(d), d[field]), xytext=at, textcoords="offset points", ha=ha, va=va, zorder=6
        )


def layout(s, aspect, room=0.0, vroom=0.0, align=None):
    # RNApuzzler coordinates, rotated to the angle at which the drawing fills a cell of the given
    # width/height best
    RNA.cvar.rna_plot_type = 4
    xy = RNA.get_xy_coordinates(s)
    x0 = [xy.get(i).X for i in range(len(s))]
    y0 = [xy.get(i).Y for i in range(len(s))]

    def rotate(a):
        c, sn = math.cos(a), math.sin(a)
        return (
            [c * x - sn * y for x, y in zip(x0, y0, strict=True)],
            [sn * x + c * y for x, y in zip(x0, y0, strict=True)],
        )

    def scale(a):
        xs, ys = rotate(a)
        return min(
            aspect / ((max(xs) - min(xs)) * (1 + room)), 1 / ((max(ys) - min(ys)) * (1 + vroom))
        )

    best = max((math.radians(d) for d in range(0, 180, 3)), key=scale)
    if align:
        # rotate so that the segment between two positions is horizontal (either direction)
        i, j = align
        a0 = -math.atan2(y0[j] - y0[i], x0[j] - x0[i])
        best = a0
    return rotate(best)


def draw_structure(fig, ax, s, motif=(), room=0.0, vroom=0.0, align=None):
    # backbone, bases colored by region, base pairs; CDS-UTR pairs in red, motif bases on top
    fig.canvas.draw()
    bb = ax.get_position()
    xs, ys = layout(
        s, bb.width * fig.get_figwidth() / (bb.height * fig.get_figheight()), room, vroom, align
    )
    ax.plot(xs, ys, color="#d9d9d9", lw=0.4, zorder=1)
    red = {k for i, j in pairs(s) if (region(i) == "cds") != (region(j) == "cds") for k in (i, j)}
    ax.scatter(
        xs,
        ys,
        s=[3.0 if i in red else 1.3 for i in range(len(s))],
        c=[COL["x"] if i in red else COL[region(i)] for i in range(len(s))],
        linewidths=0,
        zorder=2,
    )
    if motif:
        ax.scatter(
            [xs[i] for i in motif], [ys[i] for i in motif], s=6, c=COL["m"], linewidths=0, zorder=4
        )
    for i, j in pairs(s):
        x = i in red
        ax.plot(
            [xs[i], xs[j]],
            [ys[i], ys[j]],
            color=COL["x"] if x else "#cfcfcf",
            lw=1.6 if x else 0.3,
            solid_capstyle="butt",
            zorder=3 if x else 1,
        )
    w = max(xs) - min(xs)
    pad = 0.03 * w
    ax.set_xlim(min(xs) - pad - room * w, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad + vroom * (max(ys) - min(ys)))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_anchor("SW")
    return xs, ys


def motif_sites(seq, motifs):
    # every occurrence of a forbidden motif in the CDS (R = A or G): the first one's start,
    # length and spelling, and all the positions any occurrence covers
    cds, found, pos = seq[L5 : L5 + LC], [], set()
    for m in motifs:
        for hit in re.finditer("(?=" + m.replace("R", "[AG]") + ")", cds):
            found.append((L5 + hit.start(), m))
            pos |= set(range(L5 + hit.start(), L5 + hit.start() + len(m)))
    m0, m = min(found)
    return m0, len(m), seq[m0 : m0 + len(m)], pos


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


def strand_labels(xs, ys, show, cds_name=r"$\mathbf{v}$"):
    # v at the CDS strand and u_3 (or u_5) at the UTR strand of a CDS-UTR helix: (label, color,
    # anchor base) per strand, and whether the strands lie above each other (rows)
    strands = {}
    for i in show:
        strands.setdefault(region(i), []).append(i)
    names = {"cds": cds_name, "u5": r"$\mathbf{u}_5$", "u3": r"$\mathbf{u}_3$"}
    centers = [
        (sum(xs[i] for i in v) / len(v), sum(ys[i] for i in v) / len(v)) for v in strands.values()
    ]
    (xa, ya), (xb, yb) = centers
    rows = abs(ya - yb) > abs(xa - xb)

    def pick(v):
        return min(v, key=lambda i: xs[i]) if rows else max(v, key=lambda i: ys[i])

    return [(names[k], COL[k], pick(v)) for k, v in strands.items()], rows


def letter_zoom(
    ax,
    d,
    xs,
    ys,
    focus,
    color,
    bounds,
    margin=9,
    corners=None,
    name=None,
    cds_name=r"$\mathbf{v}$",
    partners=True,
):
    # the letters of the focus bases and their partners, magnified in an inset of ax
    s, seq = d["free_structure"], d["sequence"]
    show = set(focus)
    if partners:
        show |= {j for i, j in pairs(s) if i in show} | {i for i, j in pairs(s) if j in show}
    sx, sy = [xs[i] for i in show], [ys[i] for i in show]
    rx, ry = 0.5 * (max(sx) - min(sx)) + margin, 0.5 * (max(sy) - min(sy)) + margin
    cx, cy = (max(sx) + min(sx)) / 2, (max(sy) + min(sy)) / 2
    zi = ax.inset_axes(bounds)
    near = [i for i in range(len(xs)) if abs(xs[i] - cx) < rx and abs(ys[i] - cy) < ry]
    for k in range(len(xs) - 1):
        if k in near or k + 1 in near:
            zi.plot(xs[k : k + 2], ys[k : k + 2], color="#d9d9d9", lw=0.4, zorder=1)
    for i, j in pairs(s):
        if i in near and j in near:
            cross = (region(i) == "cds") != (region(j) == "cds")
            zi.plot(
                [xs[i], xs[j]],
                [ys[i], ys[j]],
                color=COL["x"] if cross else "#b0b0b0",
                lw=1.2 if cross else 0.4,
                zorder=4,
            )
    # letters are not clipped by the inset, so draw only those whose glyph fits inside it
    inside = [i for i in near if abs(xs[i] - cx) < rx - 3 and abs(ys[i] - cy) < ry - 3]
    for i in inside:
        # motif letters in the motif color; around CDS-UTR pairs the letters keep their region's
        # color and the pairs themselves are red
        c = color if i in focus and color != COL["x"] else COL[region(i)]
        letter = r"\textbf{" + seq[i] + "}" if i in focus else seq[i]
        zi.text(
            xs[i],
            ys[i],
            letter,
            ha="center",
            va="center",
            fontsize=6.0,
            color=c,
            zorder=5,
            bbox={"boxstyle": "circle,pad=0.02", "fc": "white", "ec": "none"},
        )
    # labels inside the inset: left of the rows, or above the columns, in room added there
    x0, x1, y0, y1 = cx - rx, cx + rx, cy - ry, cy + ry
    labels, rows = [], True
    if color == COL["x"]:
        labels, rows = strand_labels(xs, ys, show, cds_name)
    if name:
        labels = [(name, color, min(focus, key=lambda i: xs[i]))]
    if labels and rows:
        x0 -= (0.75 if name else 0.3) * (x1 - x0)
        for text, c, i in labels:
            zi.text(
                x0 + 0.02 * (x1 - x0), ys[i], text, ha="left", va="center", fontsize=FSS, color=c
            )
    elif labels:
        y1 += 0.3 * (y1 - y0)
        for text, c, i in labels:
            zi.text(
                xs[i], y1 - 0.03 * (y1 - y0), text, ha="center", va="top", fontsize=FSS, color=c
            )
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
    zoom_box(ax, zi, x0, x1, y0, y1, color, corners, rad=0.0, z=0.5)
