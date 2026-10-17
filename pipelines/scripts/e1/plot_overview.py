"""E1 overview figure: what the hull of the front leaves out, and what the whole front costs.

    python plot_overview.py SWEEP_JSONL WALK_DIR PIN OUT_PDF

SWEEP_JSONL holds the sweep rows; only those stamped with PIN are read. WALK_DIR holds the
walk reports `<protein>.walk.json`, read for `sweep_s`, the wall time of the pruned hull
search (the sweep rows do not carry it; see "Open items" in the README).

Panels: (a) one instance: its front with the hull vertices marked, the steps between
neighbouring designs, and each method's wall time; (b) the share of each front that a λ
reaches; (c) the steps between neighbouring designs on the hull against the whole front;
(d) wall time per method with the cells that did not finish at their budget; (e) the
speed-up of pruning the hull search, measured against the unpruned fixed-λ DP at the hull
vertices' λ values, and the cost of the whole-front walk over the pruned search.

The canvas is drawn in pixel coordinates (2450 x 960 design units) so the figure is one
fixed layout; text is sized for a full-width figure at about 6.4 pt. The palette is this
figure's own, not `derna2_analysis.figures.SLOTS`: it is the one the paper's figure uses.
"""

import json
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
import numpy as np
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

from derna2_analysis import e1_overview as ov
from derna2_analysis.io import read_jsonl

matplotlib.use("Agg")

W, H = 2450, 882
# Design units cropped off the top: the layout below is written on a 960-unit-tall canvas whose
# ink starts 39 units down and ends 897 down (the rotated method names of (a) and (d)), so 27
# off the top and a height of 882 leave 12 units of margin above and below.
TOP = 27
# Centre of the "protein length (aa)" titles of (b), (c) and (e): 12 units below their rotated
# tick labels, which end 842 down.
XLABEL_Y = 872
# Left and right margin and the gap between panels, in design units, measured on the ink;
# panel (c) takes whatever width is left over (Canvas.even_out).
SIDE, GAP, STRETCH = 12, 26, "(c)"
PX = 72 / 100 * 1.4
TK, LB, AN = 24 * PX, 26 * PX, 21 * PX
Y0, Y1 = 132.5, 713  # tops lowered 8 units to clear the descending parentheses of "(a)"-"(e)"
INSTANCE = "P47211"
PUR, PUR_FILL, PUR_DOT = "#6a4fc8", "#b6a8e4", "#7457c8"
GRY, GRY_FILL, GRY_SWATCH = "#4a4a4a", "#9a9a9a", "#7a7a7a"
METHOD_COLOR = {
    "DERNA2": "#7457c8",
    "mRNAfold": "#e9803c",
    "LinearDesign": "#3ca466",
    "LinearCDSfold": "#2ea5b9",
    "DERNA1": "#94773f",
}
LATIN_MODERN = "/usr/share/texmf/fonts/opentype/public/lm/lmroman10-{}.otf"


def setup_fonts():
    family = "serif"
    for weight in ("regular", "bold"):
        path = LATIN_MODERN.format(weight)
        if Path(path).exists():
            fm.fontManager.addfont(path)
            family = "Latin Modern Roman"
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "font.family": family,
            "mathtext.fontset": "cm",
            "axes.linewidth": 1.5,
            "xtick.major.width": 1.5,
            "ytick.major.width": 1.5,
            "xtick.major.size": 5.5,
            "ytick.major.size": 5.5,
        }
    )


class Canvas:
    """A figure addressed in design pixels, origin top left."""

    def __init__(self):
        self.fig = plt.figure(figsize=(W / 100, H / 100), dpi=100, facecolor="white")
        self.panel, self.groups = None, {}  # the panel being drawn; its axes and texts

    def fy(self, y):
        """Figure fraction of design y (measured down from the uncropped top)."""
        return 1 - (y - TOP) / H

    def axes(self, x0, y0, x1, y1):
        a = self.fig.add_axes([x0 / W, self.fy(y1), (x1 - x0) / W, (y1 - y0) / H])
        if self.panel:
            self.groups[self.panel].append(a)
        a.set_facecolor("white")
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.tick_params(labelsize=TK)
        return a

    def text(self, x, y, s, **k):
        t = self.fig.text(x / W, self.fy(y), s, **k)
        if self.panel:
            self.groups[self.panel].append(t)
        return t

    def letter(self, x, s):
        """Start panel s: its letter, and every axes and text until the next letter."""
        self.panel = s
        self.groups[s] = []
        self.text(x, 109, s, fontsize=27 * PX, fontweight="bold", ha="left", va="baseline")

    def even_out(self):
        """Slide whole panels sideways so that, measured on the drawn ink, the outer margins are
        SIDE and the gaps between panels GAP, widening the STRETCH panel to take up the rest.
        The first pass places panels by their bounding boxes, which pad text beyond its glyphs;
        the later ones measure each panel's rendered pixels alone and correct for that."""
        fig = self.fig
        groups = list(self.groups.values())

        def boxes(group, r):
            return [a.get_window_extent(r) if hasattr(a, "get_text") else a.get_tightbbox(r)
                    for a in group]

        def spans_by_box():
            fig.canvas.draw()
            r = fig.canvas.get_renderer()
            out = []
            for g in groups:
                bs = boxes(g, r)
                out.append((min(b.x0 for b in bs), max(b.x1 for b in bs),
                            min(b.y0 for b in bs), max(b.y1 for b in bs)))
            return out

        def spans_by_ink():
            """Each panel rendered alone (every other artist hidden), and its inked columns."""
            everything = list(fig.axes) + list(fig.texts) + list(fig.artists) + list(fig.lines)
            out = []
            for g in groups:
                mine = {id(a) for a in g}
                shown = [a for a in everything if a.get_visible()]
                for a in shown:
                    a.set_visible(id(a) in mine)
                fig.canvas.draw()
                img = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].min(axis=2) < 235
                cols = np.flatnonzero(img.any(axis=0))
                out.append((cols[0], cols[-1] + 1))
                for a in shown:
                    a.set_visible(True)
            return out

        def shift(spans):
            width = fig.bbox.width
            ink = sum(x1 - x0 for x0, x1 in spans)
            gap = (width - 2 * SIDE * width / W - ink) / (len(spans) - 1)
            x = SIDE * width / W
            for group, (x0, x1) in zip(groups, spans, strict=True):
                dx = (x - x0) / width
                for a in group:
                    if hasattr(a, "get_text"):
                        a.set_x(a.get_position()[0] + dx)
                    else:
                        pos = a.get_position()
                        a.set_position([pos.x0 + dx, pos.y0, pos.width, pos.height])
                x += x1 - x0 + gap
            return gap * W / width

        shift([(x0, x1) for x0, x1, _, _ in spans_by_box()])
        # Widen the STRETCH panel's axes by what is left after the margins and GAP-wide gaps;
        # its centred titles move with the axes' centre.
        # A y label placed in axes fractions moves out as its axes widen, so repeat until the
        # leftover is under a pixel.
        u = fig.bbox.width / W
        for _ in range(6):
            spans = spans_by_ink()
            ink = sum(b - a for a, b in spans)
            extra = fig.bbox.width - (2 * SIDE + (len(spans) - 1) * GAP) * u - ink
            if abs(extra) < 1:
                break
            for a in self.groups[STRETCH]:
                if hasattr(a, "get_text"):
                    if a.get_ha() == "center":
                        a.set_x(a.get_position()[0] + extra / 2 / fig.bbox.width)
                else:
                    pos = a.get_position()
                    a.set_position([pos.x0, pos.y0, pos.width + extra / fig.bbox.width, pos.height])
        gap = shift(spans_by_ink())
        fig.canvas.draw()
        return gap


def align_ylabels_left(fig, axes):
    """Move the y labels to the left edge of the leftmost one, leaving the axes in place."""
    fig.canvas.draw()
    labels = [a.yaxis.label for a in axes]
    left = min(lb.get_window_extent().x0 for lb in labels)
    for a, lb in zip(axes, labels, strict=True):
        x, y = lb.get_transform().transform(lb.get_position())
        x_axes = a.transAxes.inverted().transform((x + left - lb.get_window_extent().x0, y))[0]
        a.yaxis.set_label_coords(x_axes, 0.5)


def minor_log(ax, axis):
    loc = mt.LogLocator(base=10, subs=np.arange(2, 10), numticks=20)
    getattr(ax, axis + "axis").set_minor_locator(loc)
    getattr(ax, axis + "axis").set_minor_formatter(mt.NullFormatter())
    ax.tick_params(axis=axis, which="minor", length=3.2, width=1.2)


def group_ticks(ax):
    ax.set_xticks(range(len(ov.GROUPS)))
    ax.set_xticklabels(
        ["$\\leq$200" if i == 0 else g for i, (g, _) in enumerate(ov.GROUPS)],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=TK,
    )
    ax.tick_params(axis="x", pad=4)


def box(ax, x, v, fc, dc, w=0.32, lw=1.0):
    """Quartile box with Tukey whiskers and the median bar, vertical; log-aware."""
    log = ax.get_yscale() == "log"
    lv = np.log10(v) if log else np.asarray(v)
    f = (lambda t: 10**t) if log else (lambda t: t)
    q1, med, q3 = np.percentile(lv, [25, 50, 75])
    iqr = q3 - q1
    lo, hi = lv[lv >= q1 - 1.5 * iqr].min(), lv[lv <= q3 + 1.5 * iqr].max()
    ax.add_patch(
        Rectangle(
            (x - w / 2, f(q1)), w, f(q3) - f(q1), facecolor=fc, alpha=0.55, edgecolor=dc, lw=lw, zorder=2
        )
    )
    ax.plot([x, x], [f(lo), f(q1)], color="#222", lw=1.1, zorder=3)
    ax.plot([x, x], [f(q3), f(hi)], color="#222", lw=1.1, zorder=3)
    for yy in (lo, hi):
        ax.plot([x - w * 0.22, x + w * 0.22], [f(yy)] * 2, color="#222", lw=1.1, zorder=3)
    ax.plot([x - w / 2, x + w / 2], [f(med)] * 2, color="#111", lw=2, zorder=5, solid_capstyle="butt")


def hbox(ax, y, v, fc, dc, h=0.44):
    """The same box, horizontal, on a log x axis."""
    lv = np.log10(v)
    q1, med, q3 = np.percentile(lv, [25, 50, 75])
    iqr = q3 - q1
    lo, hi = lv[lv >= q1 - 1.5 * iqr].min(), lv[lv <= q3 + 1.5 * iqr].max()
    ax.add_patch(
        Rectangle((10**q1, y - h / 2), 10**q3 - 10**q1, h, facecolor=fc, alpha=0.6, edgecolor=dc, lw=1.0, zorder=3)
    )
    ax.plot([10**lo, 10**q1], [y, y], color="#222", lw=1.1, zorder=3)
    ax.plot([10**q3, 10**hi], [y, y], color="#222", lw=1.1, zorder=3)
    for xx in (lo, hi):
        ax.plot([10**xx] * 2, [y - 0.1, y + 0.1], color="#222", lw=1.1, zorder=3)
    ax.plot([10**med] * 2, [y - h / 2, y + h / 2], color="#111", lw=2, zorder=5, solid_capstyle="butt")


def sample(v, rng, k=150):
    return v if len(v) <= k else rng.choice(v, k, replace=False)


def dots(ax, x, v, big, color, rng, jit=0.1, s=9):
    j = rng.uniform(-jit, jit, len(v))
    ax.scatter((x + j)[~big], v[~big], s=s, color=color, alpha=0.8, lw=0, zorder=4)
    ax.scatter((x + j)[big], v[big], s=55, marker="^", color=color, edgecolor="#222", lw=0.8, zorder=6)


def legend_row(cv):
    items = list(METHOD_COLOR.items())
    items.insert(1, ("DERNA2 (hull)", GRY_SWATCH))
    r = cv.fig.canvas.get_renderer()

    def width(s):
        t = cv.text(0, 0, s, fontsize=TK)
        w = t.get_window_extent(r).width
        t.remove()
        return w

    sw, gap, y = 24, 44, 56
    widths = [sw + 10 + width(k) for k, _ in items]
    widths += [16 + 10 + width("1000 aa") + 22, 16 + 10 + width("> 1000 aa")]
    x = (W - (sum(widths) + gap * (len(widths) - 1))) / 2
    for (k, c), w in zip(items, widths, strict=False):
        cv.fig.add_artist(
            Rectangle(
                (x / W, cv.fy(y + 9)), sw / W, 18 / H, transform=cv.fig.transFigure,
                facecolor=c, alpha=0.45, edgecolor=c, lw=1.0,
            )
        )
        cv.text(x + sw + 10, y, k, fontsize=TK, va="center")
        x += w + gap
    cv.fig.add_artist(
        Line2D([(x + 8) / W], [cv.fy(y)], marker="o", markersize=8, color="#555", lw=0,
               transform=cv.fig.transFigure)
    )
    cv.text(x + 26, y, "$\\leq$1000 aa", fontsize=TK, va="center")
    x += widths[-2] + gap
    cv.fig.add_artist(
        Line2D([(x + 8) / W], [cv.fy(y)], marker="^", markersize=9, markerfacecolor="#555",
               markeredgecolor="#222", lw=0, transform=cv.fig.transFigure)
    )
    cv.text(x + 26, y, ">1000 aa", fontsize=TK, va="center")


def scatter_front(ax, f, size_hull, size_other, lw_other):
    ax.scatter(f.cai[~f.on_hull], f.mfe[~f.on_hull], s=size_other, facecolor="white",
               edgecolor=PUR_DOT, linewidths=lw_other, zorder=2)
    ax.scatter(f.cai[f.on_hull], f.mfe[f.on_hull], s=size_hull, color=GRY, lw=0, zorder=3)


def double_arrow(ax, xy, xytext):
    ax.annotate("", xy=xy, xytext=xytext, zorder=6,
                arrowprops=dict(arrowstyle="<->", color="#222", lw=1.2, shrinkA=0, shrinkB=0))


def panel_instance(cv, f, cells, rng):
    cv.letter(17, "(a)")
    scat = cv.axes(410, Y0, 710, 500)
    scatter_front(scat, f, 20, 3.0, 0.45)
    scat.set_xlim(0.715, 1.005)
    scat.set_ylim(-700, -395)
    scat.set_xticks([0.8, 0.9, 1.0])
    scat.set_xticklabels(["0.8", "0.9", "1.0"], fontsize=TK)
    scat.set_yticks([-700, -600, -500, -400])
    scat.set_yticklabels(["$-$700", "$-$600", "$-$500", "$-$400"], fontsize=TK)
    scat.set_xlabel("CAI", fontsize=LB, labelpad=2)
    scat.set_ylabel("MFE (kcal/mol)", fontsize=LB, labelpad=4)
    scat.legend(
        handles=[
            Line2D([], [], marker="o", ls="", markersize=6, color=GRY, label="on hull"),
            Line2D([], [], marker="o", ls="", markersize=6, markerfacecolor="none",
                   markeredgecolor=PUR_DOT, label="off hull"),
        ],
        frameon=False, loc="lower right", fontsize=AN * 0.8, handletextpad=0.2,
        borderaxespad=0.2, labelspacing=0.2,
    )

    # Zoom on the CAI budget region, marking the widest hull step inside it.
    ins = scat.inset_axes([0.16, 0.46, 0.62, 0.51])
    scatter_front(ins, f, 20, 3.0, 0.45)
    lo_c, hi_c = 0.912, 0.959
    ins.set_xlim(lo_c, hi_c)
    ins.set_ylim(-623, -572)
    hc, hm = f.cai[f.on_hull], f.mfe[f.on_hull]
    inside = np.where((hc[:-1] >= lo_c) & (hc[1:] <= 0.951))[0]
    k = inside[np.argmax(np.diff(hc)[inside])]
    (xl, yl), (xr, yr) = (hc[k], hm[k]), (hc[k + 1], hm[k + 1])
    ya, xb = -613.0, 0.9442
    double_arrow(ins, (xr, ya), (xl, ya))
    for x, y in ((xl, yl), (xr, yr)):
        ins.plot([x, x], [y, ya], color="#888", lw=0.9, ls=(0, (2, 2)), zorder=1)
    ins.text((xl + xr) / 2, ya - 1.2, "$\\Delta$CAI", fontsize=18 * PX, ha="center", va="top")
    double_arrow(ins, (xb, yr), (xb, yl))
    for x, y in ((xl, yl), (xr, yr)):
        ins.plot([x, xb], [y, y], color="#888", lw=0.9, ls=(0, (2, 2)), zorder=1)
    ins.text(xb + 0.0016, (yl + yr) / 2, "$\\Delta$MFE", fontsize=18 * PX, ha="center", va="top",
             rotation=90, rotation_mode="anchor")
    ins.set_xticks([])
    ins.set_yticks([])
    for sp in ins.spines.values():
        sp.set_linewidth(1.2)
        sp.set_visible(True)
    scat.indicate_inset_zoom(ins, edgecolor="0.45", linewidth=1.0, alpha=1.0)

    # Steps between neighbouring designs: CAI below the scatter, MFE left of it.
    h_cai, h_mfe = f.steps(True)
    a_cai, a_mfe = f.steps(False)
    ax = cv.axes(410, 597, 710, Y1)
    ax.set_xscale("log")
    ax.set_xlim(6e-7, 0.2)
    ax.set_ylim(1.6, -0.6)
    for y, v, fc, dc in ((0, h_cai, GRY_FILL, GRY), (1, a_cai, PUR_FILL, PUR)):
        sub = sample(v, rng)
        ax.scatter(sub, y + rng.uniform(-0.3, 0.3, len(sub)), s=6, color=dc, alpha=0.6, lw=0, zorder=1)
        hbox(ax, y, v, fc, dc)
    ax.set_yticks([0, 1])
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=5.5)
    # A tick at every decade, labelled at every other one.
    ax.set_xticks([10.0**k for k in range(-6, 0)])
    ax.set_xticklabels(["$10^{-6}$", "", "$10^{-4}$", "", "$10^{-2}$", ""], fontsize=TK)
    minor_log(ax, "x")
    ax.spines["left"].set_visible(True)
    ax.set_xlabel("$\\Delta$CAI", fontsize=LB, labelpad=2)

    ay = cv.axes(140, Y0, 240, 500)
    ay.set_yscale("log")
    ay.set_ylim(0.04, 80)
    ay.set_xlim(-0.6, 1.6)
    for x, v, fc, dc in ((0, h_mfe, GRY_FILL, GRY), (1, a_mfe, PUR_FILL, PUR)):
        sub = sample(v, rng)
        ay.scatter(x + rng.uniform(-0.3, 0.3, len(sub)), sub, s=6, color=dc, alpha=0.6, lw=0, zorder=1)
        box(ay, x, v, fc, dc, w=0.5)
    ay.set_xticks([0, 1])
    ay.set_xticklabels([])
    ay.tick_params(axis="x", length=5.5)
    ay.set_yticks([0.1, 1, 10])
    ay.set_yticklabels(["0.1", "1", "10"], fontsize=TK)
    minor_log(ay, "y")
    ay.set_ylabel("$\\Delta$MFE (kcal/mol)", fontsize=LB, labelpad=4)

    # Wall time of each method on this instance, as in d.
    bars = cv.axes(140, 535, 380, Y1)
    bars.set_xlim(-0.6, 4.6)
    bars.set_ylim(0, 11500)
    for i, (arm, name) in enumerate(ov.ARMS):
        cell = cells.get((arm, INSTANCE))
        done = ov.finished(cell)
        t = cell["wall_s"]
        bars.bar(i, t, width=0.7, color=METHOD_COLOR[name], alpha=0.75, edgecolor=METHOD_COLOR[name],
                 lw=1.0, hatch=None if done else "////", zorder=2)
        if done:
            bars.text(i, t + 250, f"{round(t):,} s", rotation=90, fontsize=AN, ha="center", va="bottom")
        else:
            # Inside the bar: a label above it would run into the scatter's axis.
            bars.text(i, t / 2, "timeout", rotation=90, fontsize=AN, ha="center", va="center",
                      bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=1.5))
    bars.set_xticks(range(len(ov.ARMS)))
    bars.set_xticklabels([n for _, n in ov.ARMS], rotation=45, ha="right", rotation_mode="anchor", fontsize=TK)
    bars.tick_params(axis="x", pad=4)
    bars.set_yticks([0, 3600, 7200])
    bars.set_yticklabels(["0", "3600", "7200"], fontsize=TK)
    bars.set_ylabel("wall time (s)", fontsize=LB, labelpad=4)
    # The two left plots keep their place; their y labels share one left edge.
    align_ylabels_left(cv.fig, [ay, bars])


def panel_reach(cv, ns, reach, rng):
    cv.letter(755, "(b)")
    ax = cv.axes(867, Y0, 1117, Y1)
    g = np.array([ov.group_of(n) for n in ns])
    for gi in range(len(ov.GROUPS)):
        m = g == gi
        box(ax, gi, reach[m], GRY_FILL, "#555", w=0.5)
        dots(ax, gi, reach[m], ns[m] > 1000, "#555", rng, jit=0.18)
    ax.set_xlim(-0.6, 4.6)
    ax.set_ylim(0, 20)
    ax.set_yticks([0, 5, 10, 15, 20])
    group_ticks(ax)
    ax.set_ylabel("front reachable by $\\lambda$ (%)", fontsize=LB, labelpad=4)
    cv.text(992, XLABEL_Y, "protein length (aa)", fontsize=LB, ha="center", va="center")


def panel_steps(cv, ns, fronts, rng):
    cv.letter(1155, "(c)")
    g = [ov.group_of(n) for n in ns]
    steps = {(hull, gi): ([], []) for hull in (True, False) for gi in range(len(ov.GROUPS))}
    for f, gi in zip(fronts, g, strict=True):
        for hull in (True, False):
            d_cai, d_mfe = f.steps(hull)
            steps[(hull, gi)][0].append(d_cai)
            steps[(hull, gi)][1].append(d_mfe)
    top = cv.axes(1267, Y0, 1557, 405)
    bottom = cv.axes(1267, 440, 1557, Y1)
    for ax in (top, bottom):
        ax.set_yscale("log")
        ax.set_xlim(-0.6, 4.6)
    all_mfe = np.concatenate([np.concatenate(steps[k][1]) for k in steps])
    for which, ax in ((1, top), (0, bottom)):
        for gi in range(len(ov.GROUPS)):
            for hull, side, fc, dc in ((True, -1, GRY_FILL, GRY), (False, 1, PUR_FILL, PUR)):
                v = np.concatenate(steps[(hull, gi)][which])
                x = gi + side * 0.2
                sub = sample(v, rng)
                ax.scatter(x + rng.uniform(-0.15, 0.15, len(sub)), sub, s=5, color=dc, alpha=0.6, lw=0, zorder=1)
                box(ax, x, v, fc, dc, w=0.3)
    top.set_ylim(0.06, all_mfe.max() * 2)
    top.set_xticks(range(len(ov.GROUPS)))
    top.set_xticklabels([])
    top.set_yticks([0.1, 1, 10, 100])
    top.set_yticklabels(["0.1", "1", "10", "100"])
    minor_log(top, "y")
    top.set_ylabel("$\\Delta$MFE (kcal/mol)", fontsize=LB, labelpad=4)
    bottom.set_ylim(6e-7, 0.2)
    group_ticks(bottom)
    bottom.tick_params(axis="y", which="both", width=1.5)
    bottom.set_ylabel("$\\Delta$CAI", fontsize=LB, labelpad=4)
    align_ylabels_left(cv.fig, [top, bottom])
    # The top y label is taller than its axes; centred a little low, it clears the panel letter.
    x, _ = top.yaxis.label.get_position()
    top.yaxis.set_label_coords(x, 0.45)
    cv.text(1412, XLABEL_Y, "protein length (aa)", fontsize=LB, ha="center", va="center")


def panel_wall(cv, ns, proteins, cells, rng):
    cv.letter(1555, "(d)")
    ax = cv.axes(1667, Y0, 1970, Y1)
    ax.set_yscale("log")
    limit = np.array([ov.budget_s(n) for n in ns])
    missed = []
    for i, (arm, name) in enumerate(ov.ARMS):
        cs = [cells.get((arm, p)) for p in proteins]
        done = np.array([ov.finished(c) for c in cs])
        wall = np.array([c["wall_s"] if c else np.nan for c in cs])
        val = np.where(done, wall, limit)
        box(ax, i, val[done], METHOD_COLOR[name], METHOD_COLOR[name], w=0.62)
        dots(ax, i, val, ns > 1000, METHOD_COLOR[name], rng, jit=0.2, s=10)
        missed.append(int((~done).sum()))
    ax.set_xlim(-0.6, 4.6)
    ax.set_ylim(12, 45000)
    # The budget lines are labelled inside the axes, 2 h above its line and 6 h below.
    # The 6 h label clears the timed-out cells' markers drawn on its line.
    for t, label, va, at in ((7200, "2 h", "bottom", 7200), (21600, "6 h", "top", 21600 / 1.2)):
        ax.axhline(t, color="#888", lw=1.0, ls=(0, (3, 3)), zorder=0)
        ax.text(0.985, at, label, transform=ax.get_yaxis_transform(), fontsize=TK, color="#777",
                ha="right", va=va)
    for i, k in enumerate(missed):
        ax.text(i, 33000, f"[{k}]", fontsize=AN, ha="center", va="center")
    ax.set_xticks(range(len(ov.ARMS)))
    ax.set_xticklabels([n for _, n in ov.ARMS], rotation=45, ha="right", rotation_mode="anchor", fontsize=TK)
    ax.set_yticks([1e2, 1e3, 1e4])
    minor_log(ax, "y")
    ax.set_ylabel("wall time (s)", fontsize=LB, labelpad=4)
    print("cells not finished per arm:", dict(zip([n for _, n in ov.ARMS], missed, strict=True)))


def panel_speedup(cv, ns, prune, walk, rng):
    cv.letter(2055, "(e)")
    g = np.array([ov.group_of(n) for n in ns])
    ax = cv.axes(2167, Y0, 2417, Y1)
    ax.set_yscale("log")
    ax.axhline(1, color="#999", lw=1.0, ls=(0, (4, 3)), zorder=0)
    for gi in range(len(ov.GROUPS)):
        m = g == gi
        for v, side, fc, dc in ((prune, -1, GRY_FILL, "#222"), (walk, 1, PUR_FILL, PUR)):
            keep = m & np.isfinite(v)
            if keep.any():
                x = gi + side * 0.2
                box(ax, x, v[keep], fc, dc, w=0.32)
                dots(ax, x, v[keep], ns[keep] > 1000, dc, rng, jit=0.08)
    ax.set_xlim(-0.6, 4.6)
    ax.set_ylim(0.5, 14)
    ax.set_yticks([0.5, 1, 2, 5, 10])
    ax.set_yticklabels(["0.5", "1", "2", "5", "10"])
    minor_log(ax, "y")
    group_ticks(ax)
    ax.set_ylabel("speed-up", fontsize=LB, labelpad=4)
    cv.text(2292, XLABEL_Y, "protein length (aa)", fontsize=LB, ha="center", va="center")
    ax.legend(
        handles=[
            Patch(facecolor=GRY_FILL, edgecolor="#222", alpha=0.7, label="C&P"),
            Patch(facecolor=PUR_FILL, edgecolor=PUR, label="C&P + ADP"),
        ],
        frameon=False, loc="center right", bbox_to_anchor=(1.0, 8.5), bbox_transform=ax.get_yaxis_transform(),
        fontsize=AN, handlelength=1.1, borderaxespad=0.1, labelspacing=0.3,
    )


def main():
    sweep, walk_dir, pin, out = sys.argv[1:]
    cells, points = ov.read_e1(read_jsonl(sweep), pin)
    # Every protein with a finished walk. A method's cell that never ran counts in (d) as not
    # finished, with no point drawn, and leaves (e) without that protein's pruning speed-up.
    proteins = sorted(p for (a, p), c in cells.items() if a == "derna2-walk" and c["exit"] == "ok")
    if INSTANCE not in proteins:
        raise SystemExit(f"{sweep}: no finished walk of {INSTANCE} stamped with pin {pin}")
    ns = np.array([cells[("derna2-walk", p)]["n_aa"] for p in proteins])
    fronts = [ov.Front(points[p]) for p in proteins]

    # The share of the front a λ reaches: the walk's own hull-vertex count over its points.
    reach = np.array(
        [100 * cells[("derna2-walk", p)]["vertices"] / cells[("derna2-walk", p)]["front_points"] for p in proteins]
    )
    # Pruning speed-up: the fixed-λ DP at the hull values over the pruned search, measured where
    # the DP finished; a protein where it stopped at its budget has no pruning speed-up.
    prune, walk = np.full(len(proteins), np.nan), np.full(len(proteins), np.nan)
    for i, p in enumerate(proteins):
        report = Path(walk_dir) / f"{p}.walk.json"
        if not report.exists():
            continue
        sweep_s = json.load(open(report)).get("sweep_s")
        _, w = ov.speedups(cells[("derna2-walk", p)], None, sweep_s)
        walk[i] = np.nan if w is None else w
        dp = ov.hull_dp_s(cells.get(("derna2-hull", p)))
        if dp and sweep_s and not dp[1]:
            prune[i] = dp[0] / sweep_s

    setup_fonts()
    rng = np.random.default_rng(1)
    cv = Canvas()
    legend_row(cv)
    panel_instance(cv, fronts[proteins.index(INSTANCE)], cells, rng)
    panel_reach(cv, ns, reach, rng)
    panel_steps(cv, ns, fronts, rng)
    panel_wall(cv, ns, proteins, cells, rng)
    panel_speedup(cv, ns, prune, walk, rng)
    print(f"gap between panels: {cv.even_out():.1f} design units")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    cv.fig.savefig(out, facecolor="white")
    cv.fig.savefig(Path(out).with_suffix(".png"), dpi=100, facecolor="white")
    print(f"{len(proteins)} proteins; pruning speed-up measured on {int(np.isfinite(prune).sum())}")


if __name__ == "__main__":
    main()
