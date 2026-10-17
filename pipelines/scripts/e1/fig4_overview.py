"""Figure 4 (E1): one front, what a weighted sum reaches, the step between neighbouring designs,
what each method costs, and what the exact front costs DERNA2 against a hull search alone.

    python fig4_overview.py STEPS_NPZ PER_PROTEIN_CSV CELLS_CSV OUT_PDF [PROTEIN]

STEPS_NPZ is fig4_steps.py's output, PER_PROTEIN_CSV fig4_per_protein.py's and CELLS_CSV
fig4_cells.py's. Hull vertices are the sweep's (walk.json), the designs a sweep over lambda
returns. PROTEIN is the protein of (a), P47211 by default; any protein of the corpus works,
since the step (a) marks, its zoom, the inset's place and every axis limit follow its data.

(a) one protein's front, a zoom on one hull step, the |dMFE| and dCAI of all steps between
neighbouring designs of its hull and of its front, and the wall time of each method on it.
(b) share of each front a weight reaches, by length group. (c) the steps of (a), pooled over the
proteins of each length group. (d) wall time per run. (e) the dichotomic search without C&P
pruning (projected) over with it, and the search with it over the complete front.

Panels sit at equal gaps, measured between their drawn extents (labels included), and every
text is checked against every other text and every mark before the file is written.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.axes._base import _TransformedBoundsLocator  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import ConnectionPatch, Patch, Rectangle  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402
from matplotlib.transforms import Bbox  # noqa: E402

FS, FSS = 7, 6.3  # pt: labels, and ticks and legend
PREAMBLE = (
    r"\usepackage[T1]{fontenc}\usepackage{mathptmx}\usepackage{textcomp}"
    r"\usepackage{amsmath}\usepackage{amssymb}"
)
# arm -> colour, in legend order: DERNA2 purple, mRNAfold orange, the UTR-unaware tools green
# and teal, as in Figure 1
PAL = {
    "derna2-walk": "#6a4bc4",
    "mrnafold-hull": "#e8772e",
    "lineardesign-hull": "#2e9e5b",
    "lcds-exact-hull": "#1f9fb4",
    "derna1-own": "#8c6d31",
}
ORDER = list(PAL)
NAMES = {
    "derna2-walk": "DERNA2",
    "mrnafold-hull": "mRNAfold",
    "lineardesign-hull": "LinearDesign",
    "lcds-exact-hull": "LinearCDSfold",
    "derna1-own": "DERNA1",
}
PLATEAUS = np.array([2 * 3600, 6 * 3600])  # s: the budgets a run was granted


def budget(cells):
    """The budget each run was granted. A stopped run sits on the plateau it was given, so read
    it from the run; 984 aa PRT-03 was given 6 h by two arms and 2 h by two more, which no length
    rule reproduces. A run that finished is held to the limit for its length."""
    near = PLATEAUS[np.abs(cells.wall_s.to_numpy()[:, None] - PLATEAUS).argmin(axis=1)]
    by_len = np.where(cells.n_aa < 1000, PLATEAUS[0], PLATEAUS[1])
    return np.where(cells.exit == "timeout", near, by_len)


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
            "xtick.minor.width": 0.4,
            "ytick.minor.width": 0.4,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "xtick.minor.size": 1.5,
            "ytick.minor.size": 1.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


PURPLE, GREY, DARK = PAL["derna2-walk"], "#8a8a8a", "#3a3a3a"
LONG = 1000
GROUPS = [0, 200, 350, 500, 1000, 10**5]
GLABELS = [r"$\le$200", "201--350", "351--500", "501--1000", r"$>$1000"]
TOL_BUDGET = 0.02
FOCUS = "P47211"  # the protein of (a); any protein of the corpus, see main()
# The inset frames one hull step with this much room around it, in points: left (the dMFE arrow
# and its label), right, below, above (the dCAI arrow and its label). The arrows go on the step's
# upper-left sides, which no front point can occupy: points before the step lie below its lower
# end, points after it right of its upper end, and points within it inside its box.
ZOOM_PAD = ((15, 4), (4, 15))
# The inset's place inside the front's axes (axes fraction), tried in order: the first that
# hides no front point and fits its labels wins.
INSET_SPOTS = [
    (0.06, 0.4, 0.6, 0.57),
    (0.05, 0.45, 0.55, 0.52),
    (0.05, 0.5, 0.5, 0.47),
    (0.05, 0.52, 0.47, 0.45),
    (0.04, 0.56, 0.44, 0.42),
    (0.42, 0.05, 0.55, 0.5),
]
BOX_IN = 0.065  # inches: every box in the figure is drawn this wide
EXP_FS = FS + 0.3  # a power-of-ten tick: its exponent is set at 70%
MIN_ZOOM = 0.08  # the zoom rectangle spans at least this share of each axis
FRONT_GAP = 0.06  # inches between (a)'s dMFE box and the front's y label
PAIR_GAP = 0.03  # inches between the two boxes of a length group in (c)
DOTS = 120  # dots drawn per box where a box pools many steps
W, H = 6.5, 2.3  # inches, the paper's text width
GAP, EDGE = 0.07, 0.03  # inches between panels' drawn extents, and at the figure's edge
LETTER_H = 0.13  # inches: the row of panel letters under the legend
LEGEND_GAP = 0.02  # inches between the shared legend and the panels below it


UNITS = {}  # (axes label, vertical box): data units per inch across the box, from a first pass


def bw(ax, vert=True, default=0.4):
    """The data width of a BOX_IN-wide box on this axes; `default` on the measuring pass."""
    u = UNITS.get((ax.get_label(), vert))
    return BOX_IN * u if u else default


def measure_units(axes):
    r = axes["a"]["front"].figure.canvas.get_renderer()
    for spec in axes.values():
        for a in spec.values():
            bb = a.get_window_extent(r)
            UNITS[(a.get_label(), True)] = np.diff(a.get_xlim())[0] / (bb.width / a.figure.dpi)
            UNITS[(a.get_label(), False)] = np.diff(a.get_ylim())[0] / (bb.height / a.figure.dpi)


def box(ax, v, x, color, w=None, vert=True, alpha=0.4, log=True):
    """A box over v. On a log axis the quartiles and the 1.5 IQR whiskers are taken on log10 v,
    so the whiskers are symmetric as drawn; a zero, which a log axis cannot show, is left out."""
    v = np.asarray(v, float)
    w = bw(ax, vert) if w is None else w
    if log:
        v = v[v > 0]
    if len(v):
        st = matplotlib.cbook.boxplot_stats(np.log10(v) if log else v)[0]
        if log:
            st = {
                k: (10**val if k in ("med", "q1", "q3", "whislo", "whishi", "mean") else val)
                for k, val in st.items()
            }
        ax.bxp(
            [st],
            positions=[x],
            widths=w,
            vert=vert,
            showfliers=False,
            patch_artist=True,
            manage_ticks=False,
            boxprops={"facecolor": color, "alpha": alpha, "lw": 0.5, "edgecolor": color},
            medianprops={"color": "#222222", "lw": 0.7, "zorder": 5},  # above the dots
            whiskerprops={"lw": 0.4},
            capprops={"lw": 0.4},
        )


def dots(ax, x, v, color, rng, long=None, half=None, s=1.0, vert=True, cap=None):
    """Jittered dots across a box (over 70% of its width by default); proteins over LONG residues
    as larger triangles."""
    v = np.asarray(v, float)
    half = 0.35 * bw(ax, vert) if half is None else half
    long = np.zeros(len(v), bool) if long is None else np.asarray(long)
    if cap and len(v) > cap:
        keep = rng.choice(len(v), cap, replace=False)
        v, long = v[keep], long[keep]
    j = x + rng.uniform(-half, half, len(v))
    for is_long in (False, True):
        m = long == is_long
        if m.any():
            xy = (j[m], v[m]) if vert else (v[m], j[m])
            ax.scatter(
                *xy,
                s=s * (3.2 if is_long else 1),
                marker="^" if is_long else "o",
                fc=color,
                ec="#222222" if is_long else "none",
                lw=0.3,
                alpha=1 if is_long else 0.75,
                zorder=4 if is_long else 3,
            )


def group_of(n):
    return np.searchsorted(GROUPS, n, side="left") - 1


def log_ticks(ax, axis, values, labels):
    getattr(ax, f"set_{axis}scale")("log")
    getattr(ax, f"set_{axis}ticks")(values)
    getattr(ax, f"set_{axis}ticklabels")(labels)
    getattr(ax, f"{axis}axis").set_minor_locator(matplotlib.ticker.NullLocator())


# ---------------------------------------------------------------------------- panel (a)


def widest_step(pts):
    """The hull step that holds the most other front designs, among steps that are neither
    flat nor steep (slope 1/3 to 3 with each objective as a share of its range on the front), so
    the inset can show both differences: (lower end, upper end)."""
    p = pts[np.argsort(pts[:, 1])]
    idx = np.flatnonzero(p[:, 2] == 1)
    lo, hi = p[idx[:-1]], p[idx[1:]]
    slope = ((hi[:, 0] - lo[:, 0]) / np.ptp(p[:, 0])) / ((hi[:, 1] - lo[:, 1]) / np.ptp(p[:, 1]))
    readable = (slope >= 1 / 3) & (slope <= 3)
    # with no readable step, the one holding the most designs all the same
    k = (
        np.argmax(np.where(readable, np.diff(idx), -1))
        if readable.any()
        else np.argmax(np.diff(idx))
    )
    return p[idx[k]], p[idx[k + 1]]


def decades(lo, hi):
    """Power-of-ten ticks covering [lo, hi], and limits padded a little past them."""
    e = np.arange(np.floor(np.log10(lo)), np.ceil(np.log10(hi)) + 1)
    lim = (lo / 1.6, hi * 1.6)
    inside = e[(10.0**e >= lim[0]) & (10.0**e <= lim[1])]
    if len(inside) < 2:  # widen to the nearest decade on each side
        lim = (
            min(lim[0], 10.0 ** np.floor(np.log10(lo))),
            max(lim[1], 10.0 ** np.ceil(np.log10(hi))),
        )
    return 10.0**e, lim


def tick_labels(ticks):
    """Plain numbers where every tick reads short (0.1 to 100), powers of ten otherwise."""
    if all(0.1 <= t <= 100 for t in ticks):
        return [f"{t:g}" for t in ticks]
    return [rf"$10^{{{int(round(np.log10(t)))}}}$" for t in ticks]


def panel_a(ax, steps, per, focus):
    pts = steps[f"{focus}|pts"]
    hull = pts[:, 2] == 1
    f, d, bars, c = ax["front"], ax["dmfe"], ax["bars"], ax["dcai"]

    def scatter(a, s):
        off = a.scatter(pts[~hull, 1], pts[~hull, 0], s=s, fc="none", ec=PURPLE, lw=0.3, zorder=2)
        on = a.scatter(pts[hull, 1], pts[hull, 0], s=s * 0.9, fc=DARK, ec="none", zorder=3)
        return on, off

    handles = scatter(f, 1.3)
    cai, mfe = pts[:, 1], pts[:, 0]
    cs, ms = np.ptp(cai), np.ptp(mfe)
    f.set_xlim(cai.min() - 0.03 * cs, cai.max() + 0.02 * cs)
    f.set_ylim(mfe.min() - 0.05 * ms, mfe.max() + 0.05 * ms)
    loc = matplotlib.ticker.MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10])
    for lim, set_ticks in ((f.get_xlim(), f.set_xticks), (f.get_ylim(), f.set_yticks)):
        t = loc.tick_values(*lim)
        set_ticks(t[(t >= lim[0]) & (t <= lim[1])])
    f.set_xlabel("CAI", labelpad=2.2)
    f.set_ylabel("MFE (kcal/mol)", labelpad=2.2)
    # the key sits above the axes, in the row beside the panel letter: no front can reach it
    f.legend(
        handles,
        ["on hull", "off hull"],
        markerscale=2.2,
        loc="lower right",
        bbox_to_anchor=(1, 1),
        ncol=2,
        frameon=False,
        fontsize=FSS,
        handletextpad=0.1,
        borderaxespad=0.15,
        columnspacing=0.8,
        handlelength=0.6,
        borderpad=0,
    )

    ins = ax["inset"]
    scatter(ins, 2.2)
    ins.set_xticks([])
    ins.set_yticks([])
    for sp in ins.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.5)
    f._step = widest_step(pts)
    f._drawn = []

    rng = np.random.default_rng(0)
    fs, hs = steps[f"{focus}|front"], steps[f"{focus}|hull"]
    for j, (v, col) in enumerate(((hs, GREY), (fs, PURPLE))):
        box(d, v[:, 0], j, col)
        dots(d, j, v[:, 0], col, rng, s=0.8, cap=2 * DOTS)
        box(c, v[:, 1], 1 - j, col, vert=False)
        dots(c, 1 - j, v[:, 1], col, rng, s=0.8, vert=False, cap=2 * DOTS)
    both = np.vstack([fs, hs])
    for axis, col_i, set_lim in ((d.yaxis, 0, d.set_ylim), (c.xaxis, 1, c.set_xlim)):
        v = both[:, col_i][both[:, col_i] > 0]
        ticks, lim = decades(v.min(), v.max())
        ticks = ticks[(ticks >= lim[0]) & (ticks <= lim[1])]  # thin only what is shown
        log_ticks(
            d if axis is d.yaxis else c,
            "y" if axis is d.yaxis else "x",
            ticks[::2] if len(ticks) > 3 else ticks,
            tick_labels(ticks[::2] if len(ticks) > 3 else ticks),
        )
        set_lim(*lim)
    c.tick_params(axis="x", labelsize=EXP_FS)
    d.set_xlim(-0.6, 1.6)
    d.set_xticks([])
    d.spines["bottom"].set_visible(False)
    d.set_ylabel(r"$\Delta$MFE (kcal/mol)", labelpad=2.2)
    c.set_ylim(-0.6, 1.6)
    c.set_yticks([])
    c.spines["left"].set_visible(False)
    c.set_xlabel(r"$\Delta$CAI", labelpad=2.2)

    r = per.loc[focus]
    limit = 7200 if r.n_aa < LONG else 21600  # the budget the sweeps were held to
    tops = []
    # the hull sweep of DERNA2 (the walk's `sweep_s`) sits next to its complete front
    order = ORDER[:1] + ["derna2-hull"] + ORDER[1:]
    for j, arm in enumerate(order):
        if arm == "derna2-hull":
            ok, h, color = True, r["dichotomic_pruned_min"] * 60, GREY
        else:
            ok = r[f"{arm}_exit"] == "ok"
            h = r[f"{arm}_wall_min"] * 60 if ok else limit
            color = PAL[arm]
        bars.bar(j, h, width=0.72, color=color, alpha=0.5, lw=0)
        tops.append(h)
        lab = f"{h:,.0f}\\,s" if ok else "timeout"
        bars.annotate(
            lab,
            (j, h),
            xytext=(0, 2),
            textcoords="offset points",
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=FSS - 0.3,
        )
    top = max(tops)
    if top >= limit:  # a run reached the budget: the spine ends at it
        ticks = [0, limit]
    else:
        loc = matplotlib.ticker.MaxNLocator(nbins=2, steps=[1, 2, 2.5, 5, 10])
        t = loc.tick_values(0, top)
        ticks = [v for v in t if v <= top] + [min(v for v in t if v >= top)]  # top tick covers it
    bars.set_ylim(0, 2.0 * top)
    bars.set_yticks(ticks)
    bars.set_yticklabels([f"{t:.0f}" for t in ticks])
    bars.spines["left"].set_bounds(0, ticks[-1])
    bars.set_xlim(-0.75, len(order) - 0.4)
    bars.set_xticks([])
    bars.spines["bottom"].set_position(("data", 0))
    bars.set_ylabel(r"wall time (s)", labelpad=2.2)


# ---------------------------------------------------------------------------- panel (b)


def panel_b(ax, steps, per):
    share = pd.Series({p: 100 * steps[f"{p}|pts"][:, 2].mean() for p in per.index})
    g = group_of(per.n_aa.to_numpy())
    long = (per.n_aa > LONG).to_numpy()
    rng = np.random.default_rng(1)
    for j in range(len(GLABELS)):
        m = g == j
        box(ax, share[m].to_numpy(), j, GREY, log=False)
        dots(ax, j, share[m].to_numpy(), DARK, rng, long[m], s=1.4)
    ax.set_ylim(0, 18)
    ax.set_yticks([0, 5, 10, 15])
    ax.set_ylabel(r"front reachable by $\lambda$ (\%)", labelpad=2.2)
    lengths(ax)


def lengths(ax, label=True):
    ax.set_xlim(-0.6, len(GLABELS) - 0.4)
    ax.set_xticks(range(len(GLABELS)))
    if label:
        ax.set_xticklabels(
            GLABELS, rotation=90, ha="right", rotation_mode="anchor", fontsize=FSS - 0.5
        )
        ax.set_xlabel("protein length (aa)", labelpad=2.2)
    else:
        ax.set_xticklabels([])


# ---------------------------------------------------------------------------- panel (c)


def panel_c(top, bot, steps, per):
    g = group_of(per.n_aa.to_numpy())
    rng = np.random.default_rng(2)
    for j in range(len(GLABELS)):
        ps = per.index[g == j]
        for k, (kind, col) in enumerate((("hull", GREY), ("front", PURPLE))):
            v = np.concatenate([steps[f"{p}|{kind}"] for p in ps])
            x = j + (k - 0.5) * (bw(top) + PAIR_GAP * bw(top) / BOX_IN)
            for ax, col_i in ((top, 0), (bot, 1)):
                box(ax, v[:, col_i], x, col)
                dots(ax, x, v[:, col_i], col, rng, s=0.5, cap=DOTS)
    log_ticks(top, "y", [0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
    top.set_ylim(0.06, 150)
    top.set_ylabel(r"$\Delta$MFE (kcal/mol)", labelpad=2.2)
    log_ticks(bot, "y", [1e-6, 1e-4, 1e-2], [r"$10^{-6}$", r"$10^{-4}$", r"$10^{-2}$"])
    bot.set_ylim(2.5e-7, 0.12)
    bot.set_ylabel(r"$\Delta$CAI", labelpad=2.2)
    bot.tick_params(axis="y", labelsize=EXP_FS)
    lengths(top, label=False)
    lengths(bot)
    top._align_with = bot  # y labels share a left edge: see align_left()


# ---------------------------------------------------------------------------- panel (d)


def panel_d(ax, cells):
    c = cells.copy()
    c["lim"] = budget(c)
    c["done"] = (c.exit == "ok") & (c.wall_s <= c.lim * (1 + TOL_BUDGET))
    c["wall"] = np.where(c.done, c.wall_s, c.lim)
    rng = np.random.default_rng(3)
    n = c.protein.nunique()
    for j, arm in enumerate(ORDER):
        g = c[c.arm == arm]
        box(ax, g[g.done].wall.to_numpy(), j, PAL[arm])
        dots(ax, j, g.wall.to_numpy(), PAL[arm], rng, (g.n_aa > LONG).to_numpy(), s=1.2)
        ax.text(j, 6.2e4, f"{n - int(g.done.sum())}", ha="center", va="bottom", fontsize=FSS - 0.3)
    for h in (2, 6):
        ax.axhline(h * 3600, color="#999999", lw=0.5, ls=(0, (2, 2)), zorder=0)
        ax.annotate(
            f"{h}\\,h",
            xy=(1, h * 3600),
            xycoords=ax.get_yaxis_transform(),
            xytext=(2, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=FSS,
            color="#666666",
            annotation_clip=False,
        )
    log_ticks(ax, "y", [1e1, 1e2, 1e3, 1e4], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"])
    ax.set_ylim(8, 1.3e5)
    ax.tick_params(axis="y", labelsize=EXP_FS)
    ax.set_xlim(-0.6, len(ORDER) - 0.4)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(
        [NAMES[a] for a in ORDER],
        rotation=52,
        ha="right",
        rotation_mode="anchor",
        fontsize=FSS - 0.5,
    )
    ax.set_ylabel("wall time (s)", labelpad=2.2)


# ---------------------------------------------------------------------------- panel (e)


def panel_e(ax, per):
    g = group_of(per.n_aa.to_numpy())
    long = (per.n_aa > LONG).to_numpy()
    # measured, no projection: the fixed-lambda DP at the hull values, over the proteins it finishes
    done = (per["derna2-hull_exit"] == "ok").to_numpy()
    prune = (per["derna2-hull_wall_min"] / per.dichotomic_pruned_min).to_numpy()
    prune = np.where(done, prune, np.nan)
    full = (per.dichotomic_pruned_min / per["derna2-walk_wall_min"]).to_numpy()
    rng = np.random.default_rng(4)
    for j in range(len(GLABELS)):
        m = g == j
        for v, col in ((prune, GREY), (full, PURPLE)):
            box(ax, v[m], j, col)
            dots(ax, j, v[m], col if col == PURPLE else DARK, rng, long[m], s=1.2)
    ax.axhline(1, color="#999999", lw=0.5, ls=(0, (2, 2)), zorder=0)
    log_ticks(ax, "y", [0.5, 1, 2, 5], ["0.5", "1", "2", "5"])
    ax.set_ylim(0.36, 9)
    ax.set_ylabel("speed-up", labelpad=2.2)
    ax.text(-0.3, 8.4, "C\\&P pruning", ha="left", va="top", fontsize=FSS, color="#111111")
    ax.text(-0.3, 0.375, "C\\&P + ADP", ha="left", va="bottom", fontsize=FSS, color=PURPLE)
    lengths(ax)


# ---------------------------------------------------------------------------- layout

# Each panel is a set of axes laid out inside the panel's own box in inches: (x, y, w, h),
# with x and w scaled with the panel's width and y and h with the common height.
PANELS = {
    "a": {
        "weight": 2.45,
        "axes": {
            "dmfe": (0.0, 0.38, 0.15, 0.62),
            "front": (0.43, 0.38, 0.57, 0.62),
            "bars": (0.0, 0.0, 0.40, 0.26),
            "dcai": (0.48, 0.0, 0.52, 0.17),
        },
    },
    "b": {"weight": 1.0, "axes": {"b": (0, 0, 1, 1)}},
    "c": {"weight": 1.45, "axes": {"top": (0, 0.53, 1, 0.47), "bot": (0, 0, 1, 0.47)}},
    "d": {"weight": 1.05, "axes": {"d": (0, 0, 1, 1)}},
    "e": {"weight": 1.0, "axes": {"e": (0, 0, 1, 1)}},
}
INSET = (0.07, 0.45, 0.52, 0.5)  # in (a)'s front axes, axes fraction


SHIFT = {}  # (panel, axes) -> left edge as a share of the panel, set by clear_front()
SIZE = {}  # (panel, axes) -> width as a share of the panel, set by clear_front()


def align_left(fig, axes):
    """Give each pair of stacked axes marked with _align_with one edge for their y labels, the
    farther-left of the two. fig.align_ylabels only aligns gridspec subplots. A y label is
    anchored at its right edge (rotated, va="bottom"), so aligning that edge aligns both sides."""
    r = fig.canvas.get_renderer()
    pairs = [
        (a, a._align_with)
        for spec in axes.values()
        for a in spec.values()
        if getattr(a, "_align_with", None) is not None
    ]
    for pair in pairs:
        for ax in pair:  # back to where matplotlib would put it
            ax.yaxis._autolabelpos = True
            ax.yaxis.label.set_transform(
                matplotlib.transforms.blended_transform_factory(
                    matplotlib.transforms.IdentityTransform(), ax.transAxes
                )
            )
    fig.canvas.draw()
    for pair in pairs:
        x = min(ax.yaxis.label.get_window_extent(r).x1 for ax in pair) / fig.dpi / W
        for ax in pair:
            ax.yaxis.set_label_coords(
                x,
                0.5,
                transform=matplotlib.transforms.blended_transform_factory(
                    fig.transFigure, ax.transAxes
                ),
            )


def clear_front(fig, axes, xs, widths):
    """Start (a)'s front and dCAI box FRONT_GAP right of the dMFE box plus the front's own y
    labels, whose width depends on the protein's energies."""
    r = fig.canvas.get_renderer()
    a, f = axes["a"], axes["a"]["front"]
    # the dMFE box gets the room across that the dCAI box gets up and down
    SIZE[("a", "dmfe")] = a["dcai"].get_window_extent(r).height / fig.dpi / widths["a"]
    over = (
        f.get_window_extent(r).x0
        - Bbox.union([f.yaxis.get_tightbbox(r), f.yaxis.label.get_window_extent(r)]).x0
    )
    left = (xs["a"] + SIZE[("a", "dmfe")] * widths["a"]) * fig.dpi + over + FRONT_GAP * fig.dpi
    share = (left / fig.dpi - xs["a"]) / widths["a"]
    for name in ("front", "dcai"):
        SHIFT[("a", name)] = share


def place(fig, axes, xs, widths, bottoms, top):
    """Every panel's axes share one top edge; each panel's bottom is its own, so that every
    panel's drawn extent reaches the figure's foot whatever its tick labels need."""
    for k, spec in PANELS.items():
        height = top - bottoms[k]
        for name, (x, y, w, h) in spec["axes"].items():
            if (k, name) in SHIFT:  # (a)'s right column, moved clear of its y labels
                x, w = SHIFT[(k, name)], 1 - SHIFT[(k, name)]
            w = SIZE.get((k, name), w)
            axes[k][name].set_position(
                (
                    (xs[k] + x * widths[k]) / W,
                    (bottoms[k] + y * height) / H,
                    w * widths[k] / W,
                    h * height / H,
                )
            )


def extents(fig, axes):
    """Each panel's drawn extent in inches, labels and annotations included."""
    r = fig.canvas.get_renderer()
    keys = [a.legend_ for spec in axes.values() for a in spec.values() if a.legend_]
    for key in keys:  # (a)'s key sits in the row of panel labels, not on a panel
        key.set_visible(False)
    out = {}
    for k in PANELS:
        # an axes' tight box keeps only the vertical extent of its x label: add the labels whole
        bb = Bbox.union(
            [
                b
                for a in axes[k].values()
                if a.get_visible()
                for b in (
                    a.get_tightbbox(r),
                    a.xaxis.label.get_window_extent(r),
                    a.yaxis.label.get_window_extent(r),
                )
            ]
        )
        out[k] = bb.transformed(fig.dpi_scale_trans.inverted())
    for key in keys:
        key.set_visible(True)
    return out


def label_rise(fig, ax):
    """Height in inches of the x label's baseline above the bottom of its axes."""
    lab = ax.xaxis.label
    r = fig.canvas.get_renderer()
    ismath = "TeX" if lab.get_usetex() else False
    descent = r.get_text_width_height_descent(lab.get_text(), lab.get_fontproperties(), ismath)[2]
    return (lab.get_window_extent(r).y0 + descent) / fig.dpi - ax.get_position().y0 * H


def layout(fig, axes):
    """Equal gaps between panels' drawn extents, no slack at the figure's edges."""
    ks = list(PANELS)
    total = sum(s["weight"] for s in PANELS.values())
    widths = {k: PANELS[k]["weight"] / total * W * 0.7 for k in ks}
    xs = {k: i * W / len(ks) for i, k in enumerate(ks)}
    r = fig.canvas.get_renderer()
    head = H - fig.legends[0].get_window_extent(r).y0 / fig.dpi + LEGEND_GAP + LETTER_H
    bottoms, top = {k: 0.4 for k in ks}, H - head - 0.1
    for _ in range(5):
        place(fig, axes, xs, widths, bottoms, top)
        clear_front(fig, axes, xs, widths)
        place(fig, axes, xs, widths, bottoms, top)
        align_left(fig, axes)
        ext = extents(fig, axes)
        over_l = {k: xs[k] - ext[k].x0 for k in ks}
        over_r = {k: ext[k].x1 - (xs[k] + widths[k]) for k in ks}
        free = W - 2 * EDGE - GAP * (len(ks) - 1) - sum(over_l.values()) - sum(over_r.values())
        widths = {k: PANELS[k]["weight"] / total * free for k in ks}
        x = EDGE
        for k in ks:
            xs[k] = x + over_l[k]
            x = xs[k] + widths[k] + over_r[k] + GAP
        above = max(ext[k].y1 - top for k in ks)
        bottoms = {k: EDGE + (bottoms[k] - ext[k].y0) for k in ks}
        shared = max(bottoms[k] for k in ks if k != "a")  # (b)-(e): one x-axis line
        bottoms.update({k: shared for k in ks if k != "a"})
        # (a)'s x label shares the baseline of (b)-(e)'s
        bottoms["a"] = shared + label_rise(fig, axes["b"]["b"]) - label_rise(fig, axes["a"]["dcai"])
        top = H - head - above
    place(fig, axes, xs, widths, bottoms, top)
    align_left(fig, axes)


def build(steps, per, cells, focus=FOCUS):
    """Two passes: the first lays the panels out and measures each axes' data units per inch,
    the second redraws with every box BOX_IN wide. Box widths do not move any label, so the
    second layout lands where the first did."""
    UNITS.clear()
    SHIFT.clear()
    SIZE.clear()
    fig, axes = draw(steps, per, cells, focus)
    measure_units(axes)
    plt.close(fig)
    fig, axes = draw(steps, per, cells, focus)
    place_inset(fig, axes["a"]["front"], axes["a"]["inset"])
    ext = extents(fig, axes)
    top = max(e.y1 for e in ext.values())
    for k in PANELS:
        fig.text(
            ext[k].x0 / W,
            (top + 0.035) / H,
            r"\textbf{(" + k + ")}",
            ha="left",
            va="bottom",
            fontsize=FS,
        )
    return fig, axes, ext


def draw(steps, per, cells, focus):
    set_style()
    fig = plt.figure(figsize=(W, H))
    axes = {
        k: {n: fig.add_axes((0, 0, 0.1, 0.1), label=f"{k}.{n}") for n in s["axes"]}
        for k, s in PANELS.items()
    }
    axes["a"]["inset"] = axes["a"]["front"].inset_axes(INSET_SPOTS[0], label="a.inset")
    panel_a(axes["a"], steps, per, focus)
    axes["a"]["dmfe"]._align_with = axes["a"]["bars"]
    panel_b(axes["b"]["b"], steps, per)
    panel_c(axes["c"]["top"], axes["c"]["bot"], steps, per)
    panel_d(axes["d"]["d"], cells)
    panel_e(axes["e"]["e"], per)
    for k in PANELS:
        for a in axes[k].values():
            a.tick_params(pad=1.2)
    legend(fig)
    layout(fig, axes)
    return fig, axes


def outer_corners(f, ins, rect):
    """The two corners (as 0/1 pairs) whose connector, from the zoom rectangle's corner to the
    same corner of the inset, crosses neither box: the outline of the two boxes together."""
    x0, x1, y0, y1 = rect
    rb = Bbox(f.transData.transform([(x0, y0), (x1, y1)]))
    ib = ins.get_window_extent(f.figure.canvas.get_renderer())
    clear = []
    for cx in (0, 1):
        for cy in (0, 1):
            p = np.array([(rb.x0, rb.x1)[cx], (rb.y0, rb.y1)[cy]])
            q = np.array([(ib.x0, ib.x1)[cx], (ib.y0, ib.y1)[cy]])
            seg = MPath(np.linspace(p, q, 50)[1:-1])
            if not (seg.intersects_bbox(rb.padded(-1)) or seg.intersects_bbox(ib.padded(-1))):
                clear.append((cx, cy))
    return clear[:2] if len(clear) >= 2 else [(0, 1), (1, 0)]


def frame_inset(fig, f, ins, min_zoom=None):
    """Zoom the inset onto (a)'s hull step with ZOOM_PAD points around it at the inset's drawn
    size, then draw the zoom rectangle, its connectors, and the two arrows with their labels."""
    for art in f._drawn:
        art.remove()
    a, b = f._step
    bb = ins.get_window_extent(fig.canvas.get_renderer())
    wp, hp = bb.width * 72 / fig.dpi, bb.height * 72 / fig.dpi
    (pl, pr), (pb, pt) = ZOOM_PAD
    dx, dy = b[1] - a[1], b[0] - a[0]
    sx, sy = dx / (wp - pl - pr), dy / (hp - pb - pt)  # data units per point
    x0, x1, y0, y1 = a[1] - pl * sx, b[1] + pr * sx, a[0] - pb * sy, b[0] + pt * sy
    for lo, hi, span, idx in (
        (x0, x1, np.diff(f.get_xlim())[0], 0),
        (y0, y1, np.diff(f.get_ylim())[0], 1),
    ):
        mz = MIN_ZOOM if min_zoom is None else min_zoom
        if hi - lo < mz * span:  # too small to see on the front: widen, step stays centred
            grow = (mz * span - (hi - lo)) / 2
            if idx == 0:
                x0, x1 = lo - grow, hi + grow
            else:
                y0, y1 = lo - grow, hi + grow
    sx, sy = (x1 - x0) / wp, (y1 - y0) / hp  # as drawn
    ins.set_xlim(x0, x1)
    ins.set_ylim(y0, y1)
    drawn = [
        f.add_patch(
            Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, lw=0.8, ec="#333333", zorder=5)
        )
    ]
    f._zoom = drawn[0]
    for cx, cy in outer_corners(f, ins, (x0, x1, y0, y1)):
        drawn.append(
            f.add_artist(
                ConnectionPatch(
                    ((x0, x1)[cx], (y0, y1)[cy]),
                    (cx, cy),
                    "data",
                    "axes fraction",
                    axesA=f,
                    axesB=ins,
                    lw=0.4,
                    color="#777777",
                    zorder=1,
                )
            )
        )
    dash = {"color": "#777777", "lw": 0.4, "ls": (0, (2, 1.5))}
    drawn += ins.plot([a[1], a[1], b[1]], [a[0], b[0], b[0]], **dash)
    arrow = {
        "arrowstyle": "<->",
        "lw": 0.6,
        "color": "#222222",
        "shrinkA": 0,
        "shrinkB": 0,
        "mutation_scale": 5,
    }
    yt = b[0] + 3.5 * sy
    # the arrow stops 1.2 pt short of the upper end, which the front leaves from; the dashed
    # guide below it spans the step exactly
    drawn.append(ins.annotate("", xy=(a[1], yt), xytext=(b[1] - 1.2 * sx, yt), arrowprops=arrow))
    drawn.append(
        ins.annotate(
            r"$\Delta$CAI",
            (b[1], yt),
            xytext=(0, 2.5),
            textcoords="offset points",
            ha="right",
            va="bottom",
            fontsize=FSS,
        )
    )
    xl = a[1] - 3.5 * sx
    drawn.append(ins.annotate("", xy=(xl, a[0]), xytext=(xl, b[0]), arrowprops=arrow))
    drawn.append(
        ins.annotate(
            r"$\Delta$MFE",
            (xl, a[0]),
            xytext=(-2.5, 0),
            textcoords="offset points",
            ha="left",
            va="bottom",
            rotation=90,
            rotation_mode="anchor",
            fontsize=FSS,
        )
    )
    f._drawn = drawn


def place_inset(fig, ax, ins):
    """The inset at the first spot in INSET_SPOTS that covers none of the front's points nor the
    rectangle it magnifies, and where its own labels clear everything (check())."""
    for spot, mz in [(sp, mz) for mz in (MIN_ZOOM, MIN_ZOOM / 2, 0) for sp in INSET_SPOTS]:
        ins.set_axes_locator(_TransformedBoundsLocator(spot, ax.transAxes))
        fig.canvas.draw()
        frame_inset(fig, ax, ins, mz)
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        path = box_path(ins.get_window_extent(r).padded(2))
        pts = np.vstack([c.get_offsets() for c in ax.collections])
        own = [
            Bbox.from_extents(x - 2, y - 2, x + 2, y + 2) for x, y in ax.transData.transform(pts)
        ]
        own.append(ax._zoom.get_window_extent(r))
        if any(path.intersects_path(box_path(m), filled=True) for m in own):
            continue
        bad = [p for p in check(fig) if "Delta$CAI'" in p or "Delta$MFE'" in p]
        if not bad:
            return
    raise SystemExit("no spot in (a) fits the inset and its labels")


def legend(fig):
    hs = [Patch(fc=PURPLE, alpha=0.5, lw=0), Patch(fc=GREY, alpha=0.5, lw=0)]
    labels = ["DERNA2", "DERNA2 (hull)"]
    for arm in ORDER[1:]:
        hs.append(Patch(fc=PAL[arm], alpha=0.5, lw=0))
        labels.append(NAMES[arm])
    hs += [
        Line2D([], [], marker="o", ls="none", ms=2.4, color="#555555"),
        Line2D([], [], marker="^", ls="none", ms=3.2, mfc="#555555", mec="#222222", mew=0.3),
    ]
    labels += [r"$\le$1000 aa", r"$>$1000 aa"]
    fig.legend(
        hs,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=len(hs),
        frameon=False,
        handlelength=0.9,
        handleheight=0.7,
        handletextpad=0.35,
        columnspacing=1.1,
        fontsize=FSS,
        borderaxespad=0.1,
    )


# ---------------------------------------------------------------------------- checks


def text_outline(t, r, air=-1.0):
    """A text's glyph rectangle as a closed path in display space, rotated as drawn. A rotated
    rectangle's bounding box shares its centre, so the outline is that centre plus the unrotated
    rectangle turned by the text's angle."""
    bb = t.get_window_extent(r)
    angle = t.get_rotation()
    t.set_rotation(0)
    flat = t.get_window_extent(r)
    t.set_rotation(angle)
    w, h = max(flat.width / 2 - air, 0.1), max(flat.height / 2 - air, 0.1)
    c, s = np.cos(np.radians(angle)), np.sin(np.radians(angle))
    corners = np.array([(-w, -h), (w, -h), (w, h), (-w, h), (-w, -h)])
    xy = corners @ np.array([[c, s], [-s, c]]) + (bb.x0 + bb.width / 2, bb.y0 + bb.height / 2)
    return MPath(xy, closed=True)


def text_boxes(fig):
    r = fig.canvas.get_renderer()
    out = []
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text().strip() and t.figure is fig:
            try:
                out.append((t, text_outline(t, r)))
            except RuntimeError:
                continue
    for ax in fig.axes:  # an in-panel key is checked whole, markers included
        if ax.legend_:
            out.append((ax.legend_, box_path(ax.legend_.get_window_extent(r))))
    return out


def box_path(b):
    return MPath(
        [(b.x0, b.y0), (b.x1, b.y0), (b.x1, b.y1), (b.x0, b.y1), (b.x0, b.y0)], closed=True
    )


def samples(xy, transform, n):
    """Pixel boxes along a polyline given in its own coordinates."""
    xy = transform.transform(np.asarray(xy, float))
    out = []
    for (x0, y0), (x1, y1) in zip(xy[:-1], xy[1:], strict=True):
        for t in np.linspace(0, 1, n):
            x, y = x0 + t * (x1 - x0), y0 + t * (y1 - y0)
            out.append(Bbox.from_extents(x - 0.5, y - 0.5, x + 0.5, y + 0.5))
    return out


def marks(fig, spines=True):
    """Display-space boxes of every drawn mark inside its axes: dots, line samples, bars, boxes."""
    r = fig.canvas.get_renderer()
    out = []
    for ax in fig.axes:
        clip = ax.get_window_extent(r)
        found = []
        for coll in ax.collections:
            off = np.ma.filled(np.ma.asarray(coll.get_offsets(), float), np.nan)
            off = off[np.isfinite(off).all(axis=1)]
            if len(off):
                rad = np.sqrt(max(coll.get_sizes(), default=1)) / 2 * fig.dpi / 72
                found += [
                    Bbox.from_extents(x - rad, y - rad, x + rad, y + rad)
                    for x, y in ax.transData.transform(off)
                ]
        for ln in ax.lines:
            xy = ln.get_transform().transform(np.column_stack(ln.get_data()))
            for (x0, y0), (x1, y1) in zip(xy[:-1], xy[1:], strict=True):
                for t in np.linspace(0, 1, 16):
                    x, y = x0 + t * (x1 - x0), y0 + t * (y1 - y0)
                    found.append(Bbox.from_extents(x - 0.5, y - 0.5, x + 0.5, y + 0.5))
        for p in ax.patches:
            if isinstance(p, ConnectionPatch):  # a line, not an area
                found += samples(p.get_path().vertices, p.get_transform(), 40)
            elif p.get_visible():
                found.append(p.get_window_extent(r))
        for sp in ax.spines.values() if spines else ():
            if sp.get_visible():
                found += samples(sp.get_path().vertices, sp.get_transform(), 60)
        out += [(ax, b) for b in found if b.overlaps(clip)]
    for p in fig.artists:
        xy = p.get_path().vertices if hasattr(p, "get_path") else None
        if xy is not None:
            xy = p.get_transform().transform(xy)
            for (x0, y0), (x1, y1) in zip(xy[:-1], xy[1:], strict=True):
                for t in np.linspace(0, 1, 40):
                    x, y = x0 + t * (x1 - x0), y0 + t * (y1 - y0)
                    out.append((None, Bbox.from_extents(x - 0.5, y - 0.5, x + 0.5, y + 0.5)))
    return out


def family(ax):
    """An axes with its insets and the axes it is inset into."""
    if ax is None:
        return set()
    out = {ax, *ax.child_axes}
    for other in ax.figure.axes:
        if ax in other.child_axes:
            out |= {other, *other.child_axes}
    return out


def name(t):
    return t.get_text() if hasattr(t, "get_text") else "in-panel key"


def check(fig):
    """Every text and in-panel key against every other, every mark (dots, lines, bars, boxes,
    frames), and every axes it does not belong to."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    tb = text_boxes(fig)
    problems = []
    for i, (t1, p1) in enumerate(tb):
        for t2, p2 in tb[i + 1 :]:
            own = [x for x in (t1, t2) if not hasattr(x, "get_text")]
            if own and {t1, t2} - {own[0]} <= set(own[0].texts):
                continue
            if p1.intersects_path(p2, filled=True):
                problems.append(f"'{name(t1)}' overlaps '{name(t2)}'")
    ms = marks(fig)
    for t, p in tb:
        ext = p.get_extents()
        near = [m for _, m in ms if ext.overlaps(m)]
        if any(p.intersects_path(box_path(m), filled=True) for m in near):
            problems.append(f"'{name(t)}' covers a mark")
        for ax in fig.axes:
            if ax not in family(t.axes) and p.intersects_path(
                box_path(ax.get_window_extent(r)), filled=True
            ):
                problems.append(f"'{name(t)}' enters another panel's axes")
    return sorted(set(problems))


def main():
    steps = np.load(sys.argv[1])
    per = pd.read_csv(sys.argv[2]).set_index("protein")
    cells = pd.read_csv(sys.argv[3])
    out = Path(sys.argv[4])
    focus = sys.argv[5] if len(sys.argv) > 5 else FOCUS
    cells = cells[cells.arm.isin(ORDER) & cells.protein.isin(per.index)]
    assert set(cells.protein) == set(per.index), "a protein of the corpus has no cell"
    assert focus in per.index, f"{focus} is not in the corpus"
    fig, axes, ext = build(steps, per, cells, focus)
    probs = check(fig)
    out.parent.mkdir(parents=True, exist_ok=True)
    # no creation date in either file, so the same inputs give the same bytes
    fig.savefig(out, metadata={"CreationDate": None})
    fig.savefig(out.with_suffix(".png"), dpi=300, metadata={"Software": None})
    ks = list(PANELS)
    gaps = [ext[b].x0 - ext[a].x1 for a, b in zip(ks[:-1], ks[1:], strict=True)]
    print(
        "gaps (in):",
        " ".join(f"{g:.3f}" for g in gaps),
        f"| left {ext[ks[0]].x0:.3f} right {W - ext[ks[-1]].x1:.3f}",
        f"| bottom {min(e.y0 for e in ext.values()):.3f}",
    )
    print(f"{len(probs)} overlap(s)")
    for p in probs:
        print("  ", p)
    if probs:  # the files stay for inspection; the rule fails
        raise SystemExit(1)


if __name__ == "__main__":
    main()
