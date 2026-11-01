"""The motif supplement plate: what a constraint can do to a protein, and what it costs.

    E2_SUPP=... python fig_supp_motif.py WIDTH HEIGHT OUT_PDF

Two rows, nine panels.

  (a)-(c)  one residue window as two columns of the codon lattice, one edge per synonymous
           pair, marked where the pair spells the pattern. The three relations a motif can
           have to a protein. (a) and (b) are the SAME motif on different residues, so the
           relation is a property of the protein and not of the motif.
  (d)      the whole-protein result, over every window and every synonymous choice.
  (e),(f)  how much of each set survives on the front as levels are added.
  (g)-(i)  what reaching each front costs against the unconstrained run.

The lattices are drawn from the standard genetic code written into supp_common, not from a
measurement, and each asserts its own verdict against the edges it drew: a wrong codon
table fails the build rather than printing a plausible picture. The jitter is seeded per
panel, so two runs over one input give one figure.
"""

import datetime
import json
import math
import os
import random
import re
import statistics as st
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import supp_common as common  # noqa: E402
import supp_gate as gate  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

WIDTH = float(sys.argv[1]) if len(sys.argv) > 1 else 7.0
HEIGHT = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
DST = sys.argv[3] if len(sys.argv) > 3 else f"{common.OUT}/fig_supp_motif.pdf"
STEM = DST[:-4] if DST.endswith(".pdf") else DST

SETC = {"I": "#2a78d6", "II": "#d1651b", "III": "#8f4fd1"}
LVLC = {"F1": SETC["I"], "F2": SETC["II"], "F3": SETC["III"]}
HIT, CLEAR = "#c0392b", "#2f6f4f"
INK, MUTED, GRID = "#0b0b0b", "#4a4742", "#e4e3df"
SETS, LEVELS, NESTS = ["I", "II", "III"], common.LEVELS, common.NESTS
# F0..F3 are the grid's internal names for the cells and mean nothing to a reader. A level
# is named by what it forbids, in the notation the main plate and the text now use. On an
# axis the sets are named cumulatively, because the full union spells to 150 px against the
# 58 px of tick pitch a cost panel has; the key below carries the unions in full.
RUNGTICK = {"F0": r"$\emptyset$", "F1": r"$\mathcal{M}_\mathrm{I}$",
            "F2": r"$\cup\,\mathcal{M}_\mathrm{II}$",
            "F3": r"$\cup\,\mathcal{M}_\mathrm{III}$"}
# The level at which each set is first forbidden; at and after it the measure is zero.
FORBID = {"I": 1, "II": 2, "III": 3}
MONO = {"family": "DejaVu Sans Mono"}
plt.rcParams.update({"font.size": 5.2, "axes.labelsize": 5.2, "xtick.labelsize": 4.8,
                     "ytick.labelsize": 4.6, "axes.edgecolor": MUTED,
                     "axes.linewidth": .5, "xtick.color": MUTED, "ytick.color": MUTED,
                     "text.color": INK, "axes.labelcolor": INK, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False,
                     # The layout is measured off the drawn figure, so the dpi it is drawn
                     # at is part of the figure's definition and not the host's to choose.
                     "figure.dpi": 200,
                     # No tight box: it crops to the ink and then adds its own pad, so a
                     # figure laid out at 6.5 x 3.0 is saved at whatever the longest tick
                     # label left over. Every artist is placed inside the canvas instead,
                     # which the gate verifies.
                     "savefig.bbox": None})


def tint(h, k):
    r, g, b = (int(h[i:i + 2], 16) for i in (1, 3, 5))
    return "#{:02x}{:02x}{:02x}".format(*(round(c + (255 - c) * k) for c in (r, g, b)))


def lattice_rows(n):
    """The y of each codon in a column, over one span whatever the column's length.

    Centres of n equal slices would give a two-codon column two thirds of the reach of a
    four-codon one, and the three lattices are one series: (a) and (b) drew in 167 px of the
    233 px (c) used, which reads as two small diagrams beside a big one.
    """
    top, bot = .875, .125
    return [(top + bot) / 2] if n == 1 else \
        [top - i * (top - bot) / (n - 1) for i in range(n)]


# Panel letters are deferred until the figure has been drawn once, because a row shares one
# baseline and that baseline is the highest panel in the row.
PENDING = []
LETFS = 6.2
# How far a panel letter stands above its row's tallest panel, in figure fractions.
TAG_RISE = .020


def tag(ax, letter, row):
    """Queue a panel letter. `row` names the band whose letters share a baseline."""
    PENDING.append((ax, letter, row))


def place_tags(fig):
    """Draw the deferred letters: one baseline per row, each at its panel's left edge.

    The left edge is everything the panel DRAWS, its axis label and tick row included, so a
    letter stands at the edge a reader takes for the panel's rather than inset over the plot
    area by whatever that panel's y label happens to be wide.

    Panels that share a gridspec column WITHIN ONE BAND share the leftmost of those edges.
    Two stacked panels have y labels of different widths -- "share" against "sites per
    design" -- and left to themselves their letters stagger by that difference. Sharing the
    edge also keeps (e) and (f) off each other: a letter set left of its own y label sits
    under the y label of the panel above, which is how (f) once came to print on (e)'s.

    The band matters. Both bands begin at the figure's left margin, so (a) shares a column
    x with (e) and (f); pooled across bands it inherited their y label's edge and stood 115
    px left of its own panel, while (b) and (c) stood on theirs -- the ragged row this rule
    exists to prevent, arriving by another route.

    Returns {axes: its letter}, which the overlap gate needs so that a letter is not read as
    intruding on the very panel it names, and the lattices need so their headings can be set
    a fixed distance from theirs.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    rows, cols = {}, {}
    for ax, letter, row in PENDING:
        col = (row[:3], round(ax.get_position().x0, 4))
        cols.setdefault(col, []).append(ax.get_tightbbox(renderer).transformed(inv).x0)
        rows.setdefault(row, []).append((ax.get_position().y1, col, ax, letter))
    left = {col: min(xs) for col, xs in cols.items()}
    placed = {}
    for group in rows.values():
        y = max(y1 for y1, _, _, _ in group)
        for _, col, ax, letter in group:
            placed[ax] = fig.text(left[col], y + TAG_RISE, letter, fontsize=LETFS,
                                  color=INK, ha="left", va="bottom", fontweight="bold")
    return placed


ADM = json.load(open(f"{common.OUT}/admits.json"))
PRES = json.load(open(f"{common.OUT}/presence.json"))["rows"]
COST = json.load(open(f"{common.OUT}/cost.json"))["rows"]

fig = plt.figure(figsize=(WIDTH, HEIGHT))
# The gap between the rows is sized by (d)'s rotated tick labels, which are the deepest
# thing on the top row by far. Under (a)-(c) that same gap is empty, so the three lattices
# are grown down into it once the geometry is known (LATTICE_DROP below) rather than left
# floating above a band of white a third of their own height.
outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.22], hspace=.73)
# The fourth column is a spacer. The three lattices are one group and (d) is a different
# kind of panel, so the gap to it is set to read as a boundary: 98 px of white against
# the 65 px between lattices, a ratio of 1.5. At the gridspec's own spacing it came to
# 26 px -- less than half the internal gap -- because (d)'s three-line y label reaches
# 140 px into whatever the spacer leaves.
top = outer[0].subgridspec(1, 5, width_ratios=[1.06, 1.06, 1.06, .46, 2.16], wspace=.16)
# (e),(f) take the width of the lower band. They carry four level groups of three boxes
# each against the cost panels' three boxes, so the content per inch is four times as
# dense, and a cost panel needs little: at the width this leaves it, its widest tick label
# is 47.6 px against 57.5 px of pitch, adjacent labels clear each other by 21.8 px, and its
# boxes are 25 px apart. The stack is four and a half times a cost panel's width: 776 px
# against 173. The gutters between the four columns do not move with this split -- wspace
# is a share of the band's AVERAGE panel width, which the ratios leave alone -- so they
# stay at 57, 59 and 59 px however the width is divided.
bot = outer[1].subgridspec(2, 4, width_ratios=[4.5, 1, 1, 1], wspace=.52, hspace=.54)

# ---- (a)-(c) the three relations a motif can have to a protein --------------------------
AAN = {"K": "Lys", "E": "Glu", "W": "Trp", "G": "Gly"}
# (window, pattern, heading, consequence, expected verdict). The verdict is asserted, so
# these four tuples are a test of the codon table as well as the figure's content.
# "Admissible" is the main plate's verb: panel (a) there counts the proteins whose coding
# sequence can spell a pattern at all. An unavoidable motif is admissible too, so (c) is a
# sub-case of (b) read strictly; the headings name what distinguishes each column and the
# caption states the containment rather than the headings pretending to a partition.
CASES = [("EE", "AAAAAA", "inadmissible", "removes no design", "none"),
         ("KK", "AAAAAA", "admissible", "removes designs", "some"),
         ("WG", "GGGG", "unavoidable", "no design avoids it", "all")]
LATTICES, MOTIF_LINES = [], []
for k, (win, pat, head, consequence, expect) in enumerate(CASES):
    ax = fig.add_subplot(top[0, k])
    LATTICES.append(ax)
    a, b = win
    left, right = common.CODONS[a], common.CODONS[b]
    rx = re.compile("".join(common.IUPAC.get(c, c) for c in pat))
    yl, yr = lattice_rows(len(left)), lattice_rows(len(right))
    marked = 0
    for i, ca in enumerate(left):
        for j, cb in enumerate(right):
            hit = bool(rx.search(ca + cb))
            marked += hit
            ax.plot([.17, .83], [yl[i], yr[j]], color=HIT if hit else CLEAR,
                    lw=1.35 if hit else .7, alpha=.95 if hit else .5,
                    zorder=2 if hit else 1, solid_capstyle="round")
    for xs, ys, cs in ((0.0, yl, left), (1.0, yr, right)):
        for y, c in zip(ys, cs, strict=True):
            ax.text(xs, y, c, ha="center", va="center", fontsize=6.0, fontdict=MONO,
                    bbox=dict(boxstyle="round,pad=.24", fc="white", ec=MUTED, lw=.45),
                    zorder=3)
    total = len(left) * len(right)
    got = "none" if marked == 0 else ("all" if marked == total else "some")
    assert got == expect, f"{win} x {pat}: {marked}/{total} is {got}, not {expect}"
    col = {"none": CLEAR, "some": INK, "all": HIT}[got]
    ax.set_title(head, fontsize=6.4, color=col, pad=6, fontweight="bold")
    # The motif line sits just INSIDE the panel. Outside, level with the heading, it shares
    # a band with the panel letter and the two came within a pixel of each other: the
    # letter is at the panel's left edge and the line is centred on a panel this narrow.
    MOTIF_LINES.append(
        ax.text(.5, .995, f"{pat} in {AAN[a]}-{AAN[b]}", transform=ax.transAxes,
                ha="center", va="top", fontsize=5.1, fontdict=MONO, color=MUTED))
    # The two readings under a lattice are offset in POINTS from one anchor, not in axes
    # fractions: the panel's height is decided after it is drawn, and in fractions the line
    # spacing would grow with it and the pair would read as two unrelated captions.
    for dy, txt, c in ((0, f"{marked} of {total} pairs spell it", col),
                       (-6.6, consequence, MUTED)):
        ax.annotate(txt, (.5, .015), xycoords="axes fraction", xytext=(0, dy),
                    textcoords="offset points", ha="center", va="top", fontsize=5.1,
                    color=c, gid=f"lattice-caption-{k}")
    ax.set_xlim(-.17, 1.17)
    # The span is the 1.14 the rows are laid out over. Where it starts decides how high the
    # lattice hangs under its heading: the two lines above it are the panel's header and
    # read as one block with the letter, so a lattice level with them reads as part of that
    # block rather than as the panel's subject. The 1.1% offset is what keeps the DRAWN
    # border of the lowest codon box 3.1 pt clear of the line of type under it -- the gate
    # measures a label's text, not the rounded box behind it, so that pair is one it cannot
    # see and it was down to 2.4 pt.
    ax.set_ylim(-.011, 1.129)
    ax.axis("off")
    tag(ax, f"({'abc'[k]})", "top")

# ---- (d) proteins that cannot spell each motif ------------------------------------------
AXD = ax = fig.add_subplot(top[0, 4])
NP = len(ADM["rows"])
order = sorted((tuple(sm) for sm in ADM["motifs"]),
               key=lambda sm: (-sum(1 for r in ADM["rows"] if not r["admits"][sm[1]]),
                               sm[1]))
cannot = [sum(1 for r in ADM["rows"] if not r["admits"][m]) for _, m in order]
ax.bar(range(len(order)), cannot, color=[SETC[s] for s, _ in order], width=.74, zorder=3)
# No value over each bar. They were set at 2.9 pt, which is half any journal's floor and so
# is not type at all, and every one of the 24 counts is printed at full size in the table
# the caption sends the reader to. Without them the bars take the panel's whole height.
ax.set_xticks(range(len(order)))
NAMES = common.names()
# 4.8 pt, the size of every other x tick row on the plate. A rotated label's footprint
# across the axis is its HEIGHT, so 24 of them are the densest type here: at 5.0 the
# closest pair cleared by 1.68 pt, at 4.8 by 2.16, and 4.6 buys nothing further.
ax.set_xticklabels([NAMES.get(m, m) for _, m in order], rotation=90, fontsize=4.8)
ax.set_ylabel(f"# proteins\ninadmissible\n(of {NP})")
ax.set_ylim(0, max(cannot) * 1.06)
ax.grid(True, axis="y", color=GRID, lw=.4)
ax.set_axisbelow(True)
ax.tick_params(axis="x", length=0, pad=1)
tag(ax, "(d)", "top")

# ---- (e),(f) how much of each set survives ----------------------------------------------
PROWS = [("share", "share", "(e)"), ("mean", "sites per\ndesign", "(f)")]
PRESAX = []
for i, (key, ylab, letter) in enumerate(PROWS):
    ax = fig.add_subplot(bot[i, 0])
    PRESAX.append(ax)
    ax.grid(True, axis="y", color=GRID, lw=.4)
    ax.set_axisbelow(True)
    for k, n in enumerate(LEVELS):
        for j, s in enumerate(SETS):
            v = [r[key][s] for r in PRES if r["level"] == n and r[key][s] is not None]
            forbidden = k >= FORBID[s]
            col = MUTED if forbidden else SETC[s]
            pos = k + (j - 1) * .30
            bp = ax.boxplot([v], positions=[pos], widths=.22, patch_artist=True,
                            showfliers=False, medianprops=dict(color=col, lw=.9),
                            whiskerprops=dict(color=col, lw=.5),
                            capprops=dict(color=col, lw=.5))
            bp["boxes"][0].set_facecolor(tint(col, .82 if forbidden else .60))
            bp["boxes"][0].set_edgecolor(col)
            bp["boxes"][0].set_linewidth(.55)
            rng = random.Random(f"{key}{s}{n}")
            ax.scatter([pos + rng.uniform(-.065, .065) for _ in v], v, s=1.3, lw=.15,
                       color=col, edgecolor="white", alpha=.85, zorder=4)
        if k:
            ax.axvline(k - .5, color=GRID, lw=.5, zorder=0)
    ax.set_xlim(-.52, len(LEVELS) - .48)
    ax.set_xticks(range(len(LEVELS)))
    # Only the lower panel names the levels: at this height two tick rows collide.
    ax.set_xticklabels([RUNGTICK[n] for n in LEVELS] if i else [""] * len(LEVELS))
    ax.set_ylabel(ylab)
    # A level that forbids a set draws it at exactly zero, and that row is the panel's
    # check on the solver. At 3.5% of the panel the grey boxes sat ON the spine and read as
    # a smudge rather than as three boxes at zero.
    ax.set_ylim(bottom=-.11 * max(ax.get_ylim()[1], 1e-9))
    ax.set_ylim(top=ax.get_ylim()[1] * 1.10)
    # Three ticks, so the medians the panel exists to show can be read off the axis; two
    # left every value between the floor and the ceiling unlabelled.
    ax.yaxis.set_major_locator(MaxNLocator(3))
    if i:
        ax.set_xlabel("motifs forbidden")
    tag(ax, letter, f"bot{i}")

# ---- (g)-(i) what reaching each front costs ---------------------------------------------
COSTP = [("wall_ratio", r"wall time, $\emptyset$ = 1", [1, 10, 100], "(g)"),
         ("rss_ratio", r"peak memory, $\emptyset$ = 1", [1, 2, 4, 8], "(h)"),
         ("size_ratio", r"front size, $\emptyset$ = 1", [.25, .5, 1, 2], "(i)")]
COSTAX = []
for k, (key, ylab, ticks, letter) in enumerate(COSTP):
    ax = fig.add_subplot(bot[:, k + 1])
    COSTAX.append(ax)
    ax.grid(True, axis="y", color=GRID, lw=.4)
    ax.set_axisbelow(True)
    counts = []
    for j, n in enumerate(NESTS):
        v = [r[key] for r in COST if r["level"] == n and r[key] is not None]
        col = LVLC[n]
        bp = ax.boxplot([v], positions=[j], widths=.56, patch_artist=True,
                        showfliers=False, medianprops=dict(color=col, lw=1.0),
                        whiskerprops=dict(color=col, lw=.55),
                        capprops=dict(color=col, lw=.55))
        bp["boxes"][0].set_facecolor(tint(col, .60))
        bp["boxes"][0].set_edgecolor(col)
        bp["boxes"][0].set_linewidth(.6)
        rng = random.Random(f"{key}{n}")
        ax.scatter([j + rng.uniform(-.16, .16) for _ in v], v, s=1.9, lw=.18, color=col,
                   edgecolor="white", alpha=.85, zorder=4)
        counts.append(len(v))
    # The n belongs to the panel: (g) and (h) hold only the cells with a timed round and
    # (i) holds every protein, so one n for the plate would be false of two panels. Per
    # LEVEL only where the levels differ -- where they do not, three copies of one number
    # say nothing three times, and at this panel width they collide.
    if len(set(counts)) == 1:
        ax.annotate(f"n={counts[0]}", (.5, 1), xycoords="axes fraction",
                    ha="center", va="top", fontsize=4.6, color=MUTED)
    else:
        for j, c in enumerate(counts):
            ax.annotate(f"n={c}", (j, 1), xycoords=("data", "axes fraction"),
                        ha="center", va="top", fontsize=4.6, color=MUTED)
    ax.axhline(1, color=MUTED, lw=.55)
    ax.set_yscale("log")
    ax.minorticks_off()
    lo, hi = ax.get_ylim()
    keep = [t for t in ticks if lo <= t <= hi]
    ax.set_yticks(keep)
    ax.set_yticklabels([f"{t:g}" for t in keep])
    # Headroom for the n= label as a share of the DECADES drawn, not as a factor on the
    # value. A factor adds a fixed number of decades, and these three panels span 3.1, 1.1
    # and 1.0 of them, so one constant bought 9% of (g) and 27% of (i) and the trio read as
    # three unrelated plots rather than one series.
    dec = math.log10(hi / lo)
    ax.set_ylim(top=10 ** (math.log10(hi) + .10 * dec))
    # And slack under the lowest point, so the rule drawn at 1 does not land on the spine
    # and read as one thick double line. (h)'s lowest ratio is 0.91.
    ax.set_ylim(bottom=10 ** (math.log10(lo) - .06 * dec))
    ax.set_xticks(range(len(NESTS)))
    ax.set_xticklabels([RUNGTICK[n] for n in NESTS])
    # Symmetric about the boxes: -.6 either side, not -.6 and len-.6. The second form
    # put the three boxes 0.1 of a unit right of the panel's midline -- 5.8 px here --
    # which showed up as the last tick label ending flush with the right spine while
    # the first stood 33 px clear of the left one. (e),(f) always had the right form.
    ax.set_xlim(-.6, len(NESTS) - 1 + .6)
    ax.set_ylabel(ylab)
    ax.set_xlabel("motifs forbidden")
    tag(ax, letter, "bot0")

# One centred key across the top, read before the panels rather than after them.
# A hue names a motif SET in
# (d)-(f) and the constraint LEVEL that newly forbids that set in (g)-(i), so each entry
# carries both readings joined by a middot -- one swatch, two meanings, which a key naming
# only the first would resolve silently in favour of the first.
# "no sets" is deliberately absent: the main plate needs it because it draws unconstrained
# points in that grey, and here no mark uses it. The grey below is the level's own set,
# which is a different statement and gets its own entry.
RUNG = [r"$\mathcal{M}_\mathrm{I}$",
        r"$\mathcal{M}_\mathrm{I}\cup\mathcal{M}_\mathrm{II}$",
        r"$\mathcal{M}_\mathrm{I}\cup\mathcal{M}_\mathrm{II}\cup\mathcal{M}_\mathrm{III}$"]
HAND = [plt.Rectangle((0, 0), 1, 1, fc=tint(SETC[s], .60), ec=SETC[s], lw=.55)
        for s in SETS]
HAND += [plt.Rectangle((0, 0), 1, 1, fc=tint(MUTED, .82), ec=MUTED, lw=.55),
         Line2D([], [], color=HIT, lw=1.0),
         Line2D([], [], color=CLEAR, lw=.55, alpha=.5)]
# The main plate's key names the levels, so this one does too. A hue still means the SET in
# (d)-(f) and the LEVEL that newly forbids it in (g)-(i); naming both in each entry doubled
# the key's width, and the caption carries the second reading instead.
LAB = list(RUNG)
LAB += ["forbidden at this level", "spells the pattern", "does not"]
assert len(LAB) == len(HAND), f"{len(LAB)} labels for {len(HAND)} handles"
fig.legend(handles=HAND, labels=LAB, loc="upper center", bbox_to_anchor=(.5, .9955),
           frameon=False, fontsize=4.6, handlelength=1.4, labelspacing=.3, borderpad=0,
           ncol=len(HAND), columnspacing=1.6)

# bottom .112 rather than .105: a text's box carries the font's descent below its ink,
# and the gate measures boxes, so the x labels' boxes came within 2.3 pt of the page
# edge where their ink was 3.4 pt clear.
fig.subplots_adjust(left=.0677, right=.990, top=.890, bottom=.126)
fig.canvas.draw()
INV = fig.transFigure.inverted()


def ink_box(artist):
    return artist.get_window_extent().transformed(INV)


def tick_inset(ax):
    """How far the y tick row reaches left of the panel, as a share of its width."""
    p = ax.get_position()
    xs = [ink_box(t).x0 for t in ax.get_yticklabels() if t.get_text().strip()]
    return (min(xs) - p.x0) / p.width


# (g)-(i) hang their y labels at one inset, taken from whichever of the three has the
# widest tick row. Left to themselves the three labels sit at three insets -- "0.25" is
# wider than "1" -- and the ink gutters between the four bottom panels came out 61, 86 and
# 55 px, so the widest white in the row fell INSIDE the cost trio while the boundary
# between it and (e),(f), which is a different claim, got less.
# The gap between the tick row and the label is a distance on the page, not a share of the
# panel. Written as a fraction of the width it shrinks when the panel does: at 216 px it
# came to 11.9 px where the gate asks 12.1, and the plate refused to save.
LAB_GAP_IN = .055
LAB_X = min(tick_inset(ax) - LAB_GAP_IN / (ax.get_position().width * WIDTH)
            for ax in COSTAX)
for ax in COSTAX:
    ax.yaxis.set_label_coords(LAB_X, .5)

# (e) and (f) are deliberately NOT given a shared label inset the way the cost trio is.
# Each sits where matplotlib puts it, a fixed distance from its own tick row, and those
# rows differ in width -- "0.0/0.5/1.0" against "0/8/16" -- so the two labels' inner edges
# stand 10 px (2.4 pt) apart. Both ways of closing that were measured and are worse. Align
# on the innermost and "share" lands 1.7 pt from its own tick labels, which the gate
# refuses. Align on the outermost, where each clears its own ticks, and the two-line label
# moves further left -- and the top band is placed from the bottom band's leftmost ink, so
# it carries (a),(b),(c) off the page. The stagger is the cheapest of the three.

# The top band starts where the bottom band's ink starts. Both are laid out from one
# gridspec margin, but the lower panels carry y labels that reach left of it and the upper
# ones carry none, so the top row began 115 px right of the figure's own left edge over a
# strip that was white. The three lattices move onto it at their own width and (d) takes
# what they leave: its 24 rotated tick labels are the tightest type on the plate, and the
# gap between (c) and (d)'s label is carried along rather than spent.
RENDERER = fig.canvas.get_renderer()
STRIP = (LATTICES[0].get_position().x0
         - min(ax.get_tightbbox(RENDERER).transformed(INV).x0 for ax in PRESAX + COSTAX))
for ax in LATTICES:
    p = ax.get_position()
    ax.set_position([p.x0 - STRIP, p.y0, p.width, p.height])
p = AXD.get_position()
AXD.set_position([p.x0 - STRIP, p.y0, p.width + STRIP, p.height])

# (a)-(c) grow down into the band between the rows. That band is sized by (d)'s rotated
# tick labels, the deepest ink on the top row; under the lattices it is white, and at this
# figure size it was a third of a lattice's own height. The floor is where the lower row's
# panel letters will stand, taken from the type size because those letters are not placed
# until the lattices have stopped moving and the letters' own left edges depend on where
# every panel ends up. Measured against the axes line instead, it put a lattice's last line
# of type on (e). Done after subplots_adjust, which re-places every gridspec axes.
LATTICE_FLOOR = (outer[1].get_position(fig).y1 + TAG_RISE
                 + LETFS * 1.1 / 72 / HEIGHT + .085 / HEIGHT)
for ax in LATTICES:
    p = ax.get_position()
    ink = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(INV)
    drop = ink.y0 - LATTICE_FLOOR
    ax.set_position([p.x0, p.y0 - drop, p.width, p.height + drop])

LETTERS = place_tags(fig)
OWNERS = {id(t): f"{id(ax)}" for ax, t in LETTERS.items()}

# The lattice headings start one fixed gap right of their own panel letter, rather than
# centred on the panel. Centred, the gap between letter and heading is set by the width of
# whichever word the heading is -- "inadmissible", "admissible", "unavoidable" ran 18.6,
# 31.8 and 23.5 px apart at 300 dpi, a spread of 70% across three panels that are otherwise
# identical, and the row read as badly spaced. Measured off the letters rather than guessed,
# because the letter's own width differs by three pixels between panels.
HEAD_GAP = .062 / WIDTH
for ax, line in zip(LATTICES, MOTIF_LINES, strict=True):
    p = ax.get_position()
    x = (ink_box(LETTERS[ax]).x1 + HEAD_GAP - p.x0) / p.width
    ax.title.set_ha("left")
    ax.title.set_x(x)
    line.set_ha("left")
    line.set_x(x)

# The plate fails rather than ships a collision, as the composite plate does.
gate.check(fig, OWNERS)
os.makedirs(os.path.dirname(os.path.abspath(STEM)), exist_ok=True)
# Byte-stable, as the repo's other figures are written: matplotlib stamps the wall clock
# into a PDF's CreationDate and its own version into a PNG's Software chunk, so two runs
# over one input gave two digests and a manifest over them could only ever fail.
fig.savefig(f"{STEM}.pdf", metadata={"CreationDate": None})
fig.savefig(f"{STEM}.png", dpi=300, metadata={"Software": None})
plt.close(fig)

REPORT = {}
for key, _, _ in PROWS:
    for s in SETS:
        REPORT[f"{key}.{s}"] = {n: round(st.median([r[key][s] for r in PRES
                                                    if r["level"] == n]), 6)
                                for n in LEVELS[:FORBID[s]]}
for key, _, _, _ in COSTP:
    REPORT[key] = {n: [round(st.median([r[key] for r in COST
                                        if r["level"] == n and r[key] is not None]), 4),
                       sum(1 for r in COST if r["level"] == n and r[key] is not None)]
                   for n in NESTS}
REPORT["cannot_spell"] = {m: v for (_, m), v in zip(order, cannot, strict=True)}
json.dump({
    "figure": os.path.basename(STEM),
    "written_utc": datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"),
    "size_in": [WIDTH, HEIGHT],
    "python": sys.version.split()[0],
    "matplotlib": matplotlib.__version__,
    # names.json decides panel (d)'s tick labels and is read straight from inputs/, so it
    # is an input of the drawing step and not of any cache.
    "inputs": {f: common.sha(f"{common.OUT}/{f}")
               for f in ("admits.json", "presence.json", "cost.json")}
    | common.digest_of(common.NAMES_JSON),
    # The gate is digested with the drawing code: it decides whether the figure is
    # written at all, so a change to it is a change to what this file can publish.
    "scripts": common.digest_of(os.path.abspath(__file__), common.__file__,
                                     gate.__file__),
    "report": REPORT,
}, open(f"{STEM}.prov.json", "w"), indent=1, sort_keys=True)
print(json.dumps(REPORT, indent=1))
print(f"wrote {STEM}.pdf  ({WIDTH} x {HEIGHT} in)")
