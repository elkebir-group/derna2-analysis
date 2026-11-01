"""Figure 2 as a single row of six panels, to the original specification.

  (a) the motifs and their categorisation, with how often each occurs across the screen
  (b) one instance's fronts, filled on the convex hull and open off it, with a zoom window
  (c)-(f) that instance's design at each nest, labelled with its MFE and CAI, with motif
      sites coloured by category

Colour means a category. In (a) and (c)-(f) it means the category *occurs* there; in (b) it
names the nest, and nest Fn is drawn in the colour of the category it newly *forbids* -- the
opposite predicate on the same hue. Nothing on a page can make that safe by itself, so the
four fronts are additionally named in the zoom, which is the only place they separate, and
the caption has to state the polarity. F0 forbids nothing and takes the one neutral.

Several choices here were measured rather than taken by eye, and each is cheap to get wrong:

  The palette. Blue / light orange / dark magenta are three hues at three distinct
  *lightnesses*, which is the only thing that survives a greyscale print; their closest pair
  is 18.2 apart in L*. The blue/purple pair they replace separated by 1.0 for a deuteranope.
  The gate is run over the full ink set and not only the three hues, and the pair it kept
  finding was category III against the axis grey at 0.17 -- twelve times closer than anything
  else in the figure. Moving the terminus marks and the CAI-floor line to body ink retired
  two *marks* and left the grey itself on every spine, every tick and the two hull legend
  swatches, so the pair was still live while this paragraph said it was not. The grey is now
  15.5 from category III and 13.7 from body ink, and the closest pair in the figure is the
  bar ground against the grid at 4.3.

  The backbone. Mid grey rather than pale, because the panel has to read as a molecule and
  the pale greys measured here put the fold below the paper. A mark is kept off the backbone
  by the white halo under it, not by contrast, so the backbone is free to be chosen for the
  fold.

  The layout. naview, with one frame shared across the four nests and each panel at its own
  best angle. Residue i moves by 0.18 of the panel between nests under naview against
  0.34-0.53 for the other fold layouts, so the four panels are the most comparable a
  fold-derived drawing gets. Its cost is stated below.

What the drawing cannot do. A 6-mer site on a 795-nt molecule in a panel this size is a path
1.0-2.4 pt long drawn with a 2.8 pt pen, so two sites within about a pen width merge into one
mark and no amount of care in the renderer separates them. The count in each caption is the
census, not the mark count; `separable_marks` measures the difference over every site at once
and the sidecar records it per panel, so the caption can state it rather than the reader
discovering it. Two of the merges on this molecule are between *different* categories, which
print as one two-coloured object.

Geometry. Nothing here is a chosen coordinate except the margin, the gaps and the two
left-hand panel widths. Every other edge is measured off a first draw: each panel's furniture
overhang on all four sides, the figure legend's height, the structure titles' height, and the
folds' own aspect after rotation. Hand-set edges were wrong in both directions when this was
tried -- (a)'s longest pattern name is 36.5 pt, not the 20 estimated, and (b)'s y-furniture
needs 25 pt, not 21 -- and a measured pass cannot drift when the type size or the data change.

Reproducibility. Every input is read from the experiment tree, never from a cached
intermediate. The sidecar records the resolved path and SHA-256 of each front and motif file,
the build line each front came from *when that line names the file actually read*, the
resolved geometry in points, and the version of every library that touches the output. The
PDF carries no creation date, so two runs on one tree are byte-identical.

Usage:
    python fig2_row.py E2_DIR OUT_PDF OUT_JSON [--cell P26436] [--cai-floor 0.87]
"""

import argparse
import hashlib
import json
import math
import platform
import statistics
import sys
import time

import matplotlib

matplotlib.use("Agg")
import fig2  # noqa: E402
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt  # noqa: E402
import RNA  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import ConnectionPatch, Patch, Rectangle  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

PANEL = ("a.", "b.", "c.", "d.", "e.", "f.", "g.")
# A variant is a set of independent features, named for convenience. `funnel` splits panel
# (a) into its three tiers; `alt` adds the fifth structure panel; the rest are one-offs.
VARIANTS = ("plain", "funnel", "price", "both", "sweep", "alt", "freq", "funnel-alt")
FUNNEL_IN = ("funnel", "both", "funnel-alt")
ALT_IN = ("alt", "funnel-alt")
FREQ_HIGH = 0.97
ALT_SLACK = 10.0
SWEEP_GRID = [0.60 + 0.025 * i for i in range(17)]
P_W = 52.0
# Just the forbidden set. `structure_titles` supplies the verb, and drops it first when a
# narrow column cannot hold the sentence.
# Commas, not plus signs. "forbids I$+$II$+$III" measures 54.2 pt against a 53.4 pt column
# at five panels; the same sentence with commas is 52.5 and fits, so every panel keeps the
# verb instead of the row falling back to a bare set name that needs the caption to decode.
FORBIDDEN = {
    "F0": "nothing",
    "F1": "I",
    "F2": "I, II",
    "F3": "I, II, III",
}

# Three hues at three lightnesses. L*: 50.4 / 69.5 / 32.2.
PALETTE = {"I": "#2a78d6", "II": "#f29254", "III": "#8c1b5e"}
# F0 forbids nothing, so it takes a neutral -- and a dark one: the lighter grey it had sat
# 6.3 L* from category I, which is the front drawn next to it.
NEUTRAL = "#636058"
BACKBONE, HALO = "#bfbeba", "#ffffff"
# BAR_BG is the lattice-admissible bar, GRID the gridlines. They were one ink 0.7 L* apart,
# so in panel (a) a gridline vanished wherever it crossed a bar it did not mean.
# The furniture grey was #4e4c47, which is 0.17 in L* from category III -- the closest pair
# in the whole ink set, by a factor of twelve over the runner-up, and the one the docstring
# claimed had been retired. Moving the terminus marks and the floor line off it retired two
# *marks*; the spines, the ticks and the two hull legend swatches still carried it, against
# 622 category-III front points and four category-III bars. This grey is 15.5 from category
# III and 13.7 from body ink.
INK, FURNITURE, BAR_BG, GRID = "#0b0b0b", "#2b2926", "#d8d7d2", "#e4e3df"
NEST_INK = {"F0": NEUTRAL} | {n: PALETTE[c] for n, c in fig2.FORBIDS.items() if c}

# Four tiers. The row had three, with a panel's heading, its numbers, its ticks and two
# legends all on one of them, so nothing told a reader which was which.
# 5.0 is the floor most journals set. The measured slack carried it: the widest title was
# 52.5 pt of 53.4 and the widest caption 47.2, so a 4% bump still clears the column, and
# `fit` catches it if a future change does not.
LETTER, BODY, HEAD, LABEL = 7.2, 6.2, 5.6, 5.0
STROKE, HALO_W = 2.8, 3.8
PLOT_TYPE = 1  # naview

# ------------------------------------------------------------------------------- geometry
FIG_W, FIG_H = 7.2, 2.4
PT_W, PT_H = 72 * FIG_W, 72 * FIG_H
MARGIN, GAP_AB, GAP_BS, S_GAP, CLEAR = 5.0, 7.0, 7.0, 4.0, 3.5
LETTER_BASE = PT_H - 10.0
A_W, B_W = 50.0, 76.0
# The category bracket and its numeral, as distances left of (a)'s frame. The numeral is set
# on its side: upright it is as wide as three bars of the panel it labels.
BRACKET_PT, ROMAN_PT = 42.0, 47.5
# Only the rotation is chosen from this; the frame is then fitted to what the rotation gives.
# A structure column is roughly half again as tall as it is wide.
ROT_ASPECT = 1.5


def rect(x, y, w, h):
    return (x / PT_W, y / PT_H, w / PT_W, h / PT_H)


def style():
    sns.set_theme(style="white", context="paper")
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "text.latex.preamble": r"\usepackage{amsmath}",
            "font.size": BODY,
            "text.color": INK,
            "axes.titlesize": HEAD,
            "axes.labelsize": BODY,
            "axes.labelcolor": INK,
            "xtick.labelsize": LABEL,
            "ytick.labelsize": LABEL,
            "axes.linewidth": 0.5,
            "axes.edgecolor": FURNITURE,
            "xtick.color": FURNITURE,
            "ytick.color": FURNITURE,
            "xtick.labelcolor": INK,
            "ytick.labelcolor": INK,
            "xtick.major.size": 1.8,
            "ytick.major.size": 1.8,
            "xtick.major.width": 0.45,
            "ytick.major.width": 0.45,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


# ------------------------------------------------------------------------------ panel (a)


def energy_model(chosen, drawn):
    """The model the fronts were walked under, and how far RNAfold's default moves it.

    Measured over the designs this variant draws. It was a literal naming "these four
    designs" while the row drew five, and the fifth's offset lay outside the stated range.
    """
    md = RNA.md()
    md.dangles = 2
    offsets = []
    for nest in drawn:
        design = chosen[nest]
        seq = design["sequence"].replace("T", "U")
        offsets.append(RNA.fold_compound(seq, md).mfe()[1] - design["mfe"])
    # Against the LADDER, not the row. The row now spans 13.3 because the alternative gives
    # up energy on purpose; the quantity the offset endangers is the thing the ladder
    # measures, which is what forbidding the three categories costs.
    rungs = [chosen[n]["mfe"] for n in drawn if n in fig2.NESTS]
    ladder = max(rungs) - min(rungs)
    verdict = "more" if -min(offsets) > ladder else "less"
    return (
        f"dangles=0, Turner 2004. Refolding these {len(drawn)} designs at RNAfold's default "
        f"d=2 gives {-max(offsets):.1f} to {-min(offsets):.1f} kcal/mol lower, {verdict} than "
        f"the {ladder:.1f} kcal/mol the F0-F3 ladder spans"
    )


def lstar(colour):
    """CIE L* of a hex colour -- what survives a greyscale print."""
    r, g, b = (int(colour[i : i + 2], 16) / 255 for i in (1, 3, 5))
    lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in (r, g, b)]
    y = 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    return 116 * y ** (1 / 3) - 16 if y > 0.008856 else 903.3 * y


def contrast_gate(groups, floor=4.0):
    """The closest pair of inks a reader has to tell apart, measured at render time.

    This lived in the module docstring twice and went stale twice -- each time because a
    feature added an ink and the paragraph was carried across it, and each time the stale
    paragraph was the one explaining why staleness here is dangerous. A prose gate is a
    claim about a palette that no longer exists; this one is recomputed from whatever the
    figure actually draws and is recorded beside the result.

    It takes GROUPS, not one palette. A pair needs separating only where a reader has to
    tell the two apart, and that is a question about the drawing, not the palette: the
    swatches of one legend row are compared with each other, and the three tiers of one bar
    are compared with each other, but a tier of a category-I row and the hue of a
    category-III row are never on the same ruler. Demanding separation over all twelve inks
    at once has no solution on the lightness axis and so could only be satisfied by giving
    up the gate.
    """
    out, worst = {}, []
    for group, inks in groups.items():
        names = sorted(inks)
        values = {k: round(lstar(v), 2) for k, v in inks.items()}
        pairs = sorted(
            (round(abs(values[a] - values[b]), 2), a, b)
            for i, a in enumerate(names)
            for b in names[i + 1 :]
        )
        out[group] = {
            "lstar": values,
            "closest": [{"delta": d, "between": [a, b]} for d, a, b in pairs[:2]],
        }
        if pairs and pairs[0][0] < floor:
            worst.append((group, *pairs[0]))
    for group, d, a, b in worst:
        print(
            f"warning: within {group}, {a} and {b} are {d:.2f} apart in L* and a reader has "
            f"to tell them apart; they print as one grey",
            file=sys.stderr,
        )
    return {"floor": floor, "groups": out}


def tint(colour, frac):
    """A category's hue lightened towards the bar ground, for the funnel's middle tier."""
    a = [int(colour[i : i + 2], 16) for i in (1, 3, 5)]
    b = [int(BAR_BG[i : i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(x + (y - x) * frac):02x}" for x, y in zip(a, b, strict=True))


def panel_frequency(ax, cat, freq):
    """How many times each pattern appears, and which way that moves with CAI.

    The presence panel answers "does this motif turn up", which is the same question for a
    pattern appearing once and one appearing three times. Here the bar is the count the
    chosen design actually spells, and the tick is the same count at a higher CAI bound: a
    tick right of its bar means the pattern gets *more* common as the design is better
    adapted, and the rows run both ways.

    The lattice's own count -- how many positions could carry the pattern -- is measured and
    recorded but not drawn. It reaches 12.7 against occurrences of at most 3.25, so sharing
    an axis with it put most of the coloured bars under two points at print size, and that
    context is already the grey bar of the presence panel.
    """
    entries, counts = cat["entries"], freq["counts"]
    y = list(range(len(entries)))[::-1]
    for row, e in zip(y, entries, strict=True):
        c = counts[e["motif"]]
        ax.barh(row, c["at_floor"], height=0.80, color=PALETTE[e["category"]], lw=0, zorder=2)
        ax.plot([c["at_high"]], [row], marker="|", ms=3.4, mew=0.7, color=INK, lw=0, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels([e["motif"].replace("T", "U") for e in entries], fontsize=LABEL)
    for tick, e in zip(ax.get_yticklabels(), entries, strict=True):
        tick.set_color(PALETTE[e["category"]])
    ax.set_ylim(-0.75, len(entries) - 0.25)
    ax.set_xlim(0, max(max(c["at_floor"], c["at_high"]) for c in counts.values()) * 1.06)
    ax.tick_params(axis="x", labelsize=LABEL, pad=1.2)
    ax.tick_params(axis="y", length=0, pad=1.2)
    ax.set_xlabel("occurrences per design", fontsize=BODY, labelpad=2.4)
    ax.grid(axis="x", color=GRID, lw=0.35)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    bracket_categories(ax, entries, y)


def bracket_categories(ax, entries, y):
    """The three category brackets down the left edge, shared by both panel (a) treatments."""
    start = 0
    for _, roman, _ in fig2.CATEGORIES:
        k = sum(1 for e in entries if e["category"] == roman)
        if not k:
            continue
        top, bot = y[start], y[start + k - 1]
        ax.annotate(
            "",
            xy=(-BRACKET_PT / A_W, bot - 0.38),
            xytext=(-BRACKET_PT / A_W, top + 0.38),
            xycoords=ax.get_yaxis_transform(),
            arrowprops={"arrowstyle": "-", "lw": 1.0, "color": PALETTE[roman]},
            annotation_clip=False,
        )
        ax.text(
            -ROMAN_PT / A_W,
            (top + bot) / 2,
            rf"\textbf{{{roman}}}",
            transform=ax.get_yaxis_transform(),
            ha="center",
            va="center",
            rotation=90,
            fontsize=BODY,
            color=PALETTE[roman],
        )
        start += k


def panel_catalog(ax, cat, census=None):
    """Every pattern, grouped by category, as a share of the screened population.

    With a census the bar becomes a funnel: the lattice admits the pattern, some design on
    the front spells it, and the *optimal* design at the bound spells it. The three tiers
    nest by construction, and the three categories narrow at very different rates -- which
    is the one comparison between the sets that needs no constrained front, and so inherits
    no front's time budget.
    """
    n = cat["population"]
    entries = cat["entries"]
    y = list(range(len(entries)))[::-1]
    for row, e in zip(y, entries, strict=True):
        hue = PALETTE[e["category"]]
        ax.barh(row, 100 * e["spellable"] / n, height=0.80, color=BAR_BG, lw=0, zorder=1)
        if census is None:
            ax.barh(row, 100 * e["occurs"] / n, height=0.80, color=hue, lw=0, zorder=2)
            continue
        ax.barh(row, 100 * e["occurs"] / n, height=0.80, color=tint(hue, 0.55), lw=0, zorder=2)
        ax.barh(
            row,
            100 * census["counts"][e["motif"]] / census["proteins"],
            height=0.80,
            color=hue,
            lw=0,
            zorder=3,
        )
    ax.set_yticks(y)
    ax.set_yticklabels([e["motif"].replace("T", "U") for e in entries], fontsize=LABEL)
    for tick, e in zip(ax.get_yticklabels(), entries, strict=True):
        tick.set_color(PALETTE[e["category"]])
    ax.set_ylim(-0.75, len(entries) - 0.25)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.tick_params(axis="x", labelsize=LABEL, pad=1.2)
    ax.tick_params(axis="y", length=0, pad=1.2)
    ax.set_xlabel(rf"\% of {n} proteins", fontsize=BODY, labelpad=2.4)
    ax.grid(axis="x", color=GRID, lw=0.35)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    bracket_categories(ax, entries, y)


def panel_sweep(ax, sweep):
    """Site density in the shipped design as the CAI bound walks up the front.

    Log scale, because the three sets differ by a factor of thirty in density and a linear
    axis prints the smallest of them as a flat line on the frame. On a log axis a given
    proportional change is the same vertical distance whatever the set's level, which is
    what the comparison is about. A curve ends where its density reaches zero.
    """
    seen = {}
    for _, roman, _ in fig2.CATEGORIES:
        xs = [row["cai"] for row in sweep if row["per_kb"][roman] > 0]
        ys = [row["per_kb"][roman] for row in sweep if row["per_kb"][roman] > 0]
        if not ys:
            continue
        ax.plot(xs, ys, color=PALETTE[roman], lw=0.9, zorder=4)
        ax.plot(xs[-1:], ys[-1:], marker="o", ms=1.7, color=PALETTE[roman], lw=0, zorder=5)
        seen[roman] = {
            "first": round(ys[0], 4),
            "last": round(ys[-1], 4),
            "ends_at_cai": xs[-1],
            "ratio": round(ys[-1] / ys[0], 4),
        }
    ax.set_yscale("log")
    # One decade of default log ticks carried a single label here; the curves span two.
    ax.set_yticks([0.03, 0.1, 0.3, 1.0, 3.0])
    ax.set_yticklabels(["0.03", "0.1", "0.3", "1", "3"])
    ax.minorticks_off()
    ax.set_xlim(min(r["cai"] for r in sweep), max(r["cai"] for r in sweep))
    ax.tick_params(labelsize=LABEL, pad=1.2)
    ax.set_xticks([0.6, 0.8, 1.0])
    ax.set_xlabel("CAI bound", fontsize=BODY, labelpad=2.4)
    ax.set_ylabel("sites per kb", fontsize=BODY, labelpad=1.6)
    ax.grid(color=GRID, lw=0.35, which="both")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return {"grid": [r["cai"] for r in sweep], "proteins": sweep[0]["proteins"], "sets": seen}


def panel_price(ax, prices):
    """What a rung cost against how many sites it had to remove, one point per cell.

    Only rows where the bound is resolved on *both* fronts the difference is taken between,
    so the populations differ by set and are stated. A point on the vertical axis is a rung
    that removed nothing: its cost is a cost against a constraint that forbade nothing, and
    it is drawn hollow so it cannot be read as a cheap constraint.
    """
    kept = [p for p in prices if p["exact"]]
    drawn = {}
    for roman in ("I", "II", "III"):
        rows = [p for p in kept if p["set"] == roman]
        bound = [p for p in rows if p["removed"] > 0]
        drawn[roman] = {"cells": len(rows), "bound": len(bound)}
        ax.scatter(
            [p["removed"] for p in rows if p["removed"] == 0],
            [p["mfe_cost"] for p in rows if p["removed"] == 0],
            s=3.0,
            facecolors="none",
            edgecolors=PALETTE[roman],
            lw=0.45,
            zorder=3,
        )
        ax.scatter(
            [p["removed"] for p in bound],
            [p["mfe_cost"] for p in bound],
            s=3.6,
            c=PALETTE[roman],
            lw=0,
            zorder=4,
        )
        if bound:
            slope = statistics.median(p["mfe_cost"] / p["removed"] for p in bound)
            x = max(p["removed"] for p in bound)
            ax.plot([0, x], [0, slope * x], color=PALETTE[roman], lw=0.5, zorder=2)
            drawn[roman]["per_site"] = round(slope, 3)
    ax.set_xlim(left=-0.4)
    ax.set_ylim(bottom=min(-0.6, ax.get_ylim()[0]))
    ax.tick_params(labelsize=LABEL, pad=1.2)
    ax.set_xlabel("sites removed", fontsize=BODY, labelpad=2.4)
    ax.set_ylabel(r"$\Delta$MFE (kcal/mol)", fontsize=BODY, labelpad=1.6)
    ax.grid(color=GRID, lw=0.35)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return drawn


# ------------------------------------------------------------------------------ panel (b)


def scatter(ax, pts, on, colour, s_on, s_off, lw, line=True):
    off = [i for i in range(len(pts)) if i not in on]
    ax.scatter(
        [pts[i]["cai"] for i in off],
        [pts[i]["mfe"] for i in off],
        s=s_off,
        facecolors="none",
        edgecolors=colour,
        lw=lw,
        zorder=3,
    )
    ins = sorted(on)
    if line:
        ax.plot(
            [pts[i]["cai"] for i in ins],
            [pts[i]["mfe"] for i in ins],
            color=colour,
            lw=0.6,
            zorder=4,
        )
    ax.scatter(
        [pts[i]["cai"] for i in ins], [pts[i]["mfe"] for i in ins], s=s_on, c=colour, lw=0, zorder=5
    )


def zoom_box(parent, child, x0, x1, y0, y1):
    """The zoomed region as a closed rectangle, with dashed ties to the inset's lower corners.

    The same treatment the overview figure uses, so a reader who has seen one reads the other.
    """
    parent.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=INK, lw=0.45, zorder=6))
    for px, cx in ((x0, 0.0), (x1, 1.0)):
        parent.add_artist(
            ConnectionPatch(
                xyA=(px, y1),
                coordsA=parent.transData,
                xyB=(cx, 0.0),
                coordsB=child.transAxes,
                color=INK,
                lw=0.45,
                ls=(0, (2.4, 1.8)),
                zorder=6,
            )
        )


def label_fronts(inset, placed, ywin, fig, renderer):
    """Name the four fronts where they separate, pushed apart so no two labels touch.

    This is the figure's answer to its own colour polarity: in (b) a hue names the nest that
    *forbids* that category, the opposite of what the same hue means everywhere else, so (b)
    has to say which front is which somewhere a reader can read it. So the separation has to
    hold: it is measured against the inset's drawn height rather than a remembered one, which
    is why this runs after the packer. A constant here enforced 3.5 pt between 4.25 pt boxes,
    and the labels it was there to separate overlapped by 0.7 pt.
    """
    height = inset.get_window_extent(renderer).height * 72 / fig.dpi
    step = 1.35 * LABEL * (ywin[1] - ywin[0]) / height
    order = sorted(placed, key=lambda n: placed[n])
    ys = [placed[n] for n in order]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + step)
    # Pushing only upwards walks the whole stack off its fronts: two of four labels then sat
    # nearer a neighbour's front than their own, which for a figure whose hue means the
    # opposite thing in this panel is the one error that cannot be allowed. Re-centring costs
    # nothing and halves every displacement.
    drift = sum(ys) / len(ys) - sum(placed.values()) / len(placed)
    ys = [y - drift for y in ys]
    lo = min(0.0, ys[0] - (ywin[0] + 0.5 * step))
    hi = max(0.0, ys[-1] - (ywin[1] - 0.5 * step))
    shift = hi + lo
    for nest, y in zip(order, ys, strict=True):
        inset.text(
            0.03,
            y - shift,
            rf"\textbf{{{nest}}}",
            transform=inset.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=LABEL,
            color=NEST_INK[nest],
            zorder=8,
            path_effects=[pe.withStroke(linewidth=1.6, foreground=HALO)],
        )


def occluded(ax, inset_rect, fronts):
    """Front points the inset hides. It is opaque, so its placement costs data unless checked.

    The empty quadrant of (b) is high MFE at low CAI, because a front rises only at the
    high-CAI end. Widening the inset eats into that quadrant, and the only way to know
    whether it has reached the fronts is to count what falls underneath.
    """
    x0, y0, w, h = inset_rect
    (xlo, xhi), (ylo, yhi) = ax.get_xlim(), ax.get_ylim()
    hidden = 0
    for _, pts in fronts.values():
        for p in pts:
            fx = (p["cai"] - xlo) / (xhi - xlo)
            fy = (p["mfe"] - ylo) / (yhi - ylo)
            if x0 <= fx <= x0 + w and y0 <= fy <= y0 + h:
                hidden += 1
    return hidden


def panel_front(ax, fronts, hulls, floor, chosen):
    """The four fronts, and the window where they separate.

    At full scale they overprint: the whole ladder costs under two percent of the free energy
    at the bound. The nests are drawn deepest-first so that F1, which coincides with F2 at
    several CAI bounds, is not buried by it.
    """
    lo = min(p["cai"] for _, pts in fronts.values() for p in pts)
    hi = max(p["cai"] for _, pts in fronts.values() for p in pts)
    for nest in reversed(fig2.NESTS):
        scatter(ax, fronts[nest][1], hulls[nest], NEST_INK[nest], 1.8, 1.3, 0.30)
    pad = 0.02 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_xticks([round(lo, 2), floor, 1.0])
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.tick_params(labelsize=LABEL, pad=1.2)
    ax.set_xlabel("CAI", fontsize=BODY, labelpad=2.4)
    ax.set_ylabel("MFE (kcal/mol)", fontsize=BODY, labelpad=1.6)
    ax.grid(color=GRID, lw=0.35)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # The selection rule, on the main axes where the designs are chosen -- body ink, not the
    # furniture grey, because it is a statement about the data and not an axis.
    ax.axvline(floor, color=INK, lw=0.45, ls=(0, (2.0, 1.4)), zorder=6)

    # Wide enough that the hull encoding is not vacuous where it is legible. At one percent
    # of the CAI span the window held 36 points and exactly one hull vertex, so "filled means
    # on the hull" was a key to a mark the reader could find once; at three percent it holds
    # 120 points and six. It widens further when an extra design has to be inside it: a mark
    # drawn in a panel the window excludes is worse than a narrower window.
    half = 0.030 * (hi - lo)
    alt = chosen.get("ALT")
    if alt is not None:
        half = max(half, 1.15 * abs(alt["cai"] - floor))
    win = (floor - half, floor + half)
    # Landscape. The fronts separate vertically over a few kcal/mol and run horizontally over
    # the whole window, so height spent past the fronts' own span is spent on nothing.
    inset_rect = (0.035, 0.560, 0.780, 0.405)
    # Above the main axes' own zorder-6 furniture. An inset defaults to 5, so the floor line
    # was drawn across it: a second dashed vertical inside the window, 6.5 pt from the real
    # one and reading at CAI 0.8719 against the inset's scale, with nothing to tell them
    # apart but which one carries the drawn-design circles.
    inset = ax.inset_axes(inset_rect, zorder=7)
    ys, placed = [], {}
    for nest in reversed(fig2.NESTS):
        pts = [p for p in fronts[nest][1] if win[0] <= p["cai"] <= win[1]]
        if not pts:
            continue
        ys += [p["mfe"] for p in pts]
        # Anchored at the window's LEFT edge. The fronts rise left to right, so the right
        # edge is where they are highest and most crowded -- and where a wider window puts
        # the most points. The left edge is the sparse end and the labels sit clear of it.
        placed[nest] = min(pts, key=lambda p: p["cai"])["mfe"]
        live = {id(fronts[nest][1][i]) for i in hulls[nest]}
        scatter(
            inset,
            pts,
            {i for i, p in enumerate(pts) if id(p) in live},
            NEST_INK[nest],
            3.2,
            2.8,
            0.30,
            line=False,
        )
    missing = []
    for nest in fig2.NESTS:
        d = chosen[nest]
        if win[0] <= d["cai"] <= win[1]:
            inset.plot(
                [d["cai"]], [d["mfe"]], marker="o", ms=3.4, mfc="none", mec=INK, mew=0.45, zorder=7
            )
        else:
            missing.append(nest)
    if alt is not None:
        # The alternative is a design the ladder's rule does not select, so it carries its own
        # mark rather than the circle that means "this panel's optimum".
        for ax_, ms in ((inset, 3.8), (ax, 2.6)):
            ax_.plot(
                [alt["cai"]],
                [alt["mfe"]],
                marker="D",
                ms=ms,
                mfc="none",
                mec=INK,
                mew=0.5,
                zorder=7,
            )
    inset.axvline(floor, color=INK, lw=0.4, ls=(0, (2.0, 1.4)), zorder=6)
    ypad = 0.12 * (max(ys) - min(ys)) or 0.5
    ywin = (min(ys) - ypad, max(ys) + ypad)
    inset.set_xlim(*win)
    inset.set_ylim(*ywin)
    inset.set_xticks([])
    inset.set_yticks([])
    for sp in inset.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.4)
        sp.set_edgecolor(INK)
    zoom_box(ax, inset, win[0], win[1], *ywin)
    # The selection rule names the dashed line, so it is set beside it, under the zoom box.
    # It used to stand on its side above the fronts; the wider inset reaches that pocket and
    # its left tie now crosses it. Below the box and right of the floor nothing is drawn: a
    # front's lowest MFE at any CAI bound is reached at the bound, so no point lies under it.
    ax.annotate(
        rf"CAI $\geq {floor:g}$",
        xy=(floor, ywin[0]),
        xytext=(1.8, -1.8),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=LABEL,
        color=INK,
    )
    return {
        "range": [lo, hi],
        "window": list(win),
        "mfe_window": list(ywin),
        "designs_outside": missing,
        "occluded": occluded(ax, inset_rect, fronts),
        "inset": inset,
        "placed": placed,
    }


# ------------------------------------------------------------------------- panels (c)-(f)


def coords(structure):
    RNA.cvar.rna_plot_type = PLOT_TYPE
    c = RNA.get_xy_coordinates(structure)
    return (
        [c.get(i).X for i in range(len(structure))],
        [c.get(i).Y for i in range(len(structure))],
    )


def separable_marks(runs, xs, ys, pt_per_unit):
    """How many marks a reader can actually tell apart, per category.

    Two site runs whose drawn paths come within a pen width of each other print as one mark,
    whatever the census says. On this molecule that is not a rounding error: a 6-mer is a path
    about a pen width long, so neighbouring sites merge unconditionally. Reporting the two
    numbers is the only honest option short of redrawing sites as leadered markers.

    The merge is run over every site at once rather than per category. Confining it to one
    category at a time measured only the kind of collision that does not change the reader's
    colour: two sites of *different* categories a pen width apart print as one two-coloured
    object, which is worse, and this molecule has two such pairs at 1.3 and 2.4 pt.
    """
    segs, owner = [], []
    for _, roman, _ in fig2.CATEGORIES:
        for lo, hi in runs[roman]:
            path = [(xs[i], ys[i]) for i in range(lo, min(hi, len(xs)))]
            if path:
                segs.append(path)
                owner.append(roman)
    parent = list(range(len(segs)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    blended = 0
    for i in range(len(segs)):
        for j in range(i + 1, len(segs)):
            if any(math.dist(p, q) * pt_per_unit < STROKE for p in segs[i] for q in segs[j]):
                parent[find(i)] = find(j)
                if owner[i] != owner[j]:
                    blended += 1
    out = {
        r: len({find(i) for i in range(len(segs)) if owner[i] == r}) for _, r, _ in fig2.CATEGORIES
    }
    return out | {"blended_pairs": blended}


def draw_structure(ax, design, runs, span, angle, pt_per_unit):
    """One nest's design. `span` is shared by every panel, so the four are one scale."""
    span_x, span_y = span
    xs, ys = fig2.rotate(*coords(design["structure"]), angle)
    segments, jumps = fig2.backbone_segments(xs, ys)
    for sx, sy in segments:
        ax.plot(sx, sy, color=BACKBONE, lw=0.9, solid_capstyle="round", zorder=2)
    for _, roman, _ in reversed(fig2.CATEGORIES):
        for lo, hi in runs[roman]:
            seg = list(range(lo, min(hi, len(xs))))
            if not seg:
                continue
            # The halo, not the backbone's lightness, is what keeps a mark legible.
            for w, colour, z in ((HALO_W, HALO, 4), (STROKE, PALETTE[roman], 5)):
                ax.plot(
                    [xs[i] for i in seg],
                    [ys[i] for i in seg],
                    color=colour,
                    lw=w,
                    solid_capstyle="round",
                    zorder=z,
                )
    term = draw_termini(ax, xs, ys, pt_per_unit)
    cx, cy = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2
    ax.set_xlim(cx - span_x / 2, cx + span_x / 2)
    ax.set_ylim(cy - span_y / 2, cy + span_y / 2)
    ax.set_aspect("equal")
    ax.axis("off")
    return jumps, term, separable_marks(runs, xs, ys, pt_per_unit)


def draw_termini(ax, xs, ys, pt_per_unit):
    """The chain's two ends, merged into one mark when they would print on top of each other.

    naview lays the two ends of a closed fold next to each other, and on this molecule they
    land 0.7-0.9 pt apart against a 1.9 pt marker: two marks and two legend entries would name
    a 5' square that the 3' square completely covers. Squares, not circles, because naview
    draws every hairpin loop as a small circle.
    """
    sep = math.dist((xs[0], ys[0]), (xs[-1], ys[-1])) * pt_per_unit
    if sep >= 3.0 * TERM_MS:
        ax.plot([xs[0]], [ys[0]], marker="s", ms=TERM_MS, mfc=HALO, mec=INK, mew=0.45, zorder=6)
        ax.plot([xs[-1]], [ys[-1]], marker="s", ms=TERM_MS, mfc=INK, mec=INK, mew=0.45, zorder=6)
        return {"mode": "split", "separation_pt": sep}
    ax.plot(
        [(xs[0] + xs[-1]) / 2],
        [(ys[0] + ys[-1]) / 2],
        marker="s",
        ms=TERM_MS,
        mfc=HALO,
        mec=INK,
        mew=0.45,
        zorder=6,
    )
    return {"mode": "merged", "separation_pt": sep}


TERM_MS = 1.9


def preserving_alternative(fronts, chosen, floor, deep, shallow):
    """The deepest nest's design whose fold is closest to the rung before it.

    The figure draws each nest's minimum-MFE design, and at the last rung that design folds
    into a different molecule. It is not the only design the constraint allows: this is the
    one, among those within a few kcal/mol of the optimum, whose fold is nearest the design
    before the rung. It obeys the same constraint -- it is simply not the cheapest way to.

    Found by filtering the recorded front, which holds only (CAI, MFE)-optimal designs, so
    the cost it reports is an upper bound: a design giving up a little on both axes would be
    dominated and therefore absent. Asking for the fold as a constraint would be exact.
    """
    before = chosen[shallow]["structure"]
    band = [
        p
        for p in fronts[deep][1]
        if p["cai"] >= floor and p["mfe"] <= chosen[deep]["mfe"] + ALT_SLACK
    ]
    keep = fig2.pair_set(before)
    return min(band, key=lambda p: len(keep ^ fig2.pair_set(p["structure"])))


def shared_span(chosen, lock=False):
    """One frame for all four nests, shaped like the folds rather than like the column.

    The rotation is chosen against an approximate column shape; the frame is then fitted to
    what the rotation actually produced. Pinning the frame to the column instead left 19 pt of
    every structure panel unfillable, because no rotation of these folds is that elongated.
    """
    # `lock` gives every panel the first one's rotation. Per-panel angles are right when the
    # panels are read one at a time, and wrong when the figure's claim is that two of the
    # folds are alike: a rotation chosen per panel can manufacture that resemblance or hide
    # it, and a reader compares shapes before reading the number underneath.
    # Chosen, not inherited. Every panel is drawn at the same scale, set by the widest fold
    # after rotation, so the angle that matters is the one minimising that width over the
    # whole set -- not the angle that happens to suit the first panel. Taking F0's best
    # angle gave a shared span of 2651 x 2962 against 2108 x 3109 here: 26% less scale and
    # 58% less drawn area, with the difference left as empty band under the frames.
    first = min(
        (math.radians(d) for d in range(0, 180)),
        key=lambda a: max(
            (lambda rx, _: max(rx) - min(rx))(*fig2.rotate(*coords(d2["structure"]), a))
            for d2 in chosen.values()
        ),
    )
    wide, tall, angles = [], [], {}
    for nest in chosen:
        raw = coords(chosen[nest]["structure"])
        angles[nest] = first if lock else fig2.best_angle(*raw, aspect=ROT_ASPECT)
        rx, ry = fig2.rotate(*raw, angles[nest])
        wide.append(max(rx) - min(rx))
        tall.append(max(ry) - min(ry))
    return (1.06 * max(wide), 1.06 * max(tall)), angles


def fit(art, forms, fig, renderer, width):
    """Set the widest form that fits, and report which one that was.

    Every form says the same thing; the later ones give up a word or take another line for
    it. Five structure columns leave 53 pt a panel, against 54 for "forbids I+II+III" and 58
    for an MFE line carrying its difference, so at this width the choice is not cosmetic --
    hand-set text overflowed its frame on two of the five panels before this existed.
    """
    for i, form in enumerate(forms):
        art.set_text(form)
        if renderer is None or art.get_window_extent(renderer).width * 72 / fig.dpi <= width:
            return i
    return len(forms) - 1


def structure_caption(
    ax, design, counts, y, base, show_cai, own_base=None, fig=None, renderer=None, width=1e9
):
    """What the drawn design cost, and the site census that says the constraint took.

    Every category is named, the measured zeros included: that a forbidden category's count
    is zero is the whole claim of the ladder, and a panel listing only the non-zero rows
    shows the claim's conclusion while hiding its subject.

    The forms are returned rather than settled here, because a panel that keeps "sites:"
    beside four neighbours that dropped the colon reads as a different kind of line. The row
    is levelled to the narrowest form any panel needed.

    CAI is printed only where it is not the bound. On this subject the four selected designs
    sit 0.00013-0.00068 above the floor, so all four print the same three digits and the line
    says nothing a reader cannot read off the dashed rule in (b). That is a MEASUREMENT, not
    a guarantee: "lowest MFE at or above the floor" picks the lowest-CAI feasible point, and
    its distance to the floor is the local front spacing, which the rule does not bound --
    over the 85-protein screen the same rule lands 0.0024 above it on one subject, and
    sweeping this subject's floor it exceeds 0.001 on about a tenth of the floors tried. The
    per-design CAI is in the sidecar either way. The fifth design is not selected that way
    and its CAI is part of why.
    """
    sites = " ".join(rf"{r}\,{counts[r]}" for _, r, _ in fig2.CATEGORIES)
    half = len(fig2.CATEGORIES) // 2 + 1
    split = (
        " ".join(rf"{r}\,{counts[r]}" for _, r, _ in fig2.CATEGORIES[:half])
        + "\n"
        + " ".join(rf"{r}\,{counts[r]}" for _, r, _ in fig2.CATEGORIES[half:])
    )
    cost = "" if design["mfe"] == base else rf" $({design['mfe'] - base:+.1f})$"
    mfe = rf"MFE ${design['mfe']:.1f}$"
    # An alternative is priced against the design it is an alternative TO. (g) against F0 is
    # +13.3, which is the ladder plus the alternative; against its own nest's optimum it is
    # +4.8, and that is what "obeys the same constraint, just not the cheapest way" costs.
    own = "" if own_base is None else rf"$+{design['mfe'] - own_base:.1f}$ on F3" + "\n"
    cai = (rf"CAI ${design['cai']:.3f}$" + "\n") if show_cai else ""
    cai = own + cai
    art = ax.text(
        0.5,
        y,
        "",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=LABEL,
        linespacing=1.45,
        color=INK,
    )
    forms = [
        f"{mfe}{cost}\n{cai}sites: {sites}",
        f"{mfe}{cost}\n{cai}sites {sites}",
        f"{mfe}\n{cost.strip()}\n{cai}sites {sites}" if cost else f"{mfe}\n{cai}sites {sites}",
        f"{mfe}\n{cost.strip()}\n{cai}{split}" if cost else f"{mfe}\n{cai}{split}",
    ]
    return art, fit(art, forms, fig, renderer, width), forms


def title_forms(forbids):
    """What a panel says it forbids, widest first.

    The sentence, then the set alone, then the set without its separators. Which form a
    panel can hold depends on how long its own set is, so they are resolved per panel and
    then levelled: a row where one panel says "forbids I" beside a neighbour saying
    "I+II+III" reads as two kinds of label rather than one scale.
    """
    # The narrow forms WRAP; they never rewrite the set. An earlier fallback stripped the
    # separators, which turned "I, II" into "III" -- the one category F2 does not forbid, in
    # a figure whose whole claim is which categories are forbidden. A title that overflows
    # is visible; a title that lies is not.
    return [f"forbids {forbids}", f"forbids\n{forbids}", forbids]


def structure_titles(ax, fig, renderer, width, frame_h, label, ink, forbids, moved):
    """The three lines over a structure panel: which nest, what it forbids, how far it moved.

    The distance sits here rather than in the caption because it is what tells (f) and (g)
    apart, and those two carry the same constraint: a reader comparing them looks at the top
    of the panel, not the bottom. Giving it to every panel that has a predecessor also puts
    the structural ladder on one line of the figure, which was previously only in the text.
    """
    # The tier is reserved whether or not this panel has a predecessor. Giving (c) one line
    # fewer dropped its whole title block 7.2 pt below the other four, so at the height where
    # they printed a distance, (c) printed a sentence, and at the height where they printed
    # the sentence, (c) printed its nest label: three tiers, three different contents.
    ax.set_title(moved or " ", fontsize=LABEL, pad=2.4, color=FURNITURE)
    base = 1 + (2.4 + LABEL * 1.5) / frame_h
    art = ax.text(
        0.5,
        base,
        "",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=HEAD,
        color=INK,
    )
    used = fit(art, title_forms(forbids), fig, renderer, width)
    ax.text(
        0.5,
        base + (HEAD * 1.45) / frame_h,
        rf"\textbf{{{label}}}",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=HEAD,
        color=ink,
    )
    return art, used


# ------------------------------------------------------------------------ measured layout


def overhang(ax, fig, renderer):
    """How far a panel's furniture reaches past its own frame, on all four sides, in points."""
    k = 72 / fig.dpi
    t, f = ax.get_tightbbox(renderer), ax.get_window_extent(renderer)
    return {
        "l": (f.x0 - t.x0) * k,
        "r": (t.x1 - f.x1) * k,
        "b": (f.y0 - t.y0) * k,
        "t": (t.y1 - f.y1) * k,
    }


def pack(fig, ax_a, ax_b, legend, renderer, ax_p=None, columns=4):
    """Place the data panels against their measured furniture; return the structures' band.

    A fifth panel, when one is drawn, is reserved out of the right-hand end before the four
    structure columns divide what is left -- so adding it narrows the folds rather than
    pushing anything off the canvas.
    """
    k = 72 / fig.dpi
    legend_top = legend.get_window_extent(renderer).y1 * k
    boxes = [overhang(ax, fig, renderer) for ax in (ax_a, ax_b) + ((ax_p,) if ax_p else ())]
    oa, ob = boxes[0], boxes[1]
    bottom = legend_top + CLEAR + max(o["b"] for o in boxes)
    top = LETTER_BASE - CLEAR - max(o["t"] for o in boxes)
    a_left = MARGIN + oa["l"]
    b_left = a_left + A_W + oa["r"] + GAP_AB + ob["l"]
    s_left = b_left + B_W + ob["r"] + GAP_BS
    reserve, p_left = 0.0, None
    if ax_p is not None:
        op = boxes[2]
        p_left = PT_W - MARGIN - op["r"] - P_W
        reserve = op["l"] + P_W + op["r"] + GAP_BS
    s_w = (PT_W - MARGIN - s_left - reserve - (columns - 1) * S_GAP) / columns
    ax_a.set_position(rect(a_left, bottom, A_W, top - bottom))
    ax_b.set_position(rect(b_left, bottom, B_W, top - bottom))
    if ax_p is not None:
        ax_p.set_position(rect(p_left, bottom, P_W, top - bottom))
    return {
        "panel_g": [p_left, P_W] if ax_p is not None else None,
        "legend_top": legend_top,
        "panel_ab_top": top,
        "panel_ab_bottom": bottom,
        "panel_a": [a_left, A_W],
        "panel_a_overhang": oa,
        "panel_a_row_pitch": None,
        "panel_b": [b_left, B_W],
        "panel_b_overhang": ob,
        "left": s_left,
        "right": PT_W - MARGIN,
        "gap": S_GAP,
        "columns": columns,
        "width": s_w,
    }


# --------------------------------------------------------------------------- provenance


def motif_digests(root):
    out = {}
    for fname, roman, _ in fig2.CATEGORIES:
        path = f"{root}/motifs/{fname}.txt"
        with open(path, "rb") as fh:
            out[roman] = {"path": path, "sha256": hashlib.sha256(fh.read()).hexdigest()}
    return out


def environment(root, motif_paths):
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "matplotlib": matplotlib.__version__,
        "seaborn": sns.__version__,
        "viennarna": RNA.__version__,
        "usetex": bool(plt.rcParams["text.usetex"]),
        "e2_dir": root,
        "motif_files": motif_paths,
    }


# ------------------------------------------------------------------------------ the figure


def legend_handles(funnel=False, drawn=fig2.NESTS, lattice=True):
    # Ordered by the panels the entries belong to: (a) and (c)-(f) first, then (b)'s three.
    handles = [Patch(facecolor=PALETTE[r], label=f"cat {r}") for _, r, _ in fig2.CATEGORIES]
    if funnel:
        # One entry showing the nesting, rather than a lone lightened swatch beside the
        # saturated ones. A tint of category I measures 0.30 in L* from category II itself,
        # so as separate legend entries two swatches two apart printed as the same grey --
        # and the inner tier, the one the panel was rebuilt to show, had no entry at all.
        # Drawn as three patches in one entry, they are read in order instead of compared.
        handles.append(
            (
                Patch(facecolor=PALETTE["I"], label="in the optimum / on the front / lattice"),
                Patch(facecolor=tint(PALETTE["I"], 0.55)),
                Patch(facecolor=BAR_BG),
            )
        )
    elif lattice:
        handles.append(Patch(facecolor=BAR_BG, label="lattice-admissible"))
    else:
        # The frequency panel draws no grey, so naming it would key a mark nothing carries.
        handles.append(
            Line2D(
                [],
                [],
                color=INK,
                marker="|",
                ms=3.4,
                mew=0.7,
                lw=0,
                label=f"at CAI $\\geq$ {FREQ_HIGH:g}",
            )
        )
    handles.append(
        Line2D(
            [],
            [],
            color=INK,
            marker="s",
            ms=TERM_MS,
            lw=0,
            mfc=HALO,
            mew=0.45,
            label="5$'$/3$'$ ends",
        )
    )
    if "ALT" in drawn:
        handles.append(
            Line2D(
                [],
                [],
                color=INK,
                marker="D",
                ms=2.4,
                lw=0,
                markerfacecolor="none",
                markeredgewidth=0.5,
                label="drawn in (g)",
            )
        )
    handles += [
        Line2D([], [], color=INK, marker="o", ms=1.7, lw=0.6, label="on hull"),
        Line2D(
            [],
            [],
            color=INK,
            marker="o",
            ms=1.7,
            lw=0,
            markerfacecolor="none",
            markeredgewidth=0.30,
            label="off hull",
        ),
        Line2D(
            [],
            [],
            color=INK,
            marker="o",
            ms=2.6,
            lw=0,
            markerfacecolor="none",
            markeredgewidth=0.45,
            label=f"drawn in (c)--({chr(ord('c') + len(drawn) - 1 - ('ALT' in drawn))})",
        ),
    ]
    return handles


def build(root, cell, floor, out_pdf, out_json, variant="plain"):
    style()
    sets = fig2.read_patterns(root)
    cat = fig2.catalog(root)
    fronts, hulls, chosen, nests = fig2.collect(root, cell, floor, sets)
    census = fig2.optimum_census(root, cat, floor, sets) if variant in FUNNEL_IN else None
    prices = fig2.set_prices(root, floor, sets) if variant in ("price", "both") else None
    sweep = fig2.front_sweep(root, cat, sets, SWEEP_GRID) if variant == "sweep" else None
    freq = (
        fig2.occurrence_frequency(root, cat, sets, floor, FREQ_HIGH) if variant == "freq" else None
    )
    drawn, alt_of = list(fig2.NESTS), None
    if variant in ALT_IN:
        deep, shallow = fig2.NESTS[-1], fig2.NESTS[-2]
        design = preserving_alternative(fronts, chosen, floor, deep, shallow)
        cover, runs, counts = fig2.scan(deep, design, sets)
        design = dict(design) | {"_runs": runs, "_cover": cover}
        chosen["ALT"] = design
        nests["ALT"] = {
            "alternative_to": deep,
            "closest_fold_to": shallow,
            "slack_searched": ALT_SLACK,
            "design": {
                "sha256": hashlib.sha256(design["sequence"].encode()).hexdigest(),
                "cai": design["cai"],
                "mfe": design["mfe"],
                "sites": counts,
            },
        }
        drawn.append("ALT")
        alt_of = (deep, shallow)

    fig = plt.figure(figsize=(FIG_W, FIG_H))
    handles = legend_handles(census is not None, drawn, freq is None)
    legend = fig.legend(
        handles=handles,
        labels=[(h[0] if isinstance(h, tuple) else h).get_label() for h in handles],
        fontsize=LABEL,
        ncol=len(handles),
        loc="lower center",
        bbox_to_anchor=(0.5, MARGIN / PT_H),
        handlelength=0.9,
        columnspacing=0.9,
        handletextpad=0.4,
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.0)},
    )
    # (b)'s own legend used to carry the last three; it overlapped the zoom rectangle and
    # printed over twelve front points, while 137 pt beside this legend sat empty.
    ax_a = fig.add_axes(rect(MARGIN, 40, A_W, 100))
    if freq is not None:
        panel_frequency(ax_a, cat, freq)
    else:
        panel_catalog(ax_a, cat, census)
    ax_b = fig.add_axes(rect(PT_W / 2, 40, B_W, 100))
    b = panel_front(ax_b, fronts, hulls, floor, chosen)
    for nest in b["designs_outside"]:
        print(f"warning: {nest}'s drawn design is outside the zoom window", file=sys.stderr)
    if b["occluded"]:
        print(f"warning: the inset covers {b['occluded']} front points", file=sys.stderr)

    # The first draw is a measurement. Nothing below is a chosen coordinate.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # The legend's ink sat 9.6 pt off the bottom against a 4.9 pt top margin, so its own
    # padding is measured and taken out. What remains is the font's descent box, which the
    # bounding box reserves and these labels do not fill: the bottom margin lands at 7.1 pt
    # rather than 4.9, and closing that last 2.2 would mean rasterising to find the ink.
    pad = legend.get_tightbbox(renderer).y0 * 72 / fig.dpi - MARGIN
    legend.set_bbox_to_anchor((0.5, (MARGIN - pad) / PT_H))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_p, price_drawn = None, None
    if prices is not None or sweep is not None:
        ax_p = fig.add_axes(rect(PT_W - MARGIN - P_W, 40, P_W, 100))
        price_drawn = panel_price(ax_p, prices) if prices else panel_sweep(ax_p, sweep)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
    geom = pack(fig, ax_a, ax_b, legend, renderer, ax_p, len(drawn))
    # The inset's height is only settled once (b) is, and the front labels are separated
    # against it rather than against a constant.
    fig.canvas.draw()
    label_fronts(b["inset"], b["placed"], b["mfe_window"], fig, fig.canvas.get_renderer())
    span, angles = shared_span(chosen, lock=variant in ALT_IN)
    s_w = geom["width"]
    frame_h = s_w * span[1] / span[0]
    pt_per_unit = s_w / span[0]

    axes, jumps, marks, termini = [ax_a, ax_b], 0, {}, {}
    tail = [ax_p] if ax_p is not None else []
    for k, nest in enumerate(drawn):
        x = geom["left"] + k * (s_w + S_GAP)
        ax = fig.add_axes(rect(x, geom["panel_ab_bottom"], s_w, frame_h))
        axes.append(ax)
        j, term, mk = draw_structure(
            ax, chosen[nest], chosen[nest]["_runs"], span, angles[nest], pt_per_unit
        )
        jumps += j
        marks[nest], termini[nest] = mk, term
        if term["mode"] == "split":
            print(
                f"warning: {nest} draws two terminus marks, and the legend has one entry "
                f"(ends are {term['separation_pt']:.2f} pt apart)",
                file=sys.stderr,
            )

    # Each panel is compared with the one before it in the ladder; the alternative is
    # compared with the rung the ladder left, which is the panel before the one beside it.
    against = {}
    for i, nest in enumerate(drawn):
        earlier = alt_of[1] if nest == "ALT" else (drawn[i - 1] if i else None)
        if earlier is not None:
            moved = len(
                fig2.pair_set(chosen[earlier]["structure"])
                ^ fig2.pair_set(chosen[nest]["structure"])
            )
            against[nest] = rf"${moved}$ bp from ({chr(ord('c') + drawn.index(earlier))})"

    # Titles before the measurement, because they ARE the measurement: the band above a
    # structure is whatever its own title block needs, and reading the overhang before the
    # titles exist reads zero -- which put the panels high enough to clip off the canvas.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    titles = {}
    for ax, nest in zip(axes[2:], drawn, strict=True):
        key = alt_of[0] if nest == "ALT" else nest
        art, used = structure_titles(
            ax,
            fig,
            renderer,
            s_w,
            frame_h,
            r"F3$'$" if nest == "ALT" else nest,
            # The nest label's colour is the row's only key to the constraint, and (g)
            # carries F3's. Giving it body ink made a reader counting hues find five nests.
            NEST_INK[alt_of[0] if nest == "ALT" else nest],
            FORBIDDEN[key],
            against.get(nest),
        )
        titles[nest] = (art, used, FORBIDDEN[key])
    level = max(used for _, used, _ in titles.values())
    for art, _, forbids in titles.values():
        art.set_text(title_forms(forbids)[level])
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    # The structures' own top, measured off those titles, and their captions centred in
    # whatever band is left between the frame and the legend.
    caption_forms = {}
    s_top = geom["panel_ab_top"] - max(overhang(ax, fig, renderer)["t"] for ax in axes[2:])
    # Hung from a fixed edge under the frame, not centred in the leftover band. Centring
    # blocks of 2 to 4 lines put the MFE line -- the one quantity the row exists to compare
    # -- on three baselines spanning 6.1 pt, and no caption line in (c) shared a baseline
    # with any line in (d)-(g).
    caption_y = -CLEAR / frame_h
    captions = {}
    for ax, nest in zip(axes[2:], drawn, strict=True):
        ax.set_position(rect(ax.get_position().x0 * PT_W, s_top - frame_h, s_w, frame_h))
        captions[nest] = structure_caption(
            ax,
            chosen[nest],
            nests[nest]["design"]["sites"],
            caption_y,
            chosen["F0"]["mfe"],
            nest == "ALT",
            chosen[alt_of[0]]["mfe"] if nest == "ALT" else None,
            fig,
            renderer,
            s_w,
        )
    step = max(used for _, used, _ in captions.values())
    for nest, (art, _, forms) in captions.items():
        art.set_text(forms[step])
        caption_forms[nest] = art.get_text()
    geom |= {"structure_top": s_top, "structure_frame": [s_w, frame_h], "caption_y": caption_y}
    geom["panel_a_row_pitch"] = (geom["panel_ab_top"] - geom["panel_ab_bottom"]) / (
        ax_a.get_ylim()[1] - ax_a.get_ylim()[0]
    )

    # One baseline for every letter, flush with each panel's own leftmost ink.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    panels = axes + tail
    for letter, ax in zip(PANEL[: len(panels)], panels, strict=True):
        tb = ax.get_tightbbox(renderer).transformed(fig.transFigure.inverted())
        fig.text(tb.x0, LETTER_BASE / PT_H, rf"\textbf{{{letter}}}", fontsize=LETTER, color=INK)
    fig.savefig(out_pdf, format="pdf", metadata={"CreationDate": None})
    plt.close(fig)

    digests = motif_digests(root)
    record = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "layout": "row",
        "figsize_in": [FIG_W, FIG_H],
        "geometry_pt": {"units": "typographic points, y up from the canvas bottom"} | geom,
        "cell": cell,
        "length_nt": len(chosen["F0"]["structure"]),
        "cai_floor": floor,
        "panel_b_cai": b["range"],
        "panel_b_zoom": b["window"],
        "panel_b_zoom_mfe": b["mfe_window"],
        "panel_b_designs_outside_zoom": b["designs_outside"],
        "panel_b_points_under_inset": b["occluded"],
        "variant": variant,
        "structure_panels": drawn,
        # The gate, recomputed from what this variant actually drew. Only inks that appear
        # together are compared: the legend's swatches, and the tiers within one bar.
        "contrast_gate": contrast_gate(
            {
                # read against each other along one legend row
                "legend": {"cat " + r: PALETTE[r] for _, r, _ in fig2.CATEGORIES}
                | {"bar ground": BAR_BG, "body ink": INK},
                # read against each other inside one structure panel
                "structure": {"cat " + r: PALETTE[r] for _, r, _ in fig2.CATEGORIES}
                | {"backbone": BACKBONE, "body ink": INK},
                # read against each other inside one axis frame
                "furniture": {"bar ground": BAR_BG, "grid": GRID, "furniture": FURNITURE},
            }
            | (
                # read against each other along one bar, per category
                {
                    f"bar {r}": {
                        "inner": PALETTE[r],
                        "middle": tint(PALETTE[r], 0.55),
                        "outer": BAR_BG,
                    }
                    for _, r, _ in fig2.CATEGORIES
                }
                if census
                else {}
            )
        ),
        # Every panel's caption, including the alternative's -- a comprehension over the four
        # nests omitted the one panel with a caption the others do not have.
        "structure_caption_lines": {k: v.count(chr(10)) + 1 for k, v in caption_forms.items()},
        "optimum_census": census,
        "occurrence_frequency": freq,
        "set_prices": prices,
        "front_sweep": sweep,
        # The fifth DATA panel's own record. In a variant whose (g) is a structure this is
        # null, and the structure's record lives under `nests.ALT` -- so the key is named for
        # what it holds rather than for a panel letter that moves between variants.
        "fifth_data_panel": price_drawn,
        "structure_caption_form": {
            n: sum(1 for _ in caption_forms[n].split(chr(10))) for n in fig2.NESTS
        },
        "structure_layout": {0: "simple", 1: "naview", 2: "circular", 3: "turtle", 4: "RNApuzzler"}[
            PLOT_TYPE
        ],
        "structure_span": list(span),
        "separable_marks": marks,
        "termini": termini,
        "backbone_jumps": jumps,
        # Measured over whatever is drawn, not a literal: the interval was written for four
        # designs and the row now draws five, whose offset falls outside it.
        "energy_model": energy_model(chosen, drawn),
        "palette": PALETTE | {"F0": NEUTRAL, "backbone": BACKBONE},
        "motif_files": digests,
        "environment": environment(root, {r: v["path"] for r, v in digests.items()}),
        "catalog": cat,
        "nests": nests,
    }
    with open(out_json, "w") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)
    return record


def main():
    ap = argparse.ArgumentParser(description="Render Figure 2 as one row of six panels.")
    ap.add_argument("e2_dir")
    ap.add_argument("out_pdf")
    ap.add_argument("out_json")
    ap.add_argument("--cell", default="P26436")
    ap.add_argument("--cai-floor", type=float, default=0.87)
    ap.add_argument("--variant", choices=VARIANTS, default="plain")
    args = ap.parse_args()
    rec = build(args.e2_dir, args.cell, args.cai_floor, args.out_pdf, args.out_json, args.variant)
    for nest in fig2.NESTS:
        n = rec["nests"][nest]
        sites = n["design"]["sites"]
        drawn = rec["separable_marks"][nest]
        merged = ", ".join(
            f"{r} {sites[r]}->{drawn[r]}" for _, r, _ in fig2.CATEGORIES if drawn[r] != sites[r]
        )
        print(
            f"{nest}: {n['points']:5d} pts  hull {n['hull_recomputed']:4d} "
            f"(flagged {n['hull_flagged']:4d})  unresolved boxes {n['open_regions']}  "
            f"stopped_by {str(n['stopped_by']):5s}  "
            f"solved this round {str(n['gaps_solved_this_round']):4s}"
            f"{'  (resumed)' if n['resumed_from'] else ''}  "
            f"pick exact at the bound {n['pick_exact_at_floor']}  "
            f"sites {sites}"
            + (f"  MARKS MERGE: {merged}" if merged else "")
            + (f"  BLENDED PAIRS: {drawn['blended_pairs']}" if drawn["blended_pairs"] else "")
        )
    env, g = rec["environment"], rec["geometry_pt"]
    print(
        f"{rec['cell']} {rec['length_nt']} nt at CAI >= {rec['cai_floor']}, "
        f"{rec['structure_layout']}, {rec['backbone_jumps']} backbone jumps, "
        f"termini {rec['termini']['F0']['mode']} "
        f"({rec['termini']['F0']['separation_pt']:.2f} pt apart), "
        f"structure frame {g['structure_frame'][0]:.1f} x {g['structure_frame'][1]:.1f} pt, "
        f"(a) row pitch {g['panel_a_row_pitch']:.2f} pt; "
        f"matplotlib {env['matplotlib']}, seaborn {env['seaborn']}, "
        f"ViennaRNA {env['viennarna']}, usetex {env['usetex']}"
    )


if __name__ == "__main__":
    main()
