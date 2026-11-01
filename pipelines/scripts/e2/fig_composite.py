"""The composite figure: one row of seven panels, one palette, one legend.

A hue does two jobs, deliberately and in one direction: blue, orange and magenta name motif
sets I, II and III wherever a site appears, and they also name the constraint LEVEL that newly
forbids that set, so orange is both "a site of set II" in the catalog and "sets I and II both
forbidden" among the fronts and the cost boxes. The key gives both readings on one line rather
than picking one. Grey is the unconstrained design and, in the cost panel, its zero. GREEN means
one thing only: a base pair the unconstrained design's fold does not have. Nothing else is green.

(a) the motifs and how far they get, at a CAI floor, over the screen catalog; (b) what each
level costs against the unconstrained design at the same cutoff, over the finished set; (c)
one instance's four fronts with an inset at the cutoff; (d)-(g) that instance's four designs.
"""
import datetime
import hashlib
import json
import math
import os
import re
import sys
import types

import matplotlib
import numpy as np
from matplotlib import patheffects

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Rectangle

# Where the cached inputs are read from and where the plate is written. Both are arguments
# so that a rule can name the files it declares: the defaults are the working directory this
# figure was developed in, and a pipeline passes its own.
#
#     fig_composite.py LAYOUT INST SUBSET BVAR HEIGHT WIDTH [DATA_DIR] [OUT_PDF]
#
# DATA_DIR holds composite-data.json, example-INST.json, motif_names.json and (for the
# two-direction (b)) dcai-at-mfe.json; OUT_PDF names the PDF, and the PNG and the provenance
# sidecar are written beside it under the same stem.
OUT = sys.argv[7] if len(sys.argv) > 7 else os.environ["E2_CACHE"]
LAYOUT = int(sys.argv[1]) if len(sys.argv) > 1 else 1
INST = sys.argv[2] if len(sys.argv) > 2 else "Q6PF06-c80-none"
# Which proteins panel (b) draws. "finished" = every protein closed on all four rungs,
# exact or bounded. "certified" = only those exact on every rung, so that each box is a cost
# and not a bound on one. "all3" narrows the certified set to the proteins whose
# unconstrained front spells all three sets somewhere, i.e. where every rung has something to
# forbid; it stays AND-ed with exactness, since testing the all-three flag alone once let the
# retired fronts back in.
SUBSET = sys.argv[3] if len(sys.argv) > 3 else "finished"
# What panel (b) draws, and in which arrangement. "strip" is the cost boxes with the
# attribution strip over them. "side" and "stack" replace that panel with BOTH directions of
# the same trade -- fix CAI and read the structure given up, fix the structure and read the CAI
# given up -- on the proteins whose unconstrained front spells all three sets, side by side or
# one above the other. The two directions need two y axes with nothing in common, which is why
# they are two axes and not two series.
BVAR = sys.argv[4] if len(sys.argv) > 4 else "strip"
# "mfe" is the stacked pair's upper panel alone: fix CAI, read the structure given up. It also
# draws two of the four designs, the unconstrained one and the deepest rung, each at twice
# the width.
assert BVAR in ("strip", "side", "stack", "mfe"), BVAR
class _Unused(types.ModuleType):
    """fig2 imports ViennaRNA for a layout helper this figure never calls; any use raises."""

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        raise RuntimeError(f"fig2 touched RNA.{name} while drawing: no ViennaRNA here.")


sys.modules.setdefault("RNA", _Unused("RNA"))
# fig2 decides which front points are filled and how every fold is oriented. In the
# analysis repo it is this file's own neighbour; outside it, the checkout is named.
# realpath, not abspath: a working copy of this figure is reached through a symlink, and
# fig2 is this file's neighbour in the checkout, not the link's.
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import fig2  # noqa: E402

D = json.load(open(f"{OUT}/composite-data.json"))
R = json.load(open(f"{OUT}/example-{INST}.json"))
# The designs the fold row draws, as indices into the example's four panels.
SEL = [0, 3] if BVAR == "mfe" else [0, 1, 2, 3]
P = [R["panels"][:4][i] for i in SEL]
NAMES = json.load(open(os.environ.get("E2_MOTIF_NAMES", f"{OUT}/motif_names.json")))

# hue = which set, and equally the level that newly forbids it; green = a pair that changed.
# P0's palette, and P0's rule: a colour means a motif category, and a constraint level is
# drawn in the colour of the category it newly forbids, so one key serves every panel. The
# unconstrained level forbids nothing and is therefore grey.
SETC = {"I": "#2a78d6", "II": "#d1651b", "III": "#8f4fd1"}
PLAIN = "#9b9a96"
LEVEL = [PLAIN, SETC["I"], SETC["II"], SETC["III"]]
GAINED = "#1f7a4d"
# What only one of (d)'s two designs has -- a base pair, a codon -- is dark, and what the two
# share is light grey. Hue is left to the motif sets, so (d) never uses one colour for both.
UNIQ, SHAREDC = "#33322f", "#d9d8d4"
# The annotation ink for "this design is drawn below". Deliberately outside the four rung
# colours and outside green, so the mark cannot be read as data on any front.
MARKC, MARKHALO = "#0b0b0b", "#ffffff"
BACK, PAIR = "#c2c1bd", "#e0dfdb"
INK, MUTED, GRID, OUTLINE = "#0b0b0b", "#4a4742", "#e4e3df", "#8e8a82"
RUNG = ["no sets", "+ I", "+ I, II", "+ I, II, III"]
# The analysis repo draws the paper's figures at FS, FSS = 7, 6.3 pt (e1/fig4_overview.py,
# which is the paper's Figure 4 and the only built renderer at the text width). Layout 3 adopts
# that scale, so every hand-set size in this file is written as a multiple of the 6.0 pt base
# and scaled once here rather than being re-chosen panel by panel.
K = (7.0 / 6.0) if LAYOUT >= 3 else 1.0
# The fold captions, at the x tick size. Layout 4 sets two lines on every panel and the
# others choose one or two by measuring; either way this is the floor, since a caption a
# notch smaller buys about a tenth of an inch of drawing and is the smallest type on a page
# whose ticks are already 6.3 pt.
CAPFS = 5.4 * K
# The panel letters. The analysis repo sets its own at FS, the figure's own label size, and
# bold -- e1/fig4_overview.py draws them as a bold "(a)" at 7 pt beside 7 pt axis labels.
# This file had them at 8.9, larger than anything else on the page, which on a half-inch
# panel is a letter as tall as the box it names.
LETFS = (6.0 if LAYOUT >= 3 else 7.6) * K
# Per-fold rotation, in radians, added to the shared fill angle. Positive is anticlockwise.
# Keyed by panel index, so {3: ...} turns the deepest rung only. This is a drawing choice and
# nothing else: a rigid rotation leaves every distance, pair and layout untouched. The two-design
# variant draws both folds in the same orientation and uses none.
TWEAK = {3: -math.radians(20)}
# (c)'s zoom window, in its axes fraction: left, bottom, width, height. Named because the
# parent's own reference line has to know where the inset stands in front of it.
# The two-design variant lets the zoom rise past the top of the panel, which is empty above it.
INSET = (.07, .48, .46, .76) if BVAR == "mfe" else (.07, .50, .46, .46)
# The paper's own serif: Figure 4 sets its type with LaTeX and mathptmx, which is Times, and
# Nimbus Roman is that face. STIXGeneral follows for any glyph Nimbus lacks.
plt.rcParams.update({"font.family": "serif", "font.serif": ["Nimbus Roman", "STIXGeneral"],
                     "mathtext.fontset": "stix"})
plt.rcParams.update({"font.size": 6.0 * K, "axes.labelsize": 6.0 * K,
                     # One tick size on both axes. At 5.4 and 5.2 the x digits set 20 px
                     # and the y digits 18 on the same panel, a tenth apart for no reason
                     # a reader can see, and (c)'s bottom-left corner shows both at once.
                     "xtick.labelsize": 5.4 * K, "ytick.labelsize": 5.4 * K,
                     "axes.edgecolor": MUTED, "axes.linewidth": .5,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
                     "axes.labelcolor": INK, "figure.dpi": 200, "pdf.fonttype": 42,
                     # No tight box: it crops to the ink and then adds its own pad, so the
                     # output was 3.04 in for a figure laid out at 3.00 and the margins were
                     # whatever the longest tick label left over. The boxes below place every
                     # artist inside the canvas, which the overlap check verifies, so the
                     # figure is saved at exactly the size it declares.
                     "savefig.bbox": None})


def tint(h, k):
    r, g, b = (int(h[i:i + 2], 16) for i in (1, 3, 5))
    return "#{:02x}{:02x}{:02x}".format(*(round(v + (255 - v) * k) for v in (r, g, b)))


TAGS = {}
PENDING = []   # (axes, letter, kind) queued by tag() until the panels have been drawn


def tag(ax, letter, dx=-.085, dy=1.035, kind="top"):
    """Panel letters.

    `kind` says which group a panel is in, explicitly. It used to be inferred from whether the
    axes sat above the figure's half-height, which silently mislabels the moment a band is
    resized: growing the fold row pushed it over the line and all four fold letters jumped onto
    the top row's baseline, on top of (b) and (c).
    """
    if LAYOUT >= 2:
        # Deferred: a letter goes at its panel's VISUAL left edge, which includes the tick
        # labels and the axis label and so is not known until the panel has been drawn once.
        PENDING.append((ax, letter, kind))
        return
    ax.annotate(letter, xy=(dx, dy), xycoords="axes fraction", fontsize=LETFS, color=INK,
                ha="left", va="bottom", fontweight="bold")


def _ink_x1(artist):
    """The right edge of an artist's drawn pixels, in device pixels.

    `get_window_extent` is the layout box, which is set from the font's ascent and descent
    and not from the glyphs: two labels can share it and still end at visibly different
    places. This reads the rendered canvas inside that box instead, so the number is what
    the page shows. The figure must have been drawn.
    """
    bb = artist.get_window_extent(artist.figure.canvas.get_renderer())
    buf = np.asarray(artist.figure.canvas.buffer_rgba())
    h, w = buf.shape[:2]
    x0, x1 = max(0, int(bb.x0) - 2), min(w, int(bb.x1) + 3)
    y0, y1 = max(0, h - int(bb.y1) - 2), min(h, h - int(bb.y0) + 3)
    cols = np.where((buf[y0:y1, x0:x1, :3].min(axis=2) < 250).any(axis=0))[0]
    return x0 + cols.max() if len(cols) else bb.x1


def _ink_x0(artist):
    """The left edge of an artist's drawn pixels, as a figure fraction. See _ink_x1."""
    bb = artist.get_window_extent(artist.figure.canvas.get_renderer())
    buf = np.asarray(artist.figure.canvas.buffer_rgba())
    h, w = buf.shape[:2]
    x0, x1 = max(0, int(bb.x0) - 2), min(w, int(bb.x1) + 3)
    y0, y1 = max(0, h - int(bb.y1) - 2), min(h, h - int(bb.y0) + 3)
    cols = np.where((buf[y0:y1, x0:x1, :3].min(axis=2) < 250).any(axis=0))[0]
    return (x0 + cols.min() if len(cols) else bb.x0) / w


def visual_x0(ax):
    """The left edge of everything the panel draws, in figure coordinates.

    `get_position` gives the data box, which on a panel with a y axis sits well to the right
    of the leftmost ink: for (b) the gap is its tick labels plus "<greek>MFE (kcal/mol)". A
    letter placed on the data box therefore floats inside the panel rather than marking where
    the panel starts, and the block as a whole reads as indented from its own label column.
    """
    bb = ax.get_tightbbox(ax.figure.canvas.get_renderer())
    return bb.transformed(ax.figure.transFigure.inverted()).x0


def place_tags():
    """Draw the deferred letters, top row on the visual left edge and folds on their own box.

    The fold row is positioned so that (d) starts at (b)'s visual left and (f) at (c)'s, so
    those letters coincide without being special-cased here.
    """
    fig.canvas.draw()
    for ax, letter, kind in PENDING:
        # A fold panel has no tick labels and no axis label, so its box IS its ink; every
        # other panel's letter goes on the ink, not on the data box.
        x = ax.get_position().x0 if kind == "fold" else visual_x0(ax)
        # In the two-design variant a panel with a y label starts at that label's ink, not at
        # the label's layout box, which is wider than the glyphs.
        if BVAR == "mfe" and kind != "fold" and ax.yaxis.label.get_text().strip():
            # Two pixels left of the ink, at the 200 dpi the layout is measured on.
            x = _ink_x0(ax.yaxis.label) - 2 / (fig.get_size_inches()[0] * fig.dpi)
        y = TAGY_K[kind]
        TAGS[id(fig.text(x, y, letter, fontsize=LETFS, color=INK,
                         ha="left", va="bottom", fontweight="bold"))] = f"{id(ax)}"


def panel_a(ax, lab):
    """Three nested shares per motif as overlaid bars, inside a full-width 100% reference box.

    The bars are overlaid, not stacked: each starts at zero and the longest is the outermost
    share. The reference box is drawn as four explicit edges, because its right-hand side used
    to be the axes spine and so was invisible.
    """
    A = D["a"]
    ent, n = A["entries"], A["n"]
    cen = A["census"][str(A["floor"])]
    # Each class opens with the union of its motifs: the proteins with some motif of the class in
    # a CDS, on the front, and in the design at the floor.
    rows_ = []
    for roman in ("I", "II", "III"):
        cl = next(c for c in A["classes"] if c["category"] == roman)
        rows_.append({"cat": roman, "label": r"$\mathcal{M}_{\mathrm{%s}}$" % roman,
                      "sp": cl["spellable"], "oc": cl["occurs"], "cn": cen["classes"][roman],
                      "union": True})
        for e in ent:
            if e["category"] == roman:
                rows_.append({"cat": roman, "union": False, "sp": e["spellable"],
                              "oc": e["occurs"], "cn": cen["counts"][e["motif"]],
                              "label": NAMES.get(e["motif"], e["motif"]).replace(
                                  ", reverse complement", " (rc)")})
    y = list(range(len(rows_)))[::-1]
    for row, e in zip(y, rows_, strict=True):
        hue = SETC[e["cat"]]
        # .76, not .84: at .84 this bar is 8% from white and its terminus against the white
        # remainder inside the frame is the weakest edge in the figure -- and for the rows
        # where the other two shares are slivers, that terminus is the quantity to read.
        ax.barh(row, 100 * e["sp"] / n, height=.80, color=tint(hue, .76), lw=0, zorder=2)
        ax.barh(row, 100 * e["oc"] / n, height=.50, color=tint(hue, .52), lw=0, zorder=3)
        # No minimum width: AAAAAA is spelled by 0 of 85 optima, and a floored bar would
        # draw that zero as indistinguishable from poly(U) run's real 3 of 85.
        ax.barh(row, 100 * e["cn"] / cen["proteins"], height=.22, color=hue, lw=0, zorder=4)
        ax.add_patch(Rectangle((0, row - .40), 100, .80, fill=False, ec=OUTLINE, lw=.5,
                               zorder=5, clip_on=False, joinstyle="miter"))
    ax.set_yticks(y)
    ax.set_yticklabels([e["label"] for e in rows_], fontsize=lab)
    for t, e in zip(ax.get_yticklabels(), rows_, strict=True):
        t.set_color(SETC[e["cat"]])
    ent = rows_
    ax.set_ylim(-.8, len(ent) - .2)
    ax.set_xlim(0, 101)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel(f"% of {n} proteins")
    ax.tick_params(axis="y", length=0, pad=3.6)   # the left spine is hidden, so the tick
    ax.grid(axis="x", color=GRID, lw=.35)         # dashes read as leader lines
    if LAYOUT == 4:
        # (a)'s two foot rows are a line apart and the fold captions' two lines are a
        # shorter line apart, so the pair cannot share both baselines once the lower one
        # is fixed. Moving the tick row down and the axis label up by the same 8 px puts
        # BOTH on the captions' lines instead of leaving the upper one 8 px off, which at
        # this size reads as an alignment that failed rather than one that was not tried.
        ax.tick_params(axis="x", pad=5.4)
        ax.xaxis.labelpad = 2.1
    ax.set_axisbelow(True)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    tag(ax, "(a)", dx=-.19)


def panel_b(ax, strip):
    """The cost of each level, over the proteins whose four fronts are all finished.

    Three populations are selectable, and they differ in what a box is allowed to claim.
    "finished" draws every protein the band owner reports as closed on all four rungs.
    "certified" keeps only those exact on every rung; on a retired front the best design at a
    floor is an upper bound, so a cost taken from one is an upper bound too, and the panel is
    then a statement about costs rather than about bounds on costs. "all3" narrows further to
    the proteins whose unconstrained front spells all three sets somewhere.

    Each is a filter on a stored flag, never a list of names, so the membership follows the
    data rather than this docstring -- and the figure's provenance sidecar records the members
    and a digest of them, since the count alone cannot report one protein leaving as another
    arrives.
    """
    B, EX = D["b_kcal"], D["b_exact"]
    cuts = sorted(B)
    A3 = D["b_all3"]
    WHY = D["b_why"]
    # "all3" stays nested inside "certified" rather than standing alone: it exists to ask
    # whether the cost looks different on the proteins that can spell every set, and that
    # question is only worth asking where the costs are exact.
    keep = {c: [i for i in range(len(EX[c]))
                if (SUBSET == "finished"
                    or (EX[c][i] and (SUBSET == "certified" or A3[c][i])))]
            for c in cuts}
    nn = len(keep[cuts[0]])
    # The boxes price only the proteins the rung BINDS on, which is what makes the middle
    # box a price rather than a mixture: at the 0.95 floor set II binds on 1 protein of 43,
    # is absent from 20 fronts entirely, and has already gone from the F1 pick on the other
    # 22. (A non-binding cell is NOT a zero here -- b_kcal is cumulative against F0 while
    # b_why is incremental, so a non-binding F2 carries F1's price. That is the reason to
    # filter rather than to average: the unfiltered middle box is the left box redrawn, the
    # same 43 numbers with one substitution at the 0.95 floor.)
    val = {c: {n: [B[c][n][i] for i in keep[c] if WHY[c][n][i] == "binds"]
               for n in ("F1", "F2", "F3")} for c in cuts}
    for i, n in enumerate(("F1", "F2", "F3")):
        pos = [k + (i - 1) * .26 for k in range(len(cuts))]
        bp = ax.boxplot([val[c][n] for c in cuts], positions=pos, widths=.22, patch_artist=True,
                        showfliers=True, medianprops=dict(color="#ffffff", lw=.9),
                        whiskerprops=dict(color=MUTED, lw=.5), capprops=dict(color=MUTED, lw=.5),
                        flierprops=dict(marker="o", ms=1.6, mfc="none", mec=MUTED, mew=.4),
                        boxprops=dict(lw=.4, edgecolor=MUTED))
        for b in bp["boxes"]:
            b.set_facecolor(LEVEL[i + 1])
            b.set_alpha(.9)
    ax.axhline(0, color=LEVEL[0], lw=.9, zorder=1)
    # The strip above: every protein attributed, so the box below cannot be read without
    # seeing how many proteins it is computed over.
    for i, n in enumerate(("F1", "F2", "F3")):
        pos = [k + (i - 1) * .26 for k in range(len(cuts))]
        tot = [len(keep[c]) for c in cuts]
        nb = [sum(1 for j in keep[c] if WHY[c][n][j] == "binds") for c in cuts]
        ng = [sum(1 for j in keep[c] if WHY[c][n][j] == "gone") for c in cuts]
        nv = [t - b_ - g for t, b_, g in zip(tot, nb, ng, strict=True)]
        strip.bar(pos, nb, width=.22, color=LEVEL[i + 1], lw=0)
        strip.bar(pos, ng, width=.22, bottom=nb, color=tint(LEVEL[i + 1], .55), lw=0)
        strip.bar(pos, nv, width=.22, bottom=[b_ + g for b_, g in zip(nb, ng, strict=True)],
                  color="#e9e7e2", lw=0)
    strip.set_xticks(range(len(cuts)))
    strip.set_xticklabels([])
    strip.set_xlim(-.45, len(cuts) - .55)
    strip.set_ylim(0, nn)
    # No ticks and no label on the strip: each bar is the full population by construction, the
    # boxes below already print that n, and at this size a tick label and an axis label here
    # both land on the panel letter.
    strip.set_yticks([])
    strip.grid(axis="y", color=GRID, lw=.3)
    strip.set_axisbelow(True)
    for sp in ("top", "right"):
        strip.spines[sp].set_visible(False)
    ax.set_xticks(range(len(cuts)))
    ax.set_xticklabels([f"{float(c):.2f}" for c in cuts])
    # boxplot's own limits leave a group's full width of air outside the end groups, which
    # was 17.3% of this panel.
    ax.set_xlim(-.45, len(cuts) - .55)
    ax.set_xlabel("CAI")
    ax.set_ylabel("\u0394MFE (kcal/mol)")
    ax.grid(axis="y", color=GRID, lw=.35)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    # The population, on the panel. (a) carries its 85 on its own axis; this is a different
    # and smaller set, and two bases side by side with one of them unstated is how one
    # measurement gets read as two that disagree.
    ax.annotate(f"n = {nn}", xy=(.015, .97), xycoords="axes fraction", fontsize=5.0 * K,
                color=MUTED, ha="left", va="top")
    # No per-group count printed here: at this size it lands on the zero line and inside the
    # boxes, and the strip above already carries the same information as a proportion.
    if LAYOUT == 2:
        ax.xaxis.labelpad = 2
    # Layout 1 puts the letter just above its axes, so it must hang off the STRIP, which is
    # now the topmost of (b)'s two boxes; layout 2 places letters from a panel's visual left
    # edge, and that edge is the boxes' y label, not the strip's bare frame.
    tag(strip if LAYOUT == 1 else ax, "(b)")


def _inw(ax):
    return ax.get_position().width * ax.figure.get_size_inches()[0]


def _inh(ax):
    return ax.get_position().height * ax.figure.get_size_inches()[1]


def panel_b2(axm, axc):
    """Both directions of the trade, over whichever population SUBSET names.

    `axm` fixes CAI at a floor and reads the structure given up; `axc` fixes the structure at a
    share of the protein's own unconstrained optimum and reads the CAI given up. The two
    quantities have no common unit, so they get two axes rather than two series on one.

    The default population is every protein closed on all four rungs: 43. "certified" keeps
    the 29 that are exact on every rung, since on a retired front the best design at a floor
    is an upper bound and so is a cost taken from it. "all3" narrows that further to the ones
    whose unconstrained front spells each of the three sets somewhere -- the filter that makes
    every rung non-vacuous, since on a protein whose front never spells set II the middle rung
    forbids nothing. It is nested inside "certified" rather than standing alone, so it draws
    9 and not the 22 that carry all three sets. The panel prints its own n, because these
    differ by a factor of five.
    """
    B, EX, A3, CELLS = D["b_kcal"], D["b_exact"], D["b_all3"], D["b_cells"]
    cuts = sorted(B)
    DC = json.load(open(f"{OUT}/dcai-at-mfe.json")) if axc is not None else None
    pcts = [f"{p:.2f}" for p in DC["pcts"]] if axc is not None else []
    keep = {c: [i for i in range(len(A3[c]))
                if (SUBSET == "finished"
                    or (EX[c][i] and (SUBSET == "certified" or A3[c][i])))]
            for c in cuts}
    members = {CELLS[cuts[0]][i] for i in keep[cuts[0]]}
    nn = len(members)
    nex = sum(1 for i in keep[cuts[0]] if EX[cuts[0]][i])
    vm = {c: {n: [B[c][n][i] for i in keep[c]] for n in ("F1", "F2", "F3")} for c in cuts}
    # Matched by accession, not by position: the two files were built by different scripts
    # over the same corpus, and a shared length is not a shared order.
    kd = {k: [i for i, cell in enumerate(DC["cells"][k]) if cell in members] for k in pcts}
    vc = {k: {n: [DC["d"][k][n][i] for i in kd[k] if DC["d"][k][n][i] is not None]
              for n in ("F1", "F2", "F3")} for k in pcts}
    # The printed n is one number over eight groups of boxes, and three things could make it
    # false without changing anything visible: a cutoff-dependent filter, a rung that cannot
    # reach a structure bound at all (dcai-at-mfe.json keeps an `infeasible` list for that,
    # and a member of it is dropped from vc and from nothing else), and the two files being
    # matched on a shared length rather than on accessions. All three are empty today. A
    # count that is true by circumstance is worth an assertion, not a comment.
    for _c in cuts:
        assert len(keep[_c]) == nn, f"(b): {len(keep[_c])} at cut {_c}, not {nn}"
    for _k in pcts:
        for _n in ("F1", "F2", "F3"):
            assert len(vm[cuts[0]][_n]) == nn and len(vc[_k][_n]) == nn, \
                f"(b): {_k} {_n} draws {len(vc[_k][_n])} of {nn}"

    rng = np.random.default_rng(11)

    def boxes(ax, data, keys):
        for i, n in enumerate(("F1", "F2", "F3")):
            pos = [k + (i - 1) * .26 for k in range(len(keys))]
            bp = ax.boxplot([data[c][n] for c in keys], positions=pos, widths=.22,
                            patch_artist=True,
                            # No fliers: the strip below draws every protein, so a flier is
                            # the same point twice, once as data and once as an exception.
                            showfliers=False,
                            medianprops=dict(color=INK, lw=.8),
                            whiskerprops=dict(color=MUTED, lw=.5),
                            capprops=dict(color=MUTED, lw=.5),
                            boxprops=dict(lw=.5), zorder=3)
            for bx in bp["boxes"]:
                # Pale fill with the rung's colour as the edge, so the strip over it shows
                # through: a solid box hides the very points it is drawn over, and a box
                # that hides its own data is a worse summary than a box alone.
                bx.set_facecolor(tint(LEVEL[i + 1], .74))
                bx.set_edgecolor(LEVEL[i + 1])
            # One point per protein behind each box. A box is five numbers; it cannot say
            # how many it was computed over, nor whether they are one cloud or two. At
            # n = 43 they fit, jittered, inside the box's own width.
            for q, c in zip(pos, keys, strict=True):
                v = data[c][n]
                if v:
                    ax.scatter(q + rng.uniform(-.080, .080, len(v)), v, s=1.1,
                               c=LEVEL[i + 1], alpha=.75, linewidths=0, zorder=4)
        ax.axhline(0, color=LEVEL[0], lw=.7, zorder=.9)
        # Room under the zero line. With the data starting at 0 it was drawn 2 px off the
        # bottom spine, so the pair read as one thick rule with no scale, and any point
        # sitting at zero was veiled by it.
        top = max(v for c in keys for n in ("F1", "F2", "F3") for v in data[c][n])
        ax.set_ylim(-.055 * top, top * 1.08)
        # The zero line IS the baseline, so the frame does not draw a second one under it.
        # With both, a 3 px grey rule at zero sat 5 px above a 2 px dark spine with white
        # between, the "0" tick label centred on the upper of the two, and a reader had to
        # work out which line the scale started from.
        ax.spines["bottom"].set_visible(False)
        ax.set_xticks(range(len(keys)))
        ax.set_xlim(-.45, len(keys) - .55)
        ax.grid(axis="y", color=GRID, lw=.35)
        ax.set_axisbelow(True)
        ax.xaxis.labelpad = 2
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    boxes(axm, vm, cuts)
    if axc is not None:
        boxes(axc, vc, pcts)
    axm.set_xticklabels([f"{float(c):.2f}" for c in cuts])
    if axc is not None:
        axc.set_xticklabels([f"{100 * float(k):.0f}%" for k in pcts])
    # Short y labels on purpose: stacked, these boxes are shorter than a long label, which
    # then reaches past the top of its own panel and into the panel letter. "against F0" moves
    # to the x label, where there is width for it.
    # On a panel half an inch tall, "\u0394MFE (kcal/mol)" is taller than its own axes: being
    # centred, it then reaches above the panel's top and prints under the panel letter. Where
    # that is so the unit moves into the corner note, which has width to spare.
    SHORT = LAYOUT == 4 and BVAR == "stack"
    # Two lines where one does not fit: "\u0394MFE (kcal/mol)" is 0.77 in of type, longer than
    # a stacked panel is tall, but stacked it is two 0.36 in lines across the label column,
    # which has the width. The corner note used to carry the unit and was asked to go.
    axm.set_ylabel("\u0394MFE\n(kcal/mol)" if SHORT else "\u0394MFE (kcal/mol)")
    if axc is not None:
        axc.set_ylabel("\u0394CAI")
    # Half-width panels cannot carry the explanatory line: at 0.9 in it is three times the
    # panel's width and lands on its neighbours and on the structures below. There the two
    # directions are named by the y axes and in the caption instead.
    if LAYOUT == 4:
        # The label sits at the right end of its own tick row rather than on a line of its
        # own. Below the ticks it costs a whole line of page -- 0.18 in -- between two
        # stacked panels that are half an inch tall each; beside them it costs nothing, and
        # at two and three characters there is room for it there.
        for ax, lab in ((axm, "CAI"), (axc, "MFE%"))[:1 if axc is None else 2]:
            ax.set_xlabel(lab)
            # On the tick row and centred, where the panel is wide enough for it: with four
            # ticks the middle of the axis falls between the second and third of them, and
            # on a 2.4 in panel that gap holds four characters with 0.13 in to spare. Under
            # the ticks the label costs a whole line of page, which is what two stacked
            # panels half an inch tall cannot pay.
            #
            # Side by side the panels are 0.99 in, the same gap is 0.13, and the label would
            # sit on the ticks either side of it -- so there it goes back under the row,
            # which that arrangement can afford because its panels are three times taller.
            if _inw(ax) >= 1.8:
                # Below the tick row, not in it. At 0.072 in the label's cap height
                # overlapped the tick band by 8 px, so "MFE%" -- which ends in a percent
                # sign, like all four of its neighbours -- sat between "85%" and "90%" and
                # the row read as five ticks. At 0.142 in it clears its own tick row by
                # the 4.3 px that (c)'s label clears (c)'s, and it stays centred in the gap.
                ax.xaxis.set_label_coords(.5, -.142 / _inh(ax))
                ax.xaxis.label.set_va("top")
            ax.xaxis.label.set_ha("center")
    else:
        axm.set_xlabel("CAI floor\nfix CAI, read the structure given up")
        axc.set_xlabel("MFE floor, share of own optimum"
                       "\nfix structure, read the CAI given up")
    if LAYOUT == 4:
        # Turning the y tick labels a quarter was tried and is wrong here, which is worth
        # recording: on their side they are 0.07 in of width instead of 0.21, but their
        # HEIGHT becomes their text width, 0.21 in each on a panel 0.47 in tall, so three of
        # them do not fit and they print through one another. Height is the scarce dimension
        # on this page, so the rotation spends the scarce one to save the plentiful one.
        # Three ticks and a short format buy the same width with no height at all.
        for ax in (axm,) if axc is None else (axm, axc):
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3, prune=None))
            ax.tick_params(axis="y", pad=1.2)
            # Both y labels at one x, so the pair reads as a column: left to itself each sits
            # outside its own tick labels, and the two tick columns are different widths, so
            # "\u0394MFE" and "\u0394CAI" came out a tenth of an inch apart.
            ax.yaxis.set_label_coords(-.292 / _inw(ax), .5)
            # Hard against its own axis. At the default pad the upper panel's tick row sat
            # 14 px under its own spine and 9 px over the next panel's frame, so the row --
            # and the "CAI" in it, which has no tick of its own beneath it -- read as
            # belonging to the panel below, which is the CAI one.
            # No tick marks, only labels. The zero rule is this panel's baseline and the
            # bottom spine is off, so a tick had nothing to hang from: four 3 x 8 px
            # dashes floating 5 px under the rule, half the length of (a)'s and attached
            # to nothing. The four categories are named by their labels.
            ax.tick_params(axis="x", pad=2.2, length=0)
    # The population, on the panel. (a) states its own share over the screen catalog and this
    # is a different set, so neither n can be read off the other.
    # The population and nothing else. The unit moved onto the y label, which carries it on
    # a second line; the exactness -- 14 of the 43 are upper bounds taken from a retired
    # front, two of them provably loose, since Q9HC38 at CAI 0.90 prices F1 at 26.0 and F2
    # at 17.7, which a nested feasible set cannot do -- belongs to the caption, which is
    # where this plate states what a bounded cost is.
    # The "mfe" variant leaves the count to the caption.
    if BVAR != "mfe":
        axm.annotate(f"n = {nn}" if LAYOUT == 4 else
                     f"n = {nn} of {len(A3[cuts[0]])}\n{nex} exact on every rung"
                     "\nboth panels are against F0",
                     xy=(.015, .97), xycoords="axes fraction", fontsize=5.0 * K,
                     color=MUTED, ha="left", va="top", linespacing=1.35)
    tag(axm, "(b)")


def _hull(front):
    """Hull membership, recomputed with the analysis repo's own definition.

    The hull is convex in the scalarized space, not in (CAI, MFE), so it is fig2's
    hull_vertices that decides it and not a local reimplementation.
    """
    import math
    pts = [{"cai_log": math.log(c), "mfe": m} for c, m in front]
    return fig2.hull_vertices(pts)


def _mark(a, size, lw):
    """Mark the four designs the fold panels draw.

    A hollow ring in the rung's own colour was tried first and read as a front point that
    happened to be larger: the fronts are already drawn as filled and hollow circles in these
    four colours, so a fifth circle in the same palette says "another point" rather than "this
    one". A square in one neutral ink does not belong to that vocabulary and so reads as an
    annotation, which is what it is -- and one colour for all four is honest, since what the
    mark means is the same on every front.

    Two strokes: a white square under the dark one, so a mark keeps its edge where it lands on
    the dense band of front points rather than dissolving into it.
    """
    for k, p in enumerate(P):
        for c, w, z in ((MARKHALO, lw + .9, 20), (MARKC, lw, 24)):
            a.plot([p["cai"]], [p["mfe"]], marker="s", ms=size, mfc="none", mec=c,
                   mew=w, ls="none", zorder=z + k, clip_on=False)


def panel_c(ax):
    """The fronts as points, drawn as P0 draws them: filled on the hull, hollow off it."""
    F = {n: sorted(R["fronts"][n]) for n in ("F0", "F1", "F2", "F3")}
    lo = min(c for f in F.values() for c, _ in f)
    hi = max(c for f in F.values() for c, _ in f)
    # P0's order, F0 to F3. The four do not coincide -- over this corpus the gap to the
    # unconstrained front is 1.0% of the MFE at the median and 8.9% at the worst, up to
    # 36.6 kcal/mol --
    # but they run close enough along most of the range that whichever is drawn last covers
    # the others there. The inset is where the figure separates them at the floor it uses.
    for k in range(4):
        n = ("F0", "F1", "F2", "F3")[k]
        pts, on = F[n], _hull(F[n])
        off = [i for i in range(len(pts)) if i not in on]
        ax.scatter([pts[i][0] for i in off], [pts[i][1] for i in off], s=1.4,
                   facecolors="none", edgecolors=LEVEL[k], lw=.28, zorder=3 + k * 3)
        ins_ = sorted(on)
        ax.plot([pts[i][0] for i in ins_], [pts[i][1] for i in ins_], color=LEVEL[k], lw=.7,
                zorder=4 + k * 3)
        ax.scatter([pts[i][0] for i in ins_], [pts[i][1] for i in ins_], s=2.6, c=LEVEL[k],
                   lw=0, zorder=5 + k * 3)
    # Stopped at the inset's foot where the inset stands in front of it. The inset is
    # opaque and covers the middle of this line, leaving a 13 px stub above it at the top
    # of the panel -- isolated, with the rule not resuming for another 140 px, so it reads
    # as a stray tick rather than as part of a dashed reference line.
    pad = .02 * (hi - lo)
    _bf = (R["bound"] - (lo - pad)) / ((hi + pad) - (lo - pad))
    ax.axvline(R["bound"], color=INK, lw=.5, ls=(0, (2.4, 1.6)), zorder=6,
               ymax=INSET[1] if INSET[0] <= _bf <= INSET[0] + INSET[2] else 1.0)
    # The four designs are marked in the INSET only. At full scale they sit at one floor
    # within 9.4 kcal/mol of each other, so four marks there overprint into a single smudge
    # on the densest part of the plot: it adds ink and resolves nothing. The window below is
    # the place where the four are separable, so it is the place to name them.
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_xlabel("CAI")
    if LAYOUT in (2, 4):
        # Centred on its own panel, which is where a panel's x label belongs. It was pushed
        # left to clear (g)'s letter, back when that letter sat on a line of its own between
        # the rows; the letters are on their panels' top edges now, so the label row below
        # (c) is clear and the offset was buying nothing but an off-centre label.
        ax.xaxis.label.set_x(.5 if LAYOUT == 4 else .38)
        # Hard against its own tick row. At labelpad 2 it sat 16 px under (c)'s ticks and
        # 16 px over (g)'s letter -- equidistant between the panel it belongs to and the
        # one below it, which at this size reads as a title line for (g).
        ax.xaxis.labelpad = 1.6 if LAYOUT == 4 else 2
    if LAYOUT == 4:
        # (b)'s tick geometry, so the two panels' label stacks land on one row. At
        # matplotlib's default -- a 15 px tick and a 15 px pad against (b)'s 0 and 4 --
        # (c)'s tick row sat where (b)'s axis label was, and the foot of the two panels
        # read as four staggered rows instead of two.
        ax.tick_params(axis="x", pad=2.2, length=0)
    ax.set_ylabel("MFE (kcal/mol)")
    ax.grid(color=GRID, lw=.4)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    # One protein, named. (a) states its population on its axis and (b) prints its own n;
    # without this, (c)-(g) read as either of those populations rather than as a single
    # instance drawn from them.
    # On white: the MFE gridline ran 4 px under the cap tops and the CAI = 1.00 gridline
    # through the last character, so the name read as struck through.
    ax.annotate(R["cell"], xy=(.985, .04), xycoords="axes fraction", fontsize=5.0 * K,
                color=MUTED, ha="right", va="bottom",
                bbox=dict(fc="#ffffff", ec="none", pad=.8))
    # Close enough over most of the range to overprint, far apart at the high-CAI end.
    # The window is 8% of the CAI span either side of the floor, not 1.2%. A front's points
    # are not spread evenly along CAI -- they bunch towards CAI 1, where synonymous choices
    # are nearly forced -- so down at the 0.80 floor the narrow window held 3, 4, 5 and 13
    # points of the four fronts and the inset was mostly empty. The same window at 8% holds
    # 25, 28, 27 and 44, which is what the inset is for.
    half = .08 * (hi - lo)
    ins = ax.inset_axes(INSET)
    for k in range(4):
        n = ("F0", "F1", "F2", "F3")[k]
        pts, on = F[n], _hull(F[n])
        near = [(i, q) for i, q in enumerate(pts) if R["bound"] - half <= q[0] <= R["bound"] + half]
        ins.scatter([q[0] for _, q in near], [q[1] for _, q in near], s=2.0,
                    c=[LEVEL[k] if i in on else "none" for i, _ in near],
                    edgecolors=LEVEL[k], lw=.25, zorder=3)
    ins.axvline(R["bound"], color=INK, lw=.4, ls=(0, (2.0, 1.4)))
    _mark(ins, 3.4, .45)
    ins.set_xlim(R["bound"] - half, R["bound"] + half)
    # The y window follows the points the x window admits, so widening one widens the other
    # and nothing drawn falls outside the frame.
    ys = [m for f in F.values() for c, m in f if abs(c - R["bound"]) <= half]
    ypad = .04 * (max(ys) - min(ys))
    ins.set_ylim(min(ys) - ypad, max(ys) + ypad)
    # Drawn the way the analysis repo cuts an inset (e1/fig4_overview.py): the window itself
    # is outlined on the parent, the inset carries a full frame, and two connectors run from
    # the rectangle's corners to the same corners of the inset -- the two whose segments
    # cross neither box, so the pair reads as one outline around both. A titled window with
    # no rectangle, which this panel used before, says in words where the inset came from
    # and leaves the reader to find it on the axes.
    ins.set_yticks([round(max(ys)), round(min(ys))])
    ins.set_xticks([R["bound"]])
    ins.set_xticklabels([f"{R['bound']:.2f}"])
    # The inset's own y labels go on its right. On the left they sit a few points from the
    # parent's y tick column and read as belonging to it -- two MFE scales, one above the
    # other, in one strip of page.
    ins.yaxis.tick_right()
    # A step below the parent's ticks, not a third smaller: at 3.8 the digits set 13 px
    # against the parent's 20, which is 4.3 pt and under any journal's floor. 5.0 is the
    # size of this figure's panel notes and still reads as another box's scale.
    # No tick marks: the inset's own are 3 px dashes drawn inside the parent's data area,
    # and the right-hand one landed on the parent's CAI = 0.90 gridline, where it reads as
    # a kink in that line rather than as the inset's scale. The labels carry it.
    ins.tick_params(labelsize=5.0 * K, length=0, pad=1.6)
    ins.set_facecolor("#ffffff")
    # The inset sits inside the parent's data area, so its tick labels are printed over the
    # parent's grid in the parent's own tick colour and size -- and they land where they
    # contradict it: "-602" fell on the main -350 gridline and the inset's "0.80" sat on the
    # main 0.85 gridline, in type a reader has no reason to tell apart. A white ground and a
    # smaller size tie each label to the box it belongs to.
    for t in ins.get_xticklabels() + ins.get_yticklabels():
        t.set_bbox(dict(fc="#ffffff", ec="none", pad=.6))
    # Above the parent's own reference lines: the dashed CAI floor was painted across the
    # inset for its full height, where it means nothing -- the inset draws its own.
    ins.set_zorder(10)
    ins.patch.set_alpha(1.0)
    for sp in ins.spines.values():
        sp.set_linewidth(.4)
        sp.set_color("#333333")
    zx0, zx1 = R["bound"] - half, R["bound"] + half
    zy0, zy1 = min(ys) - ypad, max(ys) + ypad
    # Above every front: the deepest rung's markers are drawn at zorder 13 and covered the
    # rectangle where the two meet.
    ax.add_patch(Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0, fill=False, lw=.5,
                           ec="#333333", zorder=14))
    fig.canvas.draw()
    rb = ax.transData.transform([(zx0, zy0), (zx1, zy1)])
    ib = ins.get_window_extent(fig.canvas.get_renderer())
    rbx, rby = sorted(rb[:, 0]), sorted(rb[:, 1])
    clear = []
    for cx in (0, 1):
        for cy in (0, 1):
            seg = np.linspace((rbx[cx], rby[cy]),
                              ((ib.x0, ib.x1)[cx], (ib.y0, ib.y1)[cy]), 50)[1:-1]
            inside = [((seg[:, 0] > bx0 + 1) & (seg[:, 0] < bx1 - 1)
                       & (seg[:, 1] > by0 + 1) & (seg[:, 1] < by1 - 1)).any()
                      for bx0, bx1, by0, by1 in ((rbx[0], rbx[1], rby[0], rby[1]),
                                                 (ib.x0, ib.x1, ib.y0, ib.y1))]
            if not any(inside):
                clear.append((cx, cy))
    for cx, cy in (clear[:2] if len(clear) >= 2 else [(0, 1), (1, 0)]):
        ax.add_artist(ConnectionPatch(((zx0, zx1)[cx], (zy0, zy1)[cy]), (cx, cy),
                                      "data", "axes fraction", axesA=ax, axesB=ins,
                                      lw=.4, color="#777777", zorder=2))
    tag(ax, "(c)", kind=("mid" if LAYOUT == 3 else "top"))


def _align(frames):
    """Each fold rotated onto the unconstrained one's frame, by index correspondence.

    naview picks a drawing orientation per structure, and that choice is a property of its
    own algorithm rather than of the molecule: it laid the deepest rung out portrait, 132 by
    221 units, against 226 by 119 for the other three. One shared angle therefore left that
    panel standing on end -- it read as a different kind of object -- and the union box the
    four share was 226 by 221, nearly square, while any one molecule needs about 2:1. Half the
    area inside every fold panel was paying for that mismatch.

    The four designs encode the same protein, so the drawings have the same number of points
    in the same order and a rotation between two of them is well defined. Each is rotated by
    the one that minimises the squared distance to the unconstrained drawing (Kabsch, on the
    point correspondence). It is a rigid rotation of a drawing, so it changes no geometry that
    carries meaning -- distances, pairs and the layout itself are untouched -- and it brings
    the union of the raw layouts from 226 by 221 to 234 by 164. Reflection is allowed for in
    the solve and is never chosen.

    This costs area rather than saving it, and the saving is what one would assume. Measured
    against a centre-only control through this panel's real aspect, and after the fill angle
    below has had its say on both: aligning gives a shared box of 229.9 by 161.1 and fills
    .912/.911/.913/.647 of the four panels, against 221.8 by 155.4 and .938/.937/.940/.800
    without it. So the drawings are 3.6% SMALLER aligned. What the alignment buys is that the
    four are in register -- the same helix points the same way in all four panels, which is
    the comparison the shared scale exists to support -- and the unaligned arrangement wins
    its area by rotating the deepest rung across the others, which is the thing that made it
    read as a different kind of object. The numbers above are of the raw layouts; the box the
    figure actually draws is the 229.9 by 161.1 one.
    """
    cen = [a - a.mean(axis=0) for a in frames]
    ref = cen[0]
    out = []
    for a in cen:
        u, _s, vt = np.linalg.svd(a.T @ ref)
        d = np.sign(np.linalg.det(vt.T @ u.T))
        out.append(a @ (vt.T @ np.diag([1.0, d]) @ u.T).T)
    return out


def _align_axes(frames):
    """Each fold rotated so that its long axis lies along the first fold's.

    The residue-by-residue alignment of _align fits two folds that are laid out alike. The two
    designs of (d) refold, so their layouts correspond residue to residue in name only, and the
    fit leaves one drawing turned against the other. A fold's long axis is the first principal
    axis of its coordinates. It fixes the orientation up to a half turn, and the half turn is
    the one that fits the first fold better residue by residue.
    """
    cen = [a - a.mean(axis=0) for a in frames]

    def angle(a):
        w, v = np.linalg.eigh(a.T @ a)
        return math.atan2(v[1, -1], v[0, -1])

    def rot(a, t):
        c, s_ = math.cos(t), math.sin(t)
        return a @ np.array([[c, s_], [-s_, c]])

    out = [cen[0]]
    for a in cen[1:]:
        t = angle(cen[0]) - angle(a)
        cand = [rot(a, t), rot(a, t + math.pi)]
        out.append(min(cand, key=lambda b: float(((b - cen[0]) ** 2).sum())))
    return out


def cap_parts(k):
    """The fields of fold k's caption, in order.

    Shared with the layout, which measures them before drawing: a caption wider than its
    panel runs into its neighbour's, and the overlap gate does not see that until it is 12%
    further out than the panel edge.
    """
    f0 = {tuple(q) for q in P[0]["pairs"]}
    here = {tuple(q) for q in P[k]["pairs"]}
    gained, lost = len(here - f0), len(f0 - here)
    # Green can only show a pair this design HAS. A pair it lost leaves no chord to colour,
    # so every caption carries both directions: at + I the design loses two and gains none,
    # and a single "differ" figure there described ink that is not on the panel.
    if LAYOUT == 4 and BVAR == "mfe":
        # Two designs: two lines in the top right corner of each drawing's box, naming the level
        # of M and giving the MFE.
        return [MSET[k], f"MFE {P[k]['mfe']:,.1f} kcal/mol".replace("-", "\u2212")]
    if LAYOUT == 4:
        # Two lines on every panel, the same two: the measurement on one and what the fold
        # did on the other. Four captions read as one object, so a panel whose fields wrap
        # differently from its neighbours reads as a different kind of caption rather than
        # as a longer one. Both zeroes are written -- "+0 \u22122 bp", not "\u22122 bp" --
        # because a missing term is read as a missing measurement.
        #
        # No unit here, which is the cost of the fixed two-line form: "MFE: \u2212616.0
        # kcal/mol, CAI: 0.801" is 1.46 in of type on a 1.25 in panel at any size this row
        # can carry. (c) names the same quantity on its own y axis directly above the row.
        line1 = (f"MFE: {P[k]['mfe']:,.1f}".replace("-", "\u2212")
                 + f", CAI: {P[k]['cai']:.3f}")
        # (d) is the design the other three are compared with, so its own difference is
        # zero -- written out rather than described. "the unconstrained design" was tried
        # in this slot and is prose where the others are a measurement, 334 px against
        # 124, wider than its own first line, so (d)'s caption block tapered the other way
        # from its neighbours' and the four stopped reading as one object.
        return [line1, f"+{gained} \u2212{lost} bp"]
    # The other layouts set three fields on one or two lines, so a zero term is dropped
    # rather than written: there the caption is as short as it can be made, and "same bp"
    # says what two zeroes would.
    both = " ".join(x for x in (f"+{gained}" if gained else "",
                                f"\u2212{lost}" if lost else "") if x)
    # The unit is written once, on the first panel. Four panels in a row are read as one
    # object and the other three repeat a constant; at 1.4 in wide that repetition is what
    # pushes every caption onto a second line, and a second line is paid for out of the
    # structures' own height. "bp" rather than "pairs" for the same reason.
    parts = [f"{P[k]['mfe']:,.1f}".replace("-", "\u2212") + (" kcal/mol" if k == 0 else ""),
             f"CAI {P[k]['cai']:.3f}"]
    if k:
        parts.append(f"{both} bp" if both else "same bp")
    return parts


# How much narrower than the drawing's own box a fold panel is set. The four share one
# box, so a panel's edge is where the WIDEST of the four molecules reaches, and two
# neighbours that both reach it leave almost no white between them. The layout's caption
# pre-measure has to use the same number, or it measures a caption against a panel that
# will not exist.
FOLD_SHRINK = .96
# How far the right design of (d) is drawn up and to the right in its box, at most, and the margin
# kept between its drawing and the box's top edge, in inches.
FOLD_UP, FOLD_RIGHT, FOLD_MARGIN = .0, .0, .02
FOLD_SHIFT = .14    # how far the right drawing of (d) sits left of its place in its box, in inches
CAPPAD = .02    # how far a drawing's caption stands from its box's top edge, in inches
CAPGAP = .04    # how far a caption stands above the ink beneath it, in inches
CAPBASE = [None]  # the line, in inches above the box's foot, that the first caption stands on


def fold_box(ar, slack=None):
    """The four aligned layouts, their shared scale, and the box that scale needs.

    Separated from the drawing because the LAYOUT needs it too: how much of a panel's width
    the molecules actually use decides how wide the panels should be, and that has to be
    known before the row is placed, not after it is drawn.

    Returns the rotated coordinates, the drawing's own half-extent, and the half-extent of
    the box at the panel's aspect ratio. One scale and one box for all four: per-panel limits
    let matplotlib resize each axes to its own molecule, so the deepest rung -- a wider,
    sparser fold -- was blown up relative to the others.
    """
    XY = (_align_axes if BVAR == "mfe" else _align)(
        [np.asarray(q["xy"], dtype=float) for q in P])
    # The fill angle is chosen on all four at once, not on the unconstrained one alone: the
    # row shares a box, so the angle that matters is the one that fits the union.
    allp = np.vstack(XY)
    ang = fig2.best_angle(list(allp[:, 0]), list(allp[:, 1]), aspect=ar)
    XY = [np.column_stack(fig2.rotate(list(a[:, 0]), list(a[:, 1]), ang)) for a in XY]
    # A per-panel nudge on top of the shared angle, applied BEFORE the union box is measured
    # so the shared scale accounts for it rather than letting a turned drawing overrun its
    # panel. It costs some of the register the alignment buys, which is the trade it exists
    # to make.
    for k, deg in ({} if BVAR == "mfe" else TWEAK).items():
        if k in SEL:
            k = SEL.index(k)
            XY[k] = np.column_stack(fig2.rotate(list(XY[k][:, 0]), list(XY[k][:, 1]), deg))
    lo = np.min([a.min(axis=0) for a in XY], axis=0)
    hi = np.max([a.max(axis=0) for a in XY], axis=0)
    XY = [a - (lo + hi) / 2 for a in XY]
    # The two-design variant takes its slack from the caller, which sets it so that each drawing
    # stands clear of the caption in its box's top right corner.
    sx, sy = (hi - lo) / 2 * (slack or (1.03 if BVAR == "mfe" else 1.01))
    HX, HY = (sx, sx * ar) if sy / sx < ar else (sy / ar, sy)
    return XY, sx, sy, HX, HY


def folds(axes):
    """The four designs, drawn as P0 draws them: a pale backbone, a dot per residue, chords
    for the pairs, and over the top, in that set's colour, the residues of every site the
    design still spells -- so a set's dots are present until the rung that forbids it, and
    absent from there on.

    Green marks a pair this design has and the unconstrained one does not. A pair it
    LOST leaves no chord here to colour, so the caption carries both directions and
    the key promises only what the ink can show. Drawing the lost pairs as chords was
    tried: on the deepest rung there are 315 of them and they cover the molecule.
    """
    bb = axes[0].get_position()
    ar = (bb.height * fig.get_size_inches()[1]) / (bb.width * fig.get_size_inches()[0])
    XY, sx, sy, HX, HY = fold_box(ar)
    if BVAR == "mfe":
        # Each box keeps its width, and its top right corner holds the caption. The drawings are
        # scaled down from a fit to the box only as far as it takes to keep every point of the
        # ink clear of the caption in both boxes.
        win = axes[0].get_position().width * fig.get_size_inches()[0]
        hin = axes[0].get_position().height * fig.get_size_inches()[1]
        rend_ = fig.canvas.get_renderer()
        cap = []
        for k_ in range(len(axes)):
            t_ = fig.text(0, 0, "\n".join(cap_parts(k_)), fontsize=CAPFS, linespacing=1.5)
            bb_ = t_.get_window_extent(rend_)
            cap.append((bb_.width / fig.dpi, bb_.height / fig.dpi))
            t_.remove()
        for slack in np.arange(1.01, 1.40, .005):
            XY, sx, sy, HX, HY = fold_box(ar, slack)
            clear = True
            for k_, (cw, ch) in enumerate(cap):
                ix = (XY[k_][:, 0] + HX) / (2 * HX) * win
                iy = (XY[k_][:, 1] + HY) / (2 * HY) * hin
                # The caption is centred over the drawing's ink, whose centre is not the box's,
                # and stands just above the ink beneath it. It has to fit under the box's top.
                under = iy[np.abs(ix - (ix.min() + ix.max()) / 2) < cw / 2 + .03].max()
                clear = clear and under + CAPGAP + ch <= hin - CAPPAD
            if clear:
                break
        print(f"  structures at {slack:.3f} of the fit to their box, clear of their captions")
    elif HX > sx * 1.001:
        # Height binds, so every drawing sits in the middle of its panel with the same slack
        # either side -- white that reads as nothing and pushes the row's two ends inwards.
        # Narrow each panel to the drawing's own width and spread the four over the same
        # span, so the first is flush with the row's left edge, the last with (c)'s right
        # edge, and the gaps between them are equal and real.
        x0, x1 = axes[0].get_position().x0, axes[-1].get_position().x1
        # 0.96, not 1.0: the four share one box, so a panel's edge is where the WIDEST of
        # the four molecules reaches, and two neighbours that both reach it leave almost no
        # white between them -- (f) and (g) came within 5 px against 30 for the other pairs.
        pw = axes[0].get_position().width * sx / HX * FOLD_SHRINK
        gap = ((x1 - x0) - len(axes) * pw) / (len(axes) - 1)
        for k, ax in enumerate(axes):
            q = ax.get_position()
            ax.set_position([x0 + k * (pw + gap), q.y0, pw, q.height])
        HX = sx
    # Each fold marks the base pairs the other fold does not have: the right fold in green, and
    # in the two-design variant the left fold in red, so the folds show what is unique to each
    # design and what they share.
    # A drawing is sized to fit its box, and a stroke at its tip is wider than the point it
    # ends on: clipped to the box, the tip is cut flat. The two-design variant draws unclipped.
    CLIP = BVAR != "mfe"
    for k, (ax, p) in enumerate(zip(axes, P, strict=True)):
        xy = XY[k]
        other = P[1] if (k == 0 and BVAR == "mfe") else P[0]
        f0 = {tuple(q) for q in other["pairs"]}
        mark = UNIQ if BVAR == "mfe" else GAINED
        for i, j in p["pairs"]:
            gained = (i, j) not in f0
            # A gained pair is drawn heavier than a kept one. At 0.45 pt (f)'s four gained
            # pairs were one 8 x 7 px blob -- the panel whose caption says "+4" showed a
            # speck -- and the key's handle is a solid rule, which is what the ink has to
            # look like.
            ax.plot(xy[[i, j], 0], xy[[i, j], 1], color=mark if gained else PAIR,
                    lw=.70 if gained else .35, zorder=4 if gained else 1, clip_on=CLIP)
        for bx, by in fig2.backbone_segments(list(xy[:, 0]), list(xy[:, 1]))[0]:
            ax.plot(bx, by, color=BACK, lw=.9, solid_capstyle="round", zorder=2, clip_on=CLIP)
        ax.scatter(xy[:, 0], xy[:, 1], s=.9, c=BACK, linewidths=0, zorder=3, clip_on=CLIP)
        for a_, b_, roman in p["sites"]:
            ax.scatter(xy[a_:b_, 0], xy[a_:b_, 1], s=5.0, c=SETC[roman], linewidths=0, zorder=6,
                       clip_on=CLIP)
        ax.set_xlim(-HX, HX)
        ax.set_ylim(-HY, HY)
        if BVAR == "mfe" and k == 1:
            # The right design is drawn up and to the left in its box, clear of the caption it
            # would otherwise run behind. How far is what its drawing has room for: it is set
            # on the box's top edge less a margin, and never clipped there.
            hin = ax.get_position().height * fig.get_size_inches()[1]
            win = ax.get_position().width * fig.get_size_inches()[0]
            room = (HY - xy[:, 1].max()) / (2 * HY) * hin - FOLD_MARGIN
            up_in = min(FOLD_UP, max(room, 0.0))
            up = up_in / hin * 2 * HY
            ax.set_ylim(-HY - up, HY - up)
            # Left by what its drawing has room for on that side, and no more: a drawing moved
            # past the box's edge is cut off there.
            room_r = (HX - xy[:, 0].max()) / (2 * HX) * win - FOLD_MARGIN
            right_in = min(FOLD_RIGHT, max(room_r, 0.0))
            right = right_in / win * 2 * HX
            ax.set_xlim(-HX - right, HX - right)
            print(f"  right design: {up_in:.3f} in up, {right_in:.3f} in right")
        if BVAR == "mfe":
            # The right drawing sits FOLD_SHIFT to the left of its place in the box, which takes
            # it away from the protein's name in the corner of the page. The Venn diagram is
            # set between the two drawings (see the layout).
            if k == 1:
                d = FOLD_SHIFT / (ax.get_position().width * fig.get_size_inches()[0]) * 2 * HX
                xl_, xh_ = ax.get_xlim()
                ax.set_xlim(xl_ + d, xh_ + d)
        if BVAR == "mfe":
            # What each drawing leaves between its ink and the box, in inches. Zero or less is a
            # cut-off.
            xl, xh = ax.get_xlim()
            yl, yh = ax.get_ylim()
            hin_ = ax.get_position().height * fig.get_size_inches()[1]
            win_ = ax.get_position().width * fig.get_size_inches()[0]
            gaps = {"left": (xy[:, 0].min() - xl) / (xh - xl) * win_,
                    "right": (xh - xy[:, 0].max()) / (xh - xl) * win_,
                    "bottom": (xy[:, 1].min() - yl) / (yh - yl) * hin_,
                    "top": (yh - xy[:, 1].max()) / (yh - yl) * hin_}
            print(f"  fold {k}: ink to box, in: " + ", ".join(f"{a} {v:.3f}" for a, v in gaps.items()))
            # The right drawing is meant to hang over its box's left edge, into the column the
            # layout has moved with it; every other side has to stay inside.
            checked = {a: v for a, v in gaps.items() if not (k == 1 and a == "left")}
            assert min(checked.values()) > -.03, f"fold {k} runs far past its box: {gaps}"
        ax.axis("off")
        # One line where the panel can hold it, stacked where it cannot. The rule used to be
        # a proxy -- a threshold on the panel's width in inches -- which is right only for
        # one type size: raising the caption two points made the one-line form overrun a
        # panel the proxy still called wide enough, and the caption ran into its neighbour's.
        # The test is now the thing itself, the rendered width of this caption against the
        # width of this panel.
        parts = cap_parts(k)
        if LAYOUT == 4:
            # The form is fixed, so there is nothing to choose: the two lines are the two
            # lines. The cascade below exists for the layouts whose captions are three
            # fields that may or may not fit on one line.
            # 0.046 in under the panel, in INCHES: an axes-fraction offset is a share of
            # the panel's height, so it grows with the page while the caption band below
            # it is a constant, and the second line walks off the canvas. The band is
            # CAP and the gate catches the day this stops fitting.
            # In the figure's ink, not the tick grey: these carry the per-panel numbers,
            # and in grey at 6.3 pt they were the lightest and joint-smallest primary text
            # on the page.
            if BVAR == "mfe":
                # The caption stands in the box's top right corner, clear of the drawing.
                # Centred over the drawing's ink, which is not centred in its box.
                xa, xb = ax.get_xlim()
                ya, yb = ax.get_ylim()
                wi_, hi_ = _inw(ax), _inh(ax)
                cx = (xy[:, 0].min() + xy[:, 0].max()) / 2
                cwk = cap[k][0]
                ixs = (xy[:, 0] - xa) / (xb - xa) * wi_
                iys = (xy[:, 1] - ya) / (yb - ya) * hi_
                ixc = (cx - xa) / (xb - xa) * wi_
                if k == 0:
                    # Just above the ink beneath the caption, and no higher than the box allows.
                    base = min(iys[np.abs(ixs - ixc) < cwk / 2 + .03].max() + CAPGAP,
                               hi_ - CAPPAD - cap[k][1])
                    CAPBASE[0] = base
                else:
                    # On the first caption's line, moved to the right until the drawing's ink
                    # leaves it clear.
                    base = CAPBASE[0]
                    chk = cap[k][1]
                    for step in np.arange(0.0, wi_, .01):
                        inbox = ((np.abs(ixs - (ixc + step)) < cwk / 2 + .03)
                                 & (iys > base - .03) & (iys < base + chk + .03))
                        if not inbox.any() or ixc + step + cwk / 2 > wi_ - CAPPAD:
                            break
                    ixc += step
                    print(f"  caption of fold {k}: {step:.2f} in right of its ink's centre")
                ax.annotate("\n".join(parts), xy=(ixc / wi_, base / hi_),
                            xycoords="axes fraction", fontsize=CAPFS, color=INK, ha="center",
                            va="bottom", linespacing=1.5, zorder=10,
                            path_effects=[patheffects.withStroke(linewidth=2.2,
                                                                 foreground="white")])
            else:
                ax.annotate("\n".join(parts), xy=(.5, -.046 / _inh(ax)),
                            xycoords="axes fraction", fontsize=CAPFS, color=INK,
                            ha="center", va="top", linespacing=1.5)
        else:
            # Three ways to set the same three fields, tried widest first and kept at the
            # first one that fits the panel. Every extra line is paid for out of the
            # structure's height, so the fewest lines that fit is the right answer rather
            # than always stacking.
            # Below the box, not against it: the drawings are scaled to the box, so a
            # strand that reaches the bottom edge ended 10 px from the caption -- closer
            # than the caption's own words are to each other and on the digits' own line.
            note = ax.annotate("", xy=(.5, -.046 / _inh(ax)), xycoords="axes fraction",
                               fontsize=CAPFS, color=MUTED, ha="center", va="top",
                               linespacing=1.5)
            rend = fig.canvas.get_renderer()
            for lines in ("  ".join(parts),
                          "  ".join(parts[:2]) + ("\n" + parts[2] if k else ""),
                          "\n".join(parts)):
                note.set_text(lines)
                if note.get_window_extent(rend).width <= ax.get_window_extent().width:
                    break
        # Outside the box. Inside it the letter sat over the molecule, which is what the
        # drawing is for; the fold band is placed to leave this strip clear.
        if BVAR != "mfe" or k == 0:
            tag(ax, f"({chr(100 + k)})", dx=.0, dy=1.035, kind="fold")


def middle(ax):
    """What the two folds of (d) share, between them, as a proportional Venn diagram.

    The area of each circle is the design's number of base pairs and the area of the overlap is
    the number the two designs share, so the figure shows how little the two folds have in
    common. Dark is what only one design has and light grey is what the two share, as in the
    folds.
    """
    from matplotlib.patches import Circle
    W = ax.get_position().width * fig.get_size_inches()[0]
    Hh = ax.get_position().height * fig.get_size_inches()[1]
    ax.set_xlim(0, W)
    ax.set_ylim(0, Hh)
    ax.axis("off")
    f0 = {tuple(q) for q in P[0]["pairs"]}
    f3 = {tuple(q) for q in P[1]["pairs"]}
    sets = (("base pairs", len(f0 - f3), len(f0 & f3), len(f3 - f0)),)
    # The two circles and their overlap span 3.4 radii, and the column is only as wide as the
    # gap the two folds leave.
    R1 = min(.18 if BVAR == "mfe" else .26, (W - .04) / 3.42)
    unit = math.pi * R1 ** 2 / (sets[0][1] + sets[0][2])
    TITLE, GAP, BETWEEN = .11, .02, .06
    heights = [TITLE + GAP + 2 * math.sqrt(unit * (nl + nb) / math.pi)
               for _, nl, nb, _r in sets]
    y = Hh - (Hh - sum(heights) - BETWEEN) / 2
    for (title, nl, nb, nr), hgt in zip(sets, heights, strict=True):
        ra = math.sqrt(unit * (nl + nb) / math.pi)
        rb = math.sqrt(unit * (nb + nr) / math.pi)

        def lens(d, ra=ra, rb=rb):
            if d >= ra + rb:
                return 0.0
            if d <= abs(ra - rb):
                return math.pi * min(ra, rb) ** 2
            a1 = math.acos((d * d + ra * ra - rb * rb) / (2 * d * ra))
            a2 = math.acos((d * d + rb * rb - ra * ra) / (2 * d * rb))
            return ra * ra * (a1 - math.sin(2 * a1) / 2) + rb * rb * (a2 - math.sin(2 * a2) / 2)

        lo, hi = abs(ra - rb), ra + rb
        for _ in range(60):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if lens(mid) > nb * unit else (lo, mid)
        d = (lo + hi) / 2
        cy = y - TITLE - GAP - ra
        c1x = W / 2 - (d + rb - ra) / 2
        c2x = c1x + d
        a = Circle((c1x, cy), ra, fc="#8f8e89", ec=UNIQ, lw=.5, zorder=1)
        b = Circle((c2x, cy), rb, fc="#8f8e89", ec=UNIQ, lw=.5, zorder=1)
        ax.add_patch(a)
        ax.add_patch(b)
        both = Circle((c2x, cy), rb, fc=SHAREDC, ec="none", zorder=2)
        ax.add_patch(both)
        both.set_clip_path(a)
        ax.add_patch(Circle((c1x, cy), ra, fc="none", ec=UNIQ, lw=.5, zorder=3))
        ax.add_patch(Circle((c2x, cy), rb, fc="none", ec=UNIQ, lw=.5, zorder=3))
        x_l, x_r = c1x - ra, c2x + rb
        l_in, r_in = c2x - rb, c1x + ra
        for x, n in (((x_l + l_in) / 2, nl), ((l_in + r_in) / 2, nb), ((r_in + x_r) / 2, nr)):
            ax.text(x, cy, f"{n}", ha="center", va="center", fontsize=FSS_, color=INK,
                    zorder=4)
        ax.text(W / 2, y, title, ha="center", va="top", fontsize=FSS_, color=INK)
        y -= hgt + BETWEEN


FSS_ = 5.4 * K
# How the two designs are named in (d): by the set of motifs that each forbids.
# The four levels of M in the text's notation: the empty set and the union of the sets so far.
_MI, _MII, _MIII = (r"\mathcal{M}_{\mathrm{%s}}" % r for r in ("I", "II", "III"))
MLEV = (r"\emptyset", _MI, _MI + r"\cup " + _MII, _MI + r"\cup " + _MII + r"\cup " + _MIII)
MSET = (r"$\mathcal{M}=" + MLEV[0] + "$", r"$\mathcal{M}=" + MLEV[3] + "$")


def codon_strip(ax, k):
    """One cell per codon of design k, over the whole of (d).

    A codon the two designs share is light grey and a codon that the right design changes is
    dark. In the left design's strip a codon that overlaps a site of a forbidden motif takes
    that set's colour instead.
    """
    n = P[1]["n_codons"]
    changed = set(P[1]["codons_differ"])
    site = {}
    if k == 0:
        for a, b, roman in P[0]["sites"]:
            for c in range(a // 3, (b - 1) // 3 + 1):
                site.setdefault(c, roman)
    kept, gone = SHAREDC, "#4a4742"
    row = np.array([[matplotlib.colors.to_rgb(SETC[site[c]] if c in site
                                              else gone if c in changed else kept)
                     for c in range(n)]])
    ax.imshow(row, extent=(0, n, 0, 1), aspect="auto", interpolation="nearest")
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color(OUTLINE)
        sp.set_linewidth(.4)


if LAYOUT == 1:
    LEGCOL = None
    fig = plt.figure(figsize=(12.4, 2.78))
    # Explicit boxes rather than one gridspec: a uniform gutter has to be wide enough for the
    # widest axis labels, and (c)->(d) has none on either side of it, so it was paying for (b)'s.
    TOP, BOT = .935, .205
    H_ = TOP - BOT
    FW, FGAP = .0805, .004
    panel_a(fig.add_axes([.098, BOT, .185, H_]), lab=3.9)
    # The attribution strip rides on top of (b) here too, taking its height from the panel
    # rather than from the page, so the wide-row layout keeps the same reading.
    BS = H_ * .14
    panel_b(fig.add_axes([.330, BOT, .128, H_ - BS - .030]),
            fig.add_axes([.330, BOT + H_ - BS, .128, BS]))
    panel_c(fig.add_axes([.512, BOT, .140, H_]))
    folds([fig.add_axes([.664 + k * (FW + FGAP), BOT, FW, H_]) for k in range(4)])
    LEG, LEGFS = (.5, .002), 6.2
elif LAYOUT == 2:
    # (a) runs the full height, which is what 24 stacked rows want; (b) and (c) share the top
    # right and the folds the bottom right, lying flat so they fill a wide, short box.
    # Two blocks. (b) is left-aligned with (d) and (c) with (f), so the fold pitch is set by
    # the distance between them; the letters sit on one baseline in figure coordinates.
    # The side-by-side variant needs two cost panels where there was one, so it buys the width
    # rather than halving them; the stacked variant needs the height for the same reason, and
    # (c) and the folds, which are placed in figure fractions, grow with it.
    FIGW, FIGH = {"strip": (9.4, 3.0), "side": (11.6, 3.0), "stack": (9.4, 3.9)}[BVAR]
    fig = plt.figure(figsize=(FIGW, FIGH))
    # Horizontal fractions are quoted for the 9.4 in page, so a wider page keeps the same
    # inches by scaling them. Vertical ones are not scaled: the taller page is taller
    # everywhere on purpose.
    SC = 9.4 / FIGW
    # The key is one line across the top, so the band it used to take at the bottom -- two
    # ragged rows under the folds -- goes to the panels instead.
    # The key is set two points larger than the panel type and centred, so it needs a taller
    # band than the 5.4 pt flush-left row did: 40 px against 30. The panels give up 16 px of
    # height for it rather than the key overrunning the canvas, which savefig would crop
    # silently now that the tight box is off.
    TOP = .870
    # (b) and (c) set the fold pitch, but through their VISUAL left edges rather than their
    # data boxes: a fold panel has no tick labels, so pitching the row off (b)'s data box
    # left the width of (b)'s label column unused under (b) and (c). The row now starts where
    # (b)'s ink starts, which is what lets the structures grow.
    L0, STEP = .278 * SC, .173 * SC
    # The stacked variant splits (b) in two, and two boxes with a two-line x label between
    # them do not fit in the band one box had: the top row starts lower and the fold row and
    # its letters move down with it, rather than the panels shrinking to nothing.
    # Side by side, (b)'s two-line x labels hang lower than the single-line one the strip
    # variant has, and the fold letters are directly under them, so that row starts higher.
    BROW = {"stack": .470, "side": .558}.get(BVAR, .520)
    TAGY, TAGY_F = .884, (.333 if BVAR == "stack" else .381)
    TAGY_K = {"top": TAGY, "mid": TAGY, "fold": TAGY_F}
    RIGHT = .993          # the right edge (c) and the last fold share
    FGAP = .003 * SC
    FTOP = .330 if BVAR == "stack" else .378   # just under the fold letters' baseline
    # The band under the fold row, for the caption.
    # One line of caption at CAPFS is 22.5 px; this leaves 28.8, and the save is gated by an
    # off-canvas assertion, so a caption that ever needs two lines fails loudly rather than
    # being cropped. Height is the binding dimension for the drawings -- the union of the four
    # layouts is taller, relative to the box, than it is wide -- so every pixel here is a pixel
    # off the structures.
    FBOT = .062
    # The strip under the fold row used to be a 6.8 x 0.5 in band of nothing, 11.6% of the
    # page, while (b) and (c) were a third the height of (a). The fold row sits at the
    # bottom now and the strip is spent on (b) and (c), which grow by half again.
    ax_a = fig.add_axes([.057 * SC, .120, .174 * SC, TOP - .120])
    panel_a(ax_a, lab=5.0)
    # (b) is two axes now: a thin attribution strip over the boxes. They share a left edge,
    # so the fold row's origin and the panel letter are unaffected.
    if BVAR == "strip":
        BSTRIP = .062
        ax_bs = fig.add_axes([L0, TOP - BSTRIP, .292, BSTRIP])
        ax_b = fig.add_axes([L0, BROW, .292, TOP - BROW - BSTRIP - .030])
        panel_b(ax_b, ax_bs)
    elif BVAR == "side":
        # Two panels inside the column (b) had, ending clear of where the fold pitch will put
        # (c): the fold row spans L0 to the right margin in four steps, so (c) starts one
        # third of the way across and nothing of (b)'s may reach it. The gap between the two
        # carries the right panel's tick labels and its y axis label.
        BW, BGAP = .142, .056
        ax_b = fig.add_axes([L0, BROW, BW, TOP - BROW])
        ax_b2 = fig.add_axes([L0 + BW + BGAP, BROW, BW, TOP - BROW])
        panel_b2(ax_b, ax_b2)
    else:
        # Stacked, so both keep (b)'s full width and split its height. The band between them
        # carries the upper panel's two-line x label.
        # The upper panel stops just short of TOP: its y label is centred on a short box and
        # so reaches past the box's own top, where the one-line key sits.
        BH, BMID, BCAP = .148, .108, .022
        ax_b2 = fig.add_axes([L0, BROW, .292, BH])
        ax_b = fig.add_axes([L0, BROW + BH + BMID, .292,
                             TOP - BCAP - BROW - BH - BMID])
        panel_b2(ax_b, ax_b2)
    ax_c = fig.add_axes([L0 + 2 * STEP, BROW, .341, TOP - BROW])
    panel_c(ax_c)
    # The top row has to be rendered once before its ink can be measured.
    fig.canvas.draw()
    # The fold row starts where (b)'s ink starts, but never left of (a)'s: how far (b)'s ink
    # reaches depends on how wide its tick labels are, which depends on the population drawn,
    # so the subset with the smallest costs has the narrowest labels and pushed the row left
    # into (a)'s "100". A layout that reads off the data needs a guard against the data.
    # (b)'s ink is what sets this edge; the clause on (a) is a floor that does not bind at
    # this page size (it is 0.5 in to the left). Raising the floor to widen the gap between
    # (a)'s x label and (d)'s caption was tried and is the wrong lever: it would buy 0.1 in
    # of white by giving up the alignment between the structures and the panel above them.
    L0F = max(visual_x0(ax_b),
              ax_a.get_tightbbox(fig.canvas.get_renderer())
                  .transformed(fig.transFigure.inverted()).x1 + .006)
    # Two ceilings on a fold panel's width, and the pitch decides which one binds: the pitch
    # itself, since panels must not touch, and the right margin, since the last panel must not
    # overrun the edge it shares with (c). A wider pitch relaxes the first and tightens the
    # second, so the panels are largest where the two meet -- solve for that pitch rather than
    # inheriting one and discovering which ceiling bit. Taking (c)'s ink where it happened to
    # fall gave a pitch below this and left .058 of the page, a third of a panel, unused at the
    # right margin.
    STEPF = (RIGHT - L0F + FGAP) / 4
    FWF = STEPF - FGAP
    # (c) is then placed so that its own ink starts where (f) does, which is what keeps the
    # letters aligned. Its overhang -- tick labels plus the y axis label -- does not depend on
    # where the box sits, so one measurement fixes the move.
    ax_c.set_position([L0F + 2 * STEPF + (ax_c.get_position().x0 - visual_x0(ax_c)),
                       BROW, 0, TOP - BROW])
    ax_c.set_position([ax_c.get_position().x0, BROW, RIGHT - ax_c.get_position().x0,
                       TOP - BROW])
    folds([fig.add_axes([L0F + k * STEPF, FBOT, FWF, FTOP - FBOT]) for k in range(4)])
    place_tags()
    # Clear of the canvas top. Anchored at .962 the one-line key overflowed the figure by
    # 11.7 px, and savefig's tight box then grew the output to 3.055 in.
    LEG, LEGFS = (.5, .928), 6.2
    LEGCOL = None
elif LAYOUT == 3:
    # Layout 3: the paper's text width, and the paper's type.
    #
    # 6.5 in is what the analysis repo calls the text width and what the built Figure 4 uses;
    # 7.0 is this figure's own, slightly wider, and the difference is the one the caller asked
    # for. Seven panels do not fit across it, so the plate is three bands instead of two: the
    # catalog runs down the left over the top two, the cost boxes and the fronts take the right
    # of those two, and the structures keep a full-width row of their own at the bottom. The
    # two b variants differ only in how the right of the top band is split.
    FIGH = 6.6 if BVAR == "side" else 6.8
    fig = plt.figure(figsize=(7.0, FIGH))
    # (a)'s row labels are motif names, not numbers, and they hang left of its axes: at
    # LEFT = .055 the longest ran off the page.
    LEFT, RIGHT, R0 = .085, .985, .350
    AR = .250             # (a)'s right edge, clear of (c)'s y label
    FGAP = .004
    # Each band's bottom clears the band below by more than the x label stack that hangs off
    # it -- two lines plus the tick row, which measures .060 of the page here. The letters go
    # 12 thousandths above each band, which is where the tick-free top of a panel is.
    BTOP, BBOT, CTOP, CBOT, FTOP, FBOT, LEGY = (
        (.872, .640, .565, .375, .288, .045, .912) if BVAR == "side" else
        (.890, .614, .536, .366, .280, .045, .928))
    TAGY_K = {"top": BTOP + .012, "mid": CTOP + .012, "fold": FTOP + .012}
    ax_a = fig.add_axes([LEFT, CBOT, AR - LEFT, BTOP - CBOT])
    panel_a(ax_a, lab=5.0 * K)
    if BVAR == "strip":
        raise SystemExit("layout 3 draws the two-direction (b); pass side or stack")
    if BVAR == "side":
        # The gap between the two carries the right panel's tick labels and its y label.
        BGAP = .090
        BW = (RIGHT - R0 - BGAP) / 2
        ax_b = fig.add_axes([R0, BBOT, BW, BTOP - BBOT])
        ax_b2 = fig.add_axes([R0 + BW + BGAP, BBOT, BW, BTOP - BBOT])
    else:
        BH = .104
        ax_b2 = fig.add_axes([R0, BBOT, RIGHT - R0, BH])
        ax_b = fig.add_axes([R0, BTOP - BH, RIGHT - R0, BH])
    panel_b2(ax_b, ax_b2)
    ax_c = fig.add_axes([R0, CBOT, RIGHT - R0, CTOP - CBOT])
    panel_c(ax_c)
    # The folds get the whole width here rather than being pitched off (b)'s ink: nothing
    # above them shares their left edge, so there is no alignment to preserve and the panels
    # take every pixel of the row.
    STEPF = (RIGHT - LEFT + FGAP) / 4
    folds([fig.add_axes([LEFT + k * STEPF, FBOT, STEPF - FGAP, FTOP - FBOT])
           for k in range(4)])
    place_tags()
    # Eight entries at 7 pt do not fit on one 7 in line, so the key takes two rows.
    LEG, LEGFS, LEGCOL = (.5, LEGY), 6.2 * K, None

else:
    # Layout 4 -- D50's arrangement at the paper's text width: one row of panels over one row
    # of structures, with (a) running the FULL height of the page as it does in D50 rather
    # than sitting inside the top band. Nothing is rearranged; the width comes out of (b),
    # (c) and the structures, and (a)'s 24 motif rows keep their 2.2 in.
    #
    # The vertical layout is written in INCHES and converted once, because every band except
    # the panels themselves costs a fixed height whatever the page is: the key, a panel
    # letter, an x tick row plus its label, the caption under the structures. Solving in
    # fractions and rescaling sends the slack to whichever band happens to be the largest
    # fraction, which is how the structures stayed 0.55 in tall on a page 0.2 in taller.
    FIGH = 3.0 if BVAR == "side" else (float(sys.argv[5]) if len(sys.argv) > 5 else 3.2)
    FIGW = float(sys.argv[6]) if len(sys.argv) > 6 else 7.0
    # The two-design variant has no key across the page, so the band the key took above the
    # panel letters goes: the page is that much shorter and every panel keeps its place.
    TRIM = .11 if BVAR == "mfe" else 0.0
    FIGH -= TRIM
    fig = plt.figure(figsize=(FIGW, FIGH))
    # Three fixed bands, in inches: a panel letter, the stack under (b) and (c), and the
    # caption band under the structures.
    #
    # XLAB is the band under the panel row: a tick row, an x label line, and the air that
    # keeps the row and the structures reading as two rows. (b) and (c) share it, because
    # they share a bottom -- see BROW.
    #
    # CAP is the whole caption band, measured on the page: 0.046 in under the drawing,
    # 0.215 in of two-line type at 6.3 pt, and 0.014 in of page under the descenders.
    # LET is a fold letter's box as it MEASURES, 0.095 in, not the 0.072 that was written
    # here: the difference was silently taken out of XLAB, so the band under (b) and (c)
    # was 0.023 in shallower than its own name, which is invisible while XLAB has slack
    # and is a collision the moment it is sized against a measured label stack.
    LET, XLAB, CAP = .095, .297, .275
    # With (b) a single panel the cost row needs less height, and the page gives the
    # difference to the two structures.
    FOLDMFE = 1.50
    # Under (d)'s letter, its two codon strips, and the row of insets that zoom into the motifs,
    # which are 0.66 in deep.
    STRIPBAND = .66
    # The structures' band grows with the page: a taller page is asked for to make them
    # bigger, so three fifths of whatever is added lands here and the rest in the cost
    # panels. At 3.0 in this is the band the layout was solved on.
    FOLD = (FOLDMFE if BVAR == "mfe" else 0.900) + max(0.0, FIGH - 3.0) * .6

    def _f(y):
        """Inches from the foot of the page, as a figure fraction."""
        return y / FIGH

    # The key band and the letter row, which is no longer a band. The letters used to sit on
    # a line of their own between the key and the panels: 0.16 in of page for 0.10 in of
    # type, with a strip of nothing above and below it. A letter's BASELINE is now the
    # panel's own top edge, so it rises into the gap the key's handle padding already left
    # empty, and the panels start where the letter row used to.
    # The two-design variant lifts the key and the panel letters by 0.03 in: with (b) and (c)
    # short, the letters sat on their rotated y labels.
    LIFT = .03 if BVAR == "mfe" else 0.0
    LEGY = _f(FIGH - .170 + LIFT)
    # 0.060 in between a letter's box and its panel's top edge. A top-row letter sits in
    # its panel's y tick column, and the topmost tick label reaches up into the same strip:
    # with the letter on the panel edge they were 0.9 px apart, and at 0.045 "\u2212350"
    # came within 4.2 px of "(c)" with the gate asking 4.0. Which tick lands at the top is
    # a property of the instance, so the margin has to hold for the instance that puts one
    # there -- P49755 fails at 0.045 and clears here.
    TOP = _f(FIGH - .332 + TRIM)
    TAGY = _f(FIGH - .272 + LIFT + TRIM)
    FBOT, FTOP = _f(CAP), _f(CAP + FOLD)
    # The fold letters sit below their panels' tops: a fold panel's box is its ink and the
    # drawing is centred in it, so the strip above the molecule is empty and the letter can
    # be set into it. 0.065 rather than 0.034 buys the clearance under (c)'s x label, which
    # otherwise sits the same distance from (c)'s ticks as from (g)'s letter and reads as a
    # title line for (g).
    TAGF = _f(CAP + FOLD - .065)
    BROW = _f(CAP + FOLD - .065 + LET + XLAB)
    FGAP = .004
    TAGY_K = {"top": TAGY, "mid": TAGY, "fold": TAGF}
    if BVAR == "strip":
        raise SystemExit("layout 4 draws the two-direction (b); pass side or stack")
    # Horizontal, in inches for the same reason as the vertical: the margins, the label
    # columns beside (b) and (c) and the gap between the cost panels are fixed widths, so a
    # wider page must give its extra to the PANELS and not spread it over the gutters.
    #
    # (c) takes a width of its own rather than sitting on the third fold's column. Pinning it
    # there made the fold row and (c) share a left edge, and it also fixed (c)'s width at two
    # fold pitches -- 2.4 in against 0.83 in for each cost panel, which is more than a front
    # of 1,300 points needs and less than twelve box plots do. Giving up the alignment is what
    # buys (b) its width; the two then share the page's slack in the ratio they ended at.
    # Margins are what the ink needs and nothing more: (a)'s longest motif name measures
    # 0.49 in, so 0.50 of left margin holds it with the page edge against the type, and the
    # right margin is gone. B_LAB is narrow because (b) draws three short y ticks; C_LAB is
    # wide because it is also the air between (b) and (c).
    # A_GAP is air and nothing else: B_LAB is the width (b)'s own ticks and y label need, so
    # without it (a)'s rightmost bar end and (b)'s y label are a tick row apart -- 0.035 in
    # measured, against 0.065 between (b)'s label and its own ticks -- and (a) reads as
    # abutting the block beside it.
    # M_R is half of (c)'s last x tick label and nothing else: at 0.040 "1.00" ran 2 px off
    # the canvas, where no tight box would have grown the page to hold it.
    # Both outer margins are the ink's own clearance and nothing more, but they are not
    # zero: at M_R 0.062 the right-most tick label stopped 2 px from the canvas and (a)'s
    # longest row label 7 px from the left, so any rule or trim the page adds lands on
    # them. 0.04 in each side, which is a margin and still not a gutter.
    M_L, A_W, A_GAP, B_LAB, B_GAP, C_LAB, M_R = .560, 1.150, .160, .500, .430, .610, .080
    if BVAR == "mfe":
        # The key over (a) and the page legend are gone, and (a)'s right edge no longer needs air
        # beside (b): the white that stood between them goes to the panels.
        A_GAP = .010
        # The margin holds (a)'s motif names and nothing else, which were 0.49 in wide at most.
        M_L -= .03
    FIXED = M_L + A_W + A_GAP + B_LAB + B_GAP + C_LAB + M_R
    REM = FIGW - FIXED
    BW_IN, C_IN = REM * (.980 / 3.815), REM * (1.855 / 3.815)

    def _w(x):
        """Inches from the left edge, as a figure fraction."""
        return x / FIGW

    LEFT, RIGHT = _w(M_L), 1 - _w(M_R)
    B0, BW, BGAP = _w(M_L + A_W + A_GAP + B_LAB), _w(BW_IN), _w(B_GAP)
    C0 = _w(M_L + A_W + A_GAP + B_LAB + 2 * BW_IN + B_GAP + C_LAB)
    # (a) runs to the foot of the page: its x label stack ends on the fold captions' last
    # line, so the block and the structures bottom out together and the plate has one
    # baseline rather than two.
    #
    # This is a deliberate reversal. (a) used to stop above the caption row precisely so
    # that its tick labels and the captions would not read as one run of text, and that
    # objection still holds on its own terms: at the foot of the page "100" and "MFE:" are
    # 0.14 in apart where two neighbouring captions are 0.33 in apart, so the two blocks do
    # read as two rows of type across the plate. The alignment is worth more.
    #
    # The number is an INK alignment and not the sum its parts suggest. The two blocks are
    # set at different sizes -- (a)'s x label at 7.0 pt, a caption's second line at 6.3 --
    # so equal box bottoms do not give equal baselines: the box arithmetic (a caption's
    # band bottom 0.014 plus (a)'s own 0.338 in stack) lands a pixel out, and 0.355 is what
    # puts the two inks on one line, measured at 7.3 x 3.0. It holds at every page size,
    # since the stack is type-driven, and it moves if either type size does.
    A_BOT = _f(.355)
    ax_a = fig.add_axes([LEFT, A_BOT, _w(A_W), TOP - A_BOT])
    # The motif names are tick labels and set at the tick size. At 5.0 they were the
    # annotation size -- the same as "n = 43" and the inset's scale -- which made the
    # densest and most-read text on the plate its smallest.
    panel_a(ax_a, lab=5.4 * K)
    if BVAR == "side":
        ax_b = fig.add_axes([B0, BROW, BW, TOP - BROW])
        ax_b2 = fig.add_axes([B0 + BW + BGAP, BROW, BW, TOP - BROW])
    elif BVAR == "mfe":
        ax_b = fig.add_axes([B0, BROW, 2 * BW + BGAP, TOP - BROW])
        ax_b2 = None
    else:
        # Stacked in the column the pair shares, so each is as wide as the two side-by-side
        # ones together. The band between them is one x-label stack deep, and so is the band
        # under the lower one -- which is what keeps its x label, centred over a column two
        # fold pitches wide, clear of the structures' letters.
        # The upper panel's x label sits under its own tick row, 0.142 in below the panel
        # and 0.095 in deep, centred in the gap between the second and third tick labels
        # rather than on a line of its own, so this band is that stack plus clearance --
        # 0.275 in, and at 0.195 the label printed 8.4 px inside the lower panel's box.
        BMID = _f(.275)
        # (b) and (c) share a bottom, so the lower panel's x axis and (c)'s x axis are one
        # line across the plate. They used to differ: BROW was set by what (c) needed under
        # it and (b) was placed 0.16 in lower to use the page that (b)'s own shallower label
        # stack left empty. That bought (b) its box height and cost the row its one
        # horizontal rule -- two axes an eighth of an inch apart, which reads as a miss
        # rather than as two rows. (c) comes down to (b) instead: the band under them is
        # the same 0.295 in either way, and (c) gains the 0.16.
        BH = (TOP - BROW - BMID) / 2
        ax_b2 = fig.add_axes([B0, BROW, 2 * BW + BGAP, BH])
        ax_b = fig.add_axes([B0, BROW + BH + BMID, 2 * BW + BGAP, BH])
    panel_b2(ax_b, ax_b2)
    ax_c = fig.add_axes([C0, BROW, RIGHT - C0, TOP - BROW])
    panel_c(ax_c)
    fig.canvas.draw()
    if BVAR == "stack":
        # Right-aligned on the INK, which is the only alignment a reader can see. The two
        # labels already share an anchor -- a rotated label is anchored va="bottom" under
        # rotation_mode="anchor", and a quarter turn maps "bottom" onto +x, so the anchor
        # IS the layout box's right edge whatever the line count -- and that is not where
        # the ink ends: "(kcal/mol)" carries parentheses that descend 0.2 em past the
        # baseline, so its glyphs reach 6 px further right than "\u0394CAI"'s. Aligning
        # the boxes left the two labels 1 px apart on the LEFT and 6 px apart on the
        # right, which reads as a left alignment that slipped.
        # The UPPER label moves left, not the lower one right: the lower one already sits
        # 1.5 px from its own tick column, so closing the 6 px there would put a label on
        # its own ticks. Moving the unit back instead leaves both labels the same distance
        # from the ticks they belong to.
        _dx = _ink_x1(ax_b.yaxis.label) - _ink_x1(ax_b2.yaxis.label)
        _x = ax_b.yaxis.label.get_position()[0] - _dx / ax_b.get_window_extent().width
        ax_b.yaxis.set_label_coords(_x, .5)
        fig.canvas.draw()
    # The structures still start where (b)'s ink starts and never left of (a)'s.
    L0F = max(visual_x0(ax_b),
              ax_a.get_tightbbox(fig.canvas.get_renderer())
                  .transformed(fig.transFigure.inverted()).x1 + .006)
    STEPF = (RIGHT - L0F + FGAP) / len(P)
    # The captions, against the panels they will be drawn on. The band's DEPTH is CAP, a
    # constant, because layout 4 sets two lines on every panel; what has to be checked here
    # is their WIDTH, which is known only now that the row's pitch is.
    _r = fig.canvas.get_renderer()
    # The width the captions will actually have, not the width the columns start with: the
    # row gives its drawings' side slack back as spacing, so the panels end up narrower than
    # the pitch and a caption measured against the pitch is measured against a panel that
    # will not exist.
    FTOPF = FTOP - (_f(STRIPBAND) if BVAR == "mfe" else 0)
    # The two structures run down to the foot of the page, over the band their captions stand in.
    FBOTX = _f(.03) if BVAR == "mfe" else FBOT
    # The two-design variant fixes the column between the folds and gives the rest to them.
    MIDW = .90 / FIGW
    FWID = (RIGHT - L0F - MIDW) / 2
    _pw = (STEPF - FGAP) * FIGW if BVAR != "mfe" else FWID * FIGW
    _, _sx, _, _hx, _ = fold_box((FTOPF - FBOTX) * FIGH / _pw)
    if _hx > _sx * 1.001:
        _pw *= _sx / _hx * FOLD_SHRINK
    # The band is as deep as the form the row will actually SET, which means running the
    # same cascade here that folds() runs there. Measuring only the widest form and
    # reserving two lines whenever it overruns reserved them whenever the row then found a
    # one-line form that fits -- 0.13 in of page left empty under the captions, which is a
    # tenth of the structures' own height.
    # A caption wider than its panel runs into its neighbour's, and no gate sees it.
    for _k in range(len(P)):
        # The widest line as TYPE SETS it, not as it spells: "+314 \u2212315 bp" is half the
        # characters of the line above it and 0.054 in per character against 0.047, which
        # a longest-string test gets right only by luck. This file rejects that proxy for
        # the same reason twenty lines up.
        _w_in = 0.0
        for _line in cap_parts(_k):
            _t = fig.text(0, 0, _line, fontsize=CAPFS)
            _w_in = max(_w_in, _t.get_window_extent(_r).width / fig.dpi)
            _t.remove()
        # Not an assert: `python -O` drops those, and this is the only check standing
        # between a caption and its neighbour's. The overlap gate does not see an overrun
        # until it is wide enough to touch the next caption, which is 12% further out.
        if _w_in > _pw:
            raise SystemExit(f"fold {_k}'s caption sets {_w_in:.3f} in wide on a "
                             f"{_pw:.3f} in panel: shorten it, or lower CAPFS")
    if BVAR == "mfe":
        fax = [fig.add_axes([x, FBOTX, FWID, FTOPF - FBOTX]) for x in (L0F, RIGHT - FWID)]
    else:
        fax = [fig.add_axes([L0F + k * STEPF, FBOTX, STEPF - FGAP, FTOPF - FBOTX])
               for k in range(len(P))]
    folds(fax)
    if BVAR == "mfe":
        # Between the two folds, which folds() has narrowed to their drawings.
        gl, gr = fax[0].get_position().x1, fax[1].get_position().x0
        pad = .05 / FIGW
        # Between the two drawings' ink: the left one ends 0.02 in short of its box's right edge and
        # the right one starts 0.02 in inside its box's left edge, less the shift.
        x_l = gl - .021 / FIGW
        x_r = gr + (.021 - FOLD_SHIFT) / FIGW
        ax_m = fig.add_axes([x_l + .03 / FIGW, FBOT, x_r - x_l - .06 / FIGW, FTOPF - FBOT])
        middle(ax_m)
        # The protein, in the corner of (d) as (c) names it in its own.
        _t = fig.text(RIGHT, _f(.03), R["cell"], fontsize=5.0 * K, color=MUTED, ha="right",
                      va="bottom")
        TAGS[id(_t)] = f"{id(fax[1])}"
        site_insets_pending = True
        # The strips start clear of their row labels, which sit to the left of them.
        SX0 = L0F + .80 / FIGW
        for _k, _y in ((0, .18), (1, .28)):
            codon_strip(fig.add_axes([SX0, FTOP - _f(_y), RIGHT - SX0, _f(.08)]), _k)
            fig.text(SX0 - .05 / FIGW, FTOP - _f(_y - .04), MSET[_k],
                     ha="right", va="center", fontsize=FSS_, color=INK)
    if BVAR == "mfe":
        # One inset per occurrence of a motif in the left design, in order along the CDS: that
        # stretch of the CDS in the left design, with the motif in its class's colour, and the
        # same stretch in the right design. A letter the two designs share is grey and a letter
        # that differs is dark, as in the strips above. A line runs from the occurrence's codons
        # on the strip to its inset.
        _sq = (P[0]["seq"], P[1]["seq"])
        _IU = {"R": "[AG]", "Y": "[CU]", "N": "[ACGU]", "S": "[CG]", "W": "[AU]", "K": "[GU]",
               "M": "[AC]"}
        _sites = sorted(P[0]["sites"])
        _ncod = P[1]["n_codons"]
        _slot = (RIGHT - SX0) / len(_sites)
        _hin = .17
        _ytop = FTOP - _f(.28) - _f(.09)
        for _i, (_a, _b, _roman) in enumerate(_sites):
            _lo = max(0, _a - 2)
            _hi = min(_ncod * 3, _b + 2)
            _w = (_hi - _lo) * .0405 / FIGW
            _xc = SX0 + (_i + .5) * _slot
            _ax = fig.add_axes([_xc - _w / 2, _ytop - _f(_hin), _w, _f(_hin)])
            _ax.axis("off")
            _ax.set_xlim(0, _hi - _lo)
            _ax.set_ylim(0, 2)
            for _j, _pos in enumerate(range(_lo, _hi)):
                _c0, _c1 = _sq[0][_pos], _sq[1][_pos]
                _in = _a <= _pos < _b
                _top = SETC[_roman] if _in else ("#33322f" if _c0 != _c1 else "#8f8e89")
                _ax.text(_j + .5, 1.5, _c0, ha="center", va="center", fontsize=4.2 * K,
                         family="Nimbus Mono PS", color=_top, fontweight="bold" if _in else None)
                _ax.text(_j + .5, .5, _c1, ha="center", va="center", fontsize=4.2 * K,
                         family="Nimbus Mono PS", color="#33322f" if _c0 != _c1 else "#8f8e89")
            _ax.add_patch(Rectangle((0, 0), _hi - _lo, 2, fill=False, ec=OUTLINE, lw=.3,
                                    clip_on=False))
            # The motif, named as in (a): the pattern of its class that matches these letters.
            _m = next(e["motif"] for e in D["a"]["entries"] if e["category"] == _roman
                      and re.fullmatch("".join(_IU.get(c, c) for c in e["motif"].upper().replace("T", "U")),
                                       _sq[0][_a:_b]))
            _nm = NAMES.get(_m, _m).replace(", reverse complement", " (rc)")
            _t = fig.text(_xc, _ytop - _f(_hin + .02), _nm, ha="center", va="top", fontsize=4.9 * K,
                          color=SETC[_roman])
            TAGS[id(_t)] = f"{id(_ax)}"
            _xs = SX0 + ((_a // 3 + (_b - 1) // 3) / 2 + .5) / _ncod * (RIGHT - SX0)
            fig.add_artist(plt.Line2D([_xs, _xc], [FTOP - _f(.28), _ytop], lw=.4,
                                      color=SETC[_roman], transform=fig.transFigure))
    place_tags()
    LEG, LEGFS, LEGCOL = (.5, LEGY), (FSS_ if BVAR == "mfe" else 6.2 * K), None

# One line. It used to be thirteen over two ragged rows, and seven of those
# thirteen were the same four colours keyed twice: a hue names a motif set everywhere it is
# used, and a constraint level is drawn in the hue of the set it newly forbids, so "set I"
# and "+ I" were one swatch with two captions. The funnel's three nested shares are one
# entry with three swatches rather than three entries, since the three are read together or
# not at all.
from matplotlib.legend_handler import HandlerTuple  # noqa: E402

H = [plt.Line2D([], [], color=LEVEL[0], lw=0, marker="o", ms=2.6, label="no sets")]
H += [plt.Line2D([], [], color=SETC[k], lw=0, marker="o", ms=2.6) for k in ("I", "II", "III")]
# Outlined, and the palest tint pulled off the page: at .84 it was 8% from white, which in
# the key is a swatch a reader does not see at all -- the three-tone ramp the entry exists to
# teach came across as two.
H += [tuple(plt.Rectangle((0, 0), 1, 1, fc=c, ec=MUTED, lw=.3)
            for c in (tint(SETC["III"], .76), tint(SETC["III"], .52), SETC["III"]))]
# Hull status is drawn on all four fronts, so its two handles cannot be one front's colour:
# in grey they said "this applies to the unconstrained front", and the grey filled circle
# was also the "no sets" mark, two entries apart on one line.
# (b)'s strip: three shades again, but naming a different three things from (a)'s funnel, so
# it gets its own entry rather than letting a reader carry (a)'s meaning across.
if LAYOUT == 1 or BVAR == "strip":
    H += [tuple(plt.Rectangle((0, 0), 1, 1, fc=c, ec="none")
                for c in (LEVEL[3], tint(LEVEL[3], .55), "#e9e7e2"))]
# Hull membership is the one entry the narrow pages drop. Eight entries do not fit on one
# 7 in line and a second row costs more height than the distinction is worth here: filled
# against hollow is a within-panel distinction in (c) alone, so it goes in the caption rather
# than taking a sixth and seventh of the key from every other panel.
HULLKEY = LAYOUT < 3
if HULLKEY:
    H += [plt.Line2D([], [], color=INK, marker="o", ms=1.9, lw=.7, label="on hull"),
          plt.Line2D([], [], color=INK, marker="o", ms=1.9, lw=0, mfc="none", mew=.4,
                     label="off hull")]
# (d)'s Venn diagram is the key to its base pairs, in the same colours.
if BVAR != "mfe":
    H += [plt.Line2D([], [], color=GAINED, lw=1.4, label="pair not in (d)")]
# A hue names a motif set in (a) and on the folds, and the constraint LEVEL it newly
# forbids in (b) and (c) -- so orange is "a site of set II" in one place and "set I and II
# both forbidden" in another. One swatch, two meanings, which the key used to resolve in
# favour of the first and leave the cumulative reading unstated.
# "set I   + I" read as arithmetic: the two captions on one swatch were bound no more
# tightly to each other than the swatch was to the first of them. A middot says they are two
# readings of one colour, which is what the entry means.
LAB = [RUNG[0]] + [f"set {k} \u00b7 {RUNG[i + 1]}" for i, k in enumerate(("I", "II", "III"))] \
    + ["admits / spells / optimum"] \
    + (["binds / already gone / absent"] if LAYOUT == 1 or BVAR == "strip" else []) \
    + (["on hull", "off hull"] if HULLKEY else []) \
    + ([] if BVAR == "mfe" else ["pair not in (d)"])
assert len(LAB) == len(H), f"{len(LAB)} labels for {len(H)} handles"
if BVAR == "mfe":
    # The four colours of (b) and (c) are the four sets of forbidden motifs, in (d)'s notation.
    LAB[:4] = [MSET[0]] + ["$" + m + "$" for m in MLEV[1:]]
if BVAR == "mfe":
    # The three nested shares are (a)'s own key, so they hang over (a) and not in the page's.
    _fun, _funlab = H.pop(4), LAB.pop(4)
# Centred on the page rather than flush left, which is how the overview figure keys its
# panels: a left-flush row of eight entries ends short of the right margin and reads as a
# column that ran out, while a centred one reads as a caption to the whole plate.
if BVAR == "mfe":
    # No key across the page: each colour is keyed over the panel that draws it, on the baseline
    # of the panel letters and flush with the panel's right edge. (b) draws the three sets that
    # are forbidden and (c) the unconstrained front as well.
    for _hs, _ls, _ax in ((H[1:4], LAB[1:4], ax_b), (H[:1], LAB[:1], ax_c)):
        fig.legend(handles=_hs, labels=_ls, loc="lower right", frameon=False,
                   bbox_to_anchor=(_ax.get_position().x1, TAGY), fontsize=FSS_,
                   handlelength=1.0, labelspacing=.3, borderpad=0, borderaxespad=0,
                   ncol=len(_hs), columnspacing=1.2, handletextpad=.4)
else:
    fig.legend(handles=H, labels=LAB, loc="lower center", bbox_to_anchor=LEG, frameon=False,
               # borderpad=0: the key's own border padding was 0.09 in of empty box under its
               # ink, and the overlap gate reads boxes, so that padding held the panel letters
               # -- and with them every panel -- a tenth of an inch lower than the ink needed.
               fontsize=LEGFS, handlelength=1.4, labelspacing=.3, borderpad=0,
               ncol=LEGCOL or len(H), columnspacing=1.6,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=.25)})
if BVAR == "mfe":
    # (a)'s key is one stacked bar in the order the bars show it left to right: the darkest
    # share, the middle one, the lightest, and the white that is left of the 100%. Each segment is
    # at least as wide as its label, which stands inside it. The bar's top is the top of (a)'s
    # panel letter.
    _rend = fig.canvas.get_renderer()
    KEYFS = FSS_
    _lab = ["CAI = 0.8", "front", "any", "none"]
    _col = ["#0b0b0b", "#6b6a66", "#c4c3bf", "#ffffff"]
    _w = []
    for _l in _lab:
        _t = fig.text(0, 0, _l, fontsize=KEYFS)
        _w.append(_t.get_window_extent(_rend).width / fig.dpi + .05)
        _t.remove()
    _x0 = ax_a.get_position().x0 * FIGW
    # As wide as (a)'s axes: what the labels leave over is shared out between the segments.
    _extra = (ax_a.get_position().width * FIGW - sum(_w)) / len(_w)
    _w = [_wd + _extra for _wd in _w]
    _bar = fig.add_axes([_x0 / FIGW, TAGY + (LET - .125) / FIGH, sum(_w) / FIGW, .125 / FIGH])
    _bar.set_xlim(0, sum(_w))
    _bar.set_ylim(0, 1)
    _bar.axis("off")
    _at = 0.0
    for _l, _c, _wd in zip(_lab, _col, _w, strict=True):
        _bar.add_patch(Rectangle((_at, 0), _wd, 1, fc=_c, ec=OUTLINE, lw=.4))
        _t = fig.text((_x0 + _at + _wd / 2) / FIGW, TAGY + (LET - .0625) / FIGH, _l, ha="center",
                      va="center", fontsize=KEYFS, color="#ffffff" if _c in _col[:2] else INK,
                      zorder=5)
        TAGS[id(_t)] = f"{id(_bar)}"
        _at += _wd
if BVAR == "mfe":
    # The key to (d)'s codon strips, on the panel letter's line, flush with the strips' right edge: the codons that overlap a site of
    # each set's motifs in the left design, and the codons the two designs share.
    _cod = {r: set() for r in ("I", "II", "III")}
    for _a, _b, _r in P[0]["sites"]:
        _cod[_r] |= set(range(_a // 3, (_b - 1) // 3 + 1))
    _kh = [Rectangle((0, 0), 1, 1, fc=SETC[r], ec="none") for r in ("I", "II", "III")]
    _kh.append(Rectangle((0, 0), 1, 1, fc=SHAREDC, ec=OUTLINE, lw=.4))
    _kh.append(Rectangle((0, 0), 1, 1, fc="#4a4742", ec="none"))
    _kl = [r"$\mathcal{M}_{\mathrm{%s}}$: %d codons" % (r, len(_cod[r])) for r in ("I", "II", "III")]
    _kl.append(f"shared: {P[1]['n_codons'] - len(P[1]['codons_differ'])} codons")
    _kl.append(f"different: {len(P[1]['codons_differ'])} codons")
    fig.legend(handles=_kh, labels=_kl, loc="lower right", frameon=False, ncol=5,
               bbox_to_anchor=(RIGHT, TAGF - .012), fontsize=FSS_, handlelength=.9,
               handleheight=.8, columnspacing=1.2, borderpad=0, handletextpad=.4)
def overlaps(fig, clear=4.0):
    """Every pair of drawn labels that touch, in device pixels.

    Checked rather than eyeballed: at this size a collision is a few pixels and survives a
    read-through. Text inside an axes is skipped against its own axes, and an inset against its
    parent, since those are containments rather than collisions.
    """
    fig.canvas.draw()
    items = []
    for ax in fig.axes:
        shown = [ax.xaxis.label, ax.yaxis.label, ax.title] + list(ax.texts)
        # A tick outside the view still has a label artist with a window extent, and it is not
        # drawn; counting those reports collisions that are not on the page.
        # An axes drawn with axis("off") keeps its tick-label artists, and they still report
        # themselves visible; counting them invents collisions for every fold panel.
        for t in shown:
            if t.get_visible() and t.get_text().strip():
                items.append((f"{id(ax)}", t.get_text().strip()[:28], t.get_window_extent()))
        # Tick labels carry their own owner, one per axis, so they are not compared with each
        # other -- matplotlib lays those out -- but ARE compared with this axes' own axis
        # labels. Grouping them under the axes hid the collision this test exists to catch:
        # an x label positioned by hand landed on its own tick row, and the pair was skipped
        # for sharing an owner.
        for axis, lim, who in (() if not ax.axison else
                               ((ax.xaxis, ax.get_xlim(), "tx"),
                                (ax.yaxis, ax.get_ylim(), "ty"))):
            lo, hi = sorted(lim)
            for loc, lab in zip(axis.get_ticklocs(), axis.get_ticklabels(),
                                strict=True):
                if lo - 1e-9 <= loc <= hi + 1e-9 and lab.get_visible() \
                        and lab.get_text().strip():
                    items.append((f"{id(ax)}:{who}", lab.get_text().strip()[:28],
                                  lab.get_window_extent()))
    for t in fig.texts:
        if t.get_text().strip():
            # A panel letter is drawn on the figure but belongs to its panel; left as "fig"
            # it is tested against the very box it labels and reported as intruding on it.
            items.append((TAGS.get(id(t), "fig"), t.get_text().strip()[:28],
                          t.get_window_extent()))
    if fig.legends:
        items.append(("fig", "<legend>", fig.legends[0].get_window_extent()))
    bad = []
    # Nothing may be drawn outside the canvas. With no tight box a label that overruns is
    # silently cut off rather than growing the output, so the check has to be explicit.
    W, H = (fig.get_size_inches() * fig.dpi)
    for _o, text, bb in items:
        if bb.x0 < -.5 or bb.y0 < -.5 or bb.x1 > W + .5 or bb.y1 > H + .5:
            bad.append((text, "<off the canvas>", round(bb.x0, 1), round(bb.y1, 1)))
    # A label straying into a panel it does not belong to. Comparing text against text cannot
    # see this, and it is the class that bit this figure: (b)'s "CAI" sat centred over fold (d)
    # and read as that panel's caption without touching any other label.
    boxes = {f"{id(ax)}": ax.get_window_extent() for ax in fig.axes}
    for owner, text, bb in items:
        for other, ob in boxes.items():
            if other == owner:
                continue
            dx = min(bb.x1, ob.x1) - max(bb.x0, ob.x0)
            dy = min(bb.y1, ob.y1) - max(bb.y0, ob.y0)
            # The same clearance the text-vs-text test below demands. Requiring actual
            # penetration here is what let (b)'s "CAI" sit 2.4 px off fold (d)'s box, on the
            # fold letters' own baseline, and still report clean.
            if dx > -clear and dy > -clear:
                bad.append((text, "<a panel it is not part of>", round(dx, 1), round(dy, 1)))
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            (oi, ti, bi), (oj, tj, bj) = items[i], items[j]
            if oi == oj:
                continue
            dx = min(bi.x1, bj.x1) - max(bi.x0, bj.x0)
            dy = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
            # Labels that merely come close read as touching at this size, so the test is
            # a clearance requirement rather than an intersection one.
            if dx > -clear and dy > -clear:
                bad.append((ti, tj, round(dx, 1), round(dy, 1)))
    return bad


_bad = overlaps(fig)
for ti, tj, dx, dy in _bad:
    print(f"  OVERLAP  {ti!r} x {tj!r}  by {dx} x {dy} px")
print(f"overlapping label pairs: {len(_bad)}")
# The check gates the save. Printed and not acted on, it was a number beside a figure that
# shipped either way, which is not a check.
assert not _bad, f"{len(_bad)} label collision(s); the figure was not written"

SLUG = ("fig-composite" if LAYOUT == 1 else f"fig-composite-L{LAYOUT}") \
       + ("" if SUBSET == "finished" else f"-{SUBSET}") \
       + ("" if BVAR == "strip" else f"-b{BVAR}")
# Byte-stable outputs, as the repo's other figures are written: matplotlib stamps the
# wall clock into a PDF's CreationDate and its own version into a PNG's Software chunk, so
# two runs of one script over one input produced two different digests and a manifest over
# them could only ever fail. The PNG is written at 300 dpi, the resolution the repo's
# figures are checked in at, rather than at the 200 the layout is measured on.
DST = sys.argv[8] if len(sys.argv) > 8 else f"{OUT}/{SLUG}.pdf"
STEM = DST[:-4] if DST.endswith(".pdf") else DST
os.makedirs(os.path.dirname(os.path.abspath(STEM)), exist_ok=True)
fig.savefig(f"{STEM}.pdf", metadata={"CreationDate": None})
fig.savefig(f"{STEM}.png", dpi=300, metadata={"Software": None})


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# A sidecar beside every figure: what it was drawn from and what drew it. The data file carries
# its own provenance for the measurements, so this records only the drawing step and points at
# that, rather than copying it and letting the two drift.
json.dump({
    "figure": SLUG,
    "written_utc": datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"),
    "layout": LAYOUT, "subset": SUBSET, "instance": INST, "b_variant": BVAR,
    "python": sys.version.split()[0],
    "matplotlib": matplotlib.__version__,
    "numpy": np.__version__,
    "inputs": dict({f: _sha(f"{OUT}/{f}") for f in
                    ("composite-data.json", f"example-{INST}.json")
                    + (() if BVAR == "strip" else ("dcai-at-mfe.json",))},
                   **{"motif_names.json":
                      _sha(os.environ.get("E2_MOTIF_NAMES", f"{OUT}/motif_names.json"))},
                   # The two scripts that decide what is drawn: this one, and fig2, which
                   # picks which front points are filled and orients every fold. The figure
                   # is not reproducible from the data files alone.
                   **{"fig_composite.py": _sha(os.path.abspath(__file__)),
                      "fig2.py": _sha(fig2.__file__)}),
    "data_provenance_written_utc": D.get("provenance", {}).get("written_utc"),
    # Counts alone cannot report a membership change that leaves the total fixed: one protein
    # out, another in, and every population number reads identical. The digests are over the
    # sorted accession lists, so they move whenever the members do; the lists themselves stay
    # in the data file's provenance rather than being copied here and left to drift.
    "populations": D.get("provenance", {}).get("populations"),
    "member_digests": D.get("provenance", {}).get("member_digests"),
    # The example file carries no provenance block of its own; its ViennaRNA version is the
    # one fact in it the drawing depends on, so it is copied here rather than left behind.
    "example_vienna": R.get("vienna"),
    "example_plot_type": R.get("plot_type"),
    "outputs": {f"{os.path.basename(STEM)}.{e}": _sha(f"{STEM}.{e}") for e in ("png", "pdf")},
}, open(f"{STEM}.prov.json", "w"), indent=1)
print(f"wrote {STEM}  ({fig.get_size_inches()[0]:.1f} x {fig.get_size_inches()[1]:.1f} in)")
