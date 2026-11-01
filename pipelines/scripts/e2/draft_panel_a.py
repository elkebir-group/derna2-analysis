"""Draft: panel (a) rebuilt on E1's protein set, with a full-scale outline per motif.

Three changes under consideration, drawn at three candidate thresholds so the threshold can
be chosen by looking rather than argued:

  * the population is E1's instance list rather than E2's own screen population;
  * every motif row carries an outline box running the full width, so a short bar reads as
    a small share rather than as a missing row;
  * the solid bar is the design a reader would actually ship -- the optimum at a CAI floor
    -- and the floor is a round number instead of 0.87.
"""
import collections
import json
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig2
from e3paths import E1_EXTRA, E1_SWEEP, MOTIF_NAMES, OUT, E  # noqa: F401

FLOORS = [0.80]
NAMES = json.load(open(MOTIF_NAMES))

XTRA = E1_EXTRA
# E1's instance list. The four P0DTC2_pNN values in the sweep are ladder rungs that appear
# under one arm only, not instances; the two proteins added to E1 today are not in this
# host's 10-27 snapshot and are named explicitly.
E1 = {json.loads(line)["protein"] for line in open(E1_SWEEP) if line.strip()}
E1 = {p for p in E1 if not p.startswith("P0DTC2_p")} | {"P22301", "P41159"}
cat = fig2.catalog(E)
SCREENED = set(cat["population_members"]) | {os.path.basename(f)[:-6]
                                             for f in __import__("glob").glob(f"{XTRA}/f0/*.jsonl")}
POP = sorted(SCREENED & E1)
NOT_SCREENED = sorted(E1 - SCREENED)
DROPPED = sorted(SCREENED - E1)


def rows(pattern):
    out = []
    for path in sorted(__import__("glob").glob(f"{E}/screen/{pattern}")):
        for line in open(path):
            if line.strip():
                rec = json.loads(line)
                if "motif" in rec:
                    out.append(rec)
    return out


def extra(name):
    return [json.loads(line) for line in open(f"{XTRA}/{name}.jsonl") if line.strip()]


# E2's own screen, plus the E1 instances it never covered — screened here from E1's own
# walk fronts by the driver's functions, so the rows are the same shape and the same code.
binding = rows("binding*.jsonl") + extra("binding-e2extra")
lattice = rows("lattice*.jsonl") + extra("lattice-e2extra")
occurs, spellable, category = (collections.defaultdict(set), collections.defaultdict(set), {})
for r in binding:
    category[r["motif"]] = fig2.SET_OF[r["set"]]
    if r["protein"] in POP and r["designs_carrying"] > 0:
        occurs[r["motif"]].add(r["protein"])
for r in lattice:
    if r["protein"] in POP and r["spellable"] > 0:
        spellable[r["motif"]].add(r["protein"])

entries = []
for _, roman, _ in fig2.CATEGORIES:
    pats = [m for m, c in category.items() if c == roman]
    for m in sorted(pats, key=lambda m: -len(occurs[m])):
        entries.append({"motif": m, "category": roman,
                        "occurs": len(occurs[m]), "spellable": len(spellable[m])})

# the optimum at each floor, over the same population
SETS = [(roman, [e["motif"] for e in entries if e["category"] == roman])
        for _, roman, _ in fig2.CATEGORIES]
# Each class as one pattern set: the proteins with some motif of the class in any CDS and on
# the front, which is the union of the class's motifs and not a sum of their counts.
classes = [{"category": roman,
            "occurs": len(set().union(*(occurs[m] for m in ms))),
            "spellable": len(set().union(*(spellable[m] for m in ms)))}
           for roman, ms in SETS]
census = {}
for fl in FLOORS:
    counts = {e["motif"]: 0 for e in entries}
    class_counts = {roman: 0 for roman, _ in SETS}
    seen = 0
    for cell in POP:
        path = fig2.f0_front(E, cell) or f"{XTRA}/f0/{cell}.jsonl"
        if not os.path.exists(path):
            continue
        if cell in ("P0DTC2", "PRT-04", "PRT-05", "PRT-06", "PRT-07"):
            # screen/f0-xl/ holds a 5-weight lambda grid of single optima, not a front;
            # E1's completed walk is the front, so prefer it for these.
            path = f"{XTRA}/f0/{cell}.jsonl"
        try:
            design = fig2.pick(fig2.read_front(path)[1], fl)
        except SystemExit:
            continue
        seen += 1
        spelling = design["sequence"].replace("U", "T")
        hit = {m: bool(re.search(fig2.expand(m), spelling)) for m in counts}
        for m, h in hit.items():
            counts[m] += h
        for roman, ms in SETS:
            class_counts[roman] += any(hit[m] for m in ms)
    census[fl] = {"counts": counts, "classes": class_counts, "proteins": seen}

PALETTE = {"I": "#b8651a", "II": "#8c1b5e", "III": "#2a6b4f"}
INK, MUTED, GRID, OUTLINE = "#0b0b0b", "#4a4742", "#e4e3df", "#bfbdb7"
plt.rcParams.update({
    "font.size": 6.4, "axes.labelsize": 6.2, "xtick.labelsize": 5.6, "ytick.labelsize": 5.2,
    "axes.edgecolor": MUTED, "axes.linewidth": .5, "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "axes.labelcolor": INK, "figure.dpi": 200, "savefig.bbox": "tight",
    "savefig.pad_inches": .04, "pdf.fonttype": 42})


def tint(hex_, k):
    r, g, b = (int(hex_[i:i + 2], 16) for i in (1, 3, 5))
    return "#{:02x}{:02x}{:02x}".format(*(round(v + (255 - v) * k) for v in (r, g, b)))


fig, axes = plt.subplots(1, 1, figsize=(4.6, 4.4))
axes = [axes]
n = len(POP)
y = list(range(len(entries)))[::-1]
for ax, fl in zip(axes, FLOORS, strict=True):
    cen = census[fl]
    for row, e in zip(y, entries, strict=True):
        hue = PALETTE[e["category"]]
        # the full-scale outline: a short bar now reads as a small share, not a missing row
        ax.add_patch(Rectangle((0, row - .40), 100, .80, fill=False, ec=OUTLINE, lw=.45, zorder=1))
        ax.barh(row, 100 * e["spellable"] / n, height=.80, color=tint(hue, .82), lw=0, zorder=2)
        ax.barh(row, 100 * e["occurs"] / n, height=.80, color=tint(hue, .50), lw=0, zorder=3)
        ax.barh(row, 100 * cen["counts"][e["motif"]] / cen["proteins"], height=.80,
                color=hue, lw=0, zorder=4)
    ax.set_yticks(y)
    # Name the site, not its spelling: a reader knows BsaI, not GGTCTC.
    ax.set_yticklabels([NAMES.get(e["motif"], e["motif"]).replace(", reverse complement", " (rc)")
                        for e in entries], fontsize=5.2)
    for tick, e in zip(ax.get_yticklabels(), entries, strict=True):
        tick.set_color(PALETTE[e["category"]])
    ax.set_ylim(-.75, len(entries) - .25)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel(f"% of {n} proteins", fontsize=6.0, labelpad=2.4)
    ax.grid(axis="x", color=GRID, lw=.35)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=1.2)
    ax.set_title(f"the shipped optimum is taken at CAI \u2265 {fl:.2f}",
                 fontsize=6.2, loc="left", pad=5, color=MUTED)

h = [Rectangle((0, 0), 1, 1, fill=False, ec=OUTLINE, lw=.6),
     Rectangle((0, 0), 1, 1, color=tint("#555555", .82)),
     Rectangle((0, 0), 1, 1, color=tint("#555555", .50)),
     Rectangle((0, 0), 1, 1, color="#555555")]
fig.legend(h, ["all proteins (100%)", "the lattice admits it",
               "some design on the front spells it", "the shipped optimum spells it"],
           loc="lower center", ncol=2, frameon=False, fontsize=5.4, handlelength=1.5,
           bbox_to_anchor=(.5, -.10))
fig.suptitle("Draft A — the motif funnel, named sites, CAI floor 0.80",
             fontsize=7.6, x=.02, ha="left", y=.99)
fig.savefig(f"{OUT}/draft-panel-a.png")
print("wrote draft-panel-a")
print(f"population: {n} proteins (E2 screen+extra {len(SCREENED)}, E1 {len(E1)})")
print(f"  in the screen but not in E1, so dropped ({len(DROPPED)}): {' '.join(DROPPED)}")
print(f"  in E1 but never screened ({len(NOT_SCREENED)}): {' '.join(NOT_SCREENED)}")
for fl in FLOORS:
    c = census[fl]
    tot = {r: sum(c["counts"][m] for m in pats) for r, pats in SETS}
    print(f"  floor {fl:.2f}: optimum resolved on {c['proteins']}/{n};"
          f"  set hits I {tot['I']}  II {tot['II']}  III {tot['III']}")
