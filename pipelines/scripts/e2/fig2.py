"""Figure 2: what the nested motif constraints forbid, and what forbidding them costs.

Seven panels on one subject and one CAI lower bound:

  (a) the motif catalog -- every pattern in categories I-III, how often it occurs across the
      screen, and the lattice-admissibility control beside it
  (b) the four fronts in (CAI, MFE), filled where a point is on the recomputed convex hull
  (c) the cost of each nest above the unconstrained front, with the bound swept
  (d)-(g) the chosen design at each nest, with motif sites marked by category

Colour carries one meaning throughout: a category. Nest Fn is drawn in the colour of the
category it newly forbids, so (a), (b), (c) and the marks in (d)-(g) share one key, and the
deepest nest has no colour left because it spells nothing.

Usage:
    python fig2.py E2_DIR OUT_PDF OUT_JSON [--cell P47211] [--cai-floor 0.85]
"""

import argparse
import bisect
import collections
import glob
import hashlib
import json
import math
import os
import re
import sys
import time

import matplotlib
import RNA

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

CATEGORIES = (
    ("set1-manufacturing", "I", "manufacturing"),
    ("set2-processing", "II", "processing"),
    ("set3-repeats", "III", "repeats"),
)
SET_OF = {"manufacturing": "I", "processing": "II", "repeats": "III", "initiation": "IV"}
NESTS = ("F0", "F1", "F2", "F3")
FORBIDS = {"F0": None, "F1": "I", "F2": "II", "F3": "III"}
NEST_TITLE = {
    "F0": "unconstrained",
    "F1": r"$+$ cat I",
    "F2": r"$+$ cat II",
    "F3": r"$+$ cat III",
}
PALETTE = {"I": "#2a78d6", "II": "#d1651b", "III": "#8f4fd1"}
PLAIN = "#9b9a96"
INK, INK_2, INK_3, GRID = "#0b0b0b", "#4e4c47", "#8b887f", "#e4e3df"
FS, FSS, FST = 7.0, 6.0, 5.2

# ViennaRNA's naive radial layout. RNApuzzler avoids helix crossings but lays the
# constrained nests out four times flatter than the unconstrained one off a four-base-pair
# difference, with backbone jumps over a hundred times the median step; the panels then
# differ for a reason unrelated to the figure's claim.
PLOT_TYPE = 0


def style():
    sns.set_theme(style="white", context="paper")
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "text.latex.preamble": r"\usepackage{amsmath}",
            "font.size": FS,
            "axes.titlesize": FS,
            "axes.labelsize": FS,
            "xtick.labelsize": FSS,
            "ytick.labelsize": FSS,
            "axes.linewidth": 0.5,
            "axes.edgecolor": INK_2,
            "xtick.color": INK_2,
            "ytick.color": INK_2,
            "xtick.major.size": 2.0,
            "ytick.major.size": 2.0,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


# ---------------------------------------------------------------- motif patterns and sites


def read_patterns(root):
    out = []
    for fname, roman, _ in CATEGORIES:
        pats = []
        for line in open(f"{root}/motifs/{fname}.txt"):
            pat = line.split("#")[0].strip()
            if pat:
                pats.append(pat)
        if not pats:
            raise SystemExit(f"{fname}: no patterns read")
        out.append((roman, pats))
    return out


IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[CG]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}


def expand(pattern):
    """IUPAC to a zero-width regex, so overlapping occurrences are all found.

    The whole alphabet, not only `R`. The solver reads full IUPAC; a figure that understood
    one code would enforce nothing on a pattern spelled with any other and report a census of
    zero, which reads as absence -- in the one direction the ladder's claim rests on. The
    files carry only `R` today, so this costs nothing and closes a silent failure.
    """
    unknown = sorted(set(pattern) - set(IUPAC))
    if unknown:
        raise SystemExit(f"{pattern}: not IUPAC: {unknown}")
    return "(?=" + "".join(IUPAC[c] for c in pattern) + ")"


def sites(sequence, sets):
    """Merged occurrence intervals per category, plus a position -> category map.

    Occurrences of one category that overlap are one site; a position claimed by more than
    one category is attributed to the shallowest, which is the one that forbids it first.
    """
    dna = sequence.replace("U", "T")
    runs, cover = {}, {}
    for roman, pats in sets:
        found = []
        for pat in pats:
            found += [(m.start(), m.start() + len(pat)) for m in re.finditer(expand(pat), dna)]
        merged = []
        for lo, hi in sorted(found):
            if merged and lo <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], hi)
            else:
                merged.append([lo, hi])
        runs[roman] = merged
    for roman, _ in reversed(sets):
        for lo, hi in runs[roman]:
            for i in range(lo, hi):
                cover[i] = roman
    return cover, runs


def forbidden_by(nest):
    """The categories a nest forbids: depth n forbids the first n categories."""
    return {roman for _, roman, _ in CATEGORIES[: int(nest[1:])]}


def scan(nest, design, sets):
    """Every category is scanned on every design, so a forbidden one reports a measured zero.

    A design that spells what its own nest forbids is a defect in the data or in the match,
    not something to draw; the figure refuses to render rather than show it as absence.
    """
    cover, runs = sites(design["sequence"], sets)
    violated = {r: runs[r] for r in forbidden_by(nest) if runs[r]}
    if violated:
        where = "; ".join(f"{r} at {v[0][0]}" for r, v in violated.items())
        raise SystemExit(f"{nest}: design spells a category it forbids ({where})")
    counts = {roman: len(runs[roman]) for _, roman, _ in CATEGORIES}
    cover = {i: r for i, r in cover.items() if r not in forbidden_by(nest)}
    return cover, runs, counts


# ------------------------------------------------------------------------------ front data


SETTLE_AGE = 600.0


def settled(root, cell, nest):
    """The newest round that has finished writing.

    A round counts once it has a report, or once it has gone untouched for SETTLE_AGE
    seconds. Both tests are needed and neither alone is sound: the middle band writes
    rounds continuously, so a front still being appended to would be read short; but the
    short band's driver writes no per-round report at all, so a report-only test silently
    falls back to round 0 and reads a cell's oldest and worst front as its current state.
    """
    if nest == "F0":
        for sub in ("screen/f0", "screen/f0-l"):
            path = f"{root}/{sub}/{cell}.jsonl"
            if os.path.exists(path):
                return path
        raise SystemExit(f"{cell}: no unconstrained front")
    cands = []
    for sub in ("grid", "mid"):
        base = f"{root}/{sub}/{cell}_{nest}.jsonl"
        if os.path.exists(base):
            cands.append((0, base))
        now = time.time()
        for path in glob.glob(f"{root}/{sub}/{cell}_{nest}.r*.jsonl"):
            done = os.path.exists(path.replace(".jsonl", ".report.json"))
            if done or now - os.path.getmtime(path) > SETTLE_AGE:
                cands.append((int(path.split(".r")[1].split(".")[0]), path))
    if not cands:
        raise SystemExit(f"{cell} {nest}: no settled front")
    return max(cands)[1]


def read_front(path):
    lines = open(path).read().splitlines()
    head = json.loads(lines[0])
    pts = [json.loads(ln) for ln in lines[1:] if ln]
    pts = sorted((p for p in pts if p.get("record") == "point"), key=lambda p: p["cai"])
    if not pts:
        raise SystemExit(f"{path}: no points")
    return head, pts


def hull_vertices(points):
    """Monotone chain over the lower hull in (log CAI, MFE).

    The hull is convex in the scalarized space, not in (CAI, MFE); and the stored `on_hull`
    flag records a search tolerance rather than membership, so it is recomputed here.
    """
    order = sorted(range(len(points)), key=lambda i: (points[i]["cai_log"], points[i]["mfe"]))
    chain = []
    for i in order:
        x, y = points[i]["cai_log"], points[i]["mfe"]
        while len(chain) >= 2:
            (x1, y1), (x2, y2) = chain[-2][1], chain[-1][1]
            if (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1) <= 0:
                chain.pop()
            else:
                break
        chain.append((i, (x, y)))
    return {i for i, _ in chain}


def suffix_min(points):
    """Best MFE reachable at or above each CAI, which is what a lower bound buys."""
    cai = [p["cai"] for p in points]
    best, running = [0.0] * len(points), float("inf")
    for i in range(len(points) - 1, -1, -1):
        running = min(running, points[i]["mfe"])
        best[i] = running
    return cai, best


def best_at(table, c):
    cai, best = table
    i = bisect.bisect_left(cai, c)
    return best[i] if i < len(cai) else None


def pick(points, floor):
    feasible = [p for p in points if p["cai"] >= floor]
    if not feasible:
        raise SystemExit(f"no point at CAI >= {floor}")
    return min(feasible, key=lambda p: (p["mfe"], -p["cai"], p["sequence"]))


def build_line(path):
    """The build a round ran under -- only where the record names the file actually read.

    One `.builds` file covers every round of a cell, so taking its last line credits the
    newest round's build to whichever round was read. It misfires silently and in the most
    convincing way available: the line returned carries a timestamp, a filename and a commit
    sha, so it reads as provenance. Measured: P26436's F1 was read from r1 and credited to an
    r2 that does not exist on disk. Match on the basename or say nothing.
    """
    builds = re.sub(r"\.r\d+\.jsonl$|\.jsonl$", ".builds", path)
    want = os.path.basename(path)
    if os.path.exists(builds):
        for line in reversed([ln.strip() for ln in open(builds) if ln.strip()]):
            if want in line.split():
                return line
    summary = os.path.join(os.path.dirname(path), "summary.csv")
    stem = os.path.basename(path).split(".")[0]
    if os.path.exists(summary) and "_" in stem:
        cell, nest = stem.split("_", 1)
        for row in open(summary):
            f = row.split(",")
            if len(f) > 4 and f[0] == cell and f[2] == nest:
                return f"summary.csv build={f[4]}"
    return "unrecorded"


# ------------------------------------------------------------------- the occurrence screen


def catalog(root):
    """Occurrence and lattice-admissibility per pattern, over one named population.

    The lattice screen covers a band the occurrence screen does not, so both rates are taken
    over the proteins the two screens share and the excluded band is named in the record.
    """

    def rows(pattern):
        out = []
        for path in sorted(glob.glob(f"{root}/screen/{pattern}")):
            band = os.path.basename(path).replace(".jsonl", "")
            for line in open(path):
                if line.strip():
                    rec = json.loads(line)
                    if "motif" in rec:
                        rec["band"] = band
                        out.append(rec)
        return out

    binding, lattice = rows("binding*.jsonl"), rows("lattice*.jsonl")
    if not binding or not lattice:
        raise SystemExit("screen: binding or lattice rows missing")
    population = {r["protein"] for r in binding}
    excluded = sorted({r["protein"] for r in lattice} - population)

    occurs = collections.defaultdict(set)
    spellable = collections.defaultdict(set)
    category = {}
    for r in binding:
        category[r["motif"]] = SET_OF[r["set"]]
        if r["designs_carrying"] > 0:
            occurs[r["motif"]].add(r["protein"])
    for r in lattice:
        if r["protein"] in population and r["spellable"] > 0:
            spellable[r["motif"]].add(r["protein"])

    entries = []
    for _, roman, _ in CATEGORIES:
        pats = [m for m, c in category.items() if c == roman]
        for motif in sorted(pats, key=lambda m: -len(occurs[m])):
            entries.append(
                {
                    "motif": motif,
                    "category": roman,
                    "occurs": len(occurs[motif]),
                    "spellable": len(spellable[motif]),
                }
            )
    aa = {r["protein"]: r["aa"] for r in binding}
    return {
        "population": len(population),
        # The count alone makes the exclusion unauditable: an arithmetic that balances for
        # any 19 proteins balances for the wrong 19. The members are cheap to carry.
        "population_members": sorted(population),
        "bands": sorted({r["band"] for r in binding}),
        "aa_range": [min(aa.values()), max(aa.values())],
        "excluded_from_lattice_screen": excluded,
        "entries": entries,
    }


def f0_front(root, cell):
    """The unconstrained front for one protein, wherever its band put it."""
    for band in ("f0", "f0-l", "f0-five", "f0-xl"):
        path = f"{root}/screen/{band}/{cell}.jsonl"
        if os.path.exists(path):
            return path
    return None


def optimum_census(root, cat, floor, sets):
    """Per pattern: on how many proteins the *optimal* design at the bound spells it.

    The catalogue counts a pattern as occurring when any design on a protein's front spells
    it, which is a statement about the front's reach. This is the narrower question -- does
    the design a reader would actually ship contain the site -- and the two separate by a
    factor of five between categories. It needs only the unconstrained front, so it inherits
    no constrained front's time budget.
    """
    counts = {e["motif"]: 0 for e in cat["entries"]}
    seen, missing = 0, []
    for cell in cat["population_members"]:
        path = f0_front(root, cell)
        if not path:
            missing.append(cell)
            continue
        try:
            design = pick(read_front(path)[1], floor)
        except SystemExit:
            missing.append(cell)
            continue
        spelling = design["sequence"].replace("U", "T")
        for _, patterns in sets:
            for q in patterns:
                if q in counts and re.search(expand(q), spelling):
                    counts[q] += 1
        seen += 1
    return {"counts": counts, "proteins": seen, "unavailable": missing, "floor": floor}


def front_sweep(root, cat, sets, grid):
    """Site density in the optimal design, as the CAI bound walks up the front.

    One design per protein per bound -- the minimum-MFE design at that bound, which is what
    a reader would ship -- and the density is per kilobase so proteins of different length
    pool. An MFE cutoff selects the same one-parameter family from the other end, because a
    front is a curve: sweeping either coordinate walks the same designs.
    """
    fronts = {}
    for cell in cat["population_members"]:
        path = f0_front(root, cell)
        if path:
            fronts[cell] = read_front(path)[1]
    out = []
    for bound in grid:
        per = {roman: [] for roman, _ in sets}
        for pts in fronts.values():
            if max(p["cai"] for p in pts) < bound:
                continue
            design = pick(pts, bound)
            kb = len(design["sequence"]) / 1000.0
            runs = sites(design["sequence"], sets)[1]
            for roman, _ in sets:
                per[roman].append(len(runs[roman]) / kb)
        if not per[sets[0][0]]:
            continue
        out.append(
            {
                "cai": round(bound, 4),
                "proteins": len(per[sets[0][0]]),
                "per_kb": {r: round(sum(v) / len(v), 4) for r, v in per.items()},
                "carrying": {r: sum(1 for x in v if x) / len(v) for r, v in per.items()},
            }
        )
    return out


def occurrence_frequency(root, cat, sets, floor, high):
    """How many times each pattern appears, not whether it appears.

    Three counts per pattern, all means over the screened population and all about the
    *optimal* design rather than the front's reach:

      `slots`    how many positions the codon lattice could put the pattern at, which is the
                 count the presence panel reduces to a yes/no;
      `at_floor` occurrences in the design chosen at the figure's bound;
      `at_high`  the same at a higher bound, which is what answers whether the frequency
                 moves with CAI -- and it moves in both directions, so a single bound hides
                 patterns that double and patterns that vanish.
    """
    slots = collections.defaultdict(list)
    population = set(cat["population_members"])
    for path in sorted(glob.glob(f"{root}/screen/lattice*.jsonl")):
        for line in open(path):
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("protein") in population and "motif" in row:
                slots[row["motif"]].append(row["spellable"])
    counts = {e["motif"]: {"slots": 0.0, "at_floor": 0.0, "at_high": 0.0} for e in cat["entries"]}
    seen = {floor: 0, high: 0}
    for cell in cat["population_members"]:
        path = f0_front(root, cell)
        if not path:
            continue
        points = read_front(path)[1]
        for bound, key in ((floor, "at_floor"), (high, "at_high")):
            if max(p["cai"] for p in points) < bound:
                continue
            spelling = pick(points, bound)["sequence"].replace("U", "T")
            seen[bound] += 1
            for motif in counts:
                counts[motif][key] += len(re.findall(expand(motif), spelling))
    for motif, row in counts.items():
        row["slots"] = sum(slots.get(motif, [0])) / max(1, len(slots.get(motif, [0])))
        row["at_floor"] /= max(1, seen[floor])
        row["at_high"] /= max(1, seen[high])
    return {"counts": counts, "floor": floor, "high": high, "proteins": seen}


def set_prices(root, floor, sets):
    """What each rung cost, and how many sites it had to remove, per cell.

    Two controls travel with the price. `removed` is the binding control: a rung that
    removes nothing priced a constraint that forbade nothing, and its zero is not a cheap
    constraint. `exact` says the bound was resolved on both fronts the difference is taken
    between -- without it a deeper front can read *better* than a shallower one, which is
    impossible under nesting and happens on 12 of 17 cells here.
    """
    cells = sorted(
        {
            os.path.basename(p).split("_")[0]
            for d in ("grid", "mid")
            for p in glob.glob(f"{root}/{d}/*_F3*.jsonl")
        }
    )
    new = {"F1": "I", "F2": "II", "F3": "III"}
    out = []
    for cell in cells:
        try:
            paths = {"F0": f0_front(root, cell)} | {n: settled(root, cell, n) for n in NESTS[1:]}
            if not paths["F0"]:
                continue
            fronts = {n: read_front(p) for n, p in paths.items()}
            chosen = {n: pick(fronts[n][1], floor) for n in NESTS}
        except (SystemExit, OSError):
            continue
        exact = {
            n: all(b["mfe_lo"] > chosen[n]["mfe"] for b in (fronts[n][0].get("open") or []))
            for n in NESTS
        }
        for deep, shallow in (("F1", "F0"), ("F2", "F1"), ("F3", "F2")):
            roman = new[deep]
            out.append(
                {
                    "cell": cell,
                    "set": roman,
                    "mfe_cost": round(chosen[deep]["mfe"] - chosen[shallow]["mfe"], 2),
                    "removed": len(sites(chosen[shallow]["sequence"], sets)[1][roman]),
                    "exact": exact[deep] and exact[shallow],
                    "ceiling_drop": round(
                        max(p["cai"] for p in fronts[shallow][1])
                        - max(p["cai"] for p in fronts[deep][1]),
                        6,
                    ),
                }
            )
    return out


# ------------------------------------------------------------------------- structure plots


def coords(structure):
    RNA.cvar.rna_plot_type = PLOT_TYPE
    xy = RNA.get_xy_coordinates(structure)
    return (
        [xy.get(i).X for i in range(len(structure))],
        [xy.get(i).Y for i in range(len(structure))],
    )


def rotate(xs, ys, angle):
    c, s = math.cos(angle), math.sin(angle)
    return (
        [c * x - s * y for x, y in zip(xs, ys, strict=True)],
        [s * x + c * y for x, y in zip(xs, ys, strict=True)],
    )


def best_angle(xs, ys, aspect=1.0):
    """The rotation that fills a panel of the given height-to-width ratio.

    Every nest is then drawn at one scale, so the panels stay comparable, and the scale is
    the largest the panel's real shape allows rather than the largest a square would.
    """

    def fill(a):
        rx, ry = rotate(xs, ys, a)
        return min(1 / (max(rx) - min(rx)), aspect / (max(ry) - min(ry)))

    return max((math.radians(d) for d in range(0, 180, 3)), key=fill)


def pair_list(structure):
    stack, out = [], []
    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            out.append((stack.pop(), i))
    return out


def pair_set(structure):
    """The same pairs as a set, for the symmetric difference that measures two folds apart."""
    return set(pair_list(structure))


def backbone_segments(xs, ys):
    """The backbone split wherever the layout places consecutive bases far apart.

    Drawn as one polyline, such a jump is a straight line across the whole panel that reads
    as an annotation. The naive layout makes none, but the split costs nothing and keeps the
    panel correct if the layout is ever changed.
    """
    step = [math.hypot(xs[i + 1] - xs[i], ys[i + 1] - ys[i]) for i in range(len(xs) - 1)]
    median = sorted(step)[len(step) // 2]
    cuts = {i for i, v in enumerate(step) if v > 3 * median}
    segments, start = [], 0
    for i in sorted(cuts):
        segments.append((xs[start : i + 1], ys[start : i + 1]))
        start = i + 1
    segments.append((xs[start:], ys[start:]))
    return [s for s in segments if len(s[0]) > 1], len(cuts)


def draw_structure(ax, design, runs, angle):
    xs, ys = rotate(*coords(design["structure"]), angle)
    for i, j in pair_list(design["structure"]):
        ax.plot([xs[i], xs[j]], [ys[i], ys[j]], color="#e0dfdb", lw=0.22, zorder=1)
    segments, jumps = backbone_segments(xs, ys)
    for sx, sy in segments:
        ax.plot(sx, sy, color="#c2c1bd", lw=0.5, solid_capstyle="round", zorder=2)
    for _, roman, _ in reversed(CATEGORIES):
        for lo, hi in runs[roman]:
            idx = list(range(lo, min(hi, len(xs))))
            if idx:
                ax.plot(
                    [xs[i] for i in idx],
                    [ys[i] for i in idx],
                    color=PALETTE[roman],
                    lw=1.9,
                    solid_capstyle="round",
                    zorder=5,
                )
    pad = 0.04 * max(max(xs) - min(xs), max(ys) - min(ys))
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.set_aspect("equal")
    ax.set_anchor("N")
    ax.axis("off")
    return jumps


# --------------------------------------------------------------------------------- panels


def panel_catalog(ax, cat):
    n = cat["population"]
    entries = cat["entries"]
    y = list(range(len(entries)))[::-1]
    for row, e in zip(y, entries, strict=True):
        ax.barh(row, 100 * e["spellable"] / n, height=0.74, color=GRID, lw=0, zorder=1)
        ax.barh(
            row,
            100 * e["occurs"] / n,
            height=0.74,
            color=PALETTE[e["category"]],
            lw=0,
            zorder=2,
        )
    ax.set_yticks(y)
    ax.set_yticklabels([e["motif"].replace("T", "U") for e in entries], fontsize=FST)
    for tick, e in zip(ax.get_yticklabels(), entries, strict=True):
        tick.set_color(PALETTE[e["category"]])
    ax.set_ylim(-0.9, len(entries) - 0.1)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel(rf"\% of {n} proteins", fontsize=FSS, labelpad=1.5)
    ax.tick_params(axis="y", length=0, pad=1.5)
    ax.grid(axis="x", color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    # A bracket per category, labelled with the numeral only: the category names are longer
    # than the shorter groups are tall, so spelling them here would overrun into the motifs.
    start = 0
    for _, roman, _ in CATEGORIES:
        k = sum(1 for e in entries if e["category"] == roman)
        top, bot = y[start], y[start + k - 1]
        ax.annotate(
            "",
            xy=(-0.42, bot - 0.35),
            xytext=(-0.42, top + 0.35),
            xycoords=ax.get_yaxis_transform(),
            arrowprops={"arrowstyle": "-", "lw": 1.2, "color": PALETTE[roman]},
            annotation_clip=False,
        )
        ax.text(
            -0.475,
            (top + bot) / 2,
            rf"\textbf{{{roman}}}",
            transform=ax.get_yaxis_transform(),
            ha="center",
            va="center",
            fontsize=FSS + 0.5,
            color=PALETTE[roman],
        )
        start += k
    ax.text(
        0.985,
        0.012,
        r"\textsl{bar} occurs \quad \textsl{grey} admissible",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=FST,
        color=INK_3,
    )


def panel_fronts(ax, fronts, hulls, floor):
    lo = min(p["cai"] for _, pts in fronts.values() for p in pts)
    hi = max(p["cai"] for _, pts in fronts.values() for p in pts)
    for nest in NESTS:
        pts = fronts[nest][1]
        on = hulls[nest]
        colour = PLAIN if FORBIDS[nest] is None else PALETTE[FORBIDS[nest]]
        off = [i for i in range(len(pts)) if i not in on]
        ax.scatter(
            [pts[i]["cai"] for i in off],
            [pts[i]["mfe"] for i in off],
            s=1.4,
            facecolors="none",
            edgecolors=colour,
            lw=0.28,
            zorder=3,
        )
        ins = sorted(on)
        ax.plot(
            [pts[i]["cai"] for i in ins],
            [pts[i]["mfe"] for i in ins],
            color=colour,
            lw=0.7,
            zorder=4,
        )
        ax.scatter(
            [pts[i]["cai"] for i in ins],
            [pts[i]["mfe"] for i in ins],
            s=2.6,
            c=colour,
            lw=0,
            zorder=5,
        )
    ax.axvline(floor, color=INK, lw=0.5, ls=(0, (2.4, 1.6)), zorder=6)
    ax.annotate(
        rf"CAI $\geq {floor}$",
        (floor, 0.47),
        xycoords=("data", "axes fraction"),
        xytext=(2.5, 0),
        textcoords="offset points",
        va="center",
        fontsize=FST,
        color=INK,
    )
    pad = 0.02 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_xlabel(
        rf"CAI \quad \textsl{{(axis starts at {lo:.2f}, not 0)}}", fontsize=FSS, labelpad=1.5
    )
    ax.set_ylabel("MFE (kcal/mol)", fontsize=FSS, labelpad=1.5)
    ax.grid(color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    handles = [
        Line2D([], [], color=PLAIN, marker="o", ms=1.9, lw=0.7, label="on hull"),
        Line2D(
            [],
            [],
            color=PLAIN,
            marker="o",
            ms=1.9,
            lw=0,
            markerfacecolor="none",
            markeredgewidth=0.3,
            label="off hull",
        ),
    ]
    ax.legend(
        handles=handles,
        fontsize=FST,
        loc="lower right",
        handlelength=1.2,
        borderpad=0.2,
        labelspacing=0.25,
    )

    # The four fronts differ by under one percent of the MFE, so at full scale they overprint;
    # the inset is where the figure shows that they are four fronts and not one.
    inset = ax.inset_axes((0.08, 0.52, 0.42, 0.40))
    half = 0.012 * (hi - lo)
    for nest in NESTS:
        pts = fronts[nest][1]
        near = [p for p in pts if floor - half <= p["cai"] <= floor + half]
        colour = PLAIN if FORBIDS[nest] is None else PALETTE[FORBIDS[nest]]
        on = {id(pts[i]) for i in hulls[nest]}
        inset.scatter(
            [p["cai"] for p in near],
            [p["mfe"] for p in near],
            s=2.0,
            c=[colour if id(p) in on else "none" for p in near],
            edgecolors=colour,
            lw=0.25,
            zorder=3,
        )
    inset.axvline(floor, color=INK, lw=0.4, ls=(0, (2.0, 1.4)))
    inset.set_xlim(floor - half, floor + half)
    inset.tick_params(labelsize=FST - 0.6, length=1.4, pad=1.0)
    inset.set_title(rf"CAI ${floor - half:.3f}$--${floor + half:.3f}$", fontsize=FST - 0.4, pad=1.2)
    for s in ("top", "right"):
        inset.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        inset.spines[s].set_linewidth(0.4)
    return lo, hi


def panel_sweep(ax, fronts, floor, clip):
    tables = {n: suffix_min(fronts[n][1]) for n in NESTS}
    lo = max(fronts[n][1][0]["cai"] for n in NESTS)
    hi = min(min(fronts[n][1][-1]["cai"] for n in NESTS), clip)
    grid = [lo + (hi - lo) * i / 400 for i in range(401)]
    for nest in NESTS[1:]:
        xs = [c for c in grid if best_at(tables[nest], c) is not None]
        ys = [best_at(tables[nest], c) - best_at(tables["F0"], c) for c in xs]
        ax.plot(xs, ys, color=PALETTE[FORBIDS[nest]], lw=1.0, zorder=3, label=nest)
    ax.axvline(floor, color=INK, lw=0.5, ls=(0, (2.4, 1.6)), zorder=2)
    ax.set_xlim(lo, hi)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("CAI lower bound", fontsize=FSS, labelpad=1.5)
    ax.set_ylabel(r"MFE above $F_0$ (kcal/mol)", fontsize=FSS, labelpad=1.5)
    ax.grid(color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=FST, loc="upper left", handlelength=1.1, borderpad=0.2, labelspacing=0.25)
    tail = [
        best_at(tables["F3"], c) - best_at(tables["F0"], c)
        for c in (lo + (1 - lo) * i / 2000 for i in range(2001))
        if best_at(tables["F3"], c) is not None
    ]
    full = max(tail)
    ax.text(
        0.985,
        0.03,
        rf"clipped at {hi:.2f}; reaches {full:.0f} by CAI $\to 1$",
        transform=ax.transAxes,
        ha="right",
        fontsize=FST,
        color=INK_3,
    )
    return lo, hi, full


# ------------------------------------------------------------------------------ the figure


def collect(root, cell, floor, sets):
    paths = {n: settled(root, cell, n) for n in NESTS}
    fronts = {n: read_front(paths[n]) for n in NESTS}
    hulls = {n: hull_vertices(fronts[n][1]) for n in NESTS}
    chosen = {n: pick(fronts[n][1], floor) for n in NESTS}
    nests, seen = {}, {}
    for nest in NESTS:
        head, pts = fronts[nest]
        cover, runs, counts = scan(nest, chosen[nest], sets)
        digest = hashlib.sha256(chosen[nest]["sequence"].encode()).hexdigest()
        if digest in seen:
            print(
                f"warning: {nest} emits the same design as {seen[digest]} -- category "
                f"{FORBIDS[nest]} forbids nothing at CAI >= {floor}",
                file=sys.stderr,
            )
        seen[digest] = nest
        nests[nest] = {
            "path": paths[nest],
            "front_sha256": hashlib.sha256(open(paths[nest], "rb").read()).hexdigest(),
            "build": build_line(paths[nest]),
            "points": len(pts),
            "complete": bool(head.get("complete")),
            "stopped_by": head.get("stopped_by"),
            "open_regions": len(head.get("open") or []),
            "band": "pass-grid" if "/screen/" in paths[nest] else "hull-gap",
            "stop": head.get("stop"),
            # Gaps SOLVED by this round, not gaps left -- and on a resumed round only that
            # round's increment, so it is not comparable between nests. The field that counts
            # what is unresolved is `open`. An earlier reading of this as a deficit made F3,
            # whose round solved the most of any, read as the worst front in the set.
            "gaps_solved_this_round": head.get("gaps"),
            "resumed_from": (head.get("resumed") or {}).get("from"),
            # What a reader of this figure actually needs: whether a front being cut short
            # can change the design drawn. Each open box's low corner is an exact optimum at
            # its own weight, so nothing undiscovered lies below that corner's MFE; if every
            # box floor is worse than the pick, the pick is exact at the bound whatever the
            # front's `complete` field says. Three tests that could only return False were
            # here before, including one that fired on a finished walk.
            "pick_exact_at_floor": all(
                b["mfe_lo"] > chosen[nest]["mfe"] for b in (head.get("open") or [])
            ),
            "open_box_mfe_floor": min(
                (b["mfe_lo"] for b in (head.get("open") or [])), default=None
            ),
            # Absent on every f0 header, so False would be the absence of a key rather than
            # a measurement. Recorded as a tri-state so the sidecar cannot imply one.
            "complete_recorded": "complete" in head,
            "hull_recomputed": len(hulls[nest]),
            "hull_flagged": sum(1 for p in pts if p.get("on_hull")),
            "forbids": sorted(forbidden_by(nest)),
            "design": {
                "sha256": digest,
                "cai": chosen[nest]["cai"],
                "mfe": chosen[nest]["mfe"],
                "sites": counts,
            },
        }
        # What the figure colours against what the solver actually forbade. The digests in
        # the sidecar certify the pattern files as they are at render time, which is days
        # younger than the fronts -- set3-repeats.txt was last touched three days after F3
        # was written -- so they certify what the figure read, not what the search enforced.
        # The header carries the solver's own list; comparing them costs one line.
        recorded = head.get("motifs")
        if recorded is not None:
            by_cat = dict(sets)
            want = {q for roman in forbidden_by(nest) for q in by_cat[roman]}
            got = set(recorded)
            if want != got:
                raise SystemExit(
                    f"{nest}: the front forbade a different set than this figure colours.\n"
                    f"  only in the front:  {sorted(got - want)}\n"
                    f"  only in the figure: {sorted(want - got)}"
                )
        chosen[nest]["_runs"] = runs
        chosen[nest]["_cover"] = cover
    return fronts, hulls, chosen, nests


def build(root, cell, floor, clip, out_pdf, out_json):
    style()
    sets = read_patterns(root)
    cat = catalog(root)
    fronts, hulls, chosen, nests = collect(root, cell, floor, sets)

    fig = plt.figure(figsize=(7.2, 5.9))
    top = fig.add_gridspec(
        1,
        3,
        width_ratios=[1.28, 1.30, 1.05],
        left=0.145,
        right=0.985,
        top=0.965,
        bottom=0.515,
        wspace=0.40,
    )
    # The gap between the bands carries the upper panels' axis labels and the structure
    # panels' three-line headings, which together need more room than a wspace can give.
    bottom = fig.add_gridspec(1, 4, left=0.015, right=0.985, top=0.395, bottom=0.075, wspace=0.05)

    ax_a = fig.add_subplot(top[0, 0])
    panel_catalog(ax_a, cat)
    ax_b = fig.add_subplot(top[0, 1])
    b_lo, b_hi = panel_fronts(ax_b, fronts, hulls, floor)
    ax_c = fig.add_subplot(top[0, 2])
    c_lo, c_hi, c_full = panel_sweep(ax_c, fronts, floor, clip)

    angle = best_angle(*coords(chosen["F0"]["structure"]))
    jumps, struct_axes = 0, []
    for col, nest in enumerate(NESTS):
        ax = fig.add_subplot(bottom[0, col])
        struct_axes.append(ax)
        jumps += draw_structure(ax, chosen[nest], chosen[nest]["_runs"], angle)
        counts = nests[nest]["design"]["sites"]
        live = ", ".join(f"{r} {counts[r]}" for _, r, _ in CATEGORIES if counts[r])
        mark = r" $\ast$" if nests[nest]["anytime"] else ""
        ax.set_title(
            rf"{nest}, {NEST_TITLE[nest]}{mark}"
            + "\n"
            + rf"MFE ${chosen[nest]['mfe']:.1f}$"
            + "\n"
            + (live if live else r"no site, any category"),
            fontsize=FST,
            pad=1.5,
            linespacing=1.4,
            color=INK,
        )

    for letter, ax in zip("abc", (ax_a, ax_b, ax_c), strict=True):
        box = ax.get_position()
        fig.text(box.x0 - 0.055, box.y1 + 0.012, rf"\textbf{{{letter}}}", fontsize=FS + 1.2)
    for letter, ax in zip("defg", struct_axes, strict=True):
        box = ax.get_position()
        fig.text(box.x0 + 0.004, box.y1 + 0.006, rf"\textbf{{{letter}}}", fontsize=FS + 1.2)

    handles = [Patch(facecolor=PALETTE[r], label=f"category {r}") for _, r, _ in CATEGORIES]
    handles.append(Patch(facecolor=GRID, label="lattice-admissible"))
    handles.append(Line2D([], [], color=INK_3, lw=0, marker="*", ms=4.5, label="anytime front"))
    fig.legend(
        handles=handles,
        fontsize=FST,
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.002),
        handlelength=1.1,
        columnspacing=1.4,
    )

    fig.savefig(out_pdf, format="pdf")
    plt.close(fig)

    record = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "cell": cell,
        "cai_floor": floor,
        "sweep_clip": clip,
        "categories_nested": [r for _, r, _ in CATEGORIES],
        "layout": {0: "simple", 1: "naview", 4: "RNApuzzler"}[PLOT_TYPE],
        "backbone_jumps": jumps,
        "panel_b_cai": [b_lo, b_hi],
        "panel_c_cai": [c_lo, c_hi],
        "full_ladder_cost_at_cai_1": c_full,
        "catalog": cat,
        "nests": nests,
    }
    with open(out_json, "w") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)
    return record


def main():
    ap = argparse.ArgumentParser(description="Render Figure 2 for Experiment 3.")
    ap.add_argument("e2_dir")
    ap.add_argument("out_pdf")
    ap.add_argument("out_json")
    ap.add_argument("--cell", default="P47211")
    ap.add_argument("--cai-floor", type=float, default=0.85)
    ap.add_argument("--sweep-clip", type=float, default=0.95)
    args = ap.parse_args()
    rec = build(
        args.e2_dir, args.cell, args.cai_floor, args.sweep_clip, args.out_pdf, args.out_json
    )
    for nest in NESTS:
        n = rec["nests"][nest]
        state = "complete" if n["complete"] else f"open {n['open_regions']}"
        print(
            f"{nest}: {n['points']:5d} pts  hull {n['hull_recomputed']:4d} "
            f"(flagged {n['hull_flagged']:4d})  {state:>9}  "
            f"sites {n['design']['sites']}  {n['design']['sha256'][:8]}"
        )
    print(f"population {rec['catalog']['population']}, backbone jumps {rec['backbone_jumps']}")


if __name__ == "__main__":
    main()
