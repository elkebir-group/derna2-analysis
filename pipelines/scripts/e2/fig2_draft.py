"""Figure 2 draft — the motif catalog, one instance's nested fronts, and what each nest removes.

Six panels at a fixed CAI lower bound:

  (a) the blacklist, grouped into its three categories, and how many screened proteins
      carry each entry on their unconstrained front
  (b) one instance's four nested fronts; fill separates hull vertices from off-hull points
  (c) the unconstrained design's structure, bases carrying a forbidden motif colored by
      the category that forbids it
  (d) the same under category I, (e) under I+II, (f) under I+II+III

One color per category, shared by (a), (b) and the structure panels. Panel (f) carries no
category color because no entry of any category occurs in its design — that is the claim.

Reads the E2 deposit through `data/results/e2`; writes a derived JSON beside the figure so
every plotted number is inspectable without re-reading the fronts.
"""

import json
import re
import sys
from pathlib import Path

import matplotlib
import RNA
import seaborn as sns

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# One instance carries panels (b)-(f): the only middle-band protein whose four nests are all
# on disk with a settled front. F0 comes from the occurrence screen, F1-F3 from band rounds.
INSTANCE = "P05160"
NESTS = ("F0", "F1", "F2", "F3")

# Category -> (roman label, the set file that defines it). Depth F*n* forbids the first *n*.
CATEGORIES = (
    ("manufacturing", "I", "set1-manufacturing"),
    ("processing", "II", "set2-processing"),
    ("repeats", "III", "set3-repeats"),
)

# Category colors, from the shared categorical palette's validated order. F0 is uncolored:
# it forbids nothing, so no category names it.
PALETTE = {"I": "#2a78d6", "II": "#eb6834", "III": "#1baf7a"}
NEST_COLOR = {"F0": "#52514e", "F1": PALETTE["I"], "F2": PALETTE["II"], "F3": PALETTE["III"]}
# What each nest adds, for the panel titles.
NEST_TITLE = {
    "F0": "unconstrained",
    "F1": "category I",
    "F2": "category I$+$II",
    "F3": "category I$+$II$+$III",
}

INK, INK_2, GRID, PAPER = "#0b0b0b", "#52514e", "#e4e3df", "#f2f1ee"
FS, FSS = 7.0, 6.0
PREAMBLE = r"\usepackage{amsmath}\usepackage{bm}"


def read_points(path):
    """The point records of a front file, with its header."""
    lines = Path(path).read_text().splitlines()
    header = json.loads(lines[0])
    points = [json.loads(ln) for ln in lines[1:] if ln]
    return header, [p for p in points if p.get("record") == "point"]


def read_set(motifs_dir, name):
    """The patterns of one category file, comments and blanks dropped."""
    out = []
    for line in (Path(motifs_dir) / f"{name}.txt").read_text().splitlines():
        pattern = line.split("#")[0].strip()
        if pattern:
            out.append(pattern)
    return out


def occurrences(sequence, pattern):
    """Every start at which a pattern occurs. R is A or G; the design is read as DNA."""
    dna = sequence.replace("U", "T")
    rx = "(?=" + pattern.replace("R", "[AG]") + ")"
    return [m.start() for m in re.finditer(rx, dna)]


def catalog(screen_dir):
    """Per blacklist entry: how many screened proteins' unconstrained fronts spell it.

    `occurs` is the entry being spelled by at least one design of that protein's front;
    `spellable` is the lattice admitting it anywhere at all. An entry can be spellable
    everywhere and occur nowhere, which is what makes it cheap to forbid.
    """
    rows = []
    for name in ("binding.jsonl", "binding-l.jsonl", "binding-five.jsonl"):
        path = Path(screen_dir) / name
        if path.exists():
            rows += [json.loads(ln) for ln in path.read_text().splitlines() if ln]
    by_entry, proteins = {}, set()
    for row in rows:
        proteins.add(row["protein"])
        key = (row["set"], row["motif"])
        tally = by_entry.setdefault(key, {"occurs": 0, "spellable": 0, "proteins": set()})
        tally["proteins"].add(row["protein"])
        tally["occurs"] += row["occurs"] > 0
        tally["spellable"] += row["spellable"] > 0
    out = []
    for category, roman, _ in CATEGORIES:
        for (cat, motif), tally in by_entry.items():
            if cat != category:
                continue
            n = len(tally["proteins"])
            out.append(
                {
                    "category": roman,
                    "motif": motif,
                    "proteins": n,
                    "occurs": tally["occurs"],
                    "spellable": tally["spellable"],
                    "share": tally["occurs"] / n,
                }
            )
    return out, sorted(proteins)


def nest_paths(root):
    """The front file of each nest: the screen's F0 and the latest round of F1-F3.

    A cell's round 1 writes `<cell>.jsonl` and round *N*+1 writes `<cell>.rN.jsonl`, so the
    highest `.rN` is the latest artifact and the unsuffixed file is the earliest.
    """
    paths = {"F0": Path(root) / "screen" / "f0-l" / f"{INSTANCE}.jsonl"}
    for nest in NESTS[1:]:
        rounds = sorted(
            (Path(root) / "mid").glob(f"{INSTANCE}_{nest}.r*.jsonl"),
            key=lambda p: int(p.name.split(".r")[1].split(".")[0]),
        )
        paths[nest] = rounds[-1] if rounds else Path(root) / "mid" / f"{INSTANCE}_{nest}.jsonl"
    return paths


def pick_design(points, cai_floor):
    """The lowest-MFE design at or above the CAI lower bound: the panel's fixed question."""
    feasible = [p for p in points if p["cai"] >= cai_floor]
    if not feasible:
        return None
    return min(feasible, key=lambda p: p["mfe"])


def motif_positions(sequence, sets):
    """Position -> category roman, for every base covered by an occurring forbidden entry.

    A base covered by entries of two categories takes the earliest category, so the color
    names the shallowest nest that would have removed it.
    """
    colored, counts = {}, {}
    for roman, patterns in reversed(sets):
        hits = 0
        for pattern in patterns:
            for start in occurrences(sequence, pattern):
                hits += 1
                for i in range(start, start + len(pattern)):
                    colored[i] = roman
        counts[roman] = hits
    return colored, counts


def layout(structure, aspect):
    """RNApuzzler coordinates, rotated to the angle that fills the panel's aspect best."""
    import math

    RNA.cvar.rna_plot_type = 4
    xy = RNA.get_xy_coordinates(structure)
    x0 = [xy.get(i).X for i in range(len(structure))]
    y0 = [xy.get(i).Y for i in range(len(structure))]

    def rotate(a):
        c, s = math.cos(a), math.sin(a)
        return (
            [c * x - s * y for x, y in zip(x0, y0, strict=True)],
            [s * x + c * y for x, y in zip(x0, y0, strict=True)],
        )

    def fill(a):
        xs, ys = rotate(a)
        return min(aspect / (max(xs) - min(xs)), 1 / (max(ys) - min(ys)))

    return rotate(max((math.radians(d) for d in range(0, 180, 3)), key=fill))


def pairs(structure):
    stack, out = [], []
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            out.append((stack.pop(), i))
    return out


def draw_structure(fig, ax, design, colored, title):
    """The design's MFE structure: backbone, pairs, and the motif-carrying bases on top."""
    fig.canvas.draw()
    box = ax.get_position()
    xs, ys = layout(
        design["structure"], box.width * fig.get_figwidth() / (box.height * fig.get_figheight())
    )
    for i, j in pairs(design["structure"]):
        ax.plot([xs[i], xs[j]], [ys[i], ys[j]], color="#d4d3cf", lw=0.18, zorder=1)
    ax.plot(xs, ys, color="#bdbcb8", lw=0.3, zorder=2)
    plain = [i for i in range(len(xs)) if i not in colored]
    ax.scatter([xs[i] for i in plain], [ys[i] for i in plain], s=0.35, c="#9b9a96", lw=0, zorder=3)
    for roman in ("III", "II", "I"):
        hit = [i for i, r in colored.items() if r == roman]
        if hit:
            ax.scatter(
                [xs[i] for i in hit],
                [ys[i] for i in hit],
                s=4.0,
                c=PALETTE[roman],
                lw=0,
                zorder=4,
            )
    ax.set_title(title, fontsize=FS, pad=2)
    pad = 0.03 * (max(xs) - min(xs))
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.set_aspect("equal")
    ax.set_anchor("N")
    ax.axis("off")


def draw_catalog(ax, entries, compact=False):
    """Panel (a): the blacklist by category, each entry's share of screened proteins.

    `compact` is the one-row strip, where the panel is a sixth of the page wide: the entry
    labels drop a point and the axis label is shortened, which is the whole of what fits.
    """
    ypos, labels, ticks = 0.0, [], []
    names = {roman: name for name, roman, _ in CATEGORIES}
    for roman in ("I", "II", "III"):
        ax.annotate(
            r"\textbf{" + roman + "} " + names[roman],
            (0, ypos - 0.25),
            xytext=(1, 0),
            textcoords="offset points",
            fontsize=FSS,
            color=PALETTE[roman],
            va="bottom",
        )
        ypos -= 1.05
        rows = sorted(
            (e for e in entries if e["category"] == roman), key=lambda e: -e["share"]
        )
        for entry in rows:
            ax.barh(
                ypos,
                100 * entry["share"],
                height=0.72,
                color=PALETTE[roman],
                lw=0,
                zorder=3,
            )
            labels.append(r"\texttt{" + entry["motif"] + "}")
            ticks.append(ypos)
            ypos -= 1.0
        ypos -= 0.45
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=FSS - 1.4 if compact else FSS)
    ax.set_xticks([0, 50, 100] if compact else [0, 20, 40, 60, 80, 100])
    ax.set_xlim(0, 103)
    ax.set_ylim(ypos + 0.4, 0.55)
    ax.set_xlabel(
        r"instances (\%)" if compact else r"instances carrying the motif (\%)", fontsize=FS
    )
    ax.tick_params(axis="x", labelsize=FSS)
    ax.grid(axis="x", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)



def draw_front(ax, fronts, cai_floor, chosen, compact=False):
    """Panel (b): the four nested fronts; filled marks are hull vertices, open are off-hull."""
    for nest in NESTS:
        points = fronts[nest]
        color = NEST_COLOR[nest]
        off = [p for p in points if not p["on_hull"]]
        hull = sorted((p for p in points if p["on_hull"]), key=lambda p: p["cai"])
        ax.scatter(
            [p["cai"] for p in off],
            [p["mfe"] for p in off],
            s=1.6,
            facecolors="none",
            edgecolors=color,
            lw=0.25,
            zorder=2,
        )
        ax.plot(
            [p["cai"] for p in hull],
            [p["mfe"] for p in hull],
            "-",
            color=color,
            lw=0.7,
            zorder=3,
        )
        ax.scatter(
            [p["cai"] for p in hull],
            [p["mfe"] for p in hull],
            s=3.5,
            c=color,
            lw=0,
            zorder=4,
            label=f"{nest}, {NEST_TITLE[nest]}",
        )

    ax.axvline(cai_floor, color=INK, lw=0.5, ls=(0, (2.5, 1.5)), zorder=5)
    for nest in NESTS:
        design = chosen[nest]
        ax.plot(
            design["cai"],
            design["mfe"],
            "o",
            ms=4.0,
            mfc="none",
            mec=INK,
            mew=0.6,
            zorder=6,
        )
    ax.annotate(
        rf"CAI $\geq {cai_floor:.2f}$",
        (cai_floor, ax.get_ylim()[0]),
        xytext=(6, 26),
        textcoords="offset points",
        ha="left",
        va="bottom",
        fontsize=FSS,
        color=INK,
    )
    ax.set_xlabel("CAI", fontsize=FS)
    ax.set_ylabel("MFE (kcal/mol)", fontsize=FS)
    ax.tick_params(labelsize=FSS)
    ax.grid(color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if not compact:
        # In the one-row strip the four nests are named at the end of their own curves
        # instead, because a legend box of this size can only sit on top of the data.
        ax.legend(
            fontsize=FSS,
            loc="upper left",
            bbox_to_anchor=(0.0, -0.18),
            ncol=2,
            frameon=False,
            handlelength=0.9,
            borderpad=0.0,
            labelspacing=0.25,
            columnspacing=1.0,
        )

    # The four fronts separate by a few kcal/mol out of 540, so the nesting is only legible
    # zoomed. The inset is the decade of CAI around the bound the panels are taken at.
    ins = ax.inset_axes((0.10, 0.46, 0.42, 0.40))
    lo, hi = cai_floor - 0.012, cai_floor + 0.012
    ends = {}
    for nest in NESTS:
        window = sorted(
            (p for p in fronts[nest] if lo <= p["cai"] <= hi), key=lambda p: p["cai"]
        )
        if not window:
            continue
        color = NEST_COLOR[nest]
        ins.scatter(
            [p["cai"] for p in window if not p["on_hull"]],
            [p["mfe"] for p in window if not p["on_hull"]],
            s=3.0,
            facecolors="none",
            edgecolors=color,
            lw=0.3,
        )
        hull = [p for p in window if p["on_hull"]]
        ins.plot([p["cai"] for p in hull], [p["mfe"] for p in hull], "-", color=color, lw=0.7)
        ins.scatter(
            [p["cai"] for p in hull], [p["mfe"] for p in hull], s=5.0, c=color, lw=0, zorder=4
        )
        if hull:
            ends[nest] = hull[-1]["mfe"]
    ins.axvline(cai_floor, color=INK, lw=0.5, ls=(0, (2.5, 1.5)))
    for nest in NESTS:
        ins.plot(
            chosen[nest]["cai"], chosen[nest]["mfe"], "o", ms=4.0, mfc="none", mec=INK, mew=0.6
        )
    ins.set_xlim(lo, hi + 0.009 if compact else hi)
    if compact:
        # The nests separate only here, so this is where naming them works. F1 and F2 run
        # within 0.1 kcal/mol of each other, so the labels are pushed apart to stay readable
        # while each still sits at its own curve's height.
        span = max(ends.values()) - min(ends.values())
        gap = max(0.09 * span, 0.6)
        placed = {}
        for nest in sorted(ends, key=ends.get):
            y = ends[nest]
            if placed:
                y = max(y, max(placed.values()) + gap)
            placed[nest] = y
        for nest, y in placed.items():
            ins.annotate(
                nest,
                (hi + 0.0015, y),
                fontsize=FSS - 0.5,
                color=NEST_COLOR[nest],
                va="center",
                ha="left",
            )
    ins.tick_params(labelsize=FSS - 1.2, pad=1)
    ins.grid(color=GRID, lw=0.4)
    ins.set_axisbelow(True)
    for side in ("top", "right"):
        ins.spines[side].set_visible(False)


def style():
    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "text.latex.preamble": PREAMBLE,
            "font.size": FS,
            "axes.titlesize": FS,
            "axes.labelsize": FS,
            "axes.labelcolor": INK,
            "axes.edgecolor": INK_2,
            "xtick.labelsize": FSS,
            "ytick.labelsize": FSS,
            "xtick.color": INK_2,
            "ytick.color": INK_2,
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 2.0,
            "ytick.major.size": 2.0,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "savefig.bbox": "tight",
        }
    )


# Panel order is fixed; only the shape of the grid changes. The one-row strip gives (a) and
# (b) more width than the structures, which are square and lose nothing by being narrow.
LAYOUTS = {
    "grid": {
        "figsize": (7.2, 5.6),
        "shape": (2, 3),
        "ratios": None,
        "wspace": 0.26,
        "hspace": 0.46,
    },
    "row": {
        "figsize": (14.4, 2.9),
        "shape": (1, 6),
        "ratios": [1.5, 1.35, 1, 1, 1, 1],
        "wspace": 0.34,
        "hspace": 0.0,
    },
}


def main(root, motifs_dir, out_pdf, out_json, cai_floor=0.85, layout_name="grid"):
    root = Path(root)
    entries, proteins = catalog(root / "screen")
    sets = [(roman, read_set(motifs_dir, name)) for _, roman, name in CATEGORIES]

    paths = nest_paths(root)
    fronts, headers = {}, {}
    for nest in NESTS:
        headers[nest], fronts[nest] = read_points(paths[nest])

    chosen = {nest: pick_design(fronts[nest], cai_floor) for nest in NESTS}
    missing = [n for n, d in chosen.items() if d is None]
    if missing:
        raise SystemExit(f"no design at CAI >= {cai_floor} for {', '.join(missing)}")

    colored, counts = {}, {}
    for nest in NESTS:
        # Which categories can still occur: the ones this nest does not forbid.
        forbidden = {"F0": set(), "F1": {"I"}, "F2": {"I", "II"}, "F3": {"I", "II", "III"}}[nest]
        visible = [(r, p) for r, p in sets if r not in forbidden]
        colored[nest], counts[nest] = motif_positions(chosen[nest]["sequence"], visible)

    spec = LAYOUTS[layout_name]
    compact = layout_name == "row"

    style()
    fig = plt.figure(figsize=spec["figsize"])
    rows, cols = spec["shape"]
    gs = fig.add_gridspec(
        rows, cols, width_ratios=spec["ratios"], wspace=spec["wspace"], hspace=spec["hspace"]
    )
    flat = [gs[i // cols, i % cols] for i in range(rows * cols)]
    ax_a = fig.add_subplot(flat[0])
    ax_b = fig.add_subplot(flat[1])
    draw_catalog(ax_a, entries, compact=compact)
    draw_front(ax_b, fronts, cai_floor, chosen, compact=compact)

    cells = dict(zip(NESTS, flat[2:], strict=True))
    for nest in NESTS:
        ax = fig.add_subplot(cells[nest])
        hits = sum(counts[nest].values())
        mfe = chosen[nest]["mfe"]
        tail = f"{hits} occurrences" if hits else "none occurs"
        title = f"{nest}, {NEST_TITLE[nest]}" + "\n" + rf"MFE ${mfe:.1f}$, {tail}"
        draw_structure(fig, ax, chosen[nest], colored[nest], title)

    for ax, tag in ((ax_a, "a"), (ax_b, "b")):
        ax.annotate(
            r"\textbf{" + tag + "}",
            (0, 1),
            xycoords="axes fraction",
            xytext=(-30 if compact else -26, 8),
            textcoords="offset points",
            fontsize=FS + 1.5,
            color=INK,
        )
    for nest, tag in zip(NESTS, "cdef", strict=True):
        ax = fig.axes[2 + NESTS.index(nest)]
        ax.annotate(
            r"\textbf{" + tag + "}",
            (0, 1),
            xycoords="axes fraction",
            xytext=(-6, 10),
            textcoords="offset points",
            fontsize=FS + 1.5,
            color=INK,
        )

    out_pdf = Path(out_pdf)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.with_suffix(".png"), dpi=300)
    plt.close(fig)

    Path(out_json).write_text(
        json.dumps(
            {
                "instance": INSTANCE,
                "cai_floor": cai_floor,
                "screened_proteins": len(proteins),
                "catalog": entries,
                "nests": {
                    nest: {
                        "front": str(paths[nest]),
                        "points": len(fronts[nest]),
                        "hull": sum(p["on_hull"] for p in fronts[nest]),
                        "complete": headers[nest].get("complete"),
                        "version": headers[nest].get("version"),
                        "design_mfe": chosen[nest]["mfe"],
                        "design_cai": chosen[nest]["cai"],
                        "motif_occurrences": counts[nest],
                    }
                    for nest in NESTS
                },
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main(
        *sys.argv[1:5],
        cai_floor=float(sys.argv[5]) if len(sys.argv) > 5 else 0.85,
        layout_name=sys.argv[6] if len(sys.argv) > 6 else "grid",
    )
