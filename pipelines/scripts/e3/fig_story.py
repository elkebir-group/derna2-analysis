"""The main E3 figure: what a UTR-aware design practice delivers, on one practical UTR pair and over all five.

    python fig_story.py BUNDLE E2_RESULTS HULLCHECK CAND CELL CAI_FLOOR OUT.pdf

BUNDLE      the E3 figure bundle's paper-e4/ (results_root/e3)
E2_RESULTS  E2's results tree, which holds the unconstrained fronts concatenation starts from
HULLCHECK   `derna check` of each cell's concatenation hull vertices, <cell>.out.jsonl, which carry the
            CDS-UTR pair positions the census rows lack
CAND        `derna check` of the two designs panel (e) draws besides DERNA2's, <cell>.conc.out.jsonl and
            <cell>.mrf.out.jsonl
CELL        the instance of (d) and (e), a protein on a UTR pair, and CAI_FLOOR the cutoff of its designs

(a) where the CDS-UTR pairs of concatenation's non-modular designs land, per UTR pair
(b) proteins with a modular design, and modular designs over all proteins, per method
(c) the windows of DERNA2's search by final state
(d) the fronts of the instance, with a zoom at the CAI cutoff
(e) the best design at the cutoff of concatenation, mRNAfold and DERNA2 for the instance: the codons that
    differ from DERNA2's, the CDS-UTR pairs, and the whole-molecule folds with the Kozak site

"DERNA2" is the output with window enumeration: the masked front's modular designs, which it
starts from, and the modular designs it adds. Element annotation and the 1-based `violations` reading follow
drafts/e3-plate (annotate_utr.py, plot_pairs.py).
"""
import glob
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict

import importlib.util
from pathlib import Path

import matplotlib
from matplotlib.lines import Line2D

matplotlib.use("Agg")
import matplotlib.legend_handler  # noqa: E402
import matplotlib.transforms  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e4common import refold_points, relocate, take_refold  # noqa: E402

PAIRS = ["hbb", "mrna1273", "hba", "curevac", "bnt"]
PAIR_HUES = dict(zip(PAIRS, ["#8f5f00", "#9f9f00", "#008d51", "#00c2ac", "#008ec8"]))
METHODS = ["hull", "mrnafold", "window"]
NAMES = {"hull": "concatenation", "mrnafold": "mRNAfold",
         "masked": "DERNA2 masked", "window": "DERNA2"}
# The paper's colours (pipelines/scripts/e1/fig4_overview.py): DERNA2 purple, mRNAfold orange,
# the UTR-unaware tools green; the masked front takes the lighter purple.
COL = {"hull": "#2e9e5b", "mrnafold": "#e8772e", "masked": "#b7a6e8", "window": "#6a4bc4"}
EL_NONE, EL_PYR, EL_KOZAK, EL_OTHER, EL_BOTH_EDGE = "#d8d4dd", "#1864ab", "#c2a000", "#868e96", "#1864ab"
INK, MUTED = "#3a3a3a", "#8a8a8a"
ELCLS = ["kozak", "pyr", "polyA", "are"]
CAI_FIG = 0.8


def rd(path):
    for line in open(path):
        try:
            yield json.loads(line)
        except ValueError:
            continue


def nondominated(pts):
    pts = sorted(pts, key=lambda x: (-x["cai"], x["mfe"]))
    out, best = [], math.inf
    for x in pts:
        if x["mfe"] < best - 1e-9:
            out.append(x)
            best = x["mfe"]
    return out


# ---- element annotation, as drafts/e3-plate/annotate_utr.py ----
def find_pyr(s, w=14, k=12):
    hits = [(i, i + w) for i in range(len(s) - w + 1) if sum(c in "CU" for c in s[i:i + w]) >= k]
    out = []
    for a, b in hits:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def element_masks(u5, u3):
    m5, m3 = [set() for _ in u5], [set() for _ in u3]
    for k in range(max(0, len(u5) - 6), len(u5)):
        m5[k].add("kozak")
    for m in re.finditer(r"A[AU]UAAA", u3):
        for k in range(m.start(), m.end()):
            m3[k].add("polyA")
    for t, s in ((m5, u5), (m3, u3)):
        for a, b in find_pyr(s):
            for k in range(a, b):
                t[k].add("pyr")
        for m in re.finditer(r"(?=(AUUUA))", s):
            for k in range(m.start(), min(m.start() + 5, len(t))):
                t[k].add("are")
    return m5, m3


def load(bundle, e2, hullcheck):
    cells = sorted({os.path.basename(p).split(".")[0] for p in glob.glob(f"{bundle}/seeds/*.check.*.jsonl")})
    head, masked = {}, {}
    for c in cells:
        p = f"{bundle}/seeds/{c}.check.freestop.ec083208.jsonl"
        if os.path.exists(p):
            rows = list(rd(p))
            head[c] = rows[0]
            masked[c] = [x for x in rows[1:] if x.get("record") == "point"]
        else:
            # Xinyu's note (ten_overlay.py): the non-dominated subset of the pooled fixed-stop
            # check equals the free-stop front on these cells.
            rows = list(rd(f"{bundle}/seeds/{c}.check.pooled-fixed.jsonl"))
            head[c] = rows[0]
            pts = [x for x in rows[1:] if "mfe" in x and "cai" in x]
            # the non-dominated subset is the free-stop front, which holds no modular design on
            # these cells; the modular designs of the pool are kept beside it
            front = nondominated(pts)
            seen = {(round(x["cai"], 7), round(x["mfe"], 2)) for x in front}
            masked[c] = front + [x for x in pts if x.get("gap") == 0
                                 and (round(x["cai"], 7), round(x["mfe"], 2)) not in seen]
    window, a4 = {}, {}
    # in the plate's order: a cell in two trees takes the later one's front
    for d in ("out-shed120-dd0b28fd", "out-freestop-ec083208", "out-resume-a4", "out-seedless-6e8d2e53"):
        for p in sorted(glob.glob(f"{bundle}/{d}/cells/*/a4.jsonl")):
            c = p.split("/")[-2]
            rows = list(rd(p))
            a4[c] = rows[0]
            window[c] = [x for x in rows[1:] if x.get("record") == "point" and x.get("gap") == 0]
    # Modularity of the hull vertices from the concatenation census, which covers every cell;
    # the pair positions from the refold, where it has finished. Their gaps agree exactly.
    hull, refolded = {}, set()
    for c in cells:
        p = next(q for q in (f"{bundle}/{t}/cells/{c}/a0.jsonl" for t in ("a0", "a0-all")) if os.path.exists(q))
        rows = list(rd(p))
        src = relocate(rows[0]["source"], bundle, e2)
        hv = {i for i, x in enumerate(r for r in rd(src) if r.get("record") == "point") if x.get("on_hull")}
        hull[c] = [{"index": x["index"], "cai": x["cai"], "mfe": x["constrained_kcal_mol"],
                    "gap": x["gap_kcal_mol"]} for x in rows[1:]
                   if x.get("record") == "point" and x["index"] in hv]
        if os.path.exists(f"{hullcheck}/{c}.out.jsonl"):
            v = refold_points(hullcheck, c)
            hull[c] = [take_refold(x, v[x["index"]]) for x in hull[c]]
            refolded.add(c)
    des, mf = {}, defaultdict(list)
    for p in glob.glob(f"{bundle}/a3soft/designs.*.jsonl"):
        for x in rd(p):
            if x["cell"] in head:
                des[(x["cell"], x["index"])] = x
    seen = set()
    for p in sorted(glob.glob(f"{bundle}/a3soft/census.*.jsonl")):
        for x in rd(p):
            k = (x.get("cell"), x.get("index"))
            if k in des and k not in seen and "gap_kcal_mol" in x and "guarded_kcal_mol" in x:
                seen.add(k)
                mf[k[0]].append({"cai": des[k]["cai"], "mfe": x["guarded_kcal_mol"], "gap": x["gap_kcal_mol"],
                                 "free": x["free_kcal_mol"]})
    return cells, head, masked, window, a4, hull, mf, refolded


def modular_sets(c, masked, window, hull, mf):
    key = lambda x: (round(x["cai"], 7), round(x["mfe"], 2))  # noqa: E731
    m = {key(x): x for x in masked[c] if x["gap"] == 0}
    w = {**{key(x): x for x in window.get(c, [])}, **m}
    return {"hull": [x for x in hull[c] if x["gap"] == 0], "mrnafold": [x for x in mf[c] if x["gap"] == 0],
            "masked": list(m.values()), "window": list(w.values())}


def element_census(cells, head, hull):
    out = {p: Counter() for p in PAIRS}
    share = {p: [0, 0] for p in PAIRS}
    masks = {}
    for c in cells:
        pair = c.split("__")[1]
        h = head[c]
        u5, u3 = h["utr5"]["sequence"], h["utr3"]["sequence"]
        if pair not in masks:
            masks[pair] = element_masks(u5, u3)
        m5, m3 = masks[pair]
        l5 = len(u5)
        cend = l5 + 3 * (len(h["protein"]["sequence"]) + 1)
        for x in hull[c]:
            share[pair][0] += 1
            if x["gap"] == 0:
                continue
            share[pair][1] += 1
            hit = set()
            for i, j in x.get("violations") or []:
                a, b = min(i, j) - 1, max(i, j) - 1  # 1-based, as drafts/e3-plate verified
                if a < l5 <= b < cend and m5[a]:
                    hit |= m5[a]
                elif l5 <= a < cend <= b and 0 <= b - cend < len(m3) and m3[b - cend]:
                    hit |= m3[b - cend]
            out[pair]["+".join(k for k in ELCLS if k in hit) or "none"] += 1
    return out, share


DISPLAY = {"hbb": "HBB", "mrna1273": "mRNA-1273", "hba": "HBA", "curevac": "CureVac", "bnt": "BNT"}
REGION = {"utr5": "#56b4e9", "cds": "#3a3a3a", "utr3": "#9a9a9a"}
PAIR_RED = "#d62728"
VSUB = {"concatenation": "cat", "mRNAfold": "mask", "DERNA2": "mod"}   # v_cat, v_mask, v_mod
MODC = "#555555"       # modular: the (c) key's grey
# The instance of (d) and (e): the cell and the CAI cutoff at which each method's best design is drawn; main() sets them
CELL = THRESHOLD = None
CODE = dict(zip((a + b + c for a in "UCAG" for b in "UCAG" for c in "UCAG"),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
BACK, PAIRC = "#c2c1bd", "#e0dfdb"      # Figure 5(d): a pale backbone and pale chords
SHAREDC, DIFFC, OUTLINE = "#d9d8d4", "#2f2d29", "#8e8a82"
# Panel (a): no blue, which is the 5' UTR's colour in (e)
EL = {"none": "#d4d4d4", "pyr": "#1b9e8f", "kozak": "#c2a000", "kozak+pyr": "#8a6a1f", "rest": "#a0457a"}


def style():
    """The paper's figure style, which the E1 overview defines: Times through usetex."""
    path = Path(__file__).resolve().parents[1] / "e1" / "fig4_overview.py"
    spec = importlib.util.spec_from_file_location("fig4", path)
    f4 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(f4)
    f4.set_style()
    return f4.FS, f4.FSS


def panel_elements(ax, combo, share):
    """One bar per pair over all of its concatenation designs: modular in white, the rest by the
    element the CDS-UTR pairs touch."""
    cats = [("none", "other"), ("pyr", "PRE only"), ("kozak", "Kozak only"),
            ("kozak+pyr", "Kozak and PRE"), ("rest", "poly(A) or ARE")]
    xs = np.arange(len(PAIRS))
    nonmod = np.array([100 * share[p][1] / share[p][0] for p in PAIRS])
    named = {"none", "pyr", "kozak", "kozak+pyr"}
    ax.bar(xs, 100 - nonmod, width=0.7, color=MODC, edgecolor="black", lw=0.5, zorder=3)
    handles = [matplotlib.patches.Patch(fc=MODC, ec="black", lw=0.5, label="modular")]
    for x, y in zip(xs, 100 - nonmod):
        if y >= 9:
            ax.text(x, y / 2, f"{y:.0f}", ha="center", va="center", fontsize=5.5, color="white", zorder=5)
    bottom = 100 - nonmod
    for key, lab in cats:
        tot = np.array([sum(combo[p].values()) for p in PAIRS], float)
        if key == "rest":
            v = np.array([sum(n for k, n in combo[p].items() if k not in named) for p in PAIRS])
        else:
            v = np.array([combo[p].get(key, 0) for p in PAIRS])
        h = nonmod * v / tot
        ax.bar(xs, h, bottom=bottom, width=0.7, color=EL[key], lw=0, zorder=3)
        ax.bar(xs, h, bottom=bottom, width=0.7, fill=False, edgecolor="black", lw=0.4, zorder=6)
        handles.append(matplotlib.patches.Patch(fc=EL[key], lw=0, label=lab))
        for x, y0, y in zip(xs, bottom, h):
            if y >= 9:
                ax.text(x, y0 + y / 2, f"{y:.0f}", ha="center", va="center", fontsize=5.5,
                        color="white" if key in ("pyr", "kozak+pyr", "rest") else INK, zorder=7)
        bottom = bottom + h
    ax.bar(xs, 100, width=0.7, fill=False, edgecolor="black", lw=0.6, zorder=8)
    ax.set_xticks(xs, [DISPLAY[p] for p in PAIRS], rotation=30, ha="right", rotation_mode="anchor")
    ax.set_ylim(0, 100)
    ax.set_ylabel("concatenation designs (\\%)")
    return handles


def pareto_split(cells, masked, window, sets):
    """Per method and UTR pair: (modular designs, certainly Pareto optimal, dominated).

    A modular design is certainly Pareto optimal when no point of the exact masked front beats it,
    because a modular design has no cheaper rival that the masked front does not already contain;
    it is dominated when a modular design that DERNA2 found beats it. A DERNA2 design is on the
    masked front or certified, so it is certainly Pareto optimal. The rest is undecided: no found
    design beats it and no proof says that none exists."""
    ct, et = 1e-9, 0.05

    def beaten(p, qs):
        return any(q["cai"] >= p["cai"] - ct and q["mfe"] <= p["mfe"] + et
                   and (q["cai"] > p["cai"] + ct or q["mfe"] < p["mfe"] - et) for q in qs)

    out = {m: {p: [0, 0, 0] for p in PAIRS} for m in METHODS}
    out["cell"] = {}                      # (method, cell) -> [modular designs, dominated]
    for c in cells:
        p = c.split("__")[1]
        front = [{"cai": x["cai"], "mfe": x["mfe"]} for x in masked[c]]
        known = [{"cai": x["cai"], "mfe": x["mfe"]} for x in sets[c]["window"]]
        cert = {(round(x["cai"], 7), round(x["mfe"], 2)) for x in window.get(c, [])}
        for m in METHODS:
            out["cell"][(m, c)] = [0, 0]
            for x in sets[c][m]:
                r = out[m][p]
                r[0] += 1
                out["cell"][(m, c)][0] += 1
                if (m == "window" and (round(x["cai"], 7), round(x["mfe"], 2)) in cert) or not beaten(x, front):
                    r[1] += 1
                elif beaten(x, known):
                    r[2] += 1
                    out["cell"][(m, c)][1] += 1
    return out


def panel_counts(axs, cells, sets, split):
    """Proteins with a modular design, and the modular designs, per UTR pair. In the lower panel
    each bar is drawn from the baseline twice: pale for all of a method's modular designs and solid
    for those that no modular design that DERNA2 found dominates."""
    ax_prot, ax_des = axs
    xs = np.arange(len(PAIRS))
    w = 0.26
    for i, m in enumerate(METHODS):
        prot = [sum(bool(sets[c][m]) for c in cells if c.endswith("__" + p)) for p in PAIRS]
        # proteins with a modular design no known modular design dominates
        prot_keep = [sum(split["cell"][(m, c)][0] - split["cell"][(m, c)][1] > 0
                         for c in cells if c.endswith("__" + p)) for p in PAIRS]
        tot = [split[m][p][0] for p in PAIRS]
        keep = [split[m][p][0] - split[m][p][2] for p in PAIRS]
        x = xs + (i - 1) * w
        ax_prot.bar(x, prot, width=w * 0.9, color=COL[m], alpha=0.35, lw=0, zorder=3)
        ax_prot.bar(x, prot_keep, width=w * 0.9, color=COL[m], lw=0, zorder=4)
        ax_des.bar(x, [max(t, 0.6) for t in tot], width=w * 0.9, color=COL[m], alpha=0.35, lw=0, zorder=3)
        ax_des.bar(x, [max(k, 0.6) for k in keep], width=w * 0.9, color=COL[m], lw=0, zorder=4)
        for xx, p_, t in zip(x, prot, tot):
            if p_ <= 3:
                ax_prot.text(xx, p_ + 0.4, f"{p_}", ha="center", va="bottom", fontsize=5.6, color=INK)
            if t == 0:
                ax_des.text(xx, 1.25, "0", ha="center", va="bottom", fontsize=5.6, color=INK)
    n_prot = len({c.split("__")[0] for c in cells})
    ax_prot.axhline(n_prot, color=MUTED, lw=0.4, ls=(0, (3, 2)), zorder=1)
    ax_prot.set_ylim(0, n_prot)
    ax_prot.set_yticks([0, 10, n_prot])
    ax_prot.set_ylabel("proteins")
    ax_prot.set_xticks(xs, [])
    ax_des.set_yscale("log")
    ax_des.set_ylim(0.6, 3e4)
    ax_des.set_ylabel("designs")
    ax_des.set_xticks(xs, [DISPLAY[p] for p in PAIRS], rotation=30, ha="right", rotation_mode="anchor")
    return [tuple(matplotlib.patches.Patch(fc=COL[m], lw=0) for m in METHODS),
            tuple(matplotlib.patches.Patch(fc=COL[m], alpha=0.35, lw=0) for m in METHODS)]


def concat_front(bundle, e2, c):
    """Every design of concatenation's front on one cell, with its hull flag."""
    p = next(q for q in (f"{bundle}/{t}/cells/{c}/a0.jsonl" for t in ("a0", "a0-all")) if os.path.exists(q))
    rows = list(rd(p))
    src = relocate(rows[0]["source"], bundle, e2)
    hv = {i for i, x in enumerate(r for r in rd(src) if r.get("record") == "point") if x.get("on_hull")}
    return [{"cai": x["cai"], "mfe": x["constrained_kcal_mol"], "gap": x["gap_kcal_mol"],
             "free": x["free_kcal_mol"], "hull": x["index"] in hv} for x in rows[1:] if x.get("record") == "point"]



def panel_front(ax, c, window, a4, mf, concat, sets, picks=None):
    """The designs of the three methods on the instance, with an inset at the CAI cutoff: the hull
    vertices of concatenation, the designs of mRNAfold and the modular designs of DERNA2. Colour is
    the method and square or circle is modular or not. The y axis is the MFE of the refolded mRNA."""
    lo, hi = THRESHOLD - 0.003, THRESHOLD + 0.008
    pts = {"hull": [{"cai": p["cai"], "mfe": p["free"], "gap": p["gap"]} for p in concat if p["hull"]],
           "mrnafold": [{"cai": p["cai"], "mfe": p["free"], "gap": p["gap"]} for p in mf[c]],
           "window": [{"cai": x["cai"], "mfe": x.get("free_mfe", x["mfe"]), "gap": 0} for x in sets[c]["window"]]}
    zs = {"hull": 3, "mrnafold": 6, "window": 5}

    def draw(a, size, window_only=False, ring=1.3):
        # modular designs are filled squares, non-modular designs open circles
        for m in ("hull", "mrnafold", "window"):
            for mod in (False, True):
                sel = [p for p in pts[m] if (p["gap"] == 0) == mod and (not window_only or lo <= p["cai"] <= hi)]
                if not sel:
                    continue
                # mRNAfold's modular designs are points of DERNA2's set too: draw them smaller
                # and on top, so the purple square shows round them
                sz = size * (0.45 if (m == "mrnafold" and mod) else 1.0)
                kw = dict(color=COL[m], lw=0) if mod else dict(facecolors="none", edgecolors=COL[m],
                                                               lw=0.3 if size < 3 else 0.7)
                a.scatter([p["cai"] for p in sel], [p["mfe"] for p in sel], s=sz if mod else sz * ring,
                          marker="s" if mod else "o", zorder=zs[m] + (0 if mod else -3), **kw)

    draw(ax, 2.6, ring=3.2)
    ax.set_xlabel("CAI")
    ax.set_ylabel("MFE (kcal/mol)")
    ins = ax.inset_axes([0.07, 0.37, 0.50, 0.53])
    draw(ins, 9, window_only=True)
    where = {"concatenation": (THRESHOLD - 0.0014, -587.6), "DERNA2": (THRESHOLD + 0.0038, -593.6),
             "mRNAfold": (THRESHOLD + 0.0038, -596.0)}
    for nm, x in (picks or {}).items():
        ins.scatter([x["cai"]], [x["free_mfe"]], s=30, facecolors="none", edgecolors="black", lw=0.7, zorder=7)
        ins.annotate(f"$\\mathbf{{v}}_{{\\mathrm{{{VSUB[nm]}}}}}$", (x["cai"], x["free_mfe"]), xytext=where[nm],
                     textcoords="data", ha="center", va="center", fontsize=6.3, zorder=8,
                     arrowprops=dict(arrowstyle="-", lw=0.4, color="black", shrinkA=2, shrinkB=3.5))
    ins.axvline(THRESHOLD, color=INK, lw=0.4, ls=(0, (3, 2)), zorder=1)
    sel = [p for m in pts for p in pts[m] if lo <= p["cai"] <= hi]
    ys = [p["mfe"] for p in sel]
    ins.set_xlim(lo, hi)
    ins.set_ylim(min(ys) - 0.6, max(ys) + 0.6)
    ins.yaxis.tick_right()
    ins.xaxis.tick_top()
    ins.tick_params(labelsize=5.5, length=1.5, pad=1)
    ins.set_facecolor("white")
    ins.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    for s_ in ins.spines.values():
        s_.set_visible(True)
        s_.set_linewidth(0.5)
    ax.indicate_inset_zoom(ins, edgecolor="black", lw=0.5, alpha=1.0)
    g = "#555555"
    ax.legend(handles=[Line2D([], [], ls="", marker="s", ms=3.0, color=g, label="modular"),
                       Line2D([], [], ls="", marker="o", ms=3.0, mfc="none", mec=g, mew=0.7,
                              label="non-modular")],
              loc="lower right", frameon=False, fontsize=5.4, handletextpad=0.2, borderaxespad=0.2,
              labelspacing=0.2)



# ---- panel (c): the windows of DERNA2's search ----
def window_states(h):
    """Walked windows of a cell by final state: enumerated, partial or unstarted."""
    INF = 1e9
    iv = lambda s: (-INF if s["cai_log_lo"] is None else s["cai_log_lo"], 0.0 if s["cai_log_hi"] is None else s["cai_log_hi"])  # noqa: E731
    walked, out = {}, Counter()
    for s in h["stretches"]:
        if s.get("walked"):
            k = iv(s)
            if k not in walked or s.get("round", 0) >= walked[k].get("round", 0):
                walked[k] = s
    out["enum1"] = sum(1 for s in walked.values() if s["complete"] and s["points"] > 0)
    out["enum0"] = sum(1 for s in walked.values() if s["complete"] and s["points"] == 0)
    started = [iv(s) for s in h["stretches"] if s.get("walked")]
    for o in h["open"]:
        lo, hi = (-INF if o["cai_log_lo"] is None else o["cai_log_lo"]), o["cai_log_hi"] or 0.0
        out["part" if any(min(hi, w[1]) - max(lo, w[0]) > 1e-9 for w in started) else "none"] += 1
    return out


def panel_d_windows(ax, cells, a4):
    xs = np.arange(len(PAIRS))
    st = {p: sum((window_states(a4[c]) for c in cells if c.endswith("__" + p)), Counter()) for p in PAIRS}
    bottom = np.zeros(len(PAIRS))
    handles = []
    for key, lab, col in (("enum1", "full", "#6a4bc4"), ("enum0", "empty", "#b7a6e8"),
                          ("part", "partial", "#d0d0d0"), ("none", "unstarted", "white")):
        v = np.array([st[p][key] for p in PAIRS], float)
        ax.bar(xs, v, bottom=bottom, width=0.7, color=col, edgecolor="black", lw=0.4, zorder=3)
        bottom += v
        handles.append(matplotlib.patches.Patch(fc=col, ec="black", lw=0.4, label=lab))
    ax.set_xticks(xs, [DISPLAY[p] for p in PAIRS], rotation=30, ha="right", rotation_mode="anchor")
    ax.set_ylabel("windows")
    return handles




def pairs_of(st):
    stack, out = [], []
    for k, ch in enumerate(st):
        if ch == "(":
            stack.append(k)
        elif ch == ")" and stack:
            out.append((stack.pop(), k))
    return out


def best_designs(bundle, e2, ex, c=None):
    """The best (minimum MFE) design at CAI >= THRESHOLD of concatenation, mRNAfold and DERNA2.

    Concatenation is ranked by the MFE of its CDS alone, the objective of the methods it stands
    for; mRNAfold by its MFE without CDS-UTR pairs, and so is DERNA2. Each design is shown with
    the MFE structure of its whole mRNA (`derna check`, in `ex`); `mfe_nopair` is the energy that
    panel (c) plots.
    """
    c = c or CELL
    seedp = f"{bundle}/seeds/{c}.check.freestop.ec083208.jsonl"
    if not os.path.exists(seedp):
        seedp = f"{bundle}/seeds/{c}.check.pooled-fixed.jsonl"
    seed = list(rd(seedp))
    a4p = glob.glob(f"{bundle}/*/cells/{c}/a4.jsonl")
    cert = [x for x in rd(a4p[0]) if x.get("record") == "point" and x.get("gap") == 0] if a4p else []
    d = min((x for x in cert + [y for y in seed[1:] if "mfe" in y and y.get("gap") == 0]
             if x["cai"] >= THRESHOLD), key=lambda x: x["mfe"])
    d["mfe_nopair"] = d["mfe"]
    a0 = next(q for q in (f"{bundle}/{t}/cells/{c}/a0.jsonl" for t in ("a0", "a0-all")) if os.path.exists(q))
    rows = list(rd(a0))
    src = relocate(rows[0]["source"], bundle, e2)
    hv = {i for i, x in enumerate(r for r in rd(src) if r.get("record") == "point") if x.get("on_hull")}
    conc = min((x for x in rows if x.get("record") == "point" and x["index"] in hv and x["cai"] >= THRESHOLD),
               key=lambda x: x["f0_mfe"])
    folded = [x for x in rd(f"{ex}/{c}.conc.out.jsonl") if x.get("record") == "point"][0]
    assert abs(folded["cai"] - conc["cai"]) < 1e-9
    folded["mfe_nopair"] = conc["constrained_kcal_mol"]
    des = {}
    for p in glob.glob(f"{bundle}/a3soft/designs.*.jsonl"):
        for x in rd(p):
            if x["cell"] == c:
                des[x["index"]] = x
    cen = {}
    for p in sorted(glob.glob(f"{bundle}/a3soft/census.*.jsonl")):
        for x in rd(p):
            k = x.get("index")
            if x.get("cell") == c and k in des and k not in cen and "gap_kcal_mol" in x and "guarded_kcal_mol" in x:
                cen[k] = x
    k = min((i for i in cen if des[i]["cai"] >= THRESHOLD), key=lambda i: cen[i]["guarded_kcal_mol"])
    mrf = [x for x in rd(f"{ex}/{c}.mrf.out.jsonl") if x.get("record") == "point"][0]
    assert abs(mrf["cai"] - des[k]["cai"]) < 1e-9
    mrf["mfe_nopair"] = cen[k]["guarded_kcal_mol"]
    h = seed[0]
    return {"concatenation": folded, "mRNAfold": mrf, "DERNA2": d}, \
        len(h["utr5"]["sequence"]), len(h["utr5"]["sequence"]) + len(h["protein"]["sequence"]) * 3 + 3


def differs(x, ref, l5, cend):
    """Codon indices of x's CDS that are not the codon of ref."""
    return {(c0 - l5) // 3 for c0 in range(l5, cend, 3) if x["sequence"][c0:c0 + 3] != ref["sequence"][c0:c0 + 3]}


def codon_strip(ax, x, ref, l5, cend, lab=None):
    """One cell per nucleotide of the mRNA, as Figure 5(d): the UTRs at either end in their
    colours and every codon of the CDS light grey where it is the codon of DERNA2's design and
    dark where it differs. A red tick marks each end of a CDS-UTR base pair."""
    n = len(x["sequence"])
    dif = differs(x, ref, l5, cend)
    row = []
    for k in range(n):
        if k < l5:
            row.append(matplotlib.colors.to_rgb(REGION["utr5"]))
        elif k >= cend:
            row.append(matplotlib.colors.to_rgb(REGION["utr3"]))
        else:
            row.append(matplotlib.colors.to_rgb(DIFFC if (k - l5) // 3 in dif else SHAREDC))
    ax.imshow(np.array([row]), extent=(0, n, 0, 1), aspect="auto", interpolation="nearest")
    # each CDS-UTR contact as an arc over the strip, taller the longer the contact
    t = np.linspace(0, np.pi, 60)
    for i, j in x.get("violations") or []:
        a, b = min(i, j) - 0.5, max(i, j) - 0.5
        h = 0.3 + 0.7 * min(1.0, (b - a) / 400)
        ax.plot((a + b) / 2 - (b - a) / 2 * np.cos(t), 1.0 + h * np.sin(t), color=PAIR_RED, lw=0.7,
                transform=ax.get_xaxis_transform(), clip_on=False, zorder=3)
    for mid, u in ((l5 / 2, "$\\mathbf{u}_5$"), ((cend + n) / 2, "$\\mathbf{u}_3$")):
        ax.text(mid, 0.5, u, ha="center", va="center", fontsize=6.3, color="black", zorder=4)
    if lab:
        pos = np.array([l5 + 3 * c + 1.5 for c in dif])
        half = 40                                  # nucleotides either side of the label
        centre = (l5 + cend) / 2
        for off in sorted(range(-(cend - l5) // 2 + half, (cend - l5) // 2 - half, 5), key=abs):
            if not len(pos) or np.abs(pos - (centre + off)).min() > half:
                break
        ax.text(centre + off, 0.5, f"$\\mathbf{{v}}_{{\\mathrm{{{lab}}}}}$", ha="center", va="center",
                fontsize=6.3, color="black", zorder=4)
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color(OUTLINE)
        sp.set_linewidth(0.4)
    return len(dif)


def pair_groups(viol):
    """Runs of consecutive violating pairs (i, j), (i+1, j-1), ... as one stack each."""
    groups = []
    for p in sorted(viol):
        if groups and p[0] == groups[-1][-1][0] + 1 and p[1] == groups[-1][-1][1] - 1:
            groups[-1].append(p)
        else:
            groups.append([p])
    return groups


LAB = 5.4          # letter widths left of an inset's letters for the row names, as in Figure 1


def pair_inset(ax, cx, x, base, group, l5, cend, vsub, ec="black", w=3, window=None):
    """Figure 5(d)'s two-row box for design x: the nucleotides either side of a stack of CDS-UTR
    pairs of the concatenation design (`group`), the 5' strand above and its partner below. A
    rule joins each pair x folds, red where it is a CDS-UTR pair. A letter that differs from the
    concatenation design is bold black."""
    i0, j0 = group[0]
    k = len(group)
    if window:
        # the same stretch of the 5' strand in every box: [lo, hi), the partner strand antiparallel to it
        top = list(range(*window))
        bot = [j0 + i0 - p for p in top]
        cols = range(len(top))
    else:
        cols = range(-w, k + w)
        top = [i0 + c for c in cols]
        bot = [j0 - c for c in cols]
    st, seq = x["free_structure"], x["sequence"]
    fm = {a: b for a, b in pairs_of(st)}
    fm.update({b: a for a, b in pairs_of(st)})
    ax.add_patch(matplotlib.patches.Rectangle((cx - len(cols) / 2 - 0.4 - LAB, -0.15), len(cols) + 0.8 + LAB, 2.3,
                                              fc="white", ec=ec, lw=0.8, zorder=1, clip_on=False))

    def region(p):
        return "u5" if p < l5 else "v" if p < cend else "u3"
    names = {"u5": ("$\\mathbf{u}_5$", REGION["utr5"]), "v": (f"$\\mathbf{{v}}_{{\\mathrm{{{vsub}}}}}$", INK),
             "u3": ("$\\mathbf{u}_3$", REGION["utr3"])}
    for y, p in ((1.75, top[0]), (0.25, bot[0])):
        lab, col = names[region(p)]
        ax.text(cx - len(cols) / 2 - 0.4 - LAB / 2, y, lab, ha="center", va="center", fontsize=6.3, color=col,
                zorder=3, clip_on=False)
    for c, (a, b) in enumerate(zip(top, bot)):
        px = cx - len(cols) / 2 + c + 0.5
        for y, p in ((1.75, a), (0.25, b)):
            ch = seq[p] if 0 <= p < len(seq) else " "
            changed = 0 <= p < len(seq) and seq[p] != base["sequence"][p]
            red = (a, b) in group and fm.get(a) == b
            ax.text(px, y, ch, ha="center", va="center", fontsize=5.4, family="monospace",
                    fontweight="bold" if (red or changed) else "normal",
                    color=PAIR_RED if red else (INK if changed else
                          REGION["utr5"] if region(p) == "u5" else "#6e6d69"), zorder=3, clip_on=False)
        if fm.get(a) == b:
            red = (a, b) in group
            ax.plot([px, px], [0.62, 1.38], color=PAIR_RED if red else "#9b9a96", lw=0.7 if red else 0.5,
                    zorder=2, clip_on=False)


def fold_angle(st):
    """Rotation that lays the long axis of the molecule's drawing along x."""
    import RNA
    co = RNA.naview_xy_coordinates(st)
    P = np.array([[p.X, p.Y] for p in co])
    P = P - P.mean(axis=0)
    w, v = np.linalg.eigh(np.cov(P.T))
    d = v[:, np.argmax(w)]
    return -math.atan2(d[1], d[0])


def panel_structure(ax, title, x, ref, l5, cend, centre, theta, ring_color="black"):
    """The whole fold of design x, drawn as Figure 5(d) draws its folds. Nucleotides in a codon
    that differs from DERNA2's are dark, the UTRs are in their colours and the CDS-UTR pairs are
    red; a ring marks where the concatenation design has them."""
    import RNA
    st = x["free_structure"]
    n = len(st)
    co = RNA.naview_xy_coordinates(st)
    P = np.array([[p.X, p.Y] for p in co])
    R = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    X, Y = (P @ R.T).T
    viol = {(min(i, j) - 1, max(i, j) - 1) for i, j in (x.get("violations") or [])}
    dif = differs(x, ref, l5, cend)
    for i, j in pairs_of(st):
        if (i, j) not in viol:
            ax.plot(X[[i, j]], Y[[i, j]], color=PAIRC, lw=0.4, zorder=1)
    for k in range(n - 1):
        col, lw = BACK, 0.8
        if k < l5:
            col, lw = REGION["utr5"], 1.3
        elif k >= cend:
            col, lw = REGION["utr3"], 1.3
        elif (k - l5) // 3 in dif:
            col, lw = DIFFC, 1.3
        ax.plot(X[k:k + 2], Y[k:k + 2], color=col, lw=lw, solid_capstyle="round", zorder=3)
    for i, j in viol:
        ax.plot(X[[i, j]], Y[[i, j]], color=PAIR_RED, lw=1.3, zorder=5)
    ext = max(np.ptp(X), np.ptp(Y))
    pad = 0.02 * ext
    ax.set_xlim(X.min() - pad, X.max() + pad)
    ax.set_ylim(Y.min() - pad, Y.max() + pad)
    if centre:
        fw_, fh_ = ax.figure.get_size_inches()
        bb = ax.get_position()
        per_in = min((X.max() - X.min() + 2 * pad) / (bb.width * fw_), (Y.max() - Y.min() + 2 * pad) / (bb.height * fh_))
        cx0, cy0 = np.mean([X[p] for p in centre]), np.mean([Y[p] for p in centre])
        ax.add_patch(matplotlib.patches.Circle((cx0, cy0), 0.12 * per_in, fill=False, edgecolor=ring_color, lw=0.8,
                                               zorder=6))
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=6.3, pad=2, linespacing=1.25)
    ax.set_anchor("S")
    ax.apply_aspect()
    return X, Y, ((cx0, cy0, 0.12 * per_in) if centre else None)


def main():
    global CELL, THRESHOLD
    bundle, e2, hullcheck, ex, CELL, THRESHOLD, out = sys.argv[1:8]
    THRESHOLD = float(THRESHOLD)
    fs, fss = style()
    cells, head, masked, window, a4, hull, mf, refolded = load(bundle, e2, hullcheck)
    assert len(refolded) == len(cells), f"{len(cells) - len(refolded)} cells have no refold in {hullcheck}"
    sets = {c: modular_sets(c, masked, window, hull, mf) for c in cells}
    combo, _ = element_census(sorted(refolded), head, hull)
    _, share = element_census(cells, head, hull)
    split = pareto_split(cells, masked, window, sets)
    picks, l5, cend = best_designs(bundle, e2, ex)
    conc, ref = picks["concatenation"], picks["DERNA2"]
    n = len(ref["sequence"])

    W = 7.3
    r1t, r1h = 0.30, 1.45
    ykey = r1t + r1h + 0.72 - 0.19              # baseline of the (e) letter and key
    sth, sst = 0.14, 0.22
    r2t = ykey - 0.072 + sth                    # the long arc rises sth over the strips: its top is level with the (e) letter's
    sx0, sx1 = 0.98, 7.26
    sw = sx1 - sx0
    by0 = r2t + 2 * sst + sth + 0.10           # top of the band under the strips: titles and CDS-UTR pairs
    bandh, bh = 0.0, 0.40
    fy = by0 + bandh
    fh = 1.28
    H = fy + fh + 0.08
    fig = plt.figure(figsize=(W, H))

    def box(x0, y_top, w, h):
        """Axes from inches: left edge, top edge below the figure's top, width, height."""
        return fig.add_axes([x0 / W, 1 - (y_top + h) / H, w / W, h / H])

    # row 1
    xa, wa = 0.62, 1.05
    xb, wb = 2.28, 1.54                         # ends under the DERNA2 of the key after (b)
    xc, wc = 5.85, 1.41                         # the fronts, (d); the windows, (c), are at xd
    xd, wd = 4.34, 0.98
    axa = box(xa, r1t, wa, r1h)
    axb = [box(xb, r1t, wb, 0.66), box(xb, r1t + 0.78, wb, 0.67)]
    axc = box(xc, r1t, wc, r1h)                 # the fronts
    axd = box(xd, r1t, wd, r1h)                 # the windows
    el_handles = panel_elements(axa, combo, share)
    par_handles = panel_counts(axb, cells, sets, split)
    panel_front(axc, CELL, window, a4, mf, concat_front(bundle, e2, CELL), sets, picks=picks)
    win_handles = panel_d_windows(axd, cells, a4)
    # each y label's right edge sits at a fixed distance left of its axes; the letter starts at the
    # label's left edge, one text height further left
    ycen = {}
    # headroom above the bars for the keys: each spine still ends at its last tick, and its label is
    # centred on that range
    axa.set_ylim(0, 138)
    axa.set_yticks(range(0, 101, 20))
    axa.spines["left"].set_bounds(0, 100)
    ycen[axa] = 50 / 138
    for a_ in axb:
        a_.tick_params(axis="y", length=2.0, pad=1.5)        # so that both labels can sit at one gap
    n_prot = len({c.split("__")[0] for c in cells})
    top = n_prot * 23.3 / 22                  # the same headroom over the proteins line at any count
    axb[0].set_ylim(0, top)
    axb[0].set_yticks(sorted({0, 10, *range(20, n_prot, 10), n_prot}))
    axb[0].spines["left"].set_bounds(0, n_prot)
    ycen[axb[0]] = n_prot / 2 / top
    axb[1].set_ylim(1, 1.3e4)
    axb[1].set_yticks([1, 1e2, 1e4])
    axb[1].set_yticklabels(["$10^{0}$", "$10^{2}$", "$10^{4}$"])
    axb[1].yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axb[1].spines["left"].set_bounds(1, 1e4)
    ycen[axb[1]] = 0.5 * np.log10(1e4) / np.log10(1.3e4)
    axd.tick_params(axis="y", length=2.0, pad=1.5)
    axd.set_ylim(0, 1000)
    axd.set_yticks(range(0, 1001, 200))
    axd.spines["left"].set_bounds(0, 1000)
    ycen[axd] = 0.5
    lab_gap = {"a": 0.30, "b": 0.30, "c": 0.27, "d": 0.37}
    for ax, letter in [(axa, "a"), (axb[0], "b"), (axd, "c"), (axc, "d")]:
        axes_ = axb if letter == "b" else [ax]
        # (b)'s upper label hugs its narrower ticks; the lower panel's 10^4 sets the other's gap
        gap0 = 0.255 if letter == "b" else lab_gap[letter]
        for a_ in axes_:
            pos = a_.get_position()
            g = gap0
            a_.yaxis.set_label_coords(-g / (pos.width * W), ycen.get(a_, 0.5))
        x0 = ax.get_position().x0 * W - gap0 - 0.09
        fig.text(x0 / W, 1 - (r1t - 0.03) / H, f"\\textbf{{({letter})}}", ha="left", va="bottom", fontsize=fs)
        letter_x = x0 if letter == "a" else letter_x

    # (e): strips, then the boxes at the CDS-UTR site, the changed codons, and the three folds
    axs_all, nd = [], {}
    for r, (name, x) in enumerate(picks.items()):
        axs = box(sx0, r2t + sst * r, sw, sth)
        nd[name] = codon_strip(axs, x, ref, l5, cend, VSUB[name])
        axs.text(-0.008, 0.5, name, transform=axs.transAxes,
                 ha="right", va="center", fontsize=fss)
        axs_all.append(axs)
    fig.text(letter_x / W, 1 - ykey / H, "\\textbf{(e)}", ha="left", va="baseline", fontsize=fs)
    # one pair of boxes per stack of CDS-UTR pairs of a non-modular design: the design's own
    # letters and DERNA2's at the same positions
    letter_w = 0.056                  # inches per letter in a box
    by = by0 + 0.11
    u5s, u3s = head[CELL]["utr5"]["sequence"], head[CELL]["utr3"]["sequence"]
    em5, em3 = element_masks(u5s, u3s)
    name_el = {"kozak": "Kozak", "pyr": "PRE", "polyA": "poly(A)", "are": "ARE"}
    theta = fold_angle(conc["free_structure"])
    fw = (W - 0.15) / 3
    grp = {}
    for nm in ("concatenation", "mRNAfold"):
        v = sorted((min(i, j) - 1, max(i, j) - 1) for i, j in picks[nm].get("violations") or [])
        grp[nm] = pair_groups(v)[0]
    grp["DERNA2"] = grp["mRNAfold"]               # DERNA2's letters at the Kozak site mRNAfold pairs
    base = {"concatenation": picks["concatenation"], "mRNAfold": picks["mRNAfold"], "DERNA2": picks["mRNAfold"]}
    # the columns move right until the last fold's drawing ends at the strips' right edge
    tmp = box(0.1 + 2 * fw, fy, fw * 0.98, fh)
    X_, Y_, _ = panel_structure(tmp, "", ref, ref, l5, cend, [], theta)
    shift = sx1 - tmp.transData.transform([[X_.max(), Y_[0]]])[0][0] / fig.dpi
    tmp.remove()
    for k, (name, x) in enumerate(picks.items()):
        pr = ("''", "'", "")[k]
        x_col = 0.1 + shift + k * fw
        npairs = len(x.get("violations") or [])
        l3 = f"{npairs} CDS--UTR pairs" if npairs else "no CDS--UTR pairs"
        g = grp[name]
        cls = sorted({k_ for i_, j_ in g for k_ in
                      (em5[i_] if i_ < l5 else em3[j_ - cend] if 0 <= j_ - cend < len(em3) else set())})
        el = " and ".join(name_el[k_] for k_ in cls) if cls else "no annotated element"
        # (a)'s colour for the element a contact touches, black without a contact: the ring, leaders and border
        key = ("kozak+pyr" if {"kozak", "pyr"} <= set(cls) else "kozak" if "kozak" in cls else "pyr" if "pyr" in cls
               else "rest" if cls else "none")
        mark = EL[key] if npairs else "black"
        site = sorted({p for i_, j_ in g for p in (i_, j_)})
        axp = box(x_col, fy, fw * 0.98, fh)
        X, Y, ring = panel_structure(axp, "", x, ref, l5, cend, site, theta, mark)
        drawn = axp.transData.transform(np.c_[X, Y]) / fig.dpi
        x_left, x_right = drawn[:, 0].min(), drawn[:, 0].max()
        lines = (name, f"MFE ${x['free_mfe']:.1f}$", l3)
        fig.text(x_right / W, 1 - (fy + 0.02) / H, "\n".join(lines), ha="right", va="top", fontsize=6.3,
                 linespacing=1.25)
        title_pts = np.array([(xx, yy) for xx in np.arange(x_right - 0.047 * max(len(s) for s in (name, l3)),
                                                           x_right + 1e-9, 0.04)
                              for yy in np.arange(fy + 0.02, fy + 0.40, 0.04)])
        # the inset: this design's letters at its site, in the free space nearest the ring
        win = (l5 - 6, l5 + 3)                 # the Kozak and the start codon, in every box
        assert all(win[0] <= i_ < win[1] for i_, _ in g), (name, g)
        kk = win[1] - win[0]
        bw, bh_ = (kk + 1 + LAB) * letter_w, bh
        pts = axp.transData.transform(np.c_[X, Y]) / fig.dpi
        pts = np.c_[pts[:, 0], H - pts[:, 1]]
        pts = np.vstack([pts, (pts[1:] + pts[:-1]) / 2, title_pts])
        rc = np.array(axp.transData.transform([ring[0], ring[1]])) / fig.dpi
        rc = np.array([rc[0], H - rc[1]])
        rr = ring[2] / (axp.transData.transform([1, 0])[0] - axp.transData.transform([0, 0])[0]) * fig.dpi * 0 + 0.12
        bb = axp.get_position()
        lo_x, hi_x = x_left, x_right - bw
        lo_y, hi_y = fy - 0.02, fy + fh - bh_
        best = None
        for dmin, rmin in ((0.07, 0.55), (0.05, 0.45), (0.04, 0.35), (0.03, 0.25)):
            for cx_ in np.arange(lo_x, hi_x + 1e-9, 0.04):
                for cy_ in np.arange(lo_y, hi_y + 1e-9, 0.04):
                    d = np.hypot(np.clip(pts[:, 0], cx_, cx_ + bw) - pts[:, 0],
                                 np.clip(pts[:, 1], cy_, cy_ + bh_) - pts[:, 1]).min()
                    cost = np.hypot(rc[0] - (cx_ + bw / 2), rc[1] - (cy_ + bh_ / 2))
                    if d >= dmin and cost >= rmin and (best is None or cost < best[0]):
                        best = (cost, cx_, cy_)
            if best:
                break
        _, ix, iy = best
        axb_ = box(ix, iy, bw, bh_)
        axb_.set_xlim(-LAB, kk + 1)
        axb_.set_ylim(-0.15, 2.15)
        axb_.axis("off")
        pair_inset(axb_, (kk + 1) / 2, x, base[name], g, l5, cend, VSUB[name], ec=mark, window=win)
        fig.text((ix + bw / 2) / W, 1 - (iy - 0.01) / H, el, ha="center", va="bottom", fontsize=5.6, color=MUTED)
        # dashed leaders from the ring to the inset's two corners that face it
        corners = [(ix, iy), (ix + bw, iy), (ix, iy + bh_), (ix + bw, iy + bh_)]
        order = sorted(range(4), key=lambda q: np.hypot(corners[q][0] - rc[0], corners[q][1] - rc[1]))
        for q in order[:2]:
            s = np.array(corners[q]) - rc
            s0 = rc + s / np.hypot(*s) * rr
            fig.add_artist(Line2D([s0[0] / W, corners[q][0] / W], [1 - s0[1] / H, 1 - corners[q][1] / H],
                                  transform=fig.transFigure, color=mark, lw=0.8,
                                  ls=(0, (2.5, 1.5)), zorder=4))

    # the key's rows read concatenation, mRNAfold / DERNA2: the legend fills its columns first
    mh = [matplotlib.patches.Patch(fc=COL[m], lw=0, label=NAMES[m]) for m in METHODS]
    # the method key serves (b) and (d): one row on the panel letters' line, after each
    for kx, kl in ((sx1 + 0.03, "lower right"), (2.11, "lower left")):
        fig.legend(handles=mh, loc=kl, bbox_to_anchor=(kx / W, 1 - (r1t + 0.001) / H), ncol=3, frameon=False,
                   fontsize=fss, handlelength=1.0, columnspacing=0.7, handletextpad=0.3, borderaxespad=0.0)
    # the other keys sit inside their axes, in the headroom above the bars or under the top tick
    axa.legend(handles=el_handles, loc="upper left", bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False,
               fontsize=5.7, handlelength=0.8, columnspacing=0.6, handletextpad=0.3, labelspacing=0.2,
               borderaxespad=0)
    axb[0].legend(handles=par_handles[1:], labels=["dominated"], loc="upper left", bbox_to_anchor=(0.0, n_prot / top),
                  ncol=1, frameon=False, fontsize=5.7, handlelength=1.2, columnspacing=0.5, handletextpad=0.3,
                  borderaxespad=0, handler_map={tuple: matplotlib.legend_handler.HandlerTuple(ndivide=None, pad=0)})
    # listed top to bottom in the order the bars stack, which also puts the widest label at the top, over the
    # short bars
    axd.legend(handles=win_handles[::-1], loc="upper left", ncol=1, frameon=False, fontsize=5.4, handlelength=0.7,
               columnspacing=0.5, handletextpad=0.25, labelspacing=0.15, borderaxespad=0.0)
    # the (e) key, right-aligned over the strips and on the baseline of the (e) letter
    yb = 1 - ykey / H
    x0 = sx1 - 1.62
    fig.add_artist(matplotlib.patches.Rectangle((x0 / W, yb), 0.10 / W, 0.075 / H, fc=DIFFC, lw=0,
                                                transform=fig.transFigure))
    fig.text((x0 + 0.14) / W, yb, "different codon", ha="left", va="baseline", fontsize=fss)
    x1 = x0 + 0.14 + 0.53 + 0.09
    fig.add_artist(Line2D([x1 / W, (x1 + 0.11) / W], [yb + 0.035 / H] * 2, color=PAIR_RED, lw=1.1,
                          transform=fig.transFigure))
    fig.text((x1 + 0.14) / W, yb, "CDS--UTR base pair", ha="left", va="baseline", fontsize=fss)
    # the figure's own margins, measured: 0.23 in of white at the left and 0.18 in at the top, which are cropped
    # to 0.03 in
    crop = matplotlib.transforms.Bbox([[0.20, 0.0], [W, H - 0.15]])
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, metadata={"CreationDate": None}, bbox_inches=crop, pad_inches=0)

    print("non-modular share of concatenation designs:",
          {p: f"{100 * v / n_:.0f}%" for p, (n_, v) in share.items()})
    for m in METHODS:
        print(f"{NAMES[m]:<24}", [split[m][p][0] for p in PAIRS], "dominated",
              [split[m][p][2] for p in PAIRS])
    for name, x in picks.items():
        print(f"{name:<14} CAI {x['cai']:.4f} MFE no pairs {x['mfe_nopair']} pairs {x.get('violations')}"
              f" differing codons {nd[name]}")


if __name__ == "__main__":
    main()
