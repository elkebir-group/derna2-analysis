#!/usr/bin/env python3
"""Gather every design for the plot into one JSON and draw the static PNG.

x = accumulated log-CAI, y = kcal/mol, each design at the energy its method
optimized (front 1: free MFE; front 2 and 3: the arm's guarded energy;
mRNAfold: the energy of its counterpart arm). Two panels: symmetric, and 5′
(no_clamp and unpaired); front 1 on both.

Usage: build_plot.py ROOT MRNAFOLD_COMPARE OUT_JSON OUT_PNG [OUT_HTML [NOTE]]

With OUT_HTML, the data is also embedded in `plot_template.html` for the
interactive page; NOTE is a line printed under the panels.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED = "#0b0b0b", "#52514e"
CSV = Path(__file__).resolve().parents[2] / "packages/derna-codon-dp/data/codon_usage_human.csv"
LEN5, CDS_NT = 50, 3 * 78
MRNAFOLD_MARK = {"plain": "+", "guard": "x", "unpaired5": "s", "unpaired_both": "o"}
MRNAFOLD_PANEL = {"plain": ["symmetric", "five"], "guard": ["symmetric"],
                  "unpaired5": ["five"], "unpaired_both": ["five"]}


def weights() -> dict[str, float]:
    rows = [r.split(",") for r in CSV.read_text().splitlines()[1:] if r.strip()]
    best: dict[str, float] = {}
    for aa, _, f in rows:
        best[aa] = max(best.get(aa, 0.0), float(f))
    return {codon: float(f) / best[aa] for aa, codon, f in rows}


W = weights()


def cai_log(sequence: str) -> float:
    cds = sequence[LEN5:LEN5 + CDS_NT]
    return sum(math.log(W[cds[i:i + 3]]) for i in range(0, CDS_NT, 3))


def lines(path: Path) -> list[dict]:
    return [json.loads(x) for x in path.read_text().splitlines() if x.strip()] if path.exists() else []


def front(root: Path, arm: str) -> list[dict]:
    hull = {(v["mfe_kcal_mol"], round(v["cai_log"], 6)) for v in lines(root / arm / "vertices.jsonl")}
    scored = {r["sequence"]: r for r in lines(root / arm / "front.scored.jsonl")}
    out = []
    for p in lines(root / arm / "front.jsonl"):
        s = scored.get(p["sequence"], {})
        out.append({
            "E": p["mfe_kcal_mol"], "L": p["cai_log"], "cai": p["cai"],
            "hull": (p["mfe_kcal_mol"], round(p["cai_log"], 6)) in hull,
            "free": s.get("free_mfe_kcal_mol"),
            "so": None if arm == "unguarded" else s.get(f"{arm}_so"),
        })
    return out


def front3(root: Path, arm: str) -> list[dict]:
    scored = {(r["vertex"]): r for r in lines(root / arm / "front3_designs.scored.jsonl")}
    out = []
    for r in lines(root / arm / "front3.jsonl"):
        s = scored.get(r["vertex"])
        row = {"vertex": r["vertex"], "lambda": r["lambda"], "finished": r["finished"],
               "rank": r["rank"], "walk_s": r["walk_s"],
               "vE": r["vertex_mfe_kcal_mol"], "vL": r["vertex_cai_log"]}
        if s:
            row.update({"E": s[f"{arm}_guarded_kcal_mol"], "L": cai_log(s["sequence"]),
                        "free": s["free_mfe_kcal_mol"], "so": s[f"{arm}_so"]})
        else:
            # Unfinished: drawn hollow at the vertex, labelled with the rank reached.
            row.update({"E": r["vertex_mfe_kcal_mol"], "L": r["vertex_cai_log"]})
        out.append(row)
    return out


def main() -> None:
    root, compare, out_json, out_png = map(Path, sys.argv[1:5])
    data = {
        "front1": front(root, "unguarded"),
        "front2": {a: front(root, a) for a in ("symmetric", "no_clamp", "unpaired")},
        "front3": {a: front3(root, a) for a in ("symmetric", "no_clamp", "unpaired")},
        "mrnafold": [
            # The designed structure's energy: what the arm optimized. It equals
            # our refold under the counterpart guard on every arm but
            # unpaired_both, which also keeps the 3′ UTR unpaired.
            {"arm": r["arm"], "vertex": r["vertex"], "lambda": r["lambda"],
             "E": r["designed_energy_kcal_mol"],
             "L": r["theirs_cai_log"], "free": r["free_mfe_kcal_mol"],
             "so_symmetric": r["symmetric_so"], "so_unpaired": r["unpaired_so"],
             "matches_ours": abs(r["d_obj"]) < 1e-3, "d_obj": r["d_obj"]}
            for r in lines(compare)
        ],
    }
    if len(sys.argv) > 6:
        data["note"] = sys.argv[6]
    out_json.write_text(json.dumps(data))
    if len(sys.argv) > 5:
        template = (Path(__file__).parent / "plot_template.html").read_text()
        Path(sys.argv[5]).write_text(template.replace("/*DATA*/null", json.dumps(data)))

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.2), sharey=True)
    panels = [("symmetric", ["symmetric"], "Symmetric guard (5′↔CDS, CDS↔3′ forbidden)"),
              ("five", ["no_clamp", "unpaired"], "5′ guards (no_clamp ◆, 5′ unpaired ■)")]
    shape2 = {"symmetric": "D", "no_clamp": "D", "unpaired": "s"}
    shape3 = {"symmetric": "^", "no_clamp": "^", "unpaired": "v"}
    for ax, (panel, arms, title) in zip(axes, panels):
        f1 = data["front1"]
        for hull, kw in ((False, {"facecolors": "none", "edgecolors": BLUE, "alpha": 0.6}),
                         (True, {"color": BLUE})):
            pts = [p for p in f1 if p["hull"] == hull]
            ax.scatter([p["L"] for p in pts], [p["E"] for p in pts], marker="o", s=22,
                       linewidths=1, label=f"front 1 unguarded ({'hull' if hull else 'off-hull'})", **kw)
        for arm in arms:
            f2 = data["front2"][arm]
            for hull, kw in ((False, {"facecolors": "none", "edgecolors": ORANGE, "alpha": 0.6}),
                             (True, {"color": ORANGE})):
                pts = [p for p in f2 if p["hull"] == hull]
                ax.scatter([p["L"] for p in pts], [p["E"] for p in pts], marker=shape2[arm], s=24,
                           linewidths=1, label=f"front 2 {arm} ({'hull' if hull else 'off-hull'})", **kw)
            f3 = data["front3"][arm]
            for r in f3:
                ax.plot([r["vL"], r["L"]], [r["vE"], r["E"]], color=AQUA, lw=0.8, alpha=0.7)
            for done, kw in ((True, {"color": AQUA}), (False, {"facecolors": "none", "edgecolors": AQUA})):
                pts = [r for r in f3 if r["finished"] == done]
                if pts:
                    ax.scatter([r["L"] for r in pts], [r["E"] for r in pts], marker=shape3[arm], s=60,
                               linewidths=1.4, label=f"front 3 {arm} ({'found' if done else 'unfinished'})",
                               zorder=4, **kw)
        for arm, mark in MRNAFOLD_MARK.items():
            if panel not in MRNAFOLD_PANEL[arm]:
                continue
            pts = [r for r in data["mrnafold"] if r["arm"] == arm]
            kw = {"color": INK} if mark in "+x" else {"facecolors": "none", "edgecolors": INK}
            ax.scatter([r["L"] for r in pts], [r["E"] for r in pts], marker=mark, s=50, linewidths=1,
                       label=f"mRNAfold {arm}", zorder=5, **kw)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("accumulated log-CAI  Σ ln w")
        ax.grid(color="#e6e5e1", lw=0.6)
        ax.legend(fontsize=7, loc="upper left", frameon=False)
    axes[0].set_ylabel("energy (kcal/mol): free MFE for front 1, guarded energy otherwise")
    fig.suptitle("GLPE_HUMAN n = 78 with hbb UTRs (416 nt): three fronts and mRNAfold", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)


if __name__ == "__main__":
    main()
