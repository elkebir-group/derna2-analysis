"""Peak memory of each method by protein length, in Figure 4's length groups and colours.

    python fig_memory.py SWEEP_JSONL PIN PER_PROTEIN_CSV OUT.pdf

One box per method and length group over the proteins on which the method delivered a front
within its budget, the rule of Figure 4d. A cell's peak memory is its cell row's `maxrss_kb`, else
the largest of its lambda runs' (E1-SUMMARY).
"""

import importlib.util
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).parent
spec = importlib.util.spec_from_file_location("fig4", HERE / "fig4_overview.py")
f4 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f4)

LAMBDA_KINDS = {"derna", "lineardesign", "mrnafold", "lcds", "derna1rs"}


def rank(r):
    """As fig4_cells.py: a delivered run first, then the shorter wall time."""
    return (r.get("exit") != "ok", r.get("wall_s") or 0.0)


def peaks(sweep, pin):
    """(arm, protein) -> (best cell row, peak memory in GB)."""
    cells, lam = defaultdict(list), defaultdict(float)
    for line in open(sweep):
        r = json.loads(line)
        if r.get("pin") != pin:
            continue
        key = (r["arm"], r["protein"])
        if r.get("kind") == "cell":
            cells[key].append(r)
        elif r.get("kind") in LAMBDA_KINDS and r.get("maxrss_kb"):
            lam[key] = max(lam[key], r["maxrss_kb"])
    out = {}
    for key, rs in cells.items():
        best = min(rs, key=rank)
        kb = best.get("maxrss_kb") or lam.get(key)
        if kb:
            out[key] = (best, kb / 1024**2)
    return out


def main():
    sweep, pin, per_csv, out = sys.argv[1], sys.argv[2], sys.argv[3], Path(sys.argv[4])
    per = pd.read_csv(per_csv).set_index("protein")
    corpus = set(per.index)
    rows = []
    for (arm, protein), (cell, gb) in peaks(sweep, pin).items():
        if arm in f4.ORDER and protein in corpus:
            lim = 7200 if cell["n_aa"] < f4.LONG else 21600
            done = cell["exit"] == "ok" and cell["wall_s"] <= lim * (1 + f4.TOL_BUDGET)
            rows.append((arm, protein, cell["n_aa"], gb, done))
    t = pd.DataFrame(rows, columns=["arm", "protein", "n_aa", "gb", "done"])
    t = t[t.done]
    t["g"] = f4.group_of(t.n_aa.to_numpy())

    f4.set_style()
    fig, ax = plt.subplots(figsize=(4.6, 2.3))
    rng = np.random.default_rng(4)
    k = len(f4.ORDER)
    w = 0.8 / k
    for j in range(len(f4.GLABELS)):
        for i, arm in enumerate(f4.ORDER):
            m = t[(t.arm == arm) & (t.g == j)]
            x = j + (i - (k - 1) / 2) * w
            f4.box(ax, m.gb.to_numpy(), x, f4.PAL[arm], w=w * 0.9)
            f4.dots(ax, x, m.gb.to_numpy(), f4.PAL[arm], rng, (m.n_aa > f4.LONG).to_numpy(),
                    half=w * 0.3, s=1.2)
    ax.set_yscale("log")
    ax.set_ylabel("peak memory (GB)", labelpad=2.2)
    ax.set_xticks(range(len(f4.GLABELS)))
    ax.set_xticklabels(f4.GLABELS)
    ax.set_xlim(-0.5, len(f4.GLABELS) - 0.5)
    ax.set_xlabel("protein length (aa)", labelpad=2.2)
    ax.spines[["top", "right"]].set_visible(False)
    hs = [matplotlib.patches.Patch(fc=f4.PAL[a], alpha=0.5, lw=0) for a in f4.ORDER]
    ax.legend(hs, [f4.NAMES[a] for a in f4.ORDER], loc="upper left", ncol=3, frameon=False,
              fontsize=f4.FSS, handlelength=0.9, handleheight=0.7, columnspacing=1.0)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", metadata={"CreationDate": None})
    fig.savefig(out.with_suffix(".png"), dpi=200, bbox_inches="tight", metadata={"Software": None})
    print(t.groupby("arm").gb.agg(["count", "median", "max"]).round(2))


if __name__ == "__main__":
    main()
