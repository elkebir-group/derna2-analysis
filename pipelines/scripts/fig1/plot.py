"""The data panels of the overview figure: GLPE_HUMAN (78 aa) with the hbb UTRs.

    python plot.py RESULTS_DIR OUT_PDF

a  exact CDS-only front in (log-CAI, MFE), its hull vertices marked: the only points a
   weight sweep can reach.
b  the same designs folded alone and with the UTRs attached (free fold).
c  the symmetric-guarded front (no CDS-UTR pairs) against the free fold of the same designs.
d  the guarded front and the structurally optimal front (free MFE = guarded MFE).
"""

import sys
from pathlib import Path

from derna2_analysis.figures import INK_2, SLOTS, save, style
from derna2_analysis.io import read_jsonl


def xy(rows, field):
    rows = sorted(rows, key=lambda r: r["cai_log"])
    return [r["cai_log"] for r in rows], [r[field] for r in rows]


def main():
    results, out = sys.argv[1:]
    res = Path(results)
    a = read_jsonl(res / "a_cdsonly_front.jsonl")
    b = read_jsonl(res / "b_cdsonly_with_utrs.scored.jsonl")
    cd = read_jsonl(res / "cd_symmetric_guarded_front.scored.jsonl")
    so = read_jsonl(res / "d_so_front.jsonl")
    style()
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), layout="constrained")
    (pa, pb), (pc, pd) = axes

    x, y = xy(a, "mfe_kcal_mol")
    pa.plot(x, y, color=SLOTS[0], marker="o", markersize=2, linewidth=1, label="Pareto front")
    hx, hy = xy([r for r in a if r.get("hull_vertex")], "mfe_kcal_mol")
    pa.plot(
        hx,
        hy,
        color=INK_2,
        linestyle="--",
        linewidth=0.8,
        marker="s",
        markersize=4,
        markerfacecolor="none",
        label=f"hull vertices ({len(hx)} of {len(a)})",
    )
    pa.set_title("a  CDS only", loc="left")

    for field, color, marker, label in [
        ("cds_mfe_kcal_mol", SLOTS[0], "o", "CDS alone"),
        ("free_mfe_kcal_mol", SLOTS[1], "s", "with UTRs"),
    ]:
        x, y = xy(b, field)
        pb.plot(x, y, color=color, marker=marker, markersize=2, linewidth=1, label=label)
    pb.set_title("b  the same designs, UTRs attached", loc="left")

    for field, color, marker, label in [
        ("mfe_kcal_mol", SLOTS[2], "^", "guarded MFE"),
        ("free_mfe_kcal_mol", SLOTS[1], "s", "free MFE"),
    ]:
        x, y = xy(cd, field)
        pc.plot(x, y, color=color, marker=marker, markersize=2, linewidth=1, label=label)
    pc.set_title("c  guarded front: no CDS-UTR pairs", loc="left")

    x, y = xy(cd, "mfe_kcal_mol")
    pd.plot(x, y, color=SLOTS[2], marker="^", markersize=2, linewidth=1, label="guarded front")
    x, y = xy(so, "mfe_kcal_mol")
    pd.plot(
        x,
        y,
        color=SLOTS[6],
        marker="o",
        markersize=2,
        linewidth=1,
        label=f"structurally optimal ({len(so)})",
    )
    pd.set_title("d  structurally optimal front", loc="left")

    for ax in axes.flat:
        ax.set_xlabel("log CAI (Σ ln w)")
        ax.set_ylabel("MFE (kcal/mol)")
        ax.legend(loc="upper left", fontsize=7)
    save(fig, out)


if __name__ == "__main__":
    main()
