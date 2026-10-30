"""Steps between neighbouring designs of each E2 front and of its sweep hull.

    python fig4_steps.py HULL_DIR PER_PROTEIN_CSV OUT.npz

For every protein of PER_PROTEIN_CSV: the Pareto front as the walk wrote it
(HULL_DIR/<P>.front.jsonl, energy in kcal/mol and CAI), the sweep's hull vertices
(HULL_DIR/<P>.walk.json), and the |dMFE| and dCAI between consecutive designs of each, ordered
by CAI. The front is Pareto-filtered on raw (energy, log CAI), so two designs that tie only in
the 6-decimal CAI field both stay.
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

LTOL = 1e-5


def pareto(rows):
    """(E dekacal, log CAI, CAI) sorted by E; keep a point only if log CAI strictly increases."""
    srt = sorted(set(rows), key=lambda q: (q[0], -q[1]))
    out = []
    for q in srt:
        if not out or q[1] > out[-1][1] + LTOL:
            out.append(q)
    return out


def steps(pts):
    """|dMFE| (kcal/mol) and dCAI between neighbours, pts sorted by CAI."""
    p = np.array(sorted(pts, key=lambda q: q[2]))
    return np.abs(np.diff(p[:, 0])) / 100, np.diff(p[:, 2])


def protein(hull_dir, p):
    rows = [
        (round(r["mfe_kcal_mol"] * 100), r["cai_log"], r["cai"])
        for r in map(json.loads, open(hull_dir / f"{p}.front.jsonl"))
        if "mfe_kcal_mol" in r
    ]
    front = pareto(rows)
    sweep = [
        (round(v["mfe_kcal_mol"] * 100), v["adaptation_log"])
        for v in json.load(open(hull_dir / f"{p}.walk.json"))["vertices"]
    ]
    on = np.array([any(e == q[0] and abs(lg - q[1]) < 2e-6 for e, lg in sweep) for q in front])
    assert on.sum() == len(sweep), (p, on.sum(), len(sweep))
    hull = [q for q, h in zip(front, on, strict=True) if h]
    return front, on, steps(front), steps(hull)


def main():
    hull_dir, out_path = Path(sys.argv[1]), Path(sys.argv[3])
    prots = sorted(r["protein"] for r in csv.DictReader(open(sys.argv[2])))
    out = {}
    for p in prots:
        front, on, (fm, fc), (hm, hc) = protein(hull_dir, p)
        a = np.array(front, float)
        out[f"{p}|pts"] = np.column_stack([a[:, 0] / 100, a[:, 2], on])
        out[f"{p}|front"] = np.column_stack([fm, fc])
        out[f"{p}|hull"] = np.column_stack([hm, hc])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_path, **out)
    print(len(out) // 3, "proteins")


if __name__ == "__main__":
    main()
