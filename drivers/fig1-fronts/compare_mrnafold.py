#!/usr/bin/env python3
"""Compare mRNAfold's designs with ours at the symmetric front's hull weights.

Reads `runs.scored.jsonl` (`run_mrnafold.py` output passed through
`glpe_hbb_rescore`) and the four fronts. For each run it reports the
design's free MFE, guarded energies, structural optimality under each guard,
and the forbidden pairs in mRNAfold's *designed* structure, and compares its
objective `λ·E − (1 − λ)·L` with ours, where `E` is the energy the arm's
counterpart optimizes:

- `plain` against our unguarded optimum at that weight (E = free MFE);
- `guard` against our symmetric optimum (E = guarded energy, symmetric);
- `unpaired5` against our 5′-unpaired optimum (E = guarded energy, unpaired);
- `unpaired_both` has no counterpart of ours; it is scored under the
  5′-unpaired guard for reference.

Our optimum at a weight is the minimum of the objective over the whole front,
which at a weight inside a hull vertex's interval is that vertex.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

COUNTERPART = {"plain": "unguarded", "guard": "symmetric", "unpaired5": "unpaired",
               "unpaired_both": "unpaired"}
LEN5, CDS_END = 50, 50 + 3 * 78
CLASS_NAMES = {("5", "5"): "5-5", ("5", "cds"): "5-cds", ("5", "3"): "5-3",
               ("cds", "cds"): "cds-cds", ("cds", "3"): "cds-3", ("3", "3"): "3-3"}
FORBIDS = {"symmetric": {"5-cds", "cds-3"}, "no_clamp": {"5-cds", "5-3"},
           "unpaired": {"5-5", "5-cds", "5-3"}}


def region(i: int) -> str:
    return "5" if i < LEN5 else "cds" if i < CDS_END else "3"


def classes(structure: str) -> dict[str, int]:
    stack, count = [], {}
    for j, ch in enumerate(structure):
        if ch == "(":
            stack.append(j)
        elif ch == ")":
            c = CLASS_NAMES[(region(stack.pop()), region(j))]
            count[c] = count.get(c, 0) + 1
    return count


def objective(lam: float, e: float, cai_log: float) -> float:
    return lam * e - (1.0 - lam) * cai_log


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True, type=Path)
    ap.add_argument("--fronts", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    fronts = {
        arm: [json.loads(x) for x in (args.fronts / arm / "front.jsonl").read_text().splitlines()]
        for arm in set(COUNTERPART.values())
    }
    rows = []
    for line in args.runs.read_text().splitlines():
        r = json.loads(line)
        lam, arm = r["lambda"], r["arm"]
        ours_arm = COUNTERPART[arm]
        best = min(fronts[ours_arm], key=lambda p: objective(lam, p["mfe_kcal_mol"], p["cai_log"]))
        e_theirs = r["free_mfe_kcal_mol"] if ours_arm == "unguarded" else r[f"{ours_arm}_guarded_kcal_mol"]
        designed = classes(r["structure"])
        row = {
            "arm": arm, "vertex": r["vertex"], "lambda": lam,
            "counterpart": ours_arm,
            "theirs_E": e_theirs, "theirs_cai_log": r["cai_log"],
            "theirs_obj": objective(lam, e_theirs, r["cai_log"]),
            "ours_E": best["mfe_kcal_mol"], "ours_cai_log": best["cai_log"],
            "ours_obj": objective(lam, best["mfe_kcal_mol"], best["cai_log"]),
            "same_sequence": r["sequence"] == best["sequence"],
            "free_mfe_kcal_mol": r["free_mfe_kcal_mol"],
            "designed_energy_kcal_mol": r["score_energy_kcal_mol"],
            "designed_classes": designed,
            "wall_s": r["wall_s"],
            "sequence": r["sequence"],
            "structure": r["structure"],
        }
        row["d_obj"] = row["theirs_obj"] - row["ours_obj"]
        for g, forbidden in FORBIDS.items():
            row[f"{g}_guarded_kcal_mol"] = r[f"{g}_guarded_kcal_mol"]
            row[f"{g}_so"] = r[f"{g}_so"]
            row[f"{g}_designed_forbidden"] = sum(v for c, v in designed.items() if c in forbidden)
        rows.append(row)
    with args.out.open("w") as f:
        for row in rows:
            f.write(json.dumps(row) + "\n")
    for arm in COUNTERPART:
        sub = [x for x in rows if x["arm"] == arm]
        if not sub:
            continue
        print(json.dumps({
            "arm": arm, "runs": len(sub), "counterpart": COUNTERPART[arm],
            "matches_ours": sum(abs(x["d_obj"]) < 1e-3 for x in sub),
            "same_sequence": sum(x["same_sequence"] for x in sub),
            "max_d_obj": round(max(x["d_obj"] for x in sub), 4),
            "min_d_obj": round(min(x["d_obj"] for x in sub), 4),
            "symmetric_so": sum(x["symmetric_so"] for x in sub),
            "no_clamp_so": sum(x["no_clamp_so"] for x in sub),
            "unpaired_so": sum(x["unpaired_so"] for x in sub),
            "designed_forbidden_symmetric": sum(x["symmetric_designed_forbidden"] > 0 for x in sub),
            "designed_forbidden_unpaired": sum(x["unpaired_designed_forbidden"] > 0 for x in sub),
            "wall_s_median": sorted(x["wall_s"] for x in sub)[len(sub) // 2],
        }))


if __name__ == "__main__":
    main()
