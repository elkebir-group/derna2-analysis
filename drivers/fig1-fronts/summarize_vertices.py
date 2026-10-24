#!/usr/bin/env python3
"""Summarize each arm's rescored hull vertices: how many are structurally
optimal at rank 1 under the arm's own guard, their gaps, and their weights.

Usage: summarize_vertices.py ROOT   (reads ROOT/<arm>/vertices.scored.jsonl)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ARMS = ("unguarded", "symmetric", "no_clamp", "unpaired")


def main() -> None:
    root = Path(sys.argv[1])
    for arm in ARMS:
        rows = [json.loads(x) for x in (root / arm / "vertices.scored.jsonl").read_text().splitlines()]
        out = {"arm": arm, "hull_vertices": len(rows)}
        if arm == "unguarded":
            out["free_neq_front"] = sum(abs(r["free_mfe_kcal_mol"] - r["mfe_kcal_mol"]) > 1e-6 for r in rows)
            out["so_under"] = {g: sum(r[f"{g}_so"] for r in rows) for g in ARMS[1:]}
        else:
            out["so_at_rank1"] = sum(r[f"{arm}_so"] for r in rows)
            out["guarded_neq_front"] = sum(
                abs(r[f"{arm}_guarded_kcal_mol"] - r["mfe_kcal_mol"]) > 1e-6 for r in rows
            )
            out["gaps"] = [round(r[f"{arm}_gap_kcal_mol"], 2) for r in rows]
        out["lambda_mid"] = [float(f"{r['lambda_mid']:.4g}") for r in rows]
        out["E"] = [rows[0]["mfe_kcal_mol"], rows[-1]["mfe_kcal_mol"]]
        out["L"] = [round(rows[0]["cai_log"], 3), round(rows[-1]["cai_log"], 3)]
        print(json.dumps(out))


if __name__ == "__main__":
    main()
