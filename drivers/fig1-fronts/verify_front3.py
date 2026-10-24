#!/usr/bin/env python3
"""Check each front-3 walk against an independent rescoring.

Per walk: the walk's optimum (its rank-1 objective) equals the hull vertex's
objective at the walk's weight; and for a finished walk, the rescored design
is structurally optimal under the arm's guard (guarded E = free F) and its
objective `100·(λE − (1 − λ)L)`, with L recomputed from the sequence, equals
the walk's reported objective.

Usage: verify_front3.py ROOT ARM ...
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

CSV = Path(__file__).resolve().parents[2] / "packages/derna-codon-dp/data/codon_usage_human.csv"
LEN5, CDS_NT, TOL = 50, 3 * 78, 1e-3


def weights() -> dict[str, float]:
    rows = [r.split(",") for r in CSV.read_text().splitlines()[1:] if r.strip()]
    best: dict[str, float] = {}
    for aa, _, f in rows:
        best[aa] = max(best.get(aa, 0.0), float(f))
    return {codon: float(f) / best[aa] for aa, codon, f in rows}


def main() -> None:
    root, w = Path(sys.argv[1]), weights()
    for arm in sys.argv[2:]:
        walks = [json.loads(x) for x in (root / arm / "front3.jsonl").read_text().splitlines()]
        scored_path = root / arm / "front3_designs.scored.jsonl"
        scored = {}
        if scored_path.exists():
            scored = {json.loads(x)["vertex"]: json.loads(x) for x in scored_path.read_text().splitlines()}
        bad = []
        for r in walks:
            lam = r["lambda"]
            vertex = 100 * (lam * r["vertex_mfe_kcal_mol"] - (1 - lam) * r["vertex_cai_log"])
            if abs(r["optimum"] - vertex) > TOL:
                bad.append((r["vertex"], "optimum is not the vertex", r["optimum"], vertex))
            if not r["finished"]:
                continue
            s = scored[r["vertex"]]
            cds = s["sequence"][LEN5:LEN5 + CDS_NT]
            log_cai = sum(math.log(w[cds[i:i + 3]]) for i in range(0, CDS_NT, 3))
            obj = 100 * (lam * s[f"{arm}_guarded_kcal_mol"] - (1 - lam) * log_cai)
            if not s[f"{arm}_so"]:
                bad.append((r["vertex"], "not structurally optimal", s[f"{arm}_gap_kcal_mol"]))
            if abs(obj - r["objective"]) > TOL:
                bad.append((r["vertex"], "objective mismatch", obj, r["objective"]))
        done = sum(r["finished"] for r in walks)
        print(json.dumps({"arm": arm, "walks": len(walks), "finished": done, "failures": bad}))


if __name__ == "__main__":
    main()
