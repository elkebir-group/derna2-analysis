#!/usr/bin/env python3
"""Collect the front-3 walks into one file per arm.

Reads ROOT/<arm>/front3/v<k>.json (one `lazy_kbest_probe` line each) beside
ROOT/<arm>/vertices.jsonl and writes ROOT/<arm>/front3.jsonl: per vertex the
weight, whether the walk found a structurally optimal design, its rank, its
objective and sequence (or the rank reached and the last value when it did
not), walk and verify seconds, and the stop reason; and
ROOT/<arm>/front3_designs.jsonl, the finished rows alone, for
`glpe_hbb_rescore`.

Usage: collect_front3.py ROOT [ARM ...]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def main() -> None:
    root = Path(sys.argv[1])
    arms = sys.argv[2:] or ["symmetric", "no_clamp", "unpaired"]
    for arm in arms:
        vertices = [json.loads(x) for x in (root / arm / "vertices.jsonl").read_text().splitlines()]
        rows = []
        for v in vertices:
            path = root / arm / "front3" / f"v{v['vertex']}.json"
            if not path.exists():
                continue
            r = json.loads(path.read_text().splitlines()[-1])
            finished = r["honest_sequence"] is not None
            row = {
                "arm": arm, "vertex": v["vertex"], "lambda": r["lambda"],
                "vertex_mfe_kcal_mol": v["mfe_kcal_mol"], "vertex_cai_log": v["cai_log"],
                "finished": finished,
                "rank": r["honest_rank"] if finished else r["designs"],
                "objective": r["honest_objective"] if finished else r["last_value"],
                "optimum": r["optimum"],
                "stopped": r["stopped"],
                "build_s": r["build_s"], "walk_s": r["walk_s"], "verify_s": r["verify_s"],
                "peak_rss_kb": r["peak_rss_kb"],
            }
            if finished:
                row["sequence"] = r["honest_sequence"]
            rows.append(row)
        with (root / arm / "front3.jsonl").open("w") as f:
            for row in rows:
                f.write(json.dumps(row) + "\n")
        # The finished designs alone, for the rescorer, which needs a sequence.
        with (root / arm / "front3_designs.jsonl").open("w") as f:
            for row in rows:
                if row["finished"]:
                    f.write(json.dumps(row) + "\n")
        done = sum(r["finished"] for r in rows)
        print(f"{arm}: {len(rows)}/{len(vertices)} walks, {done} finished")


if __name__ == "__main__":
    main()
