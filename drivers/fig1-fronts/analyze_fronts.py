#!/usr/bin/env python3
"""Check each front and write its hull vertices with their weight intervals.

For every arm (front 1 unguarded, front 2 per guard) this reads the walk's
`report.json` and `front.jsonl` and checks:

- no point's structure has a pair the arm's guard forbids (classes read from
  the region bounds: 5′ UTR, CDS, 3′ UTR);
- the guarded fronts differ from the unguarded one;
- the hull recomputed here in log-CAI matches the report's vertices.

It then writes `vertices.jsonl` per arm: each hull vertex with its weight
interval (tie weights with its neighbours, clipped to [LO, HI]) and the
interval's midpoint, the weight a front-3 walk runs at.

The scalarization is `λ·E − (1 − λ)·L` with `E` in kcal/mol and `L` the
accumulated log-CAI, so the tie between `a` (more stable) and `b` (more
adapted) is `ΔL / (ΔE + ΔL)`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

LO, HI = 1e-4, 0.999
CLASS_NAMES = {
    ("5", "5"): "5-5",
    ("5", "cds"): "5-cds",
    ("5", "3"): "5-3",
    ("cds", "cds"): "cds-cds",
    ("cds", "3"): "cds-3",
    ("3", "3"): "3-3",
}


def region(i: int, len5: int, cds_end: int) -> str:
    return "5" if i < len5 else "cds" if i < cds_end else "3"


def pairs(structure: str) -> list[tuple[int, int]]:
    stack, out = [], []
    for j, ch in enumerate(structure):
        if ch == "(":
            stack.append(j)
        elif ch == ")":
            out.append((stack.pop(), j))
    assert not stack, "unbalanced structure"
    return out


def classes(structure: str, len5: int, cds_end: int) -> dict[str, int]:
    count: dict[str, int] = {}
    for i, j in pairs(structure):
        c = CLASS_NAMES[(region(i, len5, cds_end), region(j, len5, cds_end))]
        count[c] = count.get(c, 0) + 1
    return count


def tie(a: dict, b: dict) -> float:
    de = b["mfe_kcal_mol"] - a["mfe_kcal_mol"]
    dl = b["cai_log"] - a["cai_log"]
    return dl / (de + dl)


def hull(points: list[dict]) -> list[dict]:
    """Points that minimize `λE − (1 − λ)L` for some λ in (LO, HI), most
    stable first, ties in both coordinates collapsed to one point."""
    distinct: dict[tuple[float, float], dict] = {}
    for p in points:
        distinct.setdefault((round(p["mfe_kcal_mol"], 2), round(p["cai_log"], 6)), p)
    pts = sorted(distinct.values(), key=lambda p: (p["mfe_kcal_mol"], -p["cai_log"]))
    # Keep only the Pareto set: strictly increasing E, strictly increasing L.
    front: list[dict] = []
    for p in pts:
        if not front or p["cai_log"] > front[-1]["cai_log"] + 1e-9:
            front.append(p)
    chain: list[dict] = []
    for p in front:
        # Pop while the last vertex is not strictly better than the chord.
        while len(chain) >= 2 and tie(chain[-2], chain[-1]) <= tie(chain[-1], p) + 1e-12:
            chain.pop()
        chain.append(p)
    # Clip to the weight range: drop vertices optimal only outside (LO, HI).
    lam_hi = [HI] + [tie(chain[k], chain[k + 1]) for k in range(len(chain) - 1)]
    lam_lo = lam_hi[1:] + [LO]
    return [
        {**v, "lambda_hi": min(h, HI), "lambda_lo": max(lo, LO)}
        for v, h, lo in zip(chain, lam_hi, lam_lo)
        if min(h, HI) > max(lo, LO)
    ]


def load(arm_dir: Path) -> tuple[dict, list[dict]]:
    report = json.loads((arm_dir / "report.json").read_text())
    front = [json.loads(x) for x in (arm_dir / "front.jsonl").read_text().splitlines() if x.strip()]
    return report, front


def check_arm(name: str, arm_dir: Path, unguarded: set | None) -> dict:
    report, front = load(arm_dir)
    len5, cds_end = report["utr5_nt"], report["utr5_nt"] + 3 * report["residues"]
    forbidden = set(report["region"])
    bad = [
        (k, c)
        for k, p in enumerate(front)
        for c in classes(p["structure"], len5, cds_end)
        if c in forbidden
    ]
    vs = hull(front)
    reported = [(round(v["mfe_kcal_mol"], 2), round(v["adaptation_log"], 4)) for v in report["vertices"]]
    ours = [(round(v["mfe_kcal_mol"], 2), round(v["cai_log"], 4)) for v in vs]
    key = {(round(p["mfe_kcal_mol"], 2), round(p["cai_log"], 6)) for p in front}
    with (arm_dir / "vertices.jsonl").open("w") as f:
        for k, v in enumerate(vs):
            row = {
                "vertex": k,
                "mfe_kcal_mol": v["mfe_kcal_mol"],
                "cai_log": v["cai_log"],
                "cai": v["cai"],
                "lambda_lo": v["lambda_lo"],
                "lambda_hi": v["lambda_hi"],
                "lambda_mid": 0.5 * (v["lambda_lo"] + v["lambda_hi"]),
                "sequence": v["sequence"],
                "structure": v["structure"],
            }
            f.write(json.dumps(row) + "\n")
    summary = {
        "arm": name,
        "region": sorted(forbidden),
        "points": len(front),
        "hull_vertices": len(vs),
        "reported_vertices": len(reported),
        "hull_matches_report": ours == reported,
        "forbidden_pairs": len(bad),
        "differs_from_unguarded": None if unguarded is None else key != unguarded,
        "shared_with_unguarded": None if unguarded is None else len(key & unguarded),
        "wall_s": report["wall_s"],
    }
    if ours != reported:
        summary["hull_diff"] = {"ours_only": sorted(set(ours) - set(reported)),
                                "report_only": sorted(set(reported) - set(ours))}
    return summary, key


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--arms", default="unguarded,symmetric,no_clamp,unpaired")
    args = ap.parse_args()
    unguarded = None
    for name in args.arms.split(","):
        if not (args.root / name / "front.jsonl").exists():
            continue
        summary, key = check_arm(name, args.root / name, unguarded)
        if name == "unguarded":
            unguarded = key
        print(json.dumps(summary))


if __name__ == "__main__":
    main()
