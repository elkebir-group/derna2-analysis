#!/usr/bin/env python3
"""Markdown tables for the README from the run outputs.

Usage: make_tables.py ROOT MRNAFOLD_COMPARE > tables.md
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ARMS = ("symmetric", "no_clamp", "unpaired")


def lines(path: Path) -> list[dict]:
    return [json.loads(x) for x in path.read_text().splitlines() if x.strip()] if path.exists() else []


def fronts(root: Path) -> str:
    out = ["| front | guard | points | hull vertices | SO points | SO hull vertices (rank 1) | wall s |",
           "|---|---|---|---|---|---|---|"]
    for arm, label in [("unguarded", "1"), *[(a, "2") for a in ARMS]]:
        report = json.loads((root / arm / "report.json").read_text())
        pts = lines(root / arm / "front.scored.jsonl")
        vs = lines(root / arm / "vertices.scored.jsonl")
        so = "–" if arm == "unguarded" else sum(p[f"{arm}_so"] for p in pts)
        so_v = "–" if arm == "unguarded" else sum(v[f"{arm}_so"] for v in vs)
        guard = ",".join(report["region"]) or "none"
        out.append(f"| {label} | {arm} (`{guard}`) | {len(pts)} | {len(vs)} | {so} | {so_v} | {report['wall_s']:.0f} |")
    return "\n".join(out)


def front3(root: Path, arm: str) -> str:
    rows = lines(root / arm / "front3.jsonl")
    scored = {r["vertex"]: r for r in lines(root / arm / "front3_designs.scored.jsonl")}
    out = [f"**{arm}**", "",
           "| vertex | λ | vertex E | found | rank | E = F | ΔE vs vertex | walk s | peak GB | start load |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        s = scored.get(r["vertex"])
        load_file = root / arm / "front3" / f"v{r['vertex']}.load"
        load = load_file.read_text().split()[0] if load_file.exists() else "–"
        e = f"{s[f'{arm}_guarded_kcal_mol']:.2f}" if s else "–"
        de = f"{s[f'{arm}_guarded_kcal_mol'] - r['vertex_mfe_kcal_mol']:+.2f}" if s else "–"
        found = "yes" if r["finished"] else f"no ({r['stopped']})"
        out.append(f"| {r['vertex']} | {r['lambda']:.4g} | {r['vertex_mfe_kcal_mol']:.2f} | {found} | {r['rank']} "
                   f"| {e} | {de} | {r['walk_s']:.1f} | {r['peak_rss_kb'] / 1e6:.1f} | {load} |")
    return "\n".join(out)


def mrnafold(compare: Path) -> str:
    rows = lines(compare)
    out = ["| arm | counterpart | runs | objective = ours | same sequence | max Δobj | SO symmetric | SO no_clamp | SO 5′ unpaired | median s |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for arm in ("plain", "guard", "unpaired5", "unpaired_both"):
        sub = [r for r in rows if r["arm"] == arm]
        walls = sorted(r["wall_s"] for r in sub)
        out.append(f"| {arm} | {sub[0]['counterpart']} | {len(sub)} | {sum(abs(r['d_obj']) < 1e-3 for r in sub)} "
                   f"| {sum(r['same_sequence'] for r in sub)} | {max(r['d_obj'] for r in sub):.3f} "
                   f"| {sum(r['symmetric_so'] for r in sub)} | {sum(r['no_clamp_so'] for r in sub)} "
                   f"| {sum(r['unpaired_so'] for r in sub)} | {walls[len(walls) // 2]:.1f} |")
    return "\n".join(out)


def main() -> None:
    root, compare = Path(sys.argv[1]), Path(sys.argv[2])
    print("## Fronts\n")
    print(fronts(root))
    print("\n## Front 3\n")
    for arm in ARMS:
        print(front3(root, arm))
        print()
    print("## mRNAfold\n")
    print(mrnafold(compare))


if __name__ == "__main__":
    main()
