"""The overview figure's structurally optimal front, assembled from its window reports.

The front was found by capped window walks between adjacent known structurally optimal
designs (seeds), plus a ceiling walk above the most adapted seed. Each report's last line
lists the designs it found; the front is the non-dominated set of seeds and finds.
"""

import json
from pathlib import Path


def last_line(path):
    with open(path) as f:
        lines = [line for line in f if line.strip()]
    return json.loads(lines[-1])


def so_front(windows_dir):
    """(front, reports): the front by increasing cai_log, and per report its relative
    path, whether the walk completed, its cai_log window and the parts of it left open."""
    root = Path(windows_dir)
    designs = {}
    for line in (root / "seed_stair.jsonl").read_text().splitlines():
        if line.strip():
            d = json.loads(line)
            designs[d["sequence"]] = (d["cai_log"], d["mfe_kcal_mol"], "seed")
    reports = []
    paths = sorted(root.glob("r*/w*.jsonl")) + sorted(root.glob("top/*.jsonl"))
    for path in sorted(paths, key=str):
        rep = last_line(path)
        rel = str(path.relative_to(root))
        reports.append(
            {
                "report": rel,
                "complete": rep["complete"],
                "lo": rep["lo"],
                "hi": rep["hi"],
                "open": [] if rep["complete"] else rep["open"],
            }
        )
        for q in rep["front"]:
            designs.setdefault(q["sequence"], (q["cai_log"], q["energy_kcal_mol"], rel))
    front, best = [], float("inf")
    ranked = sorted(designs.items(), key=lambda kv: (-kv[1][0], kv[1][1]))
    for seq, (a, e, src) in ranked:
        if e < best - 1e-9:
            front.append({"cai_log": a, "mfe_kcal_mol": e, "sequence": seq, "source": src})
            best = e
    front.reverse()
    return front, reports


def uncovered(reports, lo, hi):
    """The parts of [lo, hi] no report settled. A report settles its window less the
    intervals it lists as open: a walk stopped at its memory cap still proves the stretches
    it closed, and a narrower follow-up window covers the rest."""
    closed = []
    for r in reports:
        cuts = sorted(r["open"])
        start = r["lo"]
        for a, b in cuts:
            if a > start:
                closed.append((start, a))
            start = max(start, b)
        if start < r["hi"]:
            closed.append((start, r["hi"]))
    gaps, at = [], lo
    for a, b in sorted(closed):
        if a > at + 1e-9:
            gaps.append((at, min(a, hi)))
        at = max(at, b)
        if at >= hi:
            break
    if at < hi - 1e-9:
        gaps.append((at, hi))
    return [(a, b) for a, b in gaps if b > a + 1e-9]
