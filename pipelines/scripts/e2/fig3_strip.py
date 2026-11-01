"""The E2 constraint strip, rebuilt from the current screen and grid.

    python fig3_strip.py E2_DIR OUT_JSON

Panel (a) counts proteins, not designs: for each blacklist entry, how many of the screened
proteins the lattice admits it on, against how many it actually occurs on in a design of
that protein's unconstrained front. An entry admitted widely and never spelled forbade
nothing there.

The price is dMFE at matched adaptation, read against the unconstrained front of the same
protein. Coverage is the pair the experiment retires a cell on — `share <= 0.80 AND
hull_gap <= 0.20` per CELL-SET.md — recomputed here from each cell's open regions, which
reproduces the long band's own recorded values.

Hull vertices are computed from each front's points, not read from its `on_hull` flag:
the flag is what the search accepted at its own tolerance, and it can undercount.
"""

import csv
import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

from derna2_analysis.fronts import hull_vertices

SETS = ["manufacturing", "processing", "repeats"]
NESTS = ["F1", "F2", "F3"]
GATE = {"share": 0.80, "gap": 0.20}
INSTANCE = "P47211"


def read_front(path):
    lines = Path(path).read_text().splitlines()
    head = json.loads(lines[0])
    pts = []
    for line in lines[1:]:
        if not line.strip():
            continue
        r = json.loads(line)
        if r.get("record") == "point":
            pts.append((float(r["cai_log"]), float(r["mfe"])))
    pts.sort()
    return head, pts


def hull_xy(pts):
    """A front's hull vertices, by increasing cai_log, rounded the way the strip draws them."""
    return [[round(a, 2), e] for a, e in sorted(hull_vertices(pts))]


def coverage(head, pts):
    """(unsettled share of the energy span, largest open region as a share of it)."""
    mfes = [p[1] for p in pts]
    span = (max(mfes) - min(mfes)) if mfes else 0.0
    widths = [abs(o["mfe_hi"] - o["mfe_lo"]) for o in (head.get("open") or [])]
    if not span:
        return 1.0, 1.0, span
    return min(sum(widths) / span, 1.0), min(max(widths, default=0) / span, 1.0), span


def best_round(root, protein, nest):
    best = None
    for path in Path(root).glob(f"{protein}_{nest}*.jsonl"):
        head, pts = read_front(path)
        if best is None or len(pts) > len(best[1]):
            best = (head, pts)
    return best


def price_curve(base, pts):
    """At each adaptation the constrained front attains, how much energy it gives up."""
    xy, neg = [], 0
    for cai, mfe in pts:
        cands = [b[1] for b in base if b[0] >= cai - 1e-9]
        if not cands:
            continue
        d = round(mfe - min(cands), 2)
        if d < -1e-9:
            neg += 1
        xy.append([round(cai, 2), d])
    vals = sorted(d for _, d in xy)
    q = st.quantiles(vals, n=4) if len(vals) > 1 else [vals[0], vals[0], vals[0]]
    return {
        "xy": xy[:: max(1, len(xy) // 320)],
        "lo": vals[0],
        "q1": q[0],
        "med": st.median(vals),
        "q3": q[2],
        "hi": vals[-1],
        "n": len(vals),
        "neg": neg,
    }


def main():
    root = Path(sys.argv[1])

    rows = [
        json.loads(x)
        for x in (root / "screen" / "binding.jsonl").read_text().splitlines()
        if x.strip()
    ]
    proteins = {r["protein"] for r in rows}
    agg = defaultdict(lambda: {"occ": set(), "sp": set()})
    for r in rows:
        if r["set"] not in SETS:
            continue
        key = (r["set"], r["motif"])
        if r["designs_carrying"]:
            agg[key]["occ"].add(r["protein"])
        if r["spellable"]:
            agg[key]["sp"].add(r["protein"])
    entries = sorted(
        # The screen records DNA; the figure is about the transcript, so entries are shown
        # as RNA, which is also how the paper's blacklist is written.
        (
            {"set": s, "motif": m.replace("T", "U"), "occ": len(v["occ"]), "sp": len(v["sp"])}
            for (s, m), v in agg.items()
        ),
        key=lambda e: (SETS.index(e["set"]), -e["occ"]),
    )

    base_head, base = read_front(root / "screen" / "f0" / f"{INSTANCE}.jsonl")
    meta = [
        {
            "lab": "F0",
            "sets": "unconstrained",
            "n": len(base),
            "hull": len(hull_vertices(base)),
            "complete": base_head.get("stop") is None,
        }
    ]
    cloud = {"F0": [[round(p[0], 2), p[1]] for p in base]}
    hull = {"F0": hull_xy(base)}
    price, kept = {}, {}
    for nest in NESTS:
        got = best_round(root / "grid", INSTANCE, nest)
        if not got:
            continue
        head, pts = got
        meta.append(
            {
                "lab": nest,
                "sets": SETS[NESTS.index(nest)],
                "n": len(pts),
                "hull": len(hull_vertices(pts)),
                "complete": head.get("complete") is True,
                "build": head.get("command", [""])[0]
                if isinstance(head.get("command"), list)
                else None,
            }
        )
        hull[nest] = hull_xy(pts)
        if nest == "F3":
            cloud["F3"] = [[round(p[0], 2), p[1]] for p in pts]
        price[f"F0->{nest}"] = price_curve(base, pts)
        kept[nest] = pts

    # Each step against the one before it: this is where a set that forbids nothing on
    # this instance shows as a curve lying on zero.
    adjacent = {}
    prev_lab, prev_pts = "F0", base
    for nest in NESTS:
        if nest in kept:
            adjacent[f"{prev_lab}->{nest}"] = price_curve(prev_pts, kept[nest])
            prev_lab, prev_pts = nest, kept[nest]

    track, counts = [], defaultdict(int)
    for band, folder, header in [
        ("short", root / "grid", True),
        ("middle", root / "mid", False),
        ("long", root / "long", False),
    ]:
        summary = folder / "summary.csv"
        seen = set()
        for row in csv.reader(open(summary)):
            if not row or (header and row[0] == "protein"):
                continue
            protein, nest = row[0], row[2]
            if (protein, nest) in seen:
                continue
            seen.add((protein, nest))
            cell = folder / "cells" / f"{protein}_{nest}.jsonl"
            cand = cell if cell.exists() else None
            if cand is None:
                got = best_round(folder, protein, nest)
                if not got:
                    continue
                head, pts = got
            else:
                head, pts = read_front(cand)
            share, gap, _ = coverage(head, pts)
            status = (
                "complete"
                if head.get("complete") is True
                else "gate passed"
                if share <= GATE["share"] and gap <= GATE["gap"]
                else "budget-bound"
            )
            counts[status] += 1
            track.append(
                {
                    "band": band,
                    "protein": protein,
                    "nest": nest,
                    "share": round(share, 4),
                    "gap": round(gap, 4),
                    "status": status,
                }
            )

    json.dump(
        {
            "entries": entries,
            "nprot": len(proteins),
            "inst": INSTANCE,
            "aa": 349,
            "meta": meta,
            "cloud": cloud,
            "hull": hull,
            "price": price,
            "adjacent": adjacent,
            "gate": GATE,
            "track": track,
            "counts": dict(counts),
        },
        open(sys.argv[2], "w"),
        separators=(",", ":"),
    )

    print(f"{len(entries)} entries over {len(proteins)} proteins; {len(track)} cells")
    print("counts:", dict(counts))
    never = [e for e in entries if e["occ"] == 0]
    print("admitted but never spelled:", [(e["motif"], e["sp"]) for e in never])
    for k, v in {**price, **adjacent}.items():
        print(f"  {k}: med {v['med']:.2f}  max {v['hi']:.1f}  n {v['n']}  negatives {v['neg']}")
    print("meta:", [(m["lab"], m["n"], m["complete"]) for m in meta])


if __name__ == "__main__":
    main()
