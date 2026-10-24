#!/usr/bin/env python3
"""How many single-core cells this host runs at once without slowing each down.

Runs the same cells at each concurrency and compares each cell's wall against
its own wall when nothing else ran. A cell is the unit the grid schedules, so
the question is answered in the grid's own terms rather than in a synthetic
loop: the arms contend for memory bandwidth exactly as they will on the night.

Each level writes its own `out_dir`, so no level reads another's rows as
already done, and the levels run in ascending order, the serial baseline
first.

    python3 pilot_concurrency.py --config config.e1.tamazirt.json \
        --levels 1,4,8 --max-n 255 --out pilot/

The verdict is per level: the median and worst per-cell wall ratio against the
baseline. A level passes when its median is within `--slack` of 1.
"""

import argparse
import json
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DRIVE = HERE / "drive_e1.py"


def cells_of(path):
    """Each cell row of a sweep file, by (arm, protein)."""
    out = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("kind") == "cell":
            out[(row["arm"], row["protein"])] = row
    return out


def run_level(cfg_path, base, level, max_n, log):
    """One concurrency level: its own output root, hull then sweep."""
    cfg = json.loads(Path(cfg_path).read_text())
    cfg["concurrency"] = level
    out = base / f"c{level}"
    out.mkdir(parents=True, exist_ok=True)
    # The hull is the same work at every level and its walk is the reference
    # for the replay arms, so it is run once and copied: re-walking it would
    # time the walk under contention it will not meet in the grid's own hull
    # pass, and the lambdas must be identical across levels or the cells are
    # not the same cells.
    hull = out / "hull"
    shared = base / "hull"
    if shared.exists() and not hull.exists():
        shutil.copytree(shared, hull)
    cfg["out_dir"] = str(out)
    path = base / f"config.c{level}.json"
    path.write_text(json.dumps(cfg, indent=1) + "\n")

    for command in ("hull", "sweep"):
        argv = [sys.executable, str(DRIVE), command, "--config", str(path),
                "--phase", "single"]
        if max_n:
            argv += ["--max-n", str(max_n)]
        started = time.time()
        subprocess.run(argv, check=False)
        log(f"c{level}: {command} in {(time.time() - started) / 60:.1f} min")
    if not shared.exists():
        shutil.copytree(hull, shared)
    return cells_of(out / "sweep.jsonl")


def compare(baseline, got, level, slack, log):
    """Per-cell wall ratios of one level against the serial baseline."""
    ratios = []
    for key, row in sorted(got.items()):
        ref = baseline.get(key)
        if not ref or row.get("exit") != "ok" or ref.get("exit") != "ok":
            continue
        if not ref.get("wall_s"):
            continue
        ratios.append((row["wall_s"] / ref["wall_s"], key))
    if not ratios:
        log(f"c{level}: no comparable cell")
        return False
    walls = [r for r, _ in ratios]
    median = statistics.median(walls)
    worst, where = max(ratios)
    log(f"c{level}: {len(ratios)} cells, median {median:.2f}x, "
        f"worst {worst:.2f}x ({where[0]} {where[1]})")
    return median <= 1 + slack


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--levels", default="1,4,8")
    ap.add_argument("--max-n", type=int, default=255)
    ap.add_argument("--slack", type=float, default=0.05)
    args = ap.parse_args()
    base = Path(args.out)
    base.mkdir(parents=True, exist_ok=True)
    logfile = open(base / "pilot.log", "a")

    def log(msg):
        line = f"{time.strftime('%H:%M:%S')} {msg}"
        print(line, flush=True)
        print(line, file=logfile, flush=True)

    levels = [int(x) for x in args.levels.split(",")]
    log(f"pilot over {levels} up to {args.max_n} aa, slack {args.slack:.0%}")
    baseline, passed = None, []
    for level in levels:
        got = run_level(args.config, base, level, args.max_n, log)
        if baseline is None:
            baseline = got
            log(f"c{level}: baseline of {len(got)} cells")
            passed.append(level)
            continue
        if compare(baseline, got, level, args.slack, log):
            passed.append(level)
    log(f"within {args.slack:.0%} of the baseline: {passed}")
    log(f"concurrency to run the grid at: {max(passed)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
