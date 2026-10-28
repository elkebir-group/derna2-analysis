"""Assemble the structurally optimal front from its window reports and check it against the
deposited front.

    python so_front.py WINDOWS_DIR DEPOSITED_FRONT OUT_JSONL

Fails unless the windows settle the whole log-CAI range below the most adapted seed, and the
assembled front equals the deposited one point for point (cai_log to 1e-6, energy,
sequence). The stretch above that seed is settled by exhaustive enumeration instead
(top_exhaustive/ in the deposit), not by these windows.
"""

import sys

from derna2_analysis.fig1 import so_front, uncovered
from derna2_analysis.io import read_jsonl, write_jsonl


def key(p):
    return (round(p["cai_log"], 6), p["mfe_kcal_mol"], p["sequence"])


def main():
    windows_dir, deposited, out = sys.argv[1:]
    front, reports = so_front(windows_dir)
    windows = [r for r in reports if not r["report"].startswith("top/")]
    lo, hi = min(r["lo"] for r in windows), max(r["hi"] for r in windows)
    gaps = uncovered(windows, lo, hi)
    if gaps:
        sys.exit(f"cai_log ranges no window settled: {gaps}")
    want = {key(p) for p in read_jsonl(deposited)}
    got = {key(p) for p in front}
    if got != want:
        sys.exit(
            f"assembled front differs from {deposited}: "
            f"{len(got - want)} extra, {len(want - got)} missing"
        )
    write_jsonl(out, front)
    print(
        f"{len(front)} points from {len(reports)} window reports, equal to {deposited}; "
        f"windows settle cai_log [{lo:.4f}, {hi:.4f}]"
    )


if __name__ == "__main__":
    main()
