"""Compare a re-run front with the deposited one, energy by energy at each log-CAI level.

    python check_front.py RUN_FRONT DEPOSITED_FRONT REPORT

Equal-energy designs may be spelled differently, so the comparison is on (cai_log, energy)
and the report counts differing spellings separately. Fails on any energy difference or
any level present in one front only.
"""

import sys

from derna2_analysis.io import read_jsonl


def levels(rows):
    """Each design's energy and spelling by log-CAI level. A `derna` result names the energy
    `mfe` and opens with a header record; the deposited fronts name it `mfe_kcal_mol`."""
    return {
        round(r["cai_log"], 4): (r.get("mfe_kcal_mol", r.get("mfe")), r.get("sequence"))
        for r in rows
        if r.get("record") != "header"
    }


def main():
    run, deposited, report = sys.argv[1:]
    got, want = levels(read_jsonl(run)), levels(read_jsonl(deposited))
    shared = set(got) & set(want)
    energy = [k for k in shared if abs(got[k][0] - want[k][0]) > 1e-6]
    spelled = sum(got[k][1] != want[k][1] for k in shared)
    lines = [
        f"run {len(got)} points, deposited {len(want)}",
        f"only in run {len(set(got) - set(want))}, only deposited {len(set(want) - set(got))}",
        f"energy differs at {len(energy)} levels, spelled differently at {spelled}",
    ]
    with open(report, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if energy or set(got) ^ set(want):
        sys.exit(f"{run} does not reproduce {deposited}")


if __name__ == "__main__":
    main()
