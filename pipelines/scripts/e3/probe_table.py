"""E3's hull probe: one row per rung, and the fits the plan reads from them.

    python probe_table.py OUT_TSV OUT_FIT PROBE_DIR [PROBE_DIR ...]

A rung is one `motif_hull --arms walk --workers 1` run on a P31645 prefix. Its exit comes
from the probe's own summary.csv, its wall and peak RSS from the measurement record
(`<label>.ic.json`), and its search counts from the walk's rows: hull vertices and probes
from the `walk-total` row, folds and seconds from the per-probe rows. A rung with no summary
line is reported as `unrecorded` (still running, or lost), never as a failure.

The cost of one fold is reported three ways because the plan's finding rests on the split:
over every probe, over the λ = 10⁻⁴ critical-path probe alone, and over the ordinary probes
at λ ≥ 0.01. Walls are load-qualified lower bounds (another session builds on these cores),
so the load at launch travels with each row; peak RSS is the primary readout.
"""

import csv
import json
import math
import re
import sys
from pathlib import Path

import pandas as pd

KB_PER_GB = 1048576
LONG_ARM = (984, 1629)


def exits(probe_dir):
    """label -> exit code, from the probe's summary.csv. Round one keyed rungs by length
    (`rung`), round two by label; a stray line from a failed field capture is skipped."""
    out = {}
    path = probe_dir / "summary.csv"
    if not path.exists():
        return out
    for row in csv.DictReader(path.open()):
        label = row.get("label") or f"n{row.get('rung')}"
        if (probe_dir / f"{label}.jsonl").exists() and row.get("exit", "").isdigit():
            out[label] = int(row["exit"])
    return out


def fold_cost(probes):
    folds = sum(p["folds"] for p in probes)
    return sum(p["s"] for p in probes) / folds if folds else None


def rung(probe_dir, label, exit_code):
    lines = (probe_dir / f"{label}.jsonl").read_text().splitlines()
    rows = [json.loads(line) for line in lines if line.strip()]
    probes = [r for r in rows if r.get("arm") == "walk" and "folds" in r]
    total = next((r for r in rows if r.get("arm") == "walk-total"), {})
    ic_path = probe_dir / f"{label}.ic.json"
    ic = json.loads(ic_path.read_text()) if ic_path.exists() else {}
    status = "unrecorded" if exit_code is None else "ok" if exit_code == 0 else "limit"
    return {
        "probe": probe_dir.name,
        "label": label,
        "n": total.get("n") or int(re.sub(r"\D", "", label)),
        "status": status,
        "exit": exit_code,
        "wall_s": ic.get("wall_s") if status == "ok" else None,
        "peak_rss_gb": ic["maxrss_kb"] / KB_PER_GB if status == "ok" else None,
        "vertices": total.get("vertices"),
        "probes": total.get("probes"),
        "folds_per_probe": sum(p["folds"] for p in probes) / len(probes) if probes else None,
        "s_per_fold": fold_cost(probes),
        "s_per_fold_lambda_1e-4": fold_cost([p for p in probes if p["lambda"] == 0.0001]),
        "s_per_fold_ordinary": fold_cost([p for p in probes if p["lambda"] >= 0.01]),
    }


def loads(probe_dir):
    """label -> load average at launch, from summary.csv."""
    out = {}
    path = probe_dir / "summary.csv"
    for row in csv.DictReader(path.open()) if path.exists() else []:
        label = row.get("label") or f"n{row.get('rung')}"
        out[label] = row.get("load_at_launch")
    return out


def loglog(xs, ys):
    lx, ly = [math.log(x) for x in xs], [math.log(y) for y in ys]
    mx, my = sum(lx) / len(lx), sum(ly) / len(ly)
    b = sum((x - mx) * (y - my) for x, y in zip(lx, ly, strict=True)) / sum(
        (x - mx) ** 2 for x in lx
    )
    return b, my - b * mx


def fits(df):
    # One row per length. Where a length ran twice, the later round's run is the one launched
    # at a load comparable to its neighbours (the control), so the fit takes that one.
    done = df[df.status == "ok"].sort_values(["n", "probe"], ascending=[True, False])
    first = done.drop_duplicates("n")
    lines = [f"completed rungs: {list(first.n)}"]
    b, a = loglog(list(first.n), list(first.peak_rss_gb))
    lines.append(f"peak RSS: exponent {b:.2f}")
    for n in LONG_ARM:
        lines.append(f"    n = {n}: {math.exp(a + b * math.log(n)):.1f} GB (extrapolated)")
    steps = zip(first.itertuples(), first.iloc[1:].itertuples(), strict=False)
    lines.append(
        "wall, consecutive-rung exponents (load-qualified lower bounds): "
        + ", ".join(
            f"{p.n}->{q.n} {math.log(q.wall_s / p.wall_s) / math.log(q.n / p.n):.2f}"
            for p, q in steps
        )
    )
    censored = df[df.status != "ok"]
    for r in censored.itertuples():
        lines.append(f"not completed: {r.probe}/{r.label} (n = {r.n}), {r.status}")
    return "\n".join(lines) + "\n"


def main():
    out_tsv, out_fit, *dirs = sys.argv[1:]
    records = []
    for d in map(Path, dirs):
        codes, load = exits(d), loads(d)
        for jsonl in sorted(d.glob("*.jsonl")):
            label = jsonl.stem
            if not label.startswith(("n", "control")):
                continue
            records.append({**rung(d, label, codes.get(label)), "load": load.get(label)})
    df = pd.DataFrame(records).sort_values(["n", "probe", "label"])
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False, float_format="%.4g")
    text = fits(df)
    Path(out_fit).write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
