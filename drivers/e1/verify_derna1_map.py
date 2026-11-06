#!/usr/bin/env python3
"""Check DERNA v1's lambda map against derna2 on designs, not on arithmetic.

v1 scores `lambda*E_dcal - (1-lambda)*sum(ln w)` and derna2 scores
`lambda*E_dcal - (1-lambda)*100*sum(ln w)`, so the same point needs
`(1-l1)/l1 = 100*(1-l2)/l2`. Arithmetic on one printed objective cannot
establish that: it has to produce the same design, or at least the same
objective value, from the same table.

For each derna2 lambda it runs both tools at the mapped pair and compares the
re-scored MFE and sum(ln w). A map wrong by a transform that agrees at the
corners shows up here and nowhere else, which is why interior lambdas are the
test.

Usage: verify_derna1_map.py FASTA [lambda ...]
"""
import json
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

BIN2 = "run/paper-e1-bin/run_design_per_lambda"
BIN1 = "run/paper-e2-bin/derna1"
TABLE1 = "run/paper-e2-bin/derna1_codon_usage_human.csv"
PARAMS = ("derna2/docs/reference/_cpp/"
          "derna-main/data/InputFiles")
PIN = ["taskset", "-c", "0-26:2"]


def map_lambda(l2):
    """derna2 convex weight -> v1 convex weight."""
    return 1.0 / (1.0 + 100.0 * (1.0 - l2) / l2)


def run_derna2(fasta, lam, work):
    lam_file = work / f"l{lam}.jsonl"
    lam_file.write_text(json.dumps({"lambda": lam}) + "\n")
    out = work / f"d2-{lam}.jsonl"
    subprocess.run(
        PIN + [BIN2, "--fasta", str(fasta), "--lambdas", str(lam_file),
               "--out", str(out), "--machine", str(work / f"m{lam}.json")],
        check=True, capture_output=True, text=True, timeout=1800,
    )
    row = json.loads(out.read_text().splitlines()[0])
    return row


def run_derna1(fasta, lam, work):
    out = work / f"d1-{lam}.txt"
    proc = subprocess.run(
        PIN + [BIN1, "-i", str(fasta), "-o", str(out), "-m", "1", "-s", "2",
               "-l", repr(lam), "-c", TABLE1, "-d", PARAMS],
        check=True, capture_output=True, text=True, timeout=1800,
    )
    # The summary line: "lambda: L,O: O,cai: C,cai_s: S,mfe: M,combined: ..."
    line = [l for l in proc.stdout.splitlines() if l.startswith("lambda:")][-1]
    fields = dict(
        (k.strip(), v.strip())
        for k, v in (p.split(":", 1) for p in line.split(",") if ":" in p)
    )
    seq = ""
    for raw in out.read_text().splitlines():
        if raw.startswith("zuker cai bp") or raw.startswith("bp:"):
            continue
        if re.fullmatch(r"[ACGU]{3,}", raw.strip()):
            seq = raw.strip()
    return {
        "mfe_kcal": float(fields["mfe"]) / 100.0,
        "sum_ln_w": float(fields["cai"]),
        "cai": float(fields["cai_s"]),
        "objective": float(fields["O"]),
        "seq": seq,
    }


def main(fasta, lambdas):
    fasta = Path(fasta)
    work = Path(tempfile.mkdtemp(prefix="d1map-"))
    print(f"fasta {fasta}   work {work}")
    print(f"{'lam2':>8} {'lam1':>12} | {'derna2 mfe':>11} {'v1 mfe':>11} "
          f"{'dMFE':>7} | {'derna2 lnw':>11} {'v1 lnw':>11} {'dlnw':>7} | "
          f"{'derna2 CAI':>10} {'v1 CAI':>8}  same?")
    worst = 0.0
    for l2 in lambdas:
        l1 = map_lambda(l2)
        r2 = run_derna2(fasta, l2, work)
        r1 = run_derna1(fasta, l1, work)
        mfe2 = r2.get("mfe_kcal_mol", r2.get("mfe"))
        lnw2 = r2.get("cai_log", r2.get("sum_ln_w"))
        dm = abs(mfe2 - r1["mfe_kcal"])
        dl = abs(lnw2 - r1["sum_ln_w"])
        worst = max(worst, dm, dl)
        same = "YES" if dm < 0.05 and dl < 0.01 else "no"
        print(f"{l2:8.4f} {l1:12.8f} | {mfe2:11.2f} {r1['mfe_kcal']:11.2f} "
              f"{dm:7.2f} | {lnw2:11.4f} {r1['sum_ln_w']:11.4f} {dl:7.4f} | "
              f"{r2.get('cai', float('nan')):10.4f} {r1['cai']:8.4f}  {same}")
    print()
    print(f"worst disagreement {worst:.4f}")
    print("A matching objective at every interior lambda is what licenses the "
          "map; equal\nsequences are not required, since ties can be spelled "
          "differently.")


if __name__ == "__main__":
    args = sys.argv[1:]
    lams = [float(a) for a in args[1:]] or [0.9, 0.7, 0.5, 0.3, 0.1]
    main(args[0], lams)
