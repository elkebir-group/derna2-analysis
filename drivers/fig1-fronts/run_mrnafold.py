#!/usr/bin/env python3
"""mRNAfold at the symmetric front's hull weights, in four arms.

mRNAfold minimizes `E − λ_mf · Σ ln w(codon)` in kcal/mol over the same
Turner 2004 / dangles = 0 / loop-cap 30 model; ours minimizes
`λ·E − (1 − λ)·L` with `L = Σ ln w`, so `λ_mf = (1 − λ)/λ`. Arms:

- `plain`: no constraint;
- `guard`: `discourage_pair` on every 5′↔CDS and CDS↔3′ pair, both
  orientations (the source sets only `pair_bonus[i][j]`);
- `unpaired5`: `encourage_unpaired` over the 5′ UTR;
- `unpaired_both`: `encourage_unpaired` over both UTRs.

Each penalty is soft (+100 kcal/mol on the designed structure). One trace per
run. Writes one JSON line per (arm, vertex) with mRNAfold's sequence,
structure, score, energy recovered from the score, and `L` recomputed with our
codon weights.

`--precheck` instead runs the index-convention check on a short protein with
no flanks: plain, then each orientation of one discouraged pair the plain
fold forms, then both.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CSV = REPO / "packages/derna-codon-dp/data/codon_usage_human.csv"
FASTA = REPO / "experiments/scope-25-beam-runtime-quality/fixtures/GLPE_HUMAN_78.fasta"
UTR5 = (REPO / "experiments/scope-19-utr-support/m5/hbb5.txt").read_text().strip()
UTR3 = (REPO / "experiments/scope-19-utr-support/m5/hbb3.txt").read_text().strip()
ARMS = ("plain", "guard", "unpaired5", "unpaired_both")


def protein() -> str:
    return "".join(x.strip() for x in FASTA.read_text().splitlines() if not x.startswith(">"))


def weights() -> dict[str, float]:
    rows = [r.split(",") for r in CSV.read_text().splitlines()[1:] if r.strip()]
    best: dict[str, float] = {}
    for aa, _, f in rows:
        best[aa] = max(best.get(aa, 0.0), float(f))
    return {codon: float(f) / best[aa] for aa, codon, f in rows}


def pairs(structure: str) -> list[tuple[int, int]]:
    stack, out = [], []
    for j, ch in enumerate(structure):
        if ch == "(":
            stack.append(j)
        elif ch == ")":
            out.append((stack.pop(), j))
    return out


def run(binary: str, cfg: Path, body: list[str]) -> tuple[dict, float]:
    cfg.write_text("\n".join(body) + "\n")
    started = time.monotonic()
    out = subprocess.run([binary, str(cfg)], capture_output=True, text=True, check=True).stdout
    wall = time.monotonic() - started
    lines = [x for x in out.splitlines() if x.strip()]
    rec = {"sequence": lines[0], "structure": lines[1]}
    for x in lines[2:5]:
        k, v = x.split(":")
        rec[k.strip().lower()] = float(v)
    return rec, wall


def arm_lines(arm: str, n: int) -> list[str]:
    len5, cds_end, total = len(UTR5), len(UTR5) + 3 * n, len(UTR5) + 3 * n + len(UTR3)
    if arm == "plain":
        return []
    if arm == "guard":
        lines = []
        for i in range(len5):
            for j in range(len5, cds_end):
                lines += [f"discourage_pair {i},{j}", f"discourage_pair {j},{i}"]
        for i in range(len5, cds_end):
            for j in range(cds_end, total):
                lines += [f"discourage_pair {i},{j}", f"discourage_pair {j},{i}"]
        return lines
    if arm == "unpaired5":
        return [f"encourage_unpaired 0,{len5 - 1}"]
    if arm == "unpaired_both":
        return [f"encourage_unpaired 0,{len5 - 1}", f"encourage_unpaired {cds_end},{total - 1}"]
    raise ValueError(arm)


def precheck(args: argparse.Namespace) -> None:
    aa = protein()[:12]
    base = [f"aa_seq {aa}", "lambda 0.5", "num_subopt_traces 1", f"codon_freq_path {args.table}"]
    plain, _ = run(args.bin, args.work / "pre_plain.cfg", base)
    formed = pairs(plain["structure"])
    print(json.dumps({"case": "plain", **plain, "pairs": formed}))
    i, j = formed[0]
    for label, extra in [
        ("ij", [f"discourage_pair {i},{j}"]),
        ("ji", [f"discourage_pair {j},{i}"]),
        ("both", [f"discourage_pair {i},{j}", f"discourage_pair {j},{i}"]),
    ]:
        rec, _ = run(args.bin, args.work / f"pre_{label}.cfg", base + extra)
        print(json.dumps({"case": label, "pair": [i, j], **rec,
                          "pair_present": (i, j) in pairs(rec["structure"])}))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--table", required=True, type=Path, help="mRNAfold codon table")
    ap.add_argument("--work", required=True, type=Path)
    ap.add_argument("--vertices", type=Path, help="symmetric vertices.jsonl")
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--out", type=Path)
    ap.add_argument("--precheck", action="store_true")
    args = ap.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    if args.precheck:
        precheck(args)
        return
    aa = protein()
    n = len(aa)
    w = weights()
    vs = [json.loads(x) for x in args.vertices.read_text().splitlines() if x.strip()]
    done = set()
    if args.out.exists():
        done = {(r["arm"], r["vertex"]) for r in map(json.loads, args.out.read_text().splitlines())}
    with args.out.open("a") as sink:
        for arm in args.arms.split(","):
            extra = arm_lines(arm, n)
            for v in vs:
                if (arm, v["vertex"]) in done:
                    continue
                lam = v["lambda_mid"]
                lam_mf = (1.0 - lam) / lam
                body = [f"aa_seq {aa}", f"lambda {lam_mf!r}", "num_subopt_traces 1",
                        f"codon_freq_path {args.table}", f"utr_5p {UTR5}", f"utr_3p {UTR3}"]
                rec, wall = run(args.bin, args.work / f"{arm}_v{v['vertex']}.cfg", body + extra)
                cds = rec["sequence"][len(UTR5):len(UTR5) + 3 * n]
                cai_log = sum(math.log(w[cds[i:i + 3]]) for i in range(0, 3 * n, 3))
                row = {
                    "arm": arm, "vertex": v["vertex"], "lambda": lam, "lambda_mf": lam_mf,
                    **rec,
                    "cai_log": cai_log,
                    "score_energy_kcal_mol": rec["score"] + lam_mf * cai_log,
                    "wall_s": wall,
                }
                sink.write(json.dumps(row) + "\n")
                sink.flush()
                print(f"{arm} v{v['vertex']} λ={lam:.5f} wall={wall:.1f}s", flush=True)


if __name__ == "__main__":
    main()
