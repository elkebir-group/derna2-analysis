"""The three designs drawn in the overview figure, with their structures.

    python overview_designs.py RESULTS_DIR FOLD_SEQ WORK_DIR OUT_JSON

A weighted sum returns only hull vertices, so the two other methods are represented by the
hull vertices of their fronts: UTR masking (a proxy for mRNAfold) and concatenation (the
front without UTRs, a proxy for LinearDesign, LinearCDSfold and DERNA). The designs:
  d2  the UTR-masking hull vertex with the CAI closest to CAI_TARGET;
  d1  the concatenation hull vertex closest in CAI to d2 with another CDS;
  d3  the modular design closest in CAI to d2.
d1 and d2 keep the free structures of the deposit. The modular front was found by window
walks that report energies only, so d3 is refolded here with the pinned fold_seq.

Fails unless d3 refolds to its deposited energy with no CDS-UTR base pair and d1 and d2 are
not modular.
"""

import json
import math
import subprocess
import sys
from pathlib import Path

from derna2_analysis.io import read_jsonl

CAI_TARGET = 0.987
L5, LC = 50, 234   # hbb 5' UTR and the CDS of GLPE_HUMAN (78 aa)


def cai(r):
    return math.exp(r["cai_log"] / 78)


def cds_utr_pairs(structure):
    stack, n = [], 0
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            n += (L5 <= i < L5 + LC) != (L5 <= j < L5 + LC)
    return n


def lower_hull(pts):
    h = []
    for q in sorted(pts):
        while len(h) >= 2 and ((h[-1][0] - h[-2][0]) * (q[1] - h[-2][1])
                               - (h[-1][1] - h[-2][1]) * (q[0] - h[-2][0])) <= 0:
            h.pop()
        h.append(q)
    return set(h)


def main():
    results, fold_seq, work, out = sys.argv[1:]
    res = Path(results)
    b = read_jsonl(res / "b_cdsonly_with_utrs.scored.jsonl")
    cd = read_jsonl(res / "cd_symmetric_guarded_front.scored.jsonl")
    so = read_jsonl(res / "d_so_front.jsonl")
    hull_b = {r["sequence"] for r in read_jsonl(res / "a_cdsonly_front.jsonl")
              if r.get("hull_vertex")}
    hull_cd = lower_hull([(r["cai_log"], r["mfe_kcal_mol"]) for r in cd])
    masking = [r for r in cd if (r["cai_log"], r["mfe_kcal_mol"]) in hull_cd]
    d2 = min(masking, key=lambda r: abs(cai(r) - CAI_TARGET))
    concatenation = [r for r in b
                     if r["cds_sequence"] in hull_b and r["sequence"] != d2["sequence"]]
    d1 = min(concatenation, key=lambda r: abs(cai(r) - cai(d2)))
    d3 = min(so, key=lambda r: abs(cai(r) - cai(d2)))
    if d1["symmetric_so"] or d2["symmetric_so"]:
        sys.exit("d1 or d2 is modular")
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    (work / "d3.tsv").write_text(f"d3\t{d3['sequence']}\n")
    subprocess.run([fold_seq, "--in", work / "d3.tsv", "--out", work / "d3.jsonl"], check=True,
                   stdout=subprocess.DEVNULL)
    fold = read_jsonl(work / "d3.jsonl")[0]
    if abs(fold["mfe_kcal_mol"] - d3["mfe_kcal_mol"]) > 1e-6:
        sys.exit(f"d3 refolds to {fold['mfe_kcal_mol']}, deposited {d3['mfe_kcal_mol']}")
    if cds_utr_pairs(fold["structure"]):
        sys.exit("d3's refolded structure has a CDS-UTR base pair")
    d3 = dict(d3, free_structure=fold["structure"], free_mfe_kcal_mol=fold["mfe_kcal_mol"])
    designs = {k: dict(d, cai=cai(d)) for k, d in (("d1", d1), ("d2", d2), ("d3", d3))}
    Path(out).write_text(json.dumps(designs, indent=1) + "\n")
    for k, d in designs.items():
        print(f"{k}: CAI {d['cai']:.4f}, MFE {d['free_mfe_kcal_mol']}, "
              f"{cds_utr_pairs(d['free_structure'])} CDS-UTR pairs")


if __name__ == "__main__":
    main()
