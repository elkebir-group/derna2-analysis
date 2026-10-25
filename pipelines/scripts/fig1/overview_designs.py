"""The three designs drawn in the overview figure, with their structures.

    python overview_designs.py RESULTS_DIR MRNAFOLD_JSONL MOTIF_NAME MOTIF_SITE FOLD_SEQ WORK_DIR \
        OUT_JSON

A weighted sum returns only hull vertices. The figure compares DERNA2 with mRNAfold's designs
(CDS-UTR pairs discouraged, one per hull weight of the UTR-masking front) and with the hull
vertices of the concatenation front (the front without UTRs, a proxy for LinearDesign,
LinearCDSfold and DERNA). The forbidden motif set M is one restriction site (MOTIF_SITE, DNA or
RNA letters) and its reverse complement. The designs:
  d2  mRNAfold's design that contains a motif of M, with the CAI closest to CAI_TARGET;
  d1  the concatenation hull vertex closest in CAI to d2 with another CDS;
  d3  the modular design closest in CAI to d2.
d1 keeps the free structure of the deposit, and d2 the free structure of the UTR-masking design
with the same sequence. The modular front was found by window walks that report energies only,
so d3 is refolded here with the pinned fold_seq.

Fails unless no modular design contains a motif of M (so the deposited modular front is also the
motif-free one), d2 is also a UTR-masking design with the same free MFE, d3 refolds to its
deposited energy with no CDS-UTR base pair, and d1 and d2 are not modular.
"""

import json
import math
import subprocess
import sys
from pathlib import Path

from derna2_analysis.io import read_jsonl

CAI_TARGET = 0.87
L5, LC = 50, 234   # hbb 5' UTR and the CDS of GLPE_HUMAN (78 aa)


def cai(r):
    return math.exp(r["cai_log"] / 78)


def cds(seq):
    return seq[L5:L5 + LC]


def cds_utr_pairs(structure):
    stack, n = [], 0
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            n += (L5 <= i < L5 + LC) != (L5 <= j < L5 + LC)
    return n


def main():
    results, mrnafold, name, site, fold_seq, work, out = sys.argv[1:]
    site = site.upper().replace("T", "U")
    forbidden = {site, site[::-1].translate(str.maketrans("ACGU", "UGCA"))}

    def has(r):
        return any(m in cds(r["sequence"]) for m in forbidden)

    res = Path(results)
    b = read_jsonl(res / "b_cdsonly_with_utrs.scored.jsonl")
    cd = read_jsonl(res / "cd_symmetric_guarded_front.scored.jsonl")
    so = read_jsonl(res / "d_so_front.jsonl")
    hull_b = {r["sequence"] for r in read_jsonl(res / "a_cdsonly_front.jsonl")
              if r.get("hull_vertex")}
    mf = [dict(r, cai_log=r["theirs_cai_log"]) for r in read_jsonl(mrnafold)
          if r["arm"] == "guard"]
    if any(has(r) for r in so):
        sys.exit(f"a modular design contains {name}: the deposited front is not motif-free")
    d2 = min((r for r in mf if has(r)), key=lambda r: abs(cai(r) - CAI_TARGET))
    twin = [r for r in cd if r["sequence"] == d2["sequence"]]
    if not twin or twin[0]["free_mfe_kcal_mol"] != d2["free_mfe_kcal_mol"]:
        sys.exit("mRNAfold's design is not a UTR-masking design with the same free MFE")
    d2 = dict(d2, free_structure=twin[0]["free_structure"])
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
    designs["motif"], designs["motif_name"] = site, name
    Path(out).write_text(json.dumps(designs, indent=1) + "\n")
    for k in ("d1", "d2", "d3"):
        d = designs[k]
        print(f"{k}: CAI {d['cai']:.4f}, MFE {d['free_mfe_kcal_mol']}, "
              f"{cds_utr_pairs(d['free_structure'])} CDS-UTR pairs")


if __name__ == "__main__":
    main()
