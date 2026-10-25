"""The three designs drawn in the overview figure, with their structures.

    python overview_designs.py RESULTS_DIR MOTIFS_TXT MOTIF_NAME FOLD_SEQ WORK_DIR OUT_JSON

The forbidden motif set M is the two strands of one restriction site (MOTIF_NAME, e.g. SapI)
from MOTIFS_TXT. The designs:
  d2  the UTR-masking design that contains a motif of M, with the CAI closest to CAI_TARGET;
  d1  the concatenation design (the CDS-only front, refolded with the UTRs) closest in CAI to d2;
  d3  the modular design closest in CAI to d2.
d1 and d2 keep the free structures of the deposit. The modular front was found by window
walks that report energies only, so d3 is refolded here with the pinned fold_seq.

Fails unless: no modular design contains a motif of M (so the deposited modular front is also
the motif-free one); d3 refolds to its deposited energy with no CDS-UTR base pair; d1 and d2
are not modular.
"""

import json
import math
import subprocess
import sys
from pathlib import Path

from derna2_analysis.io import read_jsonl

CAI_TARGET = 0.886
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


def motifs(path, name):
    out = []
    for line in Path(path).read_text().splitlines():
        pattern, _, comment = line.partition("#")
        if pattern.strip() and comment.strip().split(",")[0] == name:
            out.append(pattern.strip().upper().replace("T", "U"))
    if len(out) != 2:
        sys.exit(f"{path}: want the two strands of {name}, found {out}")
    return out


def main():
    results, motif_txt, name, fold_seq, work, out = sys.argv[1:]
    res = Path(results)
    b = read_jsonl(res / "b_cdsonly_with_utrs.scored.jsonl")
    cd = read_jsonl(res / "cd_symmetric_guarded_front.scored.jsonl")
    so = read_jsonl(res / "d_so_front.jsonl")
    forbidden = motifs(motif_txt, name)

    def has(r):
        return any(m in cds(r["sequence"]) for m in forbidden)

    if any(has(r) for r in so):
        sys.exit(f"a modular design contains {name}: the deposited front is not motif-free")
    d2 = min((r for r in cd if has(r)), key=lambda r: abs(cai(r) - CAI_TARGET))
    d1 = min(b, key=lambda r: abs(cai(r) - cai(d2)))
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
    designs["motif"], designs["motif_name"] = forbidden[0], name
    Path(out).write_text(json.dumps(designs, indent=1) + "\n")
    for k in ("d1", "d2", "d3"):
        d = designs[k]
        print(f"{k}: CAI {d['cai']:.4f}, MFE {d['free_mfe_kcal_mol']}, "
              f"{cds_utr_pairs(d['free_structure'])} CDS-UTR pairs")


if __name__ == "__main__":
    main()
