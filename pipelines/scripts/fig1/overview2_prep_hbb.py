"""The hbb overview's inputs, in the layout overview2.py reads, and its three designs.

    python overview2_prep_hbb.py FIG1_RESULTS MRNAFOLD_COMPARE UTR5_TXT UTR3_TXT FOLD_SEQ \
        LAYOUT_JSON OUT_DIR

FIG1_RESULTS is the fig1 deposit (a_cdsonly_front, b_cdsonly_with_utrs.scored,
cd_symmetric_guarded_front.scored, d_so_front); MRNAFOLD_COMPARE the fig1_fronts deposit's
mrnafold/compare.jsonl. Forbidden motif: CUCGAG (XhoI). Hull vertices of the coding-only front
are computed from its points. The designs follow overview2_prep_arms.py: d2 is mRNAfold's
design (CDS-UTR pairs discouraged) that contains the motif, CAI nearest CAI_TARGET; d1 the
concatenation hull vertex nearest d2 in CAI with another sequence; d3 the structurally optimal
design nearest d2 in CAI, refolded with fold_seq because the front stores energies only.

LAYOUT_JSON (a JSON object: the UTR pair's name for the title, and any placement overview2.py
should not choose itself) is copied into the design file.
"""

import json
import math
import subprocess
import sys
from pathlib import Path

from derna2_analysis.fronts import hull_vertices
from derna2_analysis.io import read_jsonl, write_jsonl

CAI_TARGET = 0.87
MOTIF = "CUCGAG"


def cai(r):
    return math.exp(r["cai_log"] / 78)


def main():
    res, compare, u5p, u3p, fold_seq = (
        Path(sys.argv[1]),
        sys.argv[2],
        sys.argv[3],
        sys.argv[4],
        sys.argv[5],
    )
    layout, out = json.loads(sys.argv[6]), Path(sys.argv[7])
    out.mkdir(parents=True, exist_ok=True)
    u5, u3 = Path(u5p).read_text().strip(), Path(u3p).read_text().strip()
    l5, l3 = len(u5), len(u3)

    def cds(seq):
        return seq[l5 : len(seq) - l3]

    def hits(seq):
        return [MOTIF] if MOTIF in cds(seq) else []

    so = [dict(r) for r in read_jsonl(res / "d_so_front.jsonl")]
    if any(hits(r["sequence"]) for r in so):
        sys.exit("a structurally optimal design carries the motif: the front is not motif-free")
    front = read_jsonl(res / "a_cdsonly_front.jsonl")
    hv = set(hull_vertices([(r["cai_log"], r["mfe_kcal_mol"]) for r in front]))
    a0 = [dict(r, hull_vertex=(r["cai_log"], r["mfe_kcal_mol"]) in hv) for r in front]
    hull_cds = {r["sequence"] for r in a0 if r["hull_vertex"]}
    b = read_jsonl(res / "b_cdsonly_with_utrs.scored.jsonl")
    cd = read_jsonl(res / "cd_symmetric_guarded_front.scored.jsonl")
    mf = [
        dict(r, theirs_cai_log=r["theirs_cai_log"])
        for r in read_jsonl(compare)
        if r["arm"] == "guard"
    ]

    d2 = min(
        (r for r in mf if hits(r["sequence"])),
        key=lambda r: abs(math.exp(r["theirs_cai_log"] / 78) - CAI_TARGET),
    )
    twin = [r for r in cd if r["sequence"] == d2["sequence"]]
    if not twin or twin[0]["free_mfe_kcal_mol"] != d2["free_mfe_kcal_mol"]:
        sys.exit("mRNAfold's design is not a UTR-masking design with the same free MFE")
    d2 = dict(d2, cai_log=d2["theirs_cai_log"], free_structure=twin[0]["free_structure"])
    d1 = min(
        (r for r in b if r["cds_sequence"] in hull_cds and r["sequence"] != d2["sequence"]),
        key=lambda r: abs(cai(r) - cai(d2)),
    )
    d3 = min(so, key=lambda r: abs(cai(r) - cai(d2)))
    (out / "d3.tsv").write_text(f"d3\t{d3['sequence']}\n")
    subprocess.run(
        [fold_seq, "--in", out / "d3.tsv", "--out", out / "d3.jsonl"],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    fold = read_jsonl(out / "d3.jsonl")[0]
    if abs(fold["mfe_kcal_mol"] - d3["mfe_kcal_mol"]) > 1e-6:
        sys.exit(f"d3 refolds to {fold['mfe_kcal_mol']}, recorded {d3['mfe_kcal_mol']}")
    st, n, lc = [], 0, len(cds(d3["sequence"]))
    for i, c in enumerate(fold["structure"]):
        if c == "(":
            st.append(i)
        elif c == ")":
            j = st.pop()
            n += (l5 <= i < l5 + lc) != (l5 <= j < l5 + lc)
    if n:
        sys.exit("d3's refolded structure has a CDS-UTR base pair")
    d3 = dict(d3, free_structure=fold["structure"], free_mfe_kcal_mol=fold["mfe_kcal_mol"])
    if d1["symmetric_so"] or d2["symmetric_so"]:
        sys.exit("d1 or d2 is modular")

    write_jsonl(out / "a_cdsonly_front.jsonl", a0)
    write_jsonl(out / "b_cdsonly_with_utrs.scored.jsonl", b)
    write_jsonl(out / "cd_symmetric_guarded_front.scored.jsonl", cd)
    write_jsonl(out / "d_so_front.jsonl", so)
    write_jsonl(out / "mrnafold.jsonl", mf)
    (out / "u5.txt").write_text(u5 + "\n")
    (out / "u3.txt").write_text(u3 + "\n")
    designs = {
        k: dict(d, cai=cai(d), motifs_found=hits(d["sequence"]))
        for k, d in (("d1", d1), ("d2", d2), ("d3", d3))
    }
    designs["motifs"] = [MOTIF]
    designs.update(layout)
    (out / "overview_designs.json").write_text(json.dumps(designs, indent=1) + "\n")
    print(f"UTRs {l5}/{l3} nt, SO front {len(so)} pts, CDS hull {len(hull_cds)}")
    for k in ("d1", "d2", "d3"):
        d = designs[k]
        print(f"{k}: CAI {d['cai']:.4f}, MFE {d['free_mfe_kcal_mol']}, motifs {d['motifs_found']}")


if __name__ == "__main__":
    main()
