"""An arms pair's overview inputs, in the layout overview2.py reads, and its three designs.

    python overview2_prep_arms.py ARMS_DIR SO_JSONL FOLD_SEQ LAYOUT_JSON OUT_DIR

ARMS_DIR is a run3-hba arms directory. Hull vertices of the coding-only front are computed
from its points. The designs follow overview_designs.py: d2 is mRNAfold's design (CDS-UTR pairs
discouraged) that contains a forbidden motif, CAI nearest CAI_TARGET; d1 the concatenation hull
vertex nearest d2 in CAI with another sequence; d3 the A5 design nearest d2 in CAI, refolded
with fold_seq and required to reproduce its recorded free fold.

LAYOUT_JSON (a JSON object: the UTR pair's name for the title, and any placement overview2.py
should not choose itself) is copied into the design file.
"""

import json
import math
import re
import subprocess
import sys
from pathlib import Path

from derna2_analysis.fronts import hull_vertices
from derna2_analysis.io import read_jsonl, write_jsonl

CAI_TARGET = 0.87


def cai(r):
    return math.exp(r["cai_log"] / 78)


def main():
    arms, so_path, fold_seq = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    layout, out = json.loads(sys.argv[4]), Path(sys.argv[5])
    out.mkdir(parents=True, exist_ok=True)
    so_rows = read_jsonl(so_path)
    head = so_rows[0]
    u5, u3 = head["utr5"]["sequence"], head["utr3"]["sequence"]
    l5, l3 = len(u5), len(u3)
    motifs = [m.replace("T", "U") for m in head["motifs"]]
    rx = {m: re.compile("(?=" + m.replace("R", "[AG]") + ")") for m in motifs}

    def hits(seq):
        cds = seq[l5 : len(seq) - l3]
        return [m for m in motifs if rx[m].search(cds)]

    so = [dict(r, mfe_kcal_mol=r["mfe"]) for r in so_rows if r.get("record") == "point"]
    if any(r["violations"] or hits(r["sequence"]) or r["gap"] != 0 for r in so):
        sys.exit("an A5 design carries a motif or is not structurally optimal")
    front = [r for r in read_jsonl(arms / "a0" / "front.jsonl") if r.get("record") == "point"]
    hv = set(hull_vertices([(r["cai_log"], r["mfe"]) for r in front]))
    a0 = [dict(r, hull_vertex=(r["cai_log"], r["mfe"]) in hv) for r in front]
    hull_cds = {r["sequence"] for r in a0 if r["hull_vertex"]}
    b = [dict(r, cds_sequence=r["cds"]) for r in read_jsonl(arms / "a0" / "concat.scored.jsonl")]
    mf = [
        dict(r, theirs_cai_log=r["cai_log"])
        for r in read_jsonl(arms / "a3" / "mrnafold.scored.jsonl")
        if r["arm"] == "guard"
    ]
    cd = read_jsonl(arms / "a2" / "front.scored.jsonl")

    d2 = min((r for r in mf if hits(r["sequence"])), key=lambda r: abs(cai(r) - CAI_TARGET))
    d1 = min(
        (r for r in b if r["cds"] in hull_cds and r["sequence"] != d2["sequence"]),
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
    if (
        fold["structure"] != d3["free_structure"]
        or abs(fold["mfe_kcal_mol"] - d3["free_mfe"]) > 1e-6
    ):
        sys.exit(f"d3 refolds to {fold['mfe_kcal_mol']}, recorded {d3['free_mfe']}")
    d3 = dict(d3, free_mfe_kcal_mol=fold["mfe_kcal_mol"])

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
    designs["motifs"] = motifs
    designs.update(layout)
    (out / "overview_designs.json").write_text(json.dumps(designs, indent=1) + "\n")
    print(f"UTRs {l5}/{l3} nt, A5 {len(so)} pts, CDS hull {len(hull_cds)}")
    for k in ("d1", "d2", "d3"):
        d = designs[k]
        print(f"{k}: CAI {d['cai']:.4f}, MFE {d['free_mfe_kcal_mol']}, motifs {d['motifs_found']}")


if __name__ == "__main__":
    main()
