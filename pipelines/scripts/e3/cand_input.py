"""The two non-DERNA2 designs the main figure shows for one instance, as derna result files for
`derna check`.

    python cand_input.py BUNDLE E2_RESULTS CELL CAI_FLOOR OUT_PREFIX

Writes OUT_PREFIX.conc.in.jsonl and OUT_PREFIX.mrf.in.jsonl, one design each. The concatenation design is
the hull vertex of the cell's concatenation front with CAI >= CAI_FLOOR of least CDS-alone MFE
(`f0_mfe`), the objective of the methods concatenation stands for. The mRNAfold design is the one
of the cell's mRNAfold designs with CAI >= CAI_FLOOR of least MFE without CDS-UTR pairs
(`guarded_kcal_mol`). `derna check` refolds each on the whole mRNA, which gives the structure
the figure draws and the positions of its CDS-UTR pairs.
"""

import glob
import json
import sys

from e3common import a0_path, header, hull_vertices, rd


def write(path, head, x, stop):
    head = {k: head[k] for k in ("format", "format_version", "version", "protein", "utr5", "utr3", "forbid", "access",
                                 "codon_usage") if k in head}
    head.update(record="header", command="front", stop=stop)
    with open(path, "w") as f:
        f.write(json.dumps(head) + "\n")
        f.write(json.dumps({"record": "point", "index": 0, "cai": x["cai"], "cai_log": x["cai_log"], "mfe": 0.0,
                            "sequence": x["sequence"], "structure": "." * len(x["sequence"])}) + "\n")


def main():
    bundle, e2, cell, floor, prefix = sys.argv[1:6]
    floor = float(floor)
    head = header(bundle, cell)
    rows = list(rd(a0_path(bundle, cell)))
    hv = {x["index"] for x in hull_vertices(bundle, e2, cell)}
    conc = min((x for x in rows if x.get("record") == "point" and x["index"] in hv and x["cai"] >= floor),
               key=lambda x: x["f0_mfe"])
    write(f"{prefix}.conc.in.jsonl", head, conc, rows[0]["stop"])
    des = {}
    for p in glob.glob(f"{bundle}/a3soft/designs.*.jsonl"):
        for x in rd(p):
            if x["cell"] == cell:
                des[x["index"]] = x
    cen = {}
    for p in sorted(glob.glob(f"{bundle}/a3soft/census.*.jsonl")):
        for x in rd(p):
            k = x.get("index")
            if x.get("cell") == cell and k in des and k not in cen and "gap_kcal_mol" in x and "guarded_kcal_mol" in x:
                cen[k] = x
    k = min((i for i in cen if des[i]["cai"] >= floor), key=lambda i: cen[i]["guarded_kcal_mol"])
    write(f"{prefix}.mrf.in.jsonl", head, des[k], des[k]["stop"])


if __name__ == "__main__":
    main()
