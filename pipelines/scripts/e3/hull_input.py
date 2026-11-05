"""Write a derna result file of one cell's concatenation hull vertices, for `derna check`.

    python hull_input.py BUNDLE E2_RESULTS CELL OUT.jsonl

`derna check` refolds each design and reports the CDS-UTR pairs of its structure, which the
census rows do not carry. Where the census has no constrained energy (the middle band's one-fold
census, wherever the free fold forms a CDS-UTR pair), the free energy is written as `mfe`: a lower
bound on the constrained optimum, so the check computes that optimum and writes it back.
"""

import json
import sys

from e3common import PAIRS, a0_path, header, rd, relocate


def main():
    bundle, e2, cell, out = sys.argv[1:5]
    h = header(bundle, cell)
    rows = list(rd(a0_path(bundle, cell)))
    src = relocate(rows[0]["source"], bundle, e2)
    hv = {i for i, x in enumerate(r for r in rd(src) if r.get("record") == "point") if x.get("on_hull")}
    head = {k: h[k] for k in ("format", "format_version", "version", "protein", "utr5", "utr3", "forbid", "access",
                              "codon_usage") if k in h}
    if "format" not in h:
        # A seed file derived from an A4 search log (the middle band's hbb cells) carries the
        # cell's protein, UTRs and forbid but no result stamp or codon table. Those two are the
        # same on every UTR pair of a protein, so they come from a sibling cell's `derna` header.
        protein = cell.split("__")[0]
        sib = next(s for s in (header(bundle, f"{protein}__{p}") for p in PAIRS if p != cell.split("__")[1])
                   if "format" in s)
        head.update({k: sib[k] for k in ("format", "format_version", "version", "access", "codon_usage")
                     if k in sib})
    head.update(record="header", command="front", stop=rows[0]["stop"])
    with open(out, "w") as f:
        f.write(json.dumps(head) + "\n")
        for x in rows[1:]:
            if x.get("record") == "point" and x["index"] in hv:
                f.write(json.dumps({"record": "point", "index": x["index"], "cai": x["cai"], "cai_log": x["cai_log"],
                                    "mfe": x["constrained_kcal_mol"] if x["constrained_kcal_mol"] is not None
                                    else x["free_kcal_mol"], "sequence": x["sequence"],
                                    "structure": "." * len(x["sequence"]), "census_gap": x["gap_kcal_mol"]}) + "\n")


if __name__ == "__main__":
    main()
