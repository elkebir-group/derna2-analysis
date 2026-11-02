"""Where the CDS-UTR contacts of the concatenation designs fall along the mRNA, per UTR pair.

    python contacts_data.py BUNDLE E2_RESULTS HULLCHECK CELLS.txt OUT.json

A design is a non-modular hull vertex of its cell's concatenation front. Its contacts are the
CDS-UTR base pairs of its refolded structure, positions 1-based as `derna check` writes them.
A design counts once in each bin that holds a contact endpoint. The 5' UTR is binned in nt
before the start codon (3 nt), the CDS in percent of its length (1 %), the 3' UTR in nt after
the stop codon (3 nt). A pair of the 5' UTR and the CDS has its two ends in the first two, and
a pair of the CDS and the 3' UTR in the last two.
"""

import json
import sys
from pathlib import Path

import numpy as np

from e4common import PAIRS, cells_in, header, hull_vertices


def main():
    bundle, e2, hullcheck, cells_file, out = sys.argv[1:6]
    cells = [line.split()[0] for line in open(cells_file) if line.strip()]
    assert cells == cells_in(bundle), "the deposit's cells differ from inputs/e3/cells.txt"
    bc = np.linspace(0, 100, 101)
    res = {}
    for pair in PAIRS:
        cs = [c for c in cells if c.endswith("__" + pair)]
        h0 = header(bundle, cs[0])
        l5, l3 = len(h0["utr5"]["sequence"]), len(h0["utr3"]["sequence"])
        b5 = np.arange(-3 * (l5 // 3 + 1), 1, 3)
        b3 = np.arange(0, 3 * (l3 // 3 + 1) + 1, 3)
        h5, hc5, hc3, h3 = (np.zeros(len(b5) - 1), np.zeros(len(bc) - 1), np.zeros(len(bc) - 1),
                            np.zeros(len(b3) - 1))
        n = 0
        for c in cs:
            h = header(bundle, c)
            assert (len(h["utr5"]["sequence"]), len(h["utr3"]["sequence"])) == (l5, l3), c
            cend = l5 + 3 * (len(h["protein"]["sequence"]) + 1)
            for x in hull_vertices(bundle, e2, c, hullcheck):
                if x["gap"] == 0:
                    continue
                n += 1
                e5, ec5, ec3, e2_ = [], [], [], []
                for i, j in x.get("violations") or []:
                    a, b = min(i, j) - 1, max(i, j) - 1
                    if a < l5 <= b < cend:
                        e5.append(a - l5)
                        ec5.append(100 * (b - l5) / (cend - l5))
                    elif l5 <= a < cend <= b:
                        ec3.append(100 * (a - l5) / (cend - l5))
                        e2_.append(b - cend)
                    else:
                        raise ValueError(f"{c}: pair {i},{j} is not a CDS-UTR pair")
                for edges, hist, vals in ((b5, h5, e5), (bc, hc5, ec5), (bc, hc3, ec3), (b3, h3, e2_)):
                    if vals:
                        hist += np.histogram(vals, edges)[0] > 0
        res[pair] = {"cells": len(cs), "utr5_nt": l5, "utr3_nt": l3, "designs": n,
                     "utr5_edges": b5.tolist(), "utr5": h5.astype(int).tolist(),
                     "cds_edges": bc.tolist(), "cds_with_5": hc5.astype(int).tolist(),
                     "cds_with_3": hc3.astype(int).tolist(),
                     "utr3_edges": b3.tolist(), "utr3": h3.astype(int).tolist()}
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print({p: v["designs"] for p, v in res.items()})


if __name__ == "__main__":
    main()
