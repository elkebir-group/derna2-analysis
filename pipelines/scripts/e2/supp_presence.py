"""How much of each motif set survives on the front as constraint levels are added.

    E2_DEPOSIT=... E2_MOTIF_DIR=... E2_SUPP=... python supp_presence.py

Every Pareto optimal design of every finished protein is scanned for every pattern, at all
four levels, and three measures are recorded per (protein, level, set):

  share   the fraction of that level's designs spelling the set at least once
  mean    the mean number of sites per design
  chosen  the sites in the most stable design at CAI >= 0.80

A level forbids its own sets, so those cells are zero by construction. That is not padding:
it is the check that the scan agrees with what the solver emitted, and a non-zero there
would mean one of the two is wrong. The figure's subject is the sets a level does NOT
forbid.
"""

import re
import statistics as st
import sys

import supp_common as common

MOTIFS = common.motifs()
PATS = {s: [re.compile("".join(common.IUPAC.get(c, c) for c in m))
            for t, m in MOTIFS if t == s] for s in common.SETS}
FLOOR = 0.80


def measures(designs):
    hits = {s: [sum(len(p.findall(seq)) for p in PATS[s]) for _, _, seq in designs]
            for s in common.SETS}
    best = max((d for d in designs if d[0] >= FLOOR), key=lambda d: -d[1], default=None)
    return {
        "share": {s: sum(1 for h in hits[s] if h) / len(designs) for s in common.SETS},
        "mean": {s: sum(hits[s]) / len(hits[s]) for s in common.SETS},
        "median": {s: st.median(hits[s]) for s in common.SETS},
        "chosen": {s: (sum(len(p.findall(best[2])) for p in PATS[s]) if best else None)
                   for s in common.SETS},
    }


def main():
    tsv, digest = common.finished()
    _, f0path = common.sequences()
    done = sorted(c for c, v in tsv.items()
                  if all(v.get(n, {}).get("status") in ("complete", "retired")
                         for n in common.NESTS) and c in f0path)
    rows = []
    for c in done:
        for n in common.LEVELS:
            path = f0path[c] if n == "F0" else tsv[c][n]["path"]
            ds = [d for d in common.points(path) if d[2]]
            if not ds:
                continue
            rows.append(dict(
                cell=c, level=n, band=tsv[c][common.NESTS[0]]["band"], n_designs=len(ds),
                status="complete" if n == "F0" else tsv[c][n]["status"],
                **measures(ds)))
        print(c, "done", flush=True)
    common.write("presence.json", {
        "floor": FLOOR,
        "rows": rows,
        "provenance": {
            "n_proteins": len(done), "n_rows": len(rows),
            "members": done,
            "finished_set_sha256": digest,
            "inputs": common.digest_of(*common.motif_files()),
            "scripts": common.digest_of(__file__, common.__file__),
        },
    })


if __name__ == "__main__":
    sys.exit(main())
