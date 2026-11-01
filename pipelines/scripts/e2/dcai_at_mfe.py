"""Panel (b) transposed: the CAI a rung gives up at a fixed STRUCTURE quality.

Panel (b) fixes CAI and reads off MFE. This fixes MFE and reads off CAI, which is the same
trade seen from the other side and a different claim: not "what does structure cost at a given
codon quality" but "what does codon quality cost at a given structure".

A shared MFE cutoff in kcal/mol would be meaningless across a corpus running 456 to 2,157 nt,
since MFE scales with length -- the mistake D2's own panel (a) exists to reject. The cutoff is
therefore relative to each protein's own unconstrained optimum: mfe <= p * e_ideal, i.e. "at
least p of the structure this protein can reach at all". Both sides are negative, so a LARGER
p is a tighter bound and 95% is stricter than 80% -- the same direction as a CAI floor, which
is why the two panels' x axes can run the same way. The form written here before,
`mfe <= e_ideal + p*|e_ideal|`, is the same inequality with the sign lost, and it reverses
the axis.

At each cutoff and rung, the quantity is the HIGHEST CAI reachable under that structure bound.
A rung that cannot reach the bound at all is reported as infeasible, not dropped: that is the
constraint refusing the structure outright, which is a stronger statement than a large cost.
"""
import json
import os
import statistics as st
import sys

sys.argv = [sys.argv[0]]
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
exec(open(os.path.join(HERE, "pick_example3.py")).read().split("cells = sorted(")[0])
# OUT, NESTS, load and path are bound by the prefix above, which no static reader can see.
# ruff: noqa: F821

PCTS = [0.80, 0.85, 0.90, 0.95]
fin = [m for m in json.load(open(f"{OUT}/auc01_fin.json")) if m["rungs"] == NESTS]
out = {f"{p:.2f}": {n: [] for n in ("F1", "F2", "F3")} for p in PCTS}
cells = {f"{p:.2f}": [] for p in PCTS}
infeas = {f"{p:.2f}": {n: [] for n in ("F1", "F2", "F3")} for p in PCTS}
exact = {f"{p:.2f}": [] for p in PCTS}
for m in fin:
    fr = {n: load(path(m["cell"], n)) for n in NESTS}
    ideal = min(r["mfe"] for r in fr["F0"][1])
    for p in PCTS:
        k = f"{p:.2f}"
        # "at least p of the best structure this protein can reach", the mirror of a CAI
        # floor. MFE is negative, so p * ideal is LESS negative and a larger p is a TIGHTER
        # bound; F0 always satisfies it, by construction, for p <= 1.
        bound = p * ideal
        best = {}
        for n in NESTS:
            v = [r["cai"] for r in fr[n][1] if r["mfe"] <= bound]
            best[n] = max(v) if v else None
        if best["F0"] is None:
            continue                      # the floor itself is unreachable; nothing to compare
        cells[k].append(m["cell"])
        exact[k].append(bool(m["exact"]))
        for n in ("F1", "F2", "F3"):
            if best[n] is None:
                infeas[k][n].append(m["cell"])
                out[k][n].append(None)
            else:
                out[k][n].append(best["F0"] - best[n])   # CAI given up, in CAI units
json.dump({"pcts": PCTS, "d": out, "cells": cells, "exact": exact, "infeasible": infeas},
          open(f"{OUT}/dcai-at-mfe.json", "w"))
print(f"CAI given up (CAI units) at a structure floor, over {len(fin)} proteins\n")
print(f"{'MFE >=':>8} {'n':>4} | " + "  ".join(f"{n:>22}" for n in ("F1", "F2", "F3")))
for p in PCTS:
    k = f"{p:.2f}"
    row = [f"{100*p:6.0f}% {len(cells[k]):5} |"]
    for n in ("F1", "F2", "F3"):
        v = [x for x in out[k][n] if x is not None]
        nf = len(infeas[k][n])
        row.append(f"  med {st.median(v):7.4f} max {max(v):7.4f} inf {nf:2}" if v
                   else f"  {'all infeasible':>22}")
    print("".join(row))
print("\nproteins a rung cannot reach the bound on at all:")
for p in PCTS:
    k = f"{p:.2f}"
    print(f"  MFE >= {100*p:.0f}%: " + "  ".join(
        f"{n} {len(infeas[k][n])}" for n in ("F1", "F2", "F3")))
