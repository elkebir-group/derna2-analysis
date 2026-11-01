"""Cache everything the composite figure needs, so the layouts are cheap to redraw.

Panel (a) is draft_panel_a.py's computation, run once and dumped. Panel (b) is the motif cost
of each rung against F0 at the same cutoff, per protein, over the finished set.
"""
import glob
import hashlib
import json
import os
import re
import sys
import types

from e3paths import OUT  # noqa: F401


# fig2 imports ViennaRNA for a layout helper, and this interpreter does not have it. Only its
# catalog functions are wanted, so RNA -- and only RNA -- is stubbed, by a module that raises
# on any attribute access rather than returning something plausible. matplotlib is NOT stubbed
# and must be genuinely importable: the panel (a) prefix exec'd below imports and configures
# it. An earlier version of this comment described a matplotlib stub that was never written.
class _Unused(types.ModuleType):
    """A module that must not be used.

    fig2 imports ViennaRNA at module level but touches it only in its layout helper, which this
    script never calls; the catalog functions it does call are pure. Rather than fabricate a
    working stub -- which would let a future fig2 silently compute something with a fake RNA --
    any attribute access raises, so the build fails loudly if that ever stops being true.
    """

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        raise RuntimeError(f"fig2 touched RNA.{name}: this script runs without ViennaRNA, so "
                           f"the value would be fabricated. Run the data step under an "
                           f"interpreter with ViennaRNA instead.")


sys.modules.setdefault("RNA", _Unused("RNA"))

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "draft_panel_a.py")).read().split("PALETTE = {")[0]
# This half runs under the ViennaRNA interpreter, which has no matplotlib; the drawing
# half runs under the one that does. Only the computation is wanted here.
# __file__ with it: the prefix resolves its own neighbours against it, and an exec'd
# string has none unless it is given one.
ns = {"__name__": "panel_a", "__file__": os.path.join(HERE, "draft_panel_a.py")}
exec(compile(src, "draft_panel_a.py", "exec"), ns)
A = {"entries": ns["entries"], "classes": ns["classes"], "n": len(ns["POP"]),
     "census": {str(k): v for k, v in ns["census"].items()}, "floor": ns["FLOORS"][0]}
# The corpus this derives from E1's sweep is the one the supplement's panel (d) counts
# over, and that one is a committed list because deriving it needs the sweep, which is not
# in E2's deposit. Asserting them equal is what keeps a committed list from going stale
# silently: two derivations of one population agree until they do not, and the symptom is
# two counts in one paper that a reader cannot reconcile.
_CORPUS = os.environ.get("E2_CORPUS")
if _CORPUS:
    _want = json.load(open(_CORPUS))["members"]
    assert sorted(ns["POP"]) == sorted(_want), (
        f"{_CORPUS} lists {len(_want)} proteins, this step derives {len(ns['POP'])}; "
        f"only here: {sorted(set(ns['POP']) - set(_want))}, "
        f"only there: {sorted(set(_want) - set(ns['POP']))}")
# Every count is a plain int. This does not guard against the RNA stub -- that raises on
# first touch, so it can never return a value to be caught here -- it guards against a count
# arriving as a float or a string from the JSON read and serialising unnoticed.
for e in A["entries"]:
    assert all(isinstance(e[k], int) for k in ("occurs", "spellable")), e
for v in A["census"].values():
    assert isinstance(v["proteins"], int) and all(isinstance(x, int) for x in v["counts"].values())
print(f"panel (a): {len(A['entries'])} motifs over {A['n']} proteins, floor {A['floor']}")

sys.path.insert(0, HERE)
from e3common import NESTS, E, census, load, path  # noqa: E402

# 0.70 was dropped: the F0 fronts reach down to a median CAI of 0.751, so a 0.70 floor binds
# on 1 of 35 proteins and that column was the unconstrained optimum wearing a cutoff's label.
CUTS = [0.80, 0.85, 0.90, 0.95]
fin = [m for m in json.load(open(f"{OUT}/auc01_fin.json")) if m["rungs"] == NESTS]
B = {f"{c:.2f}": {n: [] for n in ("F1", "F2", "F3")} for c in CUTS}
BK = {f"{c:.2f}": {n: [] for n in ("F1", "F2", "F3")} for c in CUTS}
CELLS = {f"{c:.2f}": [] for c in CUTS}
ALL3 = {f"{c:.2f}": [] for c in CUTS}
# Proteins where each of the three sets occurs SOMEWHERE on the unconstrained front -- not
# necessarily together in one design. The subset exists to ask whether the cost looks
# different where every rung has something to forbid, and that is three independent per-rung
# questions, so the sets are tested independently. The stricter reading, "one design spells
# all three at once", was also computed and differs on exactly one protein of 42: RBD, whose
# front carries set I on all 695 designs, set II on 17 and set III on 415, with no design
# carrying both II and III. Each-set-somewhere is a recorded user decision, not a default, so
# the predicate is named in the provenance rather than left to be inferred from a count.
#
# Taken over the front, not over the design picked at one floor: at a given floor only a
# handful of proteins draw all three, so a per-floor rule would give a different population
# in each column of the panel.
THREE_PREDICATE = "each set occurs somewhere on the F0 front, not necessarily in one design"
THREE = set()
for m in fin:
    _h, _p = load(path(m["cell"], "F0"))
    _cen = [census(r["sequence"]) for r in _p]
    if all(any(c[k] for c in _cen) for k in ("I", "II", "III")):
        THREE.add(m["cell"])
print(f"all three sets on the F0 front: {len(THREE)} of {len(fin)} -> "
      + " ".join(sorted(THREE)))
EXACT = {f"{c:.2f}": [] for c in CUTS}
BINDS = {f"{c:.2f}": 0 for c in CUTS}
# Why a rung costs what it costs, per protein. A zero has three causes and a box cannot tell
# them apart: the set is absent from this protein's front entirely ("never"), it is on the
# front but the design the rung BELOW already picked does not spell it ("gone" -- the floor or
# the shallower rung displaced it first), or that pick does spell it and this rung has
# something to remove ("binds"). Only the third makes a cost a price.
STEP = {"F1": ("F0", "I"), "F2": ("F1", "II"), "F3": ("F2", "III")}
WHY = {f"{c:.2f}": {n: [] for n in STEP} for c in CUTS}
for m in fin:
    fr = {n: load(path(m["cell"], n)) for n in NESTS}
    onf0 = {s: any(census(r["sequence"])[s] for r in fr["F0"][1]) for s in ("I", "II", "III")}
    for c in CUTS:
        at = {n: [r for r in fr[n][1] if r["cai"] >= c] for n in NESTS}
        if not all(at.values()):
            continue
        pick = {n: min(v, key=lambda r: r["mfe"]) for n, v in at.items()}
        b = {n: pick[n]["mfe"] for n in NESTS}
        for n in ("F1", "F2", "F3"):
            B[f"{c:.2f}"][n].append(100 * (b[n] - b["F0"]) / abs(b["F0"]))
            BK[f"{c:.2f}"][n].append(b[n] - b["F0"])
            prev, st_ = STEP[n]
            WHY[f"{c:.2f}"][n].append(
                "never" if not onf0[st_]
                else "binds" if census(pick[prev]["sequence"])[st_]
                else "gone")
        CELLS[f"{c:.2f}"].append(m["cell"])
        ALL3[f"{c:.2f}"].append(m["cell"] in THREE)
        EXACT[f"{c:.2f}"].append(bool(m["exact"]))
        BINDS[f"{c:.2f}"] += m["cai_lo"] < c
print("why each rung costs what it does (never / already gone / binds):")
for c in CUTS:
    k = f"{c:.2f}"
    print(f"  {k}: " + "   ".join(
        f"{n} " + "/".join(str(sum(1 for x in WHY[k][n] if x == w))
                           for w in ("never", "gone", "binds")) for n in ("F1", "F2", "F3")))
print("panel (b): " + ", ".join(
    f"{k} n={len(v['F3'])} exact={sum(EXACT[k])} binds={BINDS[k]}" for k, v in B.items()))
# The dark bar is a share of cen["proteins"], the axis is labelled with len(POP); they agree
# today and the census loop has two paths that could silently part them.
for v in A["census"].values():
    assert v["proteins"] == A["n"], (v["proteins"], A["n"])


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


MEMBERS = {"panel_b_finished": sorted(m["cell"] for m in fin),
           "panel_b_certified": sorted(m["cell"] for m in fin if m["exact"]),
           "panel_b_all_three_sets": sorted(THREE)}

# Every population above is a subset of FINISHED.tsv, so a cell the record retains and the
# report never mentions is invisible to all of them: no count disagrees, no digest moves, and
# the figure's denominator is a subset of a reported set whose size it never states. This walks
# the band directories -- containers, not an enumerated cell list, so it cannot go stale as
# bands grow -- and names every cell with a round file on disk and no row in the report.
#
# An absence means the cell wants a line in the record, NOT that it wants a decision: the report
# is scoped, by its own docstring, to one host's banded rounds, so the line may equally read
# "reported in another table" or "hosted elsewhere". Reading an absence as a shortfall is reading
# it against a gate the cell was never run under. The classification lives in CELL-SET.md and the
# bands' own tables and is deliberately not duplicated here.
#
# One discrimination is cheap enough to make here, because it is the difference the blind spot
# below turns on: a `.done` marker beside the round file says the cell ran to a clean end, so an
# absence with `.done` is reported elsewhere and an absence without one never finished here.
# What remains invisible either way is a cell excluded before any run, which leaves no file at
# all and so looks exactly like a cell that does not exist.
BANDS = {"short3": "short3", "short2": "short2", "grid": "grid", "mid": "mid",
         "long": "long/cells"}
REPORTED_CELLS = set()
for _l in open(f"{E}/FINISHED.tsv"):
    if _l.startswith("#") or _l.startswith("cell\t") or not _l.strip():
        continue
    REPORTED_CELLS.add(_l.split("\t")[0])
ABSENT = {}
for _b, _d in BANDS.items():
    _on = sorted({os.path.basename(p)[:-6] for p in glob.glob(f"{E}/{_d}/*_F[1-4].jsonl")
                  if re.fullmatch(r"[A-Za-z0-9-]+_F[1-4]", os.path.basename(p)[:-6])})
    _miss = [c for c in _on if c not in REPORTED_CELLS]
    _done = [c for c in _miss if os.path.exists(f"{E}/{_d}/{c}.done")]
    ABSENT[_b] = {"on_disk": len(_on), "in_report": len(_on) - len(_miss), "absent": _miss,
                  "absent_ran_to_a_clean_end": _done}
_na = sum(len(v["absent"]) for v in ABSENT.values())
_nd = sum(len(v["absent_ran_to_a_clean_end"]) for v in ABSENT.values())
print("cells on disk and absent from the report: "
      + ", ".join(f"{b} {len(v['absent'])}/{v['on_disk']}" for b, v in ABSENT.items())
      + f" -- {_na} in all, of which {_nd} ran to a clean end and so are reported somewhere "
        f"other than this table; each wants a line in the record, not a decision")
# Provenance, to the same shape fig2 writes: every input this file was computed from, by
# checksum, plus the interpreter and the libraries. A figure is reproducible only if the thing
# it was drawn from says what it was made of.
INPUTS = {}
for f in sorted(glob.glob(f"{E}/motifs/set*.txt")):
    INPUTS[os.path.relpath(f, E)] = sha(f)
FRONTS = {}
for m in fin:
    for n in NESTS:
        q = path(m["cell"], n)
        if q:
            FRONTS[f"{m['cell']}_{n}"] = {"path": q, "sha256": sha(q)}
SCRIPTS = {os.path.basename(f): sha(f) for f in
           [__file__] + [os.path.join(HERE, n) for n in
                         ("e3common.py", "draft_panel_a.py", "fig2.py", "e3paths.py")]}
# Nothing in this file may depend on anything but the figure's own inputs, because the file is
# pinned in pins/e2_composite.outputs and the corpus underneath it is live: the bands keep
# running and the finished table keeps being rewritten. Three fields were removed for failing
# that, each of which moved while every plotted number stayed identical.
#
#   - a build time, which changes on every run and which no digest can cover;
#   - the analysis repository's own commit, which is self-referential -- the commit that would
#     make the pin match is created after the build that records it, so the field can never
#     agree with a rebuild. SCRIPTS below digests every script this file was computed by, which
#     is the code identity the commit was standing in for, and it is checkable;
#   - FINISHED.tsv's header and checksum. The table carries its own write time inside it, so
#     both move when it is regenerated over an unchanged finished set. What the figure actually
#     read is the parsed set, digested below, which moves only when that set moves.
#
# The live-corpus audit (which cells are on disk and absent from the report) is a property of
# the host at the time of the build and not of the figure, so it is written beside this file in
# composite-census.json and is deliberately not pinned.
PROV = {
    "python": sys.version.split()[0],
    "numpy": __import__("numpy").__version__,
    "finished_rows": hashlib.sha256(
        json.dumps(fin, sort_keys=True).encode()).hexdigest(),
    "census_file": "composite-census.json",
    "inputs": INPUTS, "scripts": SCRIPTS, "fronts": FRONTS,
    # A count is not a membership: one protein can leave and another arrive between two builds
    # and every number here reads unchanged, so each population names its members and carries a
    # digest of that name list. The digest is what a reader diffs; the names are what makes an
    # exclusion auditable rather than asserted.
    "populations": {"panel_a_screened": A["n"], "panel_b_finished": len(fin),
                    "panel_b_certified": sum(1 for m in fin if m["exact"]),
                    "panel_b_all_three_sets": len(THREE)},
    "all_three_predicate": THREE_PREDICATE,
    "members": MEMBERS,
    "member_digests": {k: hashlib.sha256(" ".join(v).encode()).hexdigest()[:16]
                       for k, v in MEMBERS.items()},
}
json.dump({"a": A, "b": B, "b_kcal": BK, "b_cells": CELLS, "b_exact": EXACT,
           "b_binds": BINDS, "b_all3": ALL3, "b_why": WHY, "provenance": PROV},
          open(f"{OUT}/composite-data.json", "w"))
# The census of the host the build ran on: the finished table it read, by header and checksum,
# and every cell with a round file on disk and no row in that table. Both move as the bands run,
# so this file is not pinned and nothing the figure draws is computed from it.
json.dump({"finished_tsv": {"header": open(f"{E}/FINISHED.tsv").readline().rstrip("\n"),
                            "sha256": sha(f"{E}/FINISHED.tsv")},
           "finished_rows": PROV["finished_rows"],
           "cells_absent_from_report": ABSENT},
          open(f"{OUT}/composite-census.json", "w"), indent=1)
print(f"provenance: {len(FRONTS)} front files, {len(INPUTS)} motif sets, "
      f"{len(SCRIPTS)} scripts, finished rows {PROV['finished_rows'][:8]}")
print("wrote composite-data.json and composite-census.json")
