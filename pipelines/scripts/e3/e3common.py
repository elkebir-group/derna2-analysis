"""Readers for the E3 deposit that the E3 steps share.

The deposit is the figure bundle's paper-e4/: per cell (one protein on one practical UTR pair,
named <UniProt>__<pair>) a census of designs under seeds/, the concatenation arm's points under
a0/ or a0-all/, and so on. The concatenation arm points at E2's unconstrained fronts through a
`source` path recorded by the host that ran it, which is rewritten to E2's results root here.

The recorded prefixes come from the environment (`E3_RECORDED_E2`, `E3_RECORDED_DEPOSIT`),
set by the rules from `e3.recorded_prefixes`: they name one past run's filesystem, so a
reader that spelled them out would be right about that run only.
"""

import glob
import json
import os

RECORDED_E2 = os.environ.get("E3_RECORDED_E2", "")
RECORDED_DEPOSIT = os.environ.get("E3_RECORDED_DEPOSIT", "")
PAIRS = ["hbb", "mrna1273", "hba", "curevac", "bnt"]


def rd(path):
    for line in open(path):
        try:
            yield json.loads(line)
        except ValueError:
            continue


def cells_in(bundle):
    """The cells the deposit holds a census for."""
    return sorted({os.path.basename(p).split(".")[0] for p in glob.glob(f"{bundle}/seeds/*.check.*.jsonl")})


def header(bundle, cell):
    """The cell's header row: protein, both UTRs and the run's settings."""
    for name in ("freestop.ec083208", "pooled-fixed"):
        p = f"{bundle}/seeds/{cell}.check.{name}.jsonl"
        if os.path.exists(p):
            return next(iter(rd(p)))
    raise FileNotFoundError(f"{bundle}/seeds/{cell}.check.*")


def relocate(src, bundle, e2):
    """A recorded `source` path, on this host: E2's fronts under E2's results root, and the
    deposit's own files (the middle band's one-fold fronts, a0-mid/f0/) under the deposit.

    A prefix that is not set rewrites nothing, which leaves the path as recorded rather than
    silently rewriting at the root.
    """
    for recorded, here in ((RECORDED_E2, e2), (RECORDED_DEPOSIT, bundle)):
        if recorded:
            src = src.replace(recorded, here)
    return src


def refold_points(hullcheck, cell):
    """The `derna check` rows of a cell's hull vertices, by index."""
    return {x["index"]: x for x in rd(f"{hullcheck}/{cell}.out.jsonl") if x.get("record") == "point"}


def take_refold(x, r):
    """Join a hull vertex to its refold. Where the census computed the constrained energy, the two
    gaps agree; where it did not (the middle band's one-fold census writes null whenever the free
    fold forms a CDS-UTR pair), the refold's exact energy and gap stand in for it."""
    if x["gap"] is None:
        x["mfe"], x["gap"] = r["mfe"], r["gap"]
    else:
        assert abs(r["gap"] - x["gap"]) < 0.01, (x["index"], r["gap"], x["gap"])
    x["violations"] = r.get("violations")
    return x


def a0_path(bundle, cell):
    return next(q for q in (f"{bundle}/{t}/cells/{cell}/a0.jsonl" for t in ("a0", "a0-all")) if os.path.exists(q))


def hull_vertices(bundle, e2, cell, hullcheck=None):
    """The concatenation front's hull vertices of the cell: index, CAI, constrained MFE and the
    census gap, which is 0 for a modular design. Given the refolds, each vertex also takes its
    CDS-UTR pairs, and its energy and gap where the census has none."""
    rows = list(rd(a0_path(bundle, cell)))
    src = relocate(rows[0]["source"], bundle, e2)
    hv = {i for i, x in enumerate(r for r in rd(src) if r.get("record") == "point") if x.get("on_hull")}
    out = [{"index": x["index"], "cai": x["cai"], "mfe": x["constrained_kcal_mol"], "gap": x["gap_kcal_mol"]}
           for x in rows[1:] if x.get("record") == "point" and x["index"] in hv]
    if hullcheck is None:
        return out
    r = refold_points(hullcheck, cell)
    return [take_refold(x, r[x["index"]]) for x in out]
