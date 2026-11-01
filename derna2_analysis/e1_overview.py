"""Data behind the E1 overview figure: what the hull of a front leaves out, and what the
whole front costs.

Everything is computed from the sweep rows of one pin, except the pruned hull search's wall
time (`sweep_s` of a protein's walk report), which the sweep rows do not carry.
"""

import math
from collections import defaultdict

import numpy as np

from derna2_analysis.fronts import hull_vertices, pareto

# Length groups of the figure: (label, last length in the group).
GROUPS = [("≤200", 200), ("201–350", 350), ("351–500", 500), ("501–1000", 1000), (">1000", None)]
# A cell's budget in seconds: 2 h up to 1000 aa, 6 h above.
BUDGET_S = (7200.0, 21600.0)
# A cell counts as finished when it exited ok within this factor of its budget.
BUDGET_SLACK = 1.02
# Comparator arms in the order of the figure: (arm, label).
ARMS = [
    ("derna2-walk", "DERNA2"),
    ("mrnafold-hull", "mRNAfold"),
    ("lineardesign-hull", "LinearDesign"),
    ("lcds-exact-hull", "LinearCDSfold"),
    ("derna1-own", "DERNA1"),
]


def group_of(n_aa):
    for i, (_, hi) in enumerate(GROUPS):
        if hi is None or n_aa <= hi:
            return i


def budget_s(n_aa):
    return BUDGET_S[1] if n_aa > 1000 else BUDGET_S[0]


def finished(cell):
    """Whether a cell ran to completion inside its budget; a missing cell did not."""
    if not cell or cell.get("exit") != "ok":
        return False
    return cell["wall_s"] <= BUDGET_SLACK * budget_s(cell["n_aa"])


def read_e1(rows, pin):
    """Cell rows by (arm, protein) and the walk's (CAI, MFE) front points by protein, from
    the rows of one pin. Of repeated cell rows the last replaces the earlier.

    A protein's points are the rescored ones where the walk was rescored, else the walk's
    own: its energies are on the same 0.1 kcal/mol grid, and one protein is never a mix.
    """
    cells, rescored, reported = {}, defaultdict(list), defaultdict(list)
    for r in rows:
        if r.get("pin") != pin:
            continue
        if r.get("kind") == "cell":
            cells[(r["arm"], r["protein"])] = r
        elif r.get("kind") == "front_point" and r.get("arm") == "derna2-walk":
            if r.get("cai_rescored") and r.get("mfe_rescored") is not None:
                rescored[r["protein"]].append((r["cai_rescored"], r["mfe_rescored"]))
            elif r.get("cai") and r.get("mfe") is not None:
                reported[r["protein"]].append((r["cai"], r["mfe"]))
    points = {p: rescored[p] or reported[p] for p in set(rescored) | set(reported)}
    return cells, points


class Front:
    """A protein's Pareto front by increasing CAI, and which of its points are hull vertices."""

    def __init__(self, points):
        # The front and hull are taken in (log CAI, MFE), where the hull is the set of
        # points optimal for some weight λ.
        log = [(math.log(c), e) for c, e in points]
        front = sorted(pareto(log))
        hull = set(hull_vertices(log))
        self.cai = np.array([math.exp(a) for a, _ in front])
        self.mfe = np.array([e for _, e in front])
        self.on_hull = np.array([p in hull for p in front], dtype=bool)

    def steps(self, hull_only):
        """Positive CAI steps and MFE steps (kcal/mol) between neighbouring designs."""
        keep = self.on_hull if hull_only else np.ones(len(self.cai), bool)
        d_cai, d_mfe = np.diff(self.cai[keep]), np.abs(np.diff(self.mfe[keep]))
        return d_cai[d_cai > 0], d_mfe[d_mfe > 0]


def mfe_at_cai(front, bound, hull_only):
    """The lowest MFE among designs with CAI at least `bound`, or None."""
    keep = (front.cai >= bound) & (front.on_hull if hull_only else True)
    return front.mfe[keep].min() if keep.any() else None


def speedups(walk_cell, hull_cell, sweep_s):
    """(pruning speed-up, full-walk speed-up) of one protein; the first is None where the
    unpruned search did not finish.

    Pruning: derna2's unpruned fixed-λ DP at the hull vertices' λ values (the `derna2-hull`
    cell, measured) over the pruned hull search that finds them (`sweep_s`). The first does
    not search for the λ values, so this understates the speed-up. Full walk: the pruned hull
    search over the whole-front walk, below 1 where the walk costs more.
    """
    if not sweep_s:
        return None, None
    pruning = hull_cell["wall_s"] / sweep_s if hull_cell and hull_cell.get("exit") == "ok" else None
    return pruning, sweep_s / walk_cell["wall_s"]


def hull_dp_s(hull_cell):
    """(seconds, projected) of the fixed-λ DP at every hull vertex's λ (the `derna2-hull`
    cell): its wall time where it finished; where it stopped at its budget, its time per
    finished λ times the λ values planned, flagged projected. None where no λ finished.
    """
    if not hull_cell:
        return None
    if hull_cell.get("exit") == "ok":
        return hull_cell["wall_s"], False
    done = int(str(hull_cell.get("lambdas_done") or 0).split("/")[0])
    planned = hull_cell.get("lambdas_planned")
    if not done or not planned:
        return None
    return hull_cell["wall_s"] / done * planned, True
