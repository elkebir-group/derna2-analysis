"""The single-λ benchmark (E1): one row per cell, joining the work and time passes.

The work pass carries correctness (re-scored MFE and CAI), peak RSS and instruction
counts, none of which move with host load; the time pass carries wall-clock and CPU time,
measured on a reserved socket. A cell can have several time rows (retries after
contamination); it is read from its latest clean run, else its latest run, the same rule
the driver's own summary uses.

Memory and instructions come from the work pass, except for a cell the work pass stopped by
choice (`skipped`): its work row measured a run cut short, and the time pass ran it whole.
"""

from derna2_analysis.io import read_jsonl

LADDER = "P0DTC2_p"


def latest_work(rows):
    """A deferred cell that later ran has two rows; the latest is the cell's state."""
    return {r["cell"]: r for r in rows}


def best_time(rows):
    by_cell = {}
    for r in rows:
        by_cell.setdefault((r["cell"], r["rep"]), []).append(r)
    out = {}
    for (cell, _rep), runs in by_cell.items():
        clean = [r for r in runs if r["exit"] == "ok" and not r.get("contaminated")]
        out[cell] = dict((clean or runs)[-1], clean=bool(clean), attempts=len(runs))
    return out


def cells(work_path, time_path, reference_arm="derna-exact-30"):
    """One record per cell: the work row's correctness and memory, the time row's clock,
    and the gap of its re-scored MFE to the reference arm's on the same instance and corner.
    """
    work = latest_work(read_jsonl(work_path))
    time = best_time(read_jsonl(time_path))
    ref = {
        (r["protein"], r["corner"]): r
        for r in work.values()
        if r["arm"] == reference_arm and r["exit"] == "ok"
    }
    out = []
    for cell, w in work.items():
        t = time.get(cell, {})
        base = ref.get((w["protein"], w["corner"]), {})
        # A skipped cell's work row is a terminated run; its full run is the time row.
        m = t if w["exit"] == "skipped" and t.get("exit") == "ok" else w
        gap = None
        if m.get("mfe_rescored") is not None and base.get("mfe_rescored") is not None:
            gap = m["mfe_rescored"] - base["mfe_rescored"]
        out.append(
            {
                "cell": cell,
                "arm": w["arm"],
                "variant": w["variant"],
                "protein": w["protein"],
                "n_aa": w["n_aa"],
                "corner": w["corner"],
                "ladder": w["protein"].startswith(LADDER),
                "work_exit": w["exit"],
                "time_exit": t.get("exit"),
                "wall_s": t.get("wall_s"),
                "cpu_s": (t["user_s"] + t["sys_s"]) if t.get("user_s") is not None else None,
                "time_clean": t.get("clean"),
                "time_attempts": t.get("attempts"),
                "maxrss_kb": m.get("maxrss_kb"),
                "instructions": m.get("instructions"),
                "translates": m.get("translates"),
                "mfe_rescored": m.get("mfe_rescored"),
                "cai_rescored": m.get("cai_rescored"),
                "mfe_gap_kcal_mol": gap,
            }
        )
    return out
