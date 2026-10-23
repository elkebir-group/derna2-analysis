#!/usr/bin/env python3
"""The front-enumeration grid: one cell is a method's whole lambda sweep over
one instance.

The core grid's driver (`drive.py`) measures one lambda per cell. This one
keeps its measurement, isolation and re-scoring machinery and replaces the unit
of work: a cell here runs many lambdas under one wall-clock and memory budget,
writes a row per lambda as it finishes, and stops at whichever limit it reaches
first. A cell that runs out of budget keeps every lambda it completed, which is
the anytime curve.

Commands:
  hull        derna2's chord sweep per instance -> hull vertices + the replay
              lambdas (one interior lambda per vertex, plus the corners)
  verify      each comparator at interior lambdas against derna2's answer
  plan        the cells the grid owes, against the cells it has
  sweep       the grid: every (instance, arm) cell
  summarize   fronts, recovery against the reference, anytime curves

Every lambda is a derna2 convex weight. Each arm's own parameter is derived by
`map_lambda`, never typed into a config.
"""
import argparse
import contextlib
import hashlib
import json
import math
import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import drive as e1  # noqa: E402  the core grid's driver, reused wholesale

# The sibling copy must be the one that knows about mRNAfold: an older checkout
# of it silently maps mRNAfold's lambda like LinearCDSfold's, which is wrong by
# the ratio transform and would not fail until the designs disagreed.
if not hasattr(e1, "parse_mrnafold"):
    raise SystemExit(f"{HERE / 'drive.py'} is an older copy without mRNAfold "
                     "support; use the one from chore/paper-e2-plan")

MEASURE = HERE / "measure.py"
# DERNA v1 and LinearCDSfold both leave their CAI term unscaled where derna2
# carries a factor 100 (v1's `default.cpp` even has the `*100.0` commented
# out), so the same point needs the ratio scaled by 100. Measured exactly on
# designs at five interior lambdas for each, see `verify`.
#
# The core grid never saw this: at a corner the CAI term's scale cannot change
# the answer, because the second objective only breaks ties. Every interior
# lambda is wrong without it -- LinearCDSfold handed derna2's 0.5 returns the
# MFE-corner design.
CAI_SCALE_UNSCALED = 100.0

# The derna2 binaries, whose build commit decides what a derna2 row means. The
# comparators are frozen once and never rebuilt, so they need no such stamp.
PINNED_BINS = ("walk_front", "derna", "mfe_fold")


def provenance(cfg):
    """The pin every row is stamped with.

    A re-pin is the normal case, not an exception: the walk's levers land on
    branches of their own and the arm has to be rebuilt when they merge. A row
    that does not say which build produced it cannot be separated from one that
    does afterwards, so the stamp goes on at write time and `sweep` treats a
    cell from another pin as missing rather than done.

    The commit comes from the config and the hashes from the binaries on disk,
    so a config naming a pin the binaries were not built from is still visible
    in the rows.
    """
    # The concurrency goes on every row beside the load: a wall measured
    # with three other cells on the socket is a different measurement from one
    # measured alone, and afterwards only the row can say which it was.
    stamp = {"pin": cfg.get("derna2_pin"), "concurrency": cfg.get("concurrency", 1)}
    for name in PINNED_BINS:
        path = Path(cfg["bins"][name])
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        stamp[f"sha_{name}"] = digest[:12]
    return stamp


class PinnedResults(e1.Results):
    """`Results` that stamps every row with the build that produced it."""

    def __init__(self, path, stamp):
        super().__init__(path)
        self.stamp = stamp

    def append(self, row):
        super().append({**row, **self.stamp})


def append_rescored(rows, rescore, prot, results):
    """Append `rows` in order, each with its re-scoring when it has a design.

    The folds run in parallel off the timing cores and after the tool's
    measurement has ended, so they are no part of any wall; each row is
    written as soon as it and every row before it is scored."""
    with_seq = [r for r in rows if r.get("seq")]
    scores = rescore.many([r["seq"] for r in with_seq], prot)
    for row in rows:
        if row.get("seq"):
            row.update(next(scores))
        results.append(row)


def append_walk_points(cfg, rows, rescore, prot, results):
    """A walk's front rows, re-scored now or marked for `rescore` later.

    The walk reports each point's MFE and CAI on the same model and table the
    re-scoring uses, so re-scoring it is a check rather than a conversion,
    and done in line it holds the next walk back for longer than the walk
    itself takes at the long end. With `defer_walk_rescore` the rows go out
    at once and `rescore` fills `walk_rescore.jsonl` once no walk is timing."""
    if not cfg.get("defer_walk_rescore"):
        append_rescored(rows, rescore, prot, results)
        return
    for row in rows:
        results.append({**row, "rescore_deferred": True})


def walk_rescore_path(cfg):
    return Path(cfg["out_dir"]) / "walk_rescore.jsonl"


def point_key(r):
    return (r.get("pin"), r["arm"], r["protein"], r["point_index"])


def deferred_scores(cfg):
    """{point key: scores} from `walk_rescore.jsonl`."""
    path = walk_rescore_path(cfg)
    if not path.exists():
        return {}
    out = {}
    for line in open(path):
        if line.strip():
            r = json.loads(line)
            out[point_key(r)] = {k: r[k] for k in
                                 ("translates", "cai_rescored", "mfe_rescored")}
    return out


# ---- lambda maps -------------------------------------------------------------------

def map_lambda(kind, lam):
    """A derna2 convex weight as the parameter `kind` takes."""
    if kind in ("derna", "derna1rs"):
        # The Rust port of v1 scores lambda*E_dcal - (1-lambda)*round(100*sum
        # ln w): derna2's scale, so the same weight (`verify` checks it).
        return lam
    if kind in ("lineardesign", "mrnafold"):
        return (1.0 - lam) / lam
    if kind in ("derna1", "lcds"):
        return 1.0 / (1.0 + CAI_SCALE_UNSCALED * (1.0 - lam) / lam)
    raise ValueError(f"no lambda map for kind {kind}")


def replay_lambdas(breaks, eps):
    """One lambda per hull vertex plus both corners, in the order they run.

    A vertex is optimal over an interval between two breakpoints, so its
    midpoint has one answer while a breakpoint has two tied ones and the
    recovered vertex would be a function of the tool's tie-breaking. The order
    is corners first, then bisection over the vertex index, so a cell stopped
    early has covered the whole range coarsely rather than one end finely.
    """
    edges = [1.0 - eps] + list(breaks) + [eps]
    mids = [math.sqrt(edges[i] * edges[i + 1]) for i in range(len(edges) - 1)]
    order = [0, len(mids) - 1] if len(mids) > 1 else [0]

    def bisect_order(lo, hi):
        if hi - lo <= 1:
            return
        mid = (lo + hi) // 2
        order.append(mid)
        bisect_order(lo, mid)
        bisect_order(mid, hi)

    bisect_order(0, len(mids) - 1)
    seen, out = set(), []
    for i in order:
        if i not in seen:
            seen.add(i)
            out.append(mids[i])
    return out


# ---- one lambda of one arm ---------------------------------------------------------

def command_at(cfg, arm, label, prot, lam, work, stem):
    """(argv, cwd, stdin_path, parser) for one arm at one derna2 lambda."""
    kind = arm["kind"]
    bins = cfg["bins"]
    tool_lam = map_lambda(kind, lam)
    fasta = e1.ext(stem, ".fasta")
    fasta.write_text(f">{label}\n{prot}\n")
    if kind == "derna":
        lamf = e1.ext(stem, ".lambda.jsonl")
        lamf.write_text(json.dumps({"lambda": tool_lam}) + "\n")
        out = e1.ext(stem, ".out.jsonl")
        argv = [bins["derna"], "--fasta", str(fasta),
                "--codon-usage", cfg["tables"]["derna"],
                "--lambdas", str(lamf), "--out", str(out),
                "--machine", str(e1.ext(stem, ".machine.json")),
                "--threads", str(arm.get("threads", 1))]
        argv += arm.get("args", [])
        return (argv, str(e1.REPO), None,
                lambda _s: e1.parse_derna(out, e1.ext(stem, ".machine.json")))
    if kind == "derna1rs":
        lamf = e1.ext(stem, ".lambda.jsonl")
        lamf.write_text(json.dumps({"lambda": tool_lam}) + "\n")
        out = e1.ext(stem, ".out.jsonl")
        argv = [bins["derna1rs"], "--fasta", str(fasta),
                "--codon-usage", cfg["tables"]["derna"],
                "--lambdas", str(lamf), "--out", str(out)]
        return (argv, str(e1.REPO), None,
                lambda _s: e1.parse_derna(out, e1.ext(stem, ".machine.json")))
    if kind == "derna1":
        out = e1.ext(stem, ".out.txt")
        argv = [bins["derna1"], "-i", str(fasta), "-o", str(out),
                "-m", "1", "-s", "2", "-l", repr(tool_lam),
                "-c", bins["derna1_table"], "-d", bins["derna1_params"]]
        # v1 dumps its DP tables to fixed names in the working directory
        # (`E1.bin`, `M1.bin`, ...), so two runs sharing a cwd overwrite each
        # other's tables. Every run gets its own.
        rundir = e1.ext(stem, ".run")
        rundir.mkdir(parents=True, exist_ok=True)
        return (argv, str(rundir), None,
                lambda stdout: parse_derna1(stdout, out))
    if kind == "lcds":
        out = e1.ext(stem, ".out.txt")
        argv = [bins["lcds"], "-m", arm["mode"], "-O", "DN",
                "-l", repr(tool_lam), "-c", bins["lcds_table"],
                "-o", str(out), str(fasta)]
        if "width" in arm:
            argv += ["-b", str(arm["width"])]
        return (argv, bins["lcds_dir"], None,
                lambda _s: e1.parse_lcds(out.read_text()))
    if kind == "lineardesign":
        argv = [bins["lineardesign"], repr(tool_lam), "0",
                bins["lineardesign_table"]]
        return (argv, bins["lineardesign_dir"], str(e1.ext(stem, ".stdin")),
                e1.parse_lineardesign)
    if kind == "mrnafold":
        threads = arm.get("threads", 1)
        conf = e1.ext(stem, ".cfg")
        conf.write_text(
            f"aa_seq {prot}\nlambda {tool_lam!r}\nnum_subopt_traces 1\n"
            f"codon_freq_path {e1.mrnafold_table(cfg, work)}\n"
            f"parallel {'true' if threads > 1 else 'false'}\n")
        binary = bins["mrnafold_parallel"] if threads > 1 else bins["mrnafold"]
        weights = e1.codon_weights(
            e1.read_codon_table(cfg["tables"]["canonical"]))
        return ([binary, str(conf)], str(work), None,
                lambda stdout: e1.parse_mrnafold(stdout, tool_lam, weights))
    raise ValueError(f"unknown arm kind {kind}")


def parse_derna1(stdout, out_path):
    """v1's summary line from stdout, its design from its output file.

    Its energies are dekacal and its `cai` is sum(ln w) on weights normalized
    per amino acid, so `cai_s` is the geometric mean and needs no conversion.

    The design is the `zuker cai rna:` line of the output file. The file also
    carries a `zuker rna:` line, which is the MFE-only traceback and contains
    `X` placeholders where no codon was fixed -- not a sequence, and not the
    one this arm produced.
    """
    line = [l for l in stdout.splitlines() if l.startswith("lambda:")][-1]
    f = v1_summary_fields(line)
    row = dict(mfe_reported=float(f["mfe"]) / 100.0,
               cai_reported=float(f.get("cai_s", f.get("scai"))),
               sum_ln_w_reported=float(f["cai"]))
    for raw in Path(out_path).read_text().splitlines():
        if raw.startswith("zuker cai rna:"):
            seq = raw.split(":", 1)[1].split(".size")[0].strip()
            if seq and set(seq) <= set("ACGU"):
                row["seq"] = seq
            break
    return row


def run_lambda(cfg, arm, label, prot, lam, cores, work, stem, budget_s):
    """One lambda under measure.py, capped at the cell's remaining budget."""
    argv, cwd, stdin, parse = command_at(cfg, arm, label, prot, lam,
                                         work, stem)
    if stdin:
        Path(stdin).write_text(prot + "\n")
    meas = e1.ext(stem, ".measure.json")
    wrapped = ["taskset", "-c", ",".join(map(str, cores)),
               sys.executable, str(MEASURE), str(meas),
               "--timeout", str(max(1.0, budget_s)),
               "--stdout", str(e1.ext(stem, ".stdout")),
               "--stderr", str(e1.ext(stem, ".stderr"))]
    wrapped += (["--stdin", stdin] if stdin else []) + ["--"] + capped(cfg, argv)
    subprocess.run(wrapped, cwd=cwd)
    m = json.loads(meas.read_text())
    row = {k: m[k] for k in ("exit", "wall_s", "user_s", "sys_s",
                             "maxrss_kb", "instructions", "cycles")}
    row["cores"] = cores
    if row["exit"] == "ok":
        try:
            row.update(parse(e1.ext(stem, ".stdout").read_text()))
        except (AttributeError, IndexError, KeyError, ValueError, OSError,
                StopIteration) as err:
            row.update(exit="unparsed", note=repr(err))
    return row


# ---- the sweep ---------------------------------------------------------------------

def cell_limit_s(cfg, arm, n_aa):
    """The wall limit of one arm's cell on a protein of `n_aa` residues.

    `limit_s_by_aa` is a list of `[max_aa, seconds]` steps in increasing
    `max_aa`, the last one's `max_aa` null for every longer protein; without it
    the arm's `limit_s`, else the grid's, holds at every length.
    """
    for max_aa, seconds in arm.get("limit_s_by_aa", []):
        if max_aa is None or n_aa <= max_aa:
            return seconds
    return arm.get("limit_s", cfg["limit_s"])


def walk_soft_limit_s(cfg, n_aa):
    """The walk's soft limit on a protein of `n_aa` residues, or None.

    `walk_soft_limit_s_by_aa` has `limit_s_by_aa`'s steps. A walk past it is
    logged and keeps running to the hard `limit_s`; its cell row records what
    it had found by then (`by_soft_limit`).
    """
    for max_aa, seconds in cfg.get("walk_soft_limit_s_by_aa", []):
        if max_aa is None or n_aa <= max_aa:
            return seconds
    return None


def walk_band(cfg, n_aa):
    """The walk flags `walk_flags_by_aa.json` names for `n_aa` residues: a
    list, or "pending"; no file adds none."""
    path = Path(cfg["out_dir"]) / "walk_flags_by_aa.json"
    try:
        bands = json.loads(path.read_text())["bands"]
    except (OSError, ValueError, KeyError):
        return []
    return next((f for max_aa, f in bands
                 if max_aa is None or n_aa <= max_aa), [])


def walk_band_flags(cfg, n_aa, cell, log):
    """The flags the walk adds on a protein of `n_aa` residues, read when the
    walk is about to start from `walk_flags_by_aa.json` in the output
    directory: `{"bands": [[max_aa, flags], ...]}` with `limit_s_by_aa`'s
    steps, where `flags` is a list or `"pending"`. A pending band holds its
    walks until the file names its flags, so a band can wait on an A/B
    while the others run. Without the file no flags are added."""
    said = False
    while True:
        flags = walk_band(cfg, n_aa)
        if flags != "pending":
            if said:
                log(f"{cell}: band decided, flags {flags}")
            return list(flags)
        if not said:
            log(f"{cell}: holding, the walk flags for {n_aa} aa are pending")
            said = True
        time.sleep(60)


def by_soft_limit(stream, soft_s):
    """The vertices, gaps and gap-front points a walk's stream holds by
    `soft_s` seconds of its own clock."""
    vertices, gaps, points = set(), 0, 0
    path = Path(stream)
    for ln in (path.read_text().splitlines() if path.exists() else []):
        try:
            ev = json.loads(ln)
        except json.JSONDecodeError:
            continue
        if ev.get("elapsed_s", math.inf) > soft_s:
            continue
        if ev.get("event") == "vertex":
            vertices.add(round(ev["mfe_kcal_mol"] * 100))
        elif ev.get("event") == "gap":
            gaps += 1
            points += len(ev.get("front", []))
    return dict(vertices=len(vertices), gaps=gaps, gap_front_points=points)


def sweep_cell(cfg, arm_name, arm, label, prot, lambdas, cores, out_dir,
               results, rescore, log):
    """One (instance, arm) cell: its lambdas under one wall and memory budget.

    Stops at the first limit reached and keeps every lambda already done. The
    cell's exit says which limit ended it, so a short front is never confused
    with a complete one.
    """
    work = out_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    cell = f"{arm_name}.{label}"
    limit_s = cell_limit_s(cfg, arm, len(prot))
    cap_kb = cfg["mem_cap_gb"] * 1e6
    started = time.time()
    # Re-scoring runs on the cell's thread between lambdas but is no part of
    # the arm's cost, so its time is kept off the clock the limit reads.
    rescore_s = 0.0
    done = 0
    peak_kb = 0
    verdict = "ok"
    for i, lam in enumerate(lambdas):
        spent = time.time() - started - rescore_s
        if spent >= limit_s:
            verdict = "timeout"
            break
        stem = work / f"{cell}.{i:04d}"
        row = run_lambda(cfg, arm, label, prot, lam, cores, work, stem,
                         budget_s=limit_s - spent)
        row.update(arm=arm_name, kind=arm["kind"], protein=label,
                   n_aa=len(prot), lam=lam, tool_lam=map_lambda(arm["kind"],
                                                                lam),
                   lam_index=i, ncores=len(cores),
                   elapsed_at_finish_s=time.time() - started - rescore_s,
                   load=os.getloadavg()[0], running=RUNNING.n)
        if row.get("seq"):
            scoring = time.time()
            row.update(rescore(row["seq"], prot))
            rescore_s += time.time() - scoring
        results.append(row)
        done += 1
        peak_kb = max(peak_kb, row["maxrss_kb"] or 0)
        if row["maxrss_kb"] and row["maxrss_kb"] > cap_kb:
            verdict = "oom"
            log(f"{cell}: lambda {lam:.6g} peaked at "
                f"{row['maxrss_kb'] / 1e6:.1f} GB over the "
                f"{cfg['mem_cap_gb']} GB cap")
            break
        if row["exit"] == "timeout":
            verdict = "timeout"
            break
        if row["exit"] == "oom":
            verdict = "oom"
            break
        log(f"{cell}: {done}/{len(lambdas)} lambdas, "
            f"{(time.time() - started - rescore_s) / 60:.1f} min spent")
    summary = dict(kind="cell", arm=arm_name, protein=label, n_aa=len(prot),
                   lambdas_planned=len(lambdas), lambdas_done=done,
                   exit=verdict, wall_s=time.time() - started - rescore_s,
                   rescore_s=rescore_s, maxrss_kb=peak_kb,
                   ncores=len(cores), load=os.getloadavg()[0])
    results.append(summary)
    log(f"{cell}: {verdict} after {done}/{len(lambdas)} lambdas in "
        f"{summary['wall_s'] / 60:.1f} min")
    return summary


def hull_path(cfg, label):
    return Path(cfg["out_dir"]) / "hull" / f"{label}.vertices.json"


def has_hull(cfg, label):
    """Whether a replay arm has lambdas for `label`: a complete hull, or,
    under `replay_partial_hull`, the vertices a stopped walk streamed."""
    return hull_path(cfg, label).exists() or bool(
        cfg.get("replay_partial_hull")
        and partial_hull_path(cfg, label).exists())


class SoloGate:
    """Lets a cell above `solo_above_aa` residues run with the host to itself.

    The walk has its own memory gate in the hull pass, and the replay arms
    mostly need no such rule: at the longest
    instance their peaks are projected at 19, 14 and 6 GB against the one that
    will exceed any cap this host can honour. So this is off unless a config
    asks for it, because turning it on serializes the longest cells, and those
    are the cells most likely to run the full limit: it buys memory headroom
    with days of wall.

    A cell waiting to run alone blocks cells arriving after it, or a stream of
    short cells would hold it at the gate indefinitely.
    """

    def __init__(self):
        self.cv = threading.Condition()
        self.running = 0
        self.alone = False
        self.waiting = 0

    @contextlib.contextmanager
    def hold(self, alone):
        with self.cv:
            if alone:
                self.waiting += 1
                while self.running:
                    self.cv.wait()
                self.waiting -= 1
                self.alone = True
            else:
                while self.alone or self.waiting:
                    self.cv.wait()
            self.running += 1
        try:
            yield
        finally:
            with self.cv:
                self.running -= 1
                if alone:
                    self.alone = False
                self.cv.notify_all()


class PeakModel:
    """Peak RSS per n² for each arm, read at the longest protein whose runs
    have finished, or the kind's prior until one has. A cell's memory follows a
    clean power law in n (LinearDesign: 30.9 kB·n² from 377 to 454 aa), and a
    short run's peak is mostly the binary's fixed footprint, so the longest
    observation is the one to scale from."""

    def __init__(self, cfg, rows, margin=1.1):
        self.cfg, self.margin, self.k, self.lock = cfg, margin, {}, threading.Lock()
        for r in rows:
            if r.get("arm") in cfg["arms"] and r.get("exit") == "ok" \
                    and r.get("maxrss_kb"):
                self.observe(r["arm"], r["n_aa"], r["maxrss_kb"])

    def observe(self, arm, n, rss_kb):
        if not rss_kb:
            return
        with self.lock:
            self.k[arm] = max(self.k.get(arm, (0, 0.0)), (n, rss_kb / n**2))

    def estimate_kb(self, arm_name, n):
        arm = self.cfg["arms"][arm_name]
        with self.lock:
            seen = self.k.get(arm_name)
        k = seen[1] if seen else arm.get("rss_k", e1.RSS_K.get(arm["kind"], 11.0))
        return int(self.margin * k * n**2)


class Running:
    """How many sweep cells hold a slot right now. The configured concurrency
    is how many may run; a memory gate can leave slots idle, so each λ row
    records the count it actually ran beside."""

    def __init__(self):
        self.n, self.lock = 0, threading.Lock()

    def add(self, d):
        with self.lock:
            self.n += d


RUNNING = Running()


class SweepMemGate(e1.MemGate):
    """The core grid's gate, except that a cell estimated above the whole
    budget runs alone: it waits until no other cell holds a reservation and
    the host has its estimate free, and no smaller cell is admitted while it
    waits, or a stream of them would keep it out for good. Its own cgroup cap
    still ends it as `oom` if the estimate is short, which is a measurement,
    where refusing it would leave the cell neither run nor recorded."""

    def __init__(self, budget_kb, headroom_kb, budget_file=None, cap_kb=None):
        super().__init__(budget_kb, headroom_kb)
        self.default_budget, self.budget_file = budget_kb, budget_file
        self.cap, self.draining = cap_kb, 0

    def refresh(self):
        """The budget can be changed while the sweep runs by writing a number
        of GB to `budget_file`; without the file the configured one holds."""
        try:
            self.budget = float(Path(self.budget_file).read_text()) * 1e6
        except (TypeError, OSError, ValueError):
            self.budget = self.default_budget

    def fits(self, est):
        """Whether `est` could be admitted now (caller holds `cv`)."""
        if est <= self.budget:
            return (not self.draining and self.used + est <= self.budget
                    and e1.mem_available_kb() - est >= self.headroom)
        need = min(est, self.cap) if self.cap else est
        return self.used == 0 and e1.mem_available_kb() >= need

    def try_take(self, est):
        """Admit `est` if it fits now, without waiting."""
        self.refresh()
        with self.cv:
            if self.fits(est):
                self.used += est
                return True
        return False

    def take(self, est):
        self.refresh()
        if est <= self.budget:
            with self.cv:
                while not (not self.draining
                           and self.used + est <= self.budget
                           and e1.mem_available_kb() - est >= self.headroom):
                    self.cv.wait(timeout=30)
                    self.refresh()
                self.used += est
            return True
        # Memory the other sessions on the host hold is outside the budget,
        # so an empty grid is not enough: the estimate itself must be free.
        # Above the cgroup cap the cap is what the cell can use.
        need = min(est, self.cap) if self.cap else est
        with self.cv:
            self.draining += 1
            try:
                while not (self.used == 0
                           and e1.mem_available_kb() >= need):
                    self.cv.wait(timeout=30)
            finally:
                self.draining -= 1
                self.cv.notify_all()
            self.used += est
        return True


def await_memory(cfg, cell, log):
    """Hold a cell at the gate until the host has room for it to start.

    The per-cell cap bounds one cell; it cannot bound their sum. Four cells
    each far below `mem_cap_gb` can exhaust the host between them, and at that
    point no cgroup limit fires — the kernel's own OOM killer picks a victim by
    its own reckoning, which may be a cell behaving itself, or another
    session's job on this shared host.

    So a cell is admitted only while the host still has `mem_floor_gb`
    available. This does not bound a cell that grows after it starts, which
    nothing short of an estimate could; what it does is stop the grid from
    piling a new cell onto a host already under pressure, which is how the
    collisions arise. Waiting is logged, because a grid that silently waits
    looks the same as a grid that is slow.
    """
    floor_kb = cfg.get("mem_floor_gb", 12) * 1e6
    if not floor_kb:
        return
    waited = 0.0
    while e1.mem_available_kb() < floor_kb:
        if waited == 0.0:
            log(f"{cell}: waiting for memory, "
                f"{e1.mem_available_kb() / 1e6:.0f} GB free of the "
                f"{cfg.get('mem_floor_gb', 12)} GB floor")
        time.sleep(15)
        waited += 15
    if waited:
        log(f"{cell}: admitted after {waited / 60:.1f} min at the memory gate")


def capped(cfg, argv):
    """`argv` in a memory-capped cgroup of its own, or unchanged.

    A cap the kernel enforces is the only one that holds: a cell that would
    outgrow the host is killed at `mem_cap_gb` instead of taking the host's
    other cells with it, and `measure.py` reads that SIGKILL as `oom`. The
    scope wraps the tool alone, inside the measurement, so the cell's cgroup
    holds nothing whose death would lose the measurement, and the scope's own
    cost is about 4 ms, measured.

    Swap is refused outright: a cell that swaps is no longer measuring the
    algorithm.
    """
    gb = cfg.get("mem_cap_gb")
    if not gb or not cfg.get("mem_scope", True):
        return argv
    return ["systemd-run", "--user", "--scope", "-q",
            "-p", f"MemoryMax={gb}G", "-p", "MemorySwapMax=0"] + argv


def partial_hull_path(cfg, label):
    return Path(cfg["out_dir"]) / "hull" / f"{label}.vertices.partial.json"


def hull_lambdas(cfg, label):
    """The replay lambdas for one instance, from its hull run.

    With `replay_partial_hull` set, an instance whose walk was stopped takes
    its lambdas from the vertices the walk streamed before it stopped.
    """
    path = hull_path(cfg, label)
    if not path.exists() and cfg.get("replay_partial_hull"):
        path = partial_hull_path(cfg, label)
    if not path.exists():
        raise SystemExit(f"no hull for {label}: run `hull` first ({path})")
    rec = json.loads(path.read_text())
    return replay_lambdas(rec["breakpoints"], cfg["epsilon"])


def run_sweep(cfg, args, log):
    out_dir = Path(cfg["out_dir"])
    stamp = provenance(cfg)
    results = PinnedResults(out_dir / "sweep.jsonl", stamp)
    rescore = e1.Rescorer(cfg)
    # A cell built by another pin is not done: it is a measurement of a
    # different derna2, and re-running is the only way to make the grid one
    # experiment again.
    done = {(r["arm"], r["protein"]) for r in results.rows()
            if r.get("kind") == "cell" and r.get("pin") == stamp["pin"]}
    pool = e1.CorePool(cfg["cores"])
    solo = SoloGate()
    # The sum of the running cells' estimated peaks stays within the budget,
    # and a cell starts only while the host keeps the headroom after its own
    # peak: a floor checked at admission alone let eight LinearDesign cells of
    # 7 GB each be admitted onto a 61 GB host shared with other jobs.
    peaks = PeakModel(cfg, results.rows())
    gate = SweepMemGate(cfg["mem_budget_gb"] * 1e6,
                        cfg.get("mem_headroom_gb", 6) * 1e6,
                        out_dir / "mem_budget_gb",
                        cfg["mem_cap_gb"] * 1e6)
    want = set(args.arms.split(",")) if args.arms else None

    cells = []
    for label, prot in e1.instances(cfg):
        if args.max_n and len(prot) > args.max_n:
            continue
        for name, arm in cfg["arms"].items():
            if not arm.get("enabled", True) or (want and name not in want):
                continue
            if not in_phase(args, arm.get("threads", 1)):
                continue
            if (name, label) in done:
                continue
            if arm.get("proteins") and label not in arm["proteins"]:
                continue
            if (args.walked_only and arm.get("lambdas", "hull") == "hull"
                    and not has_hull(cfg, label)):
                continue
            cells.append((name, arm, label, prot))
    # Longest first, so the run that starts last is the short one; or
    # shortest first, so most cells are in hand early and memory stays low.
    cells.sort(key=lambda c: len(c[3]) if args.shortest_first else -len(c[3]))
    if args.arm_order:
        # Whole arms in the order given (the faster comparators first), each
        # arm's cells in the order above; arms not named go last.
        rank = {a: i for i, a in enumerate(args.arm_order.split(","))}
        cells.sort(key=lambda c: rank.get(c[0], len(rank)))
    # A cell that needs the host to itself goes last, so that the grid idles
    # only at its tail while the cell waits for the memory to come free.
    alone = {tuple(c) for c in cfg.get("alone_cells", [])}
    cells.sort(key=lambda c: (c[0], c[2]) in alone)
    log(f"{len(cells)} cells to run")
    # An arm that has used its whole limit on one protein will use it again on
    # a longer one under the same or a smaller limit, so with
    # --skip-after-timeout those cells are left out; a longer protein given a
    # larger limit still runs. No row is written for a skipped cell: it stays
    # undone for a later pass to run.
    timed_out = {}
    timed_out_lock = threading.Lock()

    def one(spec):
        name, arm, label, prot = spec
        # A cell whose projected peak no host here can hold is recorded as
        # not run, with the projection, rather than started to be killed.
        not_run = cfg.get("not_run_cells", {}).get(name, {}).get(label)
        if not_run:
            log(f"{name}.{label}: not run, {not_run}")
            results.append(dict(kind="cell", arm=name, protein=label,
                                n_aa=len(prot), lambdas_planned=0,
                                lambdas_done=0, exit="not_run_mem",
                                reason=not_run, wall_s=0.0, ncores=0,
                                load=os.getloadavg()[0]))
            return None
        if args.skip_after_timeout:
            with timed_out_lock:
                first = timed_out.get(name)
            if (first is not None and len(prot) >= first[1]
                    and cell_limit_s(cfg, arm, len(prot)) <= first[2]):
                log(f"{name}.{label}: skipped, {name} timed out on {first[0]} "
                    f"({first[1]} aa, {first[2] / 3600:g} h limit)")
                return None
        est = peaks.estimate_kb(name, len(prot))
        waited = time.time()
        gate.take(est)
        if time.time() - waited > 60:
            log(f"{name}.{label}: admitted after {(time.time() - waited) / 60:.1f} "
                f"min at the memory gate (estimated {est / 1e6:.1f} GB)")
        RUNNING.add(1)
        try:
            summary = run_cell(spec)
        finally:
            RUNNING.add(-1)
            gate.give(est)
        if summary:
            peaks.observe(name, len(prot), summary.get("maxrss_kb"))
        if summary and summary.get("exit") == "timeout":
            with timed_out_lock:
                # The largest limit an arm has run out of governs: it covers
                # every cell whose limit is no larger.
                limit = cell_limit_s(cfg, arm, len(prot))
                old = timed_out.get(name)
                if (old is None or limit > old[2]
                        or (limit == old[2] and len(prot) < old[1])):
                    timed_out[name] = (label, len(prot), limit)
        return summary

    def run_cell(spec):
        name, arm, label, prot = spec
        replay = arm.get("lambdas", "hull") == "hull"
        if replay and not has_hull(cfg, label):
            # One instance whose walk timed out must not take the grid with
            # it: the cell is recorded as blocked and the other 538 go on.
            # A blocked cell is not a measurement of this arm and `summarize`
            # must not read it as a short front.
            log(f"{name}.{label}: blocked, no hull for {label}")
            results.append(dict(kind="cell", arm=name, protein=label,
                                n_aa=len(prot), lambdas_planned=0,
                                lambdas_done=0, exit="no_hull", wall_s=0.0,
                                ncores=0, load=os.getloadavg()[0]))
            return None
        lambdas = hull_lambdas(cfg, label) if replay else None
        if lambdas is None:
            return run_own_sweep(cfg, name, arm, label, prot, pool, out_dir,
                                 results, rescore, log)
        n = arm.get("threads", 1) * arm.get("workers", 1)
        above = cfg.get("solo_above_aa")
        with solo.hold(bool(above) and len(prot) >= above):
            await_memory(cfg, f"{name}.{label}", log)
            cores = pool.take(n)
            try:
                return sweep_cell(cfg, name, arm, label, prot, lambdas, cores,
                                  out_dir, results, rescore, log)
            finally:
                pool.give(cores)

    with ThreadPoolExecutor(max_workers=cfg.get("concurrency", 1)) as ex:
        list(ex.map(one, cells))


def run_own_sweep(cfg, name, arm, label, prot, pool, out_dir, results,
                  rescore, log):
    """DERNA v1 choosing its own lambdas, the one arm that pays for its own
    lambda search rather than being handed derna2's.

    `-s 3`, not `-s 4`: v1's dispatch is crossed against its own help text.
    The flag its `-h` calls "lambda sweep" (3) runs `lambda_sweep_2`, the
    adaptive search driven by the `-t`/`-p` thresholds and needing no lambda;
    the undocumented 4 runs `lambda_sweep`, a uniform grid that requires both
    `-l` and `-a` and aborts with "Invalid Value of lambda" without them.
    """
    work = out_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    stem = work / f"{name}.{label}"
    fasta = e1.ext(stem, ".fasta")
    fasta.write_text(f">{label}\n{prot}\n")
    if arm["kind"] == "derna1rs":
        return run_own_sweep_rs(cfg, name, arm, label, prot, pool, stem,
                                fasta, results, rescore, log)
    csv = e1.ext(stem, ".sweep")
    argv = [cfg["bins"]["derna1"], "-i", str(fasta),
            "-o", str(e1.ext(stem, ".out.txt")), "-m", "1", "-s", "3",
            "-O", str(csv), "-c", cfg["bins"]["derna1_table"],
            "-d", cfg["bins"]["derna1_params"],
            "-t", repr(arm.get("tau1", 0.0025)),
            "-p", repr(arm.get("tau2", 0.00075))]
    # Its own cwd: v1 writes fixed-name DP dumps beside wherever it runs.
    rundir = e1.ext(stem, ".run")
    rundir.mkdir(parents=True, exist_ok=True)
    meas = e1.ext(stem, ".measure.json")
    await_memory(cfg, f"{name}.{label}", log)
    cores = pool.take(1)
    running = RUNNING.n
    try:
        wrapped = ["taskset", "-c", ",".join(map(str, cores)),
                   sys.executable, str(MEASURE), str(meas),
                   "--timeout", str(cell_limit_s(cfg, arm, len(prot))),
                   "--stdout", str(e1.ext(stem, ".stdout")),
                   "--stderr", str(e1.ext(stem, ".stderr")), "--"] + capped(cfg, argv)
        subprocess.run(wrapped, cwd=str(rundir))
    finally:
        pool.give(cores)
    m = json.loads(meas.read_text())
    rows = read_v1_sweep(csv)
    # v1 writes its CSV through a buffered stream and flushes it at close, so a
    # run killed at the limit can leave the file empty however many lambdas it
    # finished -- and this is the arm most likely to be killed. Its stdout is
    # flushed per lambda, so the completed lambdas survive there. Those rows
    # carry no design, so they cannot be re-scored: they are kept for the
    # anytime curve on reported values and marked as such.
    if len(rows) < count_stdout_lambdas(e1.ext(stem, ".stdout")):
        from_stdout = read_v1_stdout(e1.ext(stem, ".stdout"))
        have = {round(r["lam"], 12) for r in rows}
        for r in from_stdout:
            if round(r["lam"], 12) not in have:
                r["reported_only"] = True
                rows.append(r)
        log(f"{name}.{label}: CSV held {len(have)} lambdas, stdout "
            f"{len(from_stdout)}; kept the difference as reported-only")
    for r in rows:
        r.update(arm=name, kind=arm["kind"], protein=label, n_aa=len(prot),
                 ncores=1)
    append_rescored(rows, rescore, prot, results)
    summary = dict(kind="cell", arm=name, protein=label, n_aa=len(prot),
                   lambdas_planned=None, lambdas_done=len(rows),
                   exit=m["exit"], wall_s=m["wall_s"], ncores=1,
                   maxrss_kb=m["maxrss_kb"], instructions=m["instructions"],
                   load=os.getloadavg()[0], running=running)
    results.append(summary)
    log(f"{name}.{label}: {m['exit']}, {len(rows)} lambdas of its own in "
        f"{summary['wall_s'] / 60:.1f} min")
    return summary


def run_own_sweep_rs(cfg, name, arm, label, prot, pool, stem, fasta, results,
                     rescore, log):
    """The Rust port of DERNA v1 choosing its own lambdas (`run_glpe78`).

    Its search is v1's adaptive one with v1's thresholds, read on the port's
    CAI scale, so it places more lambdas than the C++ binary does; the arm
    reports the port's search as it is. `--stream` writes each solve as it
    finishes, synced, so a cell stopped at the limit keeps every lambda it
    solved; the rows are read from the stream whether or not the run ended.
    """
    stream = e1.ext(stem, ".stream.jsonl")
    # `--cai-scale 100` places the sweep's weights where DERNA C++ places
    # them. The port quantizes CAI 100x the way C++ does and the sweep's
    # thresholds are widths in lambda, so without it the arm measures a search
    # no released tool performs -- a denser one, whose extra weights cost time
    # and buy points C++ would not have found.
    argv = [cfg["bins"]["derna1rs_sweep"], "--fasta", str(fasta),
            "--codon-usage", cfg["tables"]["derna"],
            "--cai-scale", repr(float(arm.get("cai_scale", 100))),
            "--out", str(e1.ext(stem, ".front.jsonl")),
            "--stream", str(stream)]
    meas = e1.ext(stem, ".measure.json")
    await_memory(cfg, f"{name}.{label}", log)
    cores = pool.take(1)
    running = RUNNING.n
    try:
        subprocess.run(
            ["taskset", "-c", ",".join(map(str, cores)),
             sys.executable, str(MEASURE), str(meas),
             "--timeout", str(cell_limit_s(cfg, arm, len(prot))),
             "--stdout", str(e1.ext(stem, ".stdout")),
             "--stderr", str(e1.ext(stem, ".stderr")), "--"] + capped(cfg, argv),
            cwd=str(e1.REPO))
    finally:
        pool.give(cores)
    m = json.loads(meas.read_text())
    rows = []
    for ln in (stream.read_text().splitlines() if stream.exists() else []):
        try:
            rec = json.loads(ln)
        except json.JSONDecodeError:
            continue
        rows.append(dict(lam=rec["lambda"], seq=rec["rna_seq"],
                         mfe_reported=rec["mfe_kcal_mol"],
                         cai_reported=rec["cai"], solve_s=rec.get("solve_s"),
                         elapsed_at_finish_s=rec.get("elapsed_s")))
    for r in rows:
        r.update(arm=name, kind=arm["kind"], protein=label, n_aa=len(prot),
                 ncores=1)
    append_rescored(rows, rescore, prot, results)
    summary = dict(kind="cell", arm=name, protein=label, n_aa=len(prot),
                   lambdas_planned=None, lambdas_done=len(rows),
                   exit=m["exit"], wall_s=m["wall_s"], ncores=1,
                   user_s=m.get("user_s"), sys_s=m.get("sys_s"),
                   maxrss_kb=m["maxrss_kb"], instructions=m["instructions"],
                   load=os.getloadavg()[0], running=running)
    results.append(summary)
    log(f"{name}.{label}: {m['exit']}, {len(rows)} lambdas of its own in "
        f"{summary['wall_s'] / 60:.1f} min")
    return summary


def v1_summary_fields(line):
    """`lambda: L,O: ...,cai: C,cai_s: S,mfe: M,...` as a lower-cased dict.

    Case matters: v1's fixed-lambda mode prints `cai`/`cai_s`/`mfe` while its
    sweep prints `CAI`/`sCAI`/`MFE` for the same quantities, so a parser
    written against one silently drops every row of the other.
    """
    return dict((k.strip().lower(), v.strip()) for k, v in
                (p.split(":", 1) for p in line.split(",") if ":" in p))


def count_stdout_lambdas(path):
    if not Path(path).exists():
        return 0
    return sum(1 for l in open(path, errors="replace")
               if l.startswith("lambda:"))


def read_v1_stdout(path):
    """Completed lambdas from v1's own progress output.

    One line per lambda, `lambda: L,O: ...,cai: C,cai_s: S,mfe: M,...`, with
    the energy in dekacal and `cai_s` the geometric mean. No design, so no
    re-score is possible from this source.
    """
    if not Path(path).exists():
        return []
    out = []
    for line in open(path, errors="replace"):
        if not line.startswith("lambda:"):
            continue
        f = v1_summary_fields(line)
        try:
            out.append(dict(lam=float(f["lambda"]),
                            mfe_reported=float(f["mfe"]) / 100.0,
                            cai_reported=float(f["scai"])))
        except (KeyError, ValueError):
            continue
    return out


def read_v1_sweep(csv):
    """v1's sweep CSV, whatever subset of it the run managed to write.

    A truncated final line is dropped rather than guessed at: the point of the
    partial record is that it is honest about where the run stopped.
    """
    if not Path(csv).exists():
        for alt in (Path(f"{csv}.csv"), Path(f"{csv}")):
            if alt.exists():
                csv = alt
                break
        else:
            return []
    lines = Path(csv).read_text().splitlines()
    if not lines:
        return []
    head = [h.strip() for h in lines[0].split(",")]
    out = []
    for line in lines[1:]:
        parts = line.split(",")
        if len(parts) != len(head):
            continue
        rec = dict(zip(head, parts))
        try:
            row = dict(lam=float(rec["lambda"]))
        except (KeyError, ValueError):
            continue
        for src, dst, conv in (("rc_mfe", "mfe_reported", float),
                               ("mfe", "mfe_reported", float),
                               ("rc_sCAI", "cai_reported", float),
                               ("sCAI", "cai_reported", float),
                               ("rc_seq", "seq", str),
                               ("seq", "seq", str)):
            if src in rec and dst not in row:
                try:
                    row[dst] = conv(rec[src])
                except ValueError:
                    pass
        if "mfe_reported" in row:
            row["mfe_reported"] /= 100.0
        out.append(row)
    return out


# ---- commands ----------------------------------------------------------------------

def run_hull(cfg, args, log):
    """The hull-gap walk per instance: the exact front, and the replay lambdas.

    This pass is both the experiment's reference and one of its timed arms. The
    walk is one run per instance -- it enumerates the whole front in a single
    invocation rather than one design per lambda -- so it cannot be a cell of
    the sweep, whose unit of work is a lambda. Timing it here instead of
    re-running the identical command inside `sweep` also keeps the most
    expensive arm from being paid for twice, and it has to run first in any
    case: every replay arm takes its lambdas from this pass's output.

    It is measured through the same wrapper as every lambda of every other arm,
    and writes the same rows: one per front point, plus the cell summary.
    """
    out_dir = Path(cfg["out_dir"]) / "hull"
    work = out_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    pool = e1.CorePool(cfg["cores"])
    stamp = provenance(cfg)
    results = PinnedResults(Path(cfg["out_dir"]) / "sweep.jsonl", stamp)
    rescore = e1.Rescorer(cfg)
    rows = results.rows()
    # A walk is done when this pin has its cell; with --rewalk only a cell
    # from the walk binary on disk counts, so a rebuilt walk re-walks every
    # instance. A re-walk of an instance whose hull is already present writes
    # beside it under the arm's name: the replay arms keep the hull they were
    # timed against, and the two fronts can be compared.
    done = {(r["arm"], r["protein"]) for r in rows
            if r.get("kind") == "cell" and r.get("pin") == stamp["pin"]
            and (not args.rewalk
                 or r.get("sha_walk_front") == stamp["sha_walk_front"])}
    specs = []
    for label, prot in e1.instances(cfg):
        if args.max_n and len(prot) > args.max_n:
            continue
        fasta = out_dir / f"{label}.fasta"
        fasta.write_text(f">{label}\n{prot}\n")
        for i, (arm_name, threads) in enumerate(walk_arms(cfg)):
            if args.arms and arm_name not in args.arms.split(","):
                continue
            if not in_phase(args, threads):
                continue
            present = (out_dir / f"{label}.vertices.json").exists()
            if (arm_name, label) in done and (args.rewalk or i > 0):
                log(f"{arm_name}.{label}: already run")
                continue
            if i == 0 and present and not args.rewalk:
                log(f"{label}: hull already present")
                continue
            specs.append((arm_name, threads, label, prot, fasta,
                          i == 0 and not present))
    # Listed instances first, each group shortest first when asked; a walk
    # whose band is pending goes when its band is decided.
    first = args.first.split(",") if args.first else []
    specs.sort(key=lambda sp: (sp[2] not in first,
                               len(sp[3]) if args.shortest_first else 0))
    log(f"{len(specs)} walks to run, {cfg.get('walk_concurrency', 1)} at a time")
    # Memory per walk: k·n² read at the longest walk this binary has
    # finished, else the config's prior (the binary's peak per n² near
    # 470 aa), with PeakModel's margin.
    walk_k = {}
    for r in rows:
        if (r.get("kind") == "cell" and r.get("sha_walk_front")
                == stamp["sha_walk_front"] and r.get("maxrss_kb")
                and r.get("exit") == "ok"):
            walk_k[r["arm"]] = max(walk_k.get(r["arm"], (0, 0.0)),
                                   (r["n_aa"], r["maxrss_kb"] / r["n_aa"] ** 2))
    k_lock = threading.Lock()
    gate = SweepMemGate(cfg["mem_budget_gb"] * 1e6,
                        cfg.get("mem_headroom_gb", 6) * 1e6,
                        Path(cfg["out_dir"]) / "mem_budget_gb")

    # Each free slot takes the first walk in order that the memory gate admits
    # now: a walk too large for the memory left is passed over for a smaller
    # one behind it, and taken once enough memory comes free. A walk larger
    # than the whole budget runs alone, so it is taken only when no other walk
    # is left to take, and then waits for the host to empty. A walk whose
    # flags are pending is passed over the same way until its band is written.
    queue, q_lock = list(specs), threading.Lock()
    passed, held = set(), set()

    def estimate(spec):
        with k_lock:
            k = walk_k.get(spec[0], (0, cfg.get("walk_rss_k", 14.0)))[1]
        return int(1.1 * k * len(spec[3]) ** 2)

    def take_next():
        while True:
            with q_lock:
                if not queue:
                    return None
                ready = []
                for spec in queue:
                    if walk_band(cfg, len(spec[3])) == "pending":
                        if spec[2] not in passed:
                            passed.add(spec[2])
                            log(f"{spec[0]}.{spec[2]}: passed over, the walk "
                                f"flags for {len(spec[3])} aa are pending")
                        continue
                    ready.append(spec)
                gate.refresh()
                small = [sp for sp in ready if estimate(sp) <= gate.budget]
                for spec in (small or ready):
                    est = estimate(spec)
                    if gate.try_take(est):
                        queue.remove(spec)
                        return spec, est
                    if spec[2] not in held:
                        held.add(spec[2])
                        log(f"{spec[0]}.{spec[2]}: passed over for now, "
                            f"estimated {est / 1e6:.1f} GB does not fit the "
                            f"memory left")
            with gate.cv:
                gate.cv.wait(timeout=30)

    def worker():
        while (got := take_next()) is not None:
            one(*got)

    def one(spec, est):
        arm_name, threads, label, prot, fasta, reference = spec
        n = len(prot)
        RUNNING.add(1)
        try:
            walk_cell(cfg, arm_name, threads, label, prot, fasta, out_dir,
                      work, pool, results, rescore, log, reference=reference)
        finally:
            RUNNING.add(-1)
            gate.give(est)
        last = next((r for r in reversed(results.rows())
                     if r.get("kind") == "cell" and r.get("arm") == arm_name
                     and r.get("protein") == label), None)
        if last and last.get("exit") == "ok" and last.get("maxrss_kb"):
            with k_lock:
                walk_k[arm_name] = max(walk_k.get(arm_name, (0, 0.0)),
                                       (n, last["maxrss_kb"] / n ** 2))

    with ThreadPoolExecutor(max_workers=cfg.get("walk_concurrency", 1)) as ex:
        for f in [ex.submit(worker)
                  for _ in range(cfg.get("walk_concurrency", 1))]:
            f.result()


def in_phase(args, threads):
    """Whether an arm on `threads` threads runs in this invocation's phase.

    The multi-core arms run in a phase of their own, after the single-core
    grid and one cell at a time, so a cell holding every timing core never
    shares the socket with another cell. `--phase` unset runs both.
    """
    phase = getattr(args, "phase", None)
    if phase is None:
        return True
    return (threads > 1) == (phase == "par")


def walk_arms(cfg):
    """(name, threads) of every walk arm, the reference first.

    The reference walk's hull is the one every replay arm takes its lambdas
    from; the others enumerate the same front and are timed only. A config
    without `walk_arms` has one walk, `walk_arm` on `hull_threads`.
    """
    arms = cfg.get("walk_arms") or {
        cfg.get("walk_arm", "derna2-walk"): {"threads": cfg.get("hull_threads", 1)}}
    return [(name, spec.get("threads", 1)) for name, spec in arms.items()
            if spec.get("enabled", True)]


def walk_cell(cfg, arm_name, threads, label, prot, fasta, out_dir, work, pool,
              results, rescore, log, reference):
    """One walk arm on one instance: its front rows and cell row, and, for the
    reference walk, the hull the replay arms read."""
    # The reference walk keeps the unprefixed names the replay arms and
    # `summarize` read; every other walk writes beside it under its own name.
    base = out_dir / (label if reference else f"{label}.{arm_name}")
    dest = out_dir / f"{label}.vertices.json"
    report = Path(f"{base}.walk.json")
    front = Path(f"{base}.front.jsonl")
    stream = Path(f"{base}.stream.jsonl")
    # The walk's levers come from the config so the timed walk is the one
    # the run names; `--stream` writes each vertex and gap as it is found,
    # so a walk the limit stops still reports what it found.
    flags = list(cfg.get("walk_flags", [])) + walk_band_flags(
        cfg, len(prot), f"{arm_name}.{label}", log)
    argv = [cfg["bins"]["walk_front"], "--fasta", str(fasta),
            "--n", str(len(prot)),
            "--lo", repr(cfg["epsilon"]), "--hi", repr(1 - cfg["epsilon"]),
            "--threads", str(threads)] + flags + [
            "--out", str(report), "--front", str(front),
            "--stream", str(stream)]
    stem = work / f"{arm_name}.{label}"
    meas = e1.ext(stem, ".measure.json")
    started = time.time()
    await_memory(cfg, f"{arm_name}.{label}", log)
    cores = pool.take(threads)
    running = RUNNING.n
    soft_s = walk_soft_limit_s(cfg, len(prot))
    soft_timer = None
    if soft_s is not None:
        soft_timer = threading.Timer(soft_s, log, args=(
            f"{arm_name}.{label}: past the soft limit of {soft_s / 3600:.1f} h; "
            f"still running to the hard limit (reported, not stopped)",))
        soft_timer.daemon = True
        soft_timer.start()
    try:
        subprocess.run(
            ["taskset", "-c", ",".join(map(str, cores)),
             sys.executable, str(MEASURE), str(meas),
             "--timeout", str(cfg["limit_s"]),
             "--stdout", str(e1.ext(stem, ".stdout")),
             "--stderr", str(e1.ext(stem, ".stderr")), "--"] + capped(cfg, argv),
            cwd=str(e1.REPO))
    finally:
        pool.give(cores)
        if soft_timer is not None:
            soft_timer.cancel()
    m = json.loads(meas.read_text())
    soft = {} if soft_s is None else dict(
        soft_limit_s=soft_s, over_soft_limit=m["wall_s"] > soft_s,
        by_soft_limit=by_soft_limit(stream, soft_s))
    soft.update(walk_flags=flags, running=running, reference=reference)
    wall_min = (time.time() - started) / 60
    log(f"{arm_name}.{label}: walk {m['exit']} in {wall_min:.1f} min")
    if m["exit"] != "ok":
        # No report means no certified hull: record what the stream holds
        # as a partial cell, write its vertices beside (not as) the
        # `.vertices.json`, and let the replay arms refuse the instance
        # unless the config opts into the partial hull's lambdas.
        vertices, breaks, gaps = partial_from(stream, len(prot))
        rows = [dict(kind="front_point", arm=arm_name, protein=label,
                     n_aa=len(prot), point_index=i, ncores=len(cores),
                     partial=True, found_s=gap.get("elapsed_s"),
                     mfe=pt.get("mfe_kcal_mol"), cai=pt.get("cai"),
                     seq=pt.get("sequence"))
                for i, (gap, pt) in enumerate(
                    (gap, pt) for gap in gaps for pt in gap.get("front", []))]
        append_walk_points(cfg, rows, rescore, prot, results)
        npoints = len(rows)
        results.append(dict(kind="cell", arm=arm_name, protein=label,
                            n_aa=len(prot), lambdas_planned=0,
                            lambdas_done=npoints, exit=m["exit"],
                            wall_s=m["wall_s"], user_s=m.get("user_s"),
                            sys_s=m.get("sys_s"),
                            maxrss_kb=m.get("maxrss_kb"),
                            ncores=len(cores), partial=True,
                            vertices=len(vertices), gaps=len(gaps),
                            front_points=npoints,
                            load=os.getloadavg()[0], **soft))
        if vertices and reference:
            partial_hull_path(cfg, label).write_text(json.dumps(
                dict(instance=label, n_aa=len(prot), complete=False,
                     vertices=vertices, breakpoints=breaks,
                     replay_lambdas=replay_lambdas(breaks, cfg["epsilon"]),
                     gaps_solved=len(gaps), walk_exit=m["exit"]),
                indent=1))
        return
    points = [json.loads(ln) for ln in
              front.read_text().splitlines() if ln.strip()]
    rows = [dict(kind="front_point", arm=arm_name, protein=label,
                 n_aa=len(prot), point_index=i, ncores=len(cores),
                 mfe=pt.get("mfe_kcal_mol"), cai=pt.get("cai"),
                 seq=pt.get("sequence"))
            for i, pt in enumerate(points)]
    append_walk_points(cfg, rows, rescore, prot, results)
    vertices, breaks = vertices_from(report, len(prot))
    walk = json.loads(report.read_text())
    # One run, so the whole front shares the cell's wall and peak: the arm
    # has no per-lambda cost to report and `lambdas_done` would be a
    # category error. The point count is what it produced.
    results.append(dict(kind="cell", arm=arm_name, protein=label,
                        n_aa=len(prot), lambdas_planned=0,
                        lambdas_done=len(points), exit=m["exit"],
                        wall_s=m["wall_s"], user_s=m["user_s"],
                        sys_s=m["sys_s"], maxrss_kb=m["maxrss_kb"],
                        instructions=m["instructions"],
                        cycles=m["cycles"], ncores=len(cores),
                        vertices=len(vertices), front_points=len(points),
                        load=os.getloadavg()[0], **soft))
    if not reference:
        return
    dest.write_text(json.dumps(
        dict(instance=label, n_aa=len(prot), vertices=vertices,
             breakpoints=breaks,
             replay_lambdas=replay_lambdas(breaks, cfg["epsilon"]),
             walk_wall_s=walk.get("wall_s"),
             walk_peak_rss_kb=walk.get("peak_rss_kb"),
             walk_points=walk.get("points")), indent=1))




def vertices_from(report, n_aa):
    """Hull vertices and the breakpoint lambdas between them.

    The walk's report lists vertices by energy and adaptation, and carries the
    lambdas in its gap records: a gap's `lambda_star` is the weight at which
    its two vertices tie, which is exactly a breakpoint of the hull. V
    vertices therefore come with V-1 breakpoints, and the interval each vertex
    is optimal over sits between two of them.

    `adaptation_log` is sum(ln w), so the geometric-mean CAI is
    exp(adaptation_log / n).
    """
    rec = json.loads(Path(report).read_text())
    if "vertices" not in rec or "gaps" not in rec:
        raise SystemExit(f"{report} has no vertices/gaps; keys: {sorted(rec)}")
    vertices = []
    for v in rec["vertices"]:
        log_w = v.get("adaptation_log")
        vertices.append(dict(
            mfe=v.get("mfe_kcal_mol"), sum_ln_w=log_w,
            cai=math.exp(log_w / n_aa) if log_w is not None else None))
    vertices.sort(key=lambda v: v["mfe"])
    breaks = sorted((g["lambda_star"] for g in rec["gaps"]
                     if g.get("lambda_star") is not None), reverse=True)
    if len(breaks) != len(vertices) - 1:
        # Not fatal, but it means the lambda intervals cannot all be named, so
        # say so rather than silently emitting fewer replay lambdas.
        print(f"warning: {len(vertices)} vertices against {len(breaks)} "
              f"breakpoints in {report}")
    return vertices, breaks


def partial_from(stream, n_aa):
    """What a stopped walk's `--stream` file holds: its vertices, the tie
    weights between consecutive ones, and its solved gaps.

    Every streamed vertex is an exact optimum at some weight, so it is a hull
    vertex; consecutive ones need not be adjacent yet, so the weights between
    them are breakpoints of the hull found so far, not of the final one. A
    line cut short by the kill is skipped.
    """
    seen, vertices, gaps = set(), [], []
    path = Path(stream)
    lines = path.read_text().splitlines() if path.exists() else []
    for ln in lines:
        try:
            ev = json.loads(ln)
        except json.JSONDecodeError:
            continue
        if ev.get("event") == "gap":
            gaps.append(ev)
        elif ev.get("event") == "vertex":
            key = round(ev["mfe_kcal_mol"] * 100)
            if key not in seen:
                seen.add(key)
                log_w = ev["adaptation_log"]
                vertices.append(dict(mfe=ev["mfe_kcal_mol"], sum_ln_w=log_w,
                                     cai=math.exp(log_w / n_aa),
                                     found_s=ev.get("elapsed_s")))
    vertices.sort(key=lambda v: v["mfe"])
    # The walk's own tie weight (`derna_lib::tie_weight`): energies in
    # dekacalories, adaptation scaled by 100.
    breaks = []
    for a, b in zip(vertices, vertices[1:]):
        de = round(b["mfe"] * 100) - round(a["mfe"] * 100)
        dl = 100.0 * (b["sum_ln_w"] - a["sum_ln_w"])
        if de > 0 and dl > 0:
            breaks.append(dl / (de + dl))
    return vertices, sorted(breaks, reverse=True), gaps


def verify(cfg, args, log):
    """Each comparator at interior lambdas against derna2's own answer.

    The corners cannot establish a map: a transform wrong in the middle can
    agree at both ends. Objectives are compared, not sequences, since ties can
    be spelled differently.
    """
    out_dir = Path(cfg["out_dir"]) / "verify"
    out_dir.mkdir(parents=True, exist_ok=True)
    rescore = e1.Rescorer(cfg)
    label, prot = min(e1.instances(cfg), key=lambda i: len(i[1]))
    lams = [float(x) for x in (args.lambdas or "0.9,0.7,0.5,0.3,0.1").split(",")]
    log(f"verify on {label} ({len(prot)} aa) at lambdas {lams}")
    ref = {}
    bad = []
    for name, arm in cfg["arms"].items():
        if not arm.get("enabled", True) or arm.get("lambdas") == "own":
            continue
        for lam in lams:
            stem = out_dir / f"{name}.{lam}"
            row = run_lambda(cfg, arm, label, prot, lam, [cfg["cores"][0]],
                             out_dir, stem, budget_s=cfg.get("verify_s", 1800))
            if row["exit"] != "ok" or not row.get("seq"):
                log(f"  {name:20} lambda {lam:<7} {row['exit']}")
                bad.append((name, lam, row["exit"]))
                continue
            sc = rescore(row["seq"], prot)
            obj = (lam * sc["mfe_rescored"] * 100
                   - (1 - lam) * 100 * math.log(sc["cai_rescored"])
                   * len(prot))
            if arm["kind"] == "derna":
                ref[lam] = obj
            delta = obj - ref.get(lam, obj)
            flag = "" if abs(delta) < 1.0 else "  <-- MISMATCH"
            log(f"  {name:20} lambda {lam:<7} mfe {sc['mfe_rescored']:9.2f} "
                f"cai {sc['cai_rescored']:.4f}  objective {obj:12.1f} "
                f"delta {delta:+9.1f}{flag}")
            if flag:
                bad.append((name, lam, f"objective off by {delta:+.1f}"))
    if bad:
        log("")
        log("FAILURES — these arms are not optimizing the same problem:")
        for name, lam, why in bad:
            log(f"  {name} at lambda {lam}: {why}")
        return 1
    log("")
    log("every arm matches derna2's objective at every lambda tried")
    return 0


def plan(cfg, args, log):
    """The cells the grid owes, against the cells it has.

    A short output file is indistinguishable from a small grid, so the count is
    stated before the run and checked after it rather than inferred from what
    arrived. Every cell is named, so a missing one can be named too.
    """
    stamp = provenance(cfg)
    rows = e1.Results(Path(cfg["out_dir"]) / "sweep.jsonl").rows()
    have = {(r["arm"], r["protein"]): r for r in rows
            if r.get("kind") == "cell" and r.get("pin") == stamp["pin"]}
    other = [r for r in rows if r.get("kind") == "cell"
             and r.get("pin") != stamp["pin"]]
    walks = [name for name, _ in walk_arms(cfg)]
    arms = walks + [a for a in cfg["arms"]
                    if cfg["arms"][a].get("enabled", True)]
    want, missing, hull_missing = [], [], []
    for label, prot in e1.instances(cfg):
        if args.max_n and len(prot) > args.max_n:
            continue
        hull = (Path(cfg["out_dir"]) / "hull" / f"{label}.vertices.json")
        for arm in arms:
            want.append((arm, label))
            if (arm, label) not in have:
                missing.append((arm, label, len(prot)))
                # A replay arm with no hull cannot even be planned: its lambda
                # list comes from the hull run, so it is blocked, not pending.
                if (arm not in walks
                        and cfg["arms"].get(arm, {}).get("lambdas", "hull")
                        == "hull" and not has_hull(cfg, label)):
                    hull_missing.append((arm, label))
    present = set(have) & set(want)
    log(f"{len(want)} cells planned: {len(arms)} arms over "
        f"{len(want) // len(arms)} instances")
    log(f"pin {stamp['pin']}, walk {stamp['sha_walk_front']}")
    log(f"{len(present)} present, {len(missing)} missing, "
        f"{len(hull_missing)} of those blocked on a hull that has not run")
    if other:
        pins = sorted({r.get("pin") for r in other})
        log(f"{len(other)} cells present from another pin ({pins}) — they "
            f"measure a different derna2 and count as missing here")
    stale = sorted(set(have) - set(want))
    if stale:
        log(f"{len(stale)} cells present that this config does not plan "
            f"(a config change, or another grid writing here): "
            f"{stale[:5]}{' ...' if len(stale) > 5 else ''}")
    by_exit = {}
    for key in present:
        by_exit.setdefault(have[key]["exit"], []).append(key)
    for verdict, keys in sorted(by_exit.items()):
        log(f"  {verdict}: {len(keys)}")
    if args.arms:
        want_arms = set(args.arms.split(","))
        for arm, label, n in missing:
            if arm in want_arms:
                log(f"  missing {arm} {label} ({n} aa)")
    return 0


def summarize(cfg, args, log):
    """Fronts, recovery against the exact hull, and the anytime curves.

    Reads one pin's rows only. A file accumulated across re-pins holds rows
    from several builds, and pooling them silently doubles every count and
    merges fronts from binaries that are not the same program -- so the filter
    is applied here rather than left to whoever reads the output. `--pin`
    selects another one, or `all` to pool deliberately.
    """
    rows = e1.Results(Path(cfg["out_dir"]) / "sweep.jsonl").rows()
    want = args.pin or cfg.get("derna2_pin")
    seen = sorted({r.get("pin") for r in rows}, key=lambda p: (p is None, p))
    if want != "all":
        kept = [r for r in rows if r.get("pin") == want]
        if len(kept) != len(rows):
            log(f"pin {want}: {len(rows) - len(kept)} of {len(rows)} rows are "
                f"from another build and are excluded ({seen})")
        rows = kept
    else:
        log(f"POOLING {len(seen)} pins deliberately ({seen}): counts below "
            f"are sums over builds, not one experiment")
    cells = [r for r in rows if r.get("kind") == "cell"]
    points = [r for r in rows if r.get("kind") != "cell"]
    log(f"{len(cells)} cells, {len(points)} lambda runs, pin {want}")
    late = deferred_scores(cfg)
    waiting = 0
    for r in points:
        if r.get("rescore_deferred"):
            if point_key(r) in late:
                r.update(late[point_key(r)])
            else:
                waiting += 1
    if waiting:
        log(f"{waiting} walk points not yet re-scored by `rescore`; they are "
            f"left out of the fronts below")
    if not rows:
        log("nothing to summarize at this pin")
        return 0
    log("")
    blocked = [c for c in cells if c["exit"] == "no_hull"]
    if blocked:
        log(f"{len(blocked)} cells blocked with no hull — not short fronts, "
            f"no measurement at all: "
            f"{sorted({c['protein'] for c in blocked})}")
    log(f"{'arm':22} {'protein':14} {'n':>5} {'lambdas':>12} {'exit':>9} "
        f"{'wall':>9}")
    for c in sorted(cells, key=lambda c: (c["n_aa"], c["arm"])):
        planned = c["lambdas_planned"]
        got = f"{c['lambdas_done']}/{planned}" if planned else \
            str(c["lambdas_done"])
        log(f"{c['arm']:22} {c['protein']:14} {c['n_aa']:5d} {got:>12} "
            f"{c['exit']:>9} {c['wall_s'] / 60:8.1f}m")
    by_arm = {}
    for r in points:
        if r.get("mfe_rescored") is None or r.get("cai_rescored") is None:
            continue
        by_arm.setdefault((r["arm"], r["protein"]), []).append(
            (r["mfe_rescored"], r["cai_rescored"]))
    log("")
    log("Points and hull-vertex recovery per cell, on re-scored values")
    log(f"  {'arm':22} {'protein':14} {'runs':>5} {'pareto':>7} "
        f"{'vertices':>10}")
    hulls = {}
    for key in sorted(by_arm):
        arm, protein = key
        pts = pareto(by_arm[key])
        if protein not in hulls:
            hulls[protein] = load_hull(cfg, protein)
        verts = hulls[protein]
        if verts:
            hit = sum(1 for v in verts if any(
                abs(v["mfe"] - m) <= 0.05 and abs(v["cai"] - c) <= 1e-3
                for m, c in pts))
            rec = f"{hit}/{len(verts)}"
        else:
            rec = "no hull"
        log(f"  {arm:22} {protein:14} {len(by_arm[key]):5d} {len(pts):7d} "
            f"{rec:>10}")


def load_hull(cfg, protein):
    """The exact hull vertices for one instance, if its hull has been run."""
    path = Path(cfg["out_dir"]) / "hull" / f"{protein}.vertices.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())["vertices"]


def pareto(points):
    """Non-dominated (mfe, cai): lower mfe better, higher cai better.

    Walking in increasing mfe, every earlier point has a better energy, so the
    current one survives only by beating all of them on CAI. Earlier points
    are never removed: a lower mfe at a lower CAI is not dominated by a higher
    mfe at a higher CAI, it is the other end of the same trade-off.
    """
    out, best = [], -math.inf
    for mfe, cai in sorted(set(points)):
        if cai > best:
            out.append((mfe, cai))
            best = cai
    return out


def walk_timing(cfg):
    """Whether a walk pass is under way, read from the grid's log: its last
    `=== hull` line has no matching `done`."""
    path = Path(cfg["out_dir"]) / "grid.log"
    if not path.exists():
        return False
    last = [ln for ln in path.read_text().splitlines() if "=== hull" in ln]
    return bool(last) and " done in " not in last[-1]


def rescore_walks(cfg, args, log):
    """Re-score every deferred walk point into `walk_rescore.jsonl`.

    Runs beside the grid, never beside a walk: it waits while a walk pass is
    under way, so no walk is timed next to its folds. Each score is synced as
    it is written, and a re-run skips the points already scored. With
    `--follow` it keeps going until the grid has finished and nothing is left.
    """
    rescore = e1.Rescorer(cfg)
    out = e1.Results(walk_rescore_path(cfg))
    prots = dict(e1.instances(cfg))
    chunk = cfg.get("rescore_chunk", 200)
    while True:
        done = set(deferred_scores(cfg))
        todo = [r for r in e1.Results(Path(cfg["out_dir"]) / "sweep.jsonl").rows()
                if r.get("rescore_deferred") and r.get("seq")
                and point_key(r) not in done]
        if todo:
            log(f"{len(todo)} walk points to re-score")
        for i in range(0, len(todo), chunk):
            while walk_timing(cfg):
                time.sleep(60)
            part = todo[i:i + chunk]
            by_prot = {}
            for r in part:
                by_prot.setdefault(r["protein"], []).append(r)
            for protein, rs in by_prot.items():
                scores = rescore.many([r["seq"] for r in rs], prots[protein])
                for r in rs:
                    out.append({"pin": r.get("pin"), "arm": r["arm"],
                                "protein": protein,
                                "point_index": r["point_index"], **next(scores)})
        grid = Path(cfg["out_dir"]) / "grid.log"
        finished = grid.exists() and "grid finished" in grid.read_text()
        if not args.follow or (finished and not todo):
            log("walk re-scoring caught up")
            return 0
        time.sleep(300)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=["hull", "verify", "plan",
                                        "sweep", "summarize", "lambdas",
                                        "rescore"])
    ap.add_argument("--follow", action="store_true",
                    help="rescore: keep going until the grid has finished")
    ap.add_argument("--config", required=True)
    ap.add_argument("--arms")
    ap.add_argument("--lambdas")
    ap.add_argument("--max-n", type=int)
    ap.add_argument("--shortest-first", action="store_true",
                    help="sweep, hull: run the shortest proteins' cells first")
    ap.add_argument("--rewalk", action="store_true",
                    help="hull: a walk counts as run only from the walk binary "
                         "on disk; a re-walk of a present hull writes beside it")
    ap.add_argument("--first",
                    help="hull: comma-separated instances to walk before the rest")
    ap.add_argument("--arm-order",
                    help="sweep: run whole arms in this comma-separated order")
    ap.add_argument("--walked-only", action="store_true",
                    help="sweep: leave out replay cells whose walk has not run, "
                         "rather than record them as blocked")
    ap.add_argument("--skip-after-timeout", action="store_true",
                    help="sweep: once an arm times out, leave out its cells on "
                         "proteins at least that long whose limit is no larger "
                         "(no row is written)")
    ap.add_argument("--phase", choices=["single", "par"],
                    help="run only the single-core arms, or only the "
                         "multi-core ones")
    ap.add_argument("--pin", help="summarize this pin's rows; 'all' pools "
                                  "every build, which is not one experiment")
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    # `drive.py` derives the repository root from its own location, which is
    # wrong whenever this driver runs from a copy outside the tree (a pinned
    # worktree, a scratch checkout). The config may say where it is.
    if cfg.get("repo"):
        e1.REPO = Path(cfg["repo"])
    cfg = e1.resolve_paths(cfg)
    out_dir = Path(cfg["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    logfile = open(out_dir / f"{args.command}.log", "a")

    def log(msg):
        line = f"{time.strftime('%H:%M:%S')} {msg}"
        print(line, flush=True)
        print(line, file=logfile, flush=True)

    if args.command == "hull":
        return run_hull(cfg, args, log)
    if args.command == "verify":
        return verify(cfg, args, log)
    if args.command == "plan":
        return plan(cfg, args, log)
    if args.command == "sweep":
        return run_sweep(cfg, args, log)
    if args.command == "rescore":
        return rescore_walks(cfg, args, log)
    if args.command == "lambdas":
        for label, prot in e1.instances(cfg):
            if args.max_n and len(prot) > args.max_n:
                continue
            lams = hull_lambdas(cfg, label)
            log(f"{label} ({len(prot)} aa): {len(lams)} lambdas, "
                f"{lams[0]:.6g} .. {min(lams):.6g}")
        return 0
    return summarize(cfg, args, log)


if __name__ == "__main__":
    sys.exit(main() or 0)
