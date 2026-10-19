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
    stamp = {"pin": cfg.get("derna2_pin")}
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
    wrapped += (["--stdin", stdin] if stdin else []) + ["--"] + argv
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
    limit_s = arm.get("limit_s", cfg["limit_s"])
    cap_kb = cfg["mem_cap_gb"] * 1e6
    started = time.time()
    done = 0
    verdict = "ok"
    for i, lam in enumerate(lambdas):
        spent = time.time() - started
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
                   elapsed_at_finish_s=time.time() - started,
                   load=os.getloadavg()[0])
        if row.get("seq"):
            row.update(rescore(row["seq"], prot))
        results.append(row)
        done += 1
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
            f"{(time.time() - started) / 60:.1f} min spent")
    summary = dict(kind="cell", arm=arm_name, protein=label, n_aa=len(prot),
                   lambdas_planned=len(lambdas), lambdas_done=done,
                   exit=verdict, wall_s=time.time() - started,
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
            cells.append((name, arm, label, prot))
    # Longest first, so the run that starts last is the short one.
    cells.sort(key=lambda c: -len(c[3]))
    log(f"{len(cells)} cells to run")

    def one(spec):
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
    started = time.time()
    cores = pool.take(1)
    try:
        wrapped = ["taskset", "-c", ",".join(map(str, cores)),
                   sys.executable, str(MEASURE), str(meas),
                   "--timeout", str(arm.get("limit_s", cfg["limit_s"])),
                   "--stdout", str(e1.ext(stem, ".stdout")),
                   "--stderr", str(e1.ext(stem, ".stderr")), "--"] + argv
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
        if r.get("seq"):
            r.update(rescore(r["seq"], prot))
        results.append(r)
    summary = dict(kind="cell", arm=name, protein=label, n_aa=len(prot),
                   lambdas_planned=None, lambdas_done=len(rows),
                   exit=m["exit"], wall_s=time.time() - started, ncores=1,
                   maxrss_kb=m["maxrss_kb"], instructions=m["instructions"],
                   load=os.getloadavg()[0])
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
    argv = [cfg["bins"]["derna1rs_sweep"], "--fasta", str(fasta),
            "--codon-usage", cfg["tables"]["derna"],
            "--out", str(e1.ext(stem, ".front.jsonl")),
            "--stream", str(stream)]
    meas = e1.ext(stem, ".measure.json")
    started = time.time()
    cores = pool.take(1)
    try:
        subprocess.run(
            ["taskset", "-c", ",".join(map(str, cores)),
             sys.executable, str(MEASURE), str(meas),
             "--timeout", str(arm.get("limit_s", cfg["limit_s"])),
             "--stdout", str(e1.ext(stem, ".stdout")),
             "--stderr", str(e1.ext(stem, ".stderr")), "--"] + argv,
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
        if r.get("seq"):
            r.update(rescore(r["seq"], prot))
        results.append(r)
    summary = dict(kind="cell", arm=name, protein=label, n_aa=len(prot),
                   lambdas_planned=None, lambdas_done=len(rows),
                   exit=m["exit"], wall_s=time.time() - started, ncores=1,
                   user_s=m.get("user_s"), sys_s=m.get("sys_s"),
                   maxrss_kb=m["maxrss_kb"], instructions=m["instructions"],
                   load=os.getloadavg()[0])
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
    done = {(r["arm"], r["protein"]) for r in results.rows()
            if r.get("kind") == "cell" and r.get("pin") == stamp["pin"]}
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
            if i == 0 and (out_dir / f"{label}.vertices.json").exists():
                log(f"{label}: hull already present")
                continue
            if i > 0 and (arm_name, label) in done:
                log(f"{arm_name}.{label}: already run")
                continue
            walk_cell(cfg, arm_name, threads, label, prot, fasta, out_dir,
                      work, pool, results, rescore, log, reference=i == 0)


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
    argv = [cfg["bins"]["walk_front"], "--fasta", str(fasta),
            "--n", str(len(prot)),
            "--lo", repr(cfg["epsilon"]), "--hi", repr(1 - cfg["epsilon"]),
            "--threads", str(threads)] + list(cfg.get("walk_flags", [])) + [
            "--out", str(report), "--front", str(front),
            "--stream", str(stream)]
    stem = work / f"{arm_name}.{label}"
    meas = e1.ext(stem, ".measure.json")
    started = time.time()
    cores = pool.take(threads)
    try:
        subprocess.run(
            ["taskset", "-c", ",".join(map(str, cores)),
             sys.executable, str(MEASURE), str(meas),
             "--timeout", str(cfg["limit_s"]),
             "--stdout", str(e1.ext(stem, ".stdout")),
             "--stderr", str(e1.ext(stem, ".stderr")), "--"] + argv,
            cwd=str(e1.REPO))
    finally:
        pool.give(cores)
    m = json.loads(meas.read_text())
    wall_min = (time.time() - started) / 60
    log(f"{arm_name}.{label}: walk {m['exit']} in {wall_min:.1f} min")
    if m["exit"] != "ok":
        # No report means no certified hull: record what the stream holds
        # as a partial cell, write its vertices beside (not as) the
        # `.vertices.json`, and let the replay arms refuse the instance
        # unless the config opts into the partial hull's lambdas.
        vertices, breaks, gaps = partial_from(stream, len(prot))
        npoints = 0
        for gap in gaps:
            for pt in gap.get("front", []):
                row = dict(kind="front_point", arm=arm_name, protein=label,
                           n_aa=len(prot), point_index=npoints,
                           ncores=len(cores), partial=True,
                           found_s=gap.get("elapsed_s"),
                           mfe=pt.get("mfe_kcal_mol"), cai=pt.get("cai"),
                           seq=pt.get("sequence"))
                if row["seq"]:
                    row.update(rescore(row["seq"], prot))
                results.append(row)
                npoints += 1
        results.append(dict(kind="cell", arm=arm_name, protein=label,
                            n_aa=len(prot), lambdas_planned=0,
                            lambdas_done=npoints, exit=m["exit"],
                            wall_s=m["wall_s"], user_s=m.get("user_s"),
                            sys_s=m.get("sys_s"),
                            maxrss_kb=m.get("maxrss_kb"),
                            ncores=len(cores), partial=True,
                            vertices=len(vertices), gaps=len(gaps),
                            front_points=npoints,
                            load=os.getloadavg()[0]))
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
    for i, pt in enumerate(points):
        row = dict(kind="front_point", arm=arm_name, protein=label,
                   n_aa=len(prot), point_index=i, ncores=len(cores),
                   mfe=pt.get("mfe_kcal_mol"), cai=pt.get("cai"),
                   seq=pt.get("sequence"))
        if row["seq"]:
            row.update(rescore(row["seq"], prot))
        results.append(row)
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
                        load=os.getloadavg()[0]))
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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=["hull", "verify", "plan",
                                        "sweep", "summarize", "lambdas"])
    ap.add_argument("--config", required=True)
    ap.add_argument("--arms")
    ap.add_argument("--lambdas")
    ap.add_argument("--max-n", type=int)
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
