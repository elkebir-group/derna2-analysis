#!/usr/bin/env python3
"""The UTR-aware grid: one cell is one protein behind one UTR pair, in four arms.

The question is what it costs to take the UTRs seriously. A cell runs the same
molecule four ways and reports, for each, the share of its points whose own
free fold respects the constraint:

  A1  unconstrained   the UTRs inside the DP, no pair forbidden
  A2  constrained     the UTRs inside the DP, every forbidden pair rejected
  A3  mRNAfold        the comparator: UTRs in its energy model, nothing
                      rejecting a CDS-UTR pair, run at A2's hull weights
  A4  optimal         the structurally optimal front, partial under a budget

A1 and A3 are not fronts of structurally optimal designs and are not claimed
to be. What they contribute is that share.

Commands:
  plan       the cells the grid owes, against the cells it has
  run        the grid: every (instance, pair) cell, arms in dependency order
  score      the re-scoring of timed output, off every wall
  summarize  per-pair shares, the gap distribution, the A2 -> A4 price

Timing and scoring are separate commands on purpose. Every fold that is not
the measured search runs under `score`, on other cores, after the measurement
has ended, so no arm's wall carries the cost of reading it.
"""
import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
MEASURE = HERE / "measure.py"

# `systemd-run --user` needs a session bus, which a host reached over ssh
# without one does not have. There every capped run would fail in under a
# second, so `mem_cap_gb` refuses to start rather than emptying a queue into
# crash rows.
HAVE_SYSTEMD_RUN = subprocess.run(
    ["systemd-run", "--user", "--scope", "-q", "true"],
    capture_output=True, text=True).returncode == 0

# The binaries whose build decides what a row means. `derna` must be at or
# past the merge that made `--front-from` exist and the default constraint the
# two classes that reach the coding sequence (ADR 0044); the comparator and the
# folders are frozen once.
PINNED_BINS = ("derna", "mfe_fold", "mrnafold")

THREE_LETTER = dict(A="Ala", R="Arg", N="Asn", D="Asp", C="Cys", Q="Gln", E="Glu", G="Gly", H="His",
                    I="Ile", L="Leu", K="Lys", M="Met", F="Phe", P="Pro", S="Ser", T="Thr", W="Trp",
                    Y="Tyr", V="Val")

ARMS = ("A1", "A2", "A3", "A4")


# ---- configuration and provenance --------------------------------------------------

def load_cfg(path):
    cfg = json.loads(Path(path).read_text())
    cfg.setdefault("forbid", "5-cds,cds-3")
    cfg.setdefault("stop", "UAA")
    cfg.setdefault("epsilon", 2.0 ** -23)
    cfg.setdefault("threads", 1)
    cfg.setdefault("score_threads", 8)
    cfg.setdefault("budget_s", 43200)
    # `budget_s` is A4's `--max-secs` on the structurally optimal search;
    # `walk_budget_s` is the walk's `--max-secs` on A1 and A2. Each run then
    # returns the part of its front it found, and the outer kill sits `emit_s`
    # past its budget so that writing the front is never what gets killed.
    # `a3_budget_s` bounds a cell's whole run of comparator weights.
    cfg.setdefault("walk_budget_s", cfg["budget_s"])
    cfg.setdefault("a3_budget_s", cfg["budget_s"])
    cfg.setdefault("emit_s", 3600)
    cfg.setdefault("mrnafold_budget_s", cfg["budget_s"])
    cfg.setdefault("a3_max_lambdas", 24)
    cfg.setdefault("a4_threads", cfg["threads"])
    # `arm_builds` lets one arm run on another build than the grid's:
    # {"A4": {"derna2_pin": ..., "bins": {"derna": ...}}}. Its rows carry that
    # build's stamp, so a table can hold an arm whose path a later build changed
    # beside arms whose path it did not, each row saying which build ran it.
    cfg.setdefault("arm_builds", {})
    return cfg


def arm_cfg(cfg, arm):
    """The config an arm runs under: the grid's, with the arm's build over it."""
    build = cfg["arm_builds"].get(arm)
    if not build:
        return cfg
    return {**cfg, "derna2_pin": build.get("derna2_pin", cfg.get("derna2_pin")),
            "bins": {**cfg["bins"], **build.get("bins", {})}}


def provenance(cfg):
    """The pin every row is stamped with.

    The constraint goes on the stamp beside the build, because `derna`'s
    default moved from three crossing classes to two on 2026-10-19 (ADR 0044):
    afterwards a row that records only its build cannot say what constraint it
    ran under, and one that records only its constraint cannot say whether the
    build honoured it. Both, or neither is readable.
    """
    stamp = {"pin": cfg.get("derna2_pin"), "forbid": cfg["forbid"], "stop": cfg["stop"]}
    for name in PINNED_BINS:
        digest = hashlib.sha256(Path(cfg["bins"][name]).read_bytes()).hexdigest()
        stamp[f"sha_{name}"] = digest[:12]
    return stamp


class Results:
    """Append-only JSONL, synced per row, keyed for resume.

    A row is a finished run. A host that dies must not take back what a run
    already paid for, which at twelve hours a cell is the whole grid.
    """

    def __init__(self, path, stamp, arm_stamps=None):
        self.path, self.stamp = Path(path), stamp
        self.arm_stamps = arm_stamps or {}

    def stamp_for(self, arm):
        """The stamp a row of `arm` carries: its own build's, or the grid's."""
        return self.arm_stamps.get(arm, self.stamp)

    def rows(self):
        if not self.path.exists():
            return []
        return [json.loads(line) for line in self.path.open() if line.strip()]

    def append(self, row):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a") as f:
            f.write(json.dumps({**row, **self.stamp_for(row.get("arm"))}) + "\n")
            f.flush()
            os.fsync(f.fileno())


def done_keys(results):
    """The (cell, arm, index) keys already recorded under this pin.

    A row from another pin than its arm's is not a row for this grid: it counts
    as missing, so a re-pin re-runs rather than mixing two builds in one arm.

    An A1 row without `flanks_designed` designed the coding sequence alone,
    which this grid no longer runs, so it counts as missing too.
    """
    return {(r["cell"], r["arm"], r.get("index")) for r in results.rows()
            if all(r.get(k) == v for k, v in results.stamp_for(r["arm"]).items())
            and (r["arm"] != "A1" or r.get("flanks_designed"))}


# ---- inputs ------------------------------------------------------------------------

def read_fasta(path):
    seq = "".join(x.strip() for x in Path(path).read_text().splitlines() if not x.startswith(">"))
    return seq.rstrip("*")


def read_utrs(path):
    """{record name: RNA sequence} from the pinned `utrs.fasta`.

    The pin is written as DNA because that is how its source documents print
    it. Everything downstream here is RNA: `derna` would read the `T` as `U`
    itself, but mRNAfold's config wants `U`, and one
    spelling through the whole driver means no arm can differ from another by a
    conversion that happened in only one place.
    """
    out, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
        elif line.strip() and name:
            out[name] += line.strip().upper().replace("T", "U")
    return out


def read_codon_table(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("amino_acid") or not line.strip():
            continue
        aa, codon, frac = line.strip().split(",")
        rows.append((aa, codon.replace("T", "U"), float(frac)))
    return rows


def codon_weights(rows):
    """Relative adaptiveness: each codon's fraction over its amino acid's largest."""
    best = {}
    for aa, _, f in rows:
        best[aa] = max(best.get(aa, 0.0), f)
    return {codon: f / best[aa] for aa, codon, f in rows}


def cells(cfg):
    """Every (instance, UTR pair) cell, longest instance first.

    Longest-first inside the grid, because the cells that do not finish are the
    ones the experiment is about: their partiality is the measurement, and it
    should arrive before the budget runs out rather than after.

    The grid is not a cross product. An instance may name the pairs it runs
    against, and each cell carries the panel it belongs to, because two panels
    answering different questions are reported separately and never pooled: one
    sweeps length at a single pair, the other sweeps the pair on a few
    instances. The panel is a property of the cell rather than of the instance,
    which is what keeps the two from overlapping -- an instance in both panels
    runs its cell at the length grid's pair once, and that one cell serves as
    the pair panel's column for that pair too. An instance that names no pairs
    runs against all of them.
    """
    utrs = read_utrs(cfg["utrs_fasta"])
    by_key = {pair["key"]: pair for pair in cfg["pairs"]}
    axis = cfg.get("length_grid_pair")
    out = []
    for inst in cfg["instances"]:
        prot = read_fasta(inst["fasta"])
        for key in inst.get("pairs") or list(by_key):
            if key not in by_key:
                raise SystemExit(f"{inst['id']}: no UTR pair {key} in this config")
            pair = by_key[key]
            out.append({"cell": f"{inst['id']}__{key}", "instance": inst["id"],
                        "pair": key, "protein": prot, "n_aa": len(prot),
                        "panel": "length" if key == axis else "pairs",
                        "band": inst.get("band", "-"),
                        "utr5": utrs[pair["utr5"]], "utr3": utrs[pair["utr3"]],
                        "utr5_name": pair["utr5"], "utr3_name": pair["utr3"]})
    out.sort(key=lambda c: (-c["n_aa"], c["cell"]))
    return out


def cell_dir(cfg, cell):
    path = Path(cfg["out_dir"]) / "cells" / cell["cell"]
    path.mkdir(parents=True, exist_ok=True)
    return path


def residues(cfg, cell):
    """The length of the coding region in residues, the stop codon included.

    The constraint counts a fixed stop codon as coding sequence, so a window over the
    coding region has to cover it: a pair from the stop into the 3' UTR is `cds-3` and
    forbidden, and a window three nucleotides short would class it `3-3` and allow it.
    """
    return cell["n_aa"] + (1 if cfg["stop"].lower() != "none" else 0)


def flank_files(cfg, cell):
    """The pair's two sequences as FASTA beside the cell, and their paths.

    `derna` takes a literal too, but a 295-nucleotide literal in an argv that
    `measure.py` records is unreadable afterwards, and the file is what the row
    can name.
    """
    out = []
    for end in ("utr5", "utr3"):
        path = cell_dir(cfg, cell) / f"{end}.fa"
        path.write_text(f">{cell[end + '_name']}\n{cell[end]}\n")
        out.append(str(path))
    return out


# ---- the hull weights A3 replays ---------------------------------------------------

def breakpoints(vertices):
    """The weights at which consecutive hull vertices tie, highest first.

    `derna front --report` lists the hull's vertices as points, not as the
    weights that reach them, so the weights are recovered from the points. The
    objective is `lambda*E - (1 - lambda)*L`, so two points tie at
    `r / (1 + r)` with `r = dL / dE`, and a vertex is optimal over the interval
    between its two neighbouring ties.
    """
    ordered = sorted(vertices, key=lambda v: v["mfe"])
    out = []
    for a, b in zip(ordered, ordered[1:]):
        d_e, d_l = a["mfe"] - b["mfe"], a["adaptation_log"] - b["adaptation_log"]
        if d_e == 0:
            continue
        r = d_l / d_e
        if r > 0:
            out.append(r / (1.0 + r))
    return sorted(out, reverse=True)


def replay_lambdas(breaks, eps):
    """One weight per hull vertex plus both corners, in the order they run.

    A vertex owns an interval, so its midpoint has one answer while a
    breakpoint has two tied ones and the recovered vertex would be a function
    of the comparator's tie-breaking. The order is corners first, then
    bisection over the vertex index, so a cell stopped early has covered the
    whole range coarsely rather than one end of it finely.
    """
    edges = [1.0 - eps] + list(breaks) + [eps]
    mids = [math.sqrt(edges[i] * edges[i + 1]) for i in range(len(edges) - 1)]
    order = [0, len(mids) - 1] if len(mids) > 1 else [0]

    def bisect(lo, hi):
        if hi - lo <= 1:
            return
        mid = (lo + hi) // 2
        order.append(mid)
        bisect(lo, mid)
        bisect(mid, hi)

    bisect(0, len(mids) - 1)
    seen, out = set(), []
    for i in order:
        if i not in seen:
            seen.add(i)
            out.append(mids[i])
    return out


def hull_lambdas(cfg, cell):
    """A2's hull weights for this cell, or `[]` when A2 has left no report.

    A2 stopped by its budget still writes the part of the front it walked, so
    this is the weights that actually completed rather than the weights the
    plan hoped for.
    """
    path = cell_dir(cfg, cell) / "a2.report.json"
    if not path.exists():
        return []
    report = json.loads(path.read_text())
    return replay_lambdas(breakpoints(report.get("vertices", [])), cfg["epsilon"])


# ---- running one measured command --------------------------------------------------

def capped(cfg, argv):
    """`argv` in a memory-capped cgroup of its own, or unchanged.

    A cap the kernel enforces is the only one that holds: a cell that would
    outgrow the host is killed at `mem_cap_gb` instead of taking the host's
    other work with it, and `measure.py` reads that SIGKILL as `oom`. Swap is
    refused outright, since a cell that swaps is no longer measuring the
    algorithm.
    """
    gb = cfg.get("mem_cap_gb")
    if not gb:
        return argv
    if not HAVE_SYSTEMD_RUN:
        raise RuntimeError(
            "mem_cap_gb is set but `systemd-run --user` has no bus on this host: "
            "every run would fail at once. Use rss_cap_gb, which the binary enforces "
            "in process, or clear mem_cap_gb.")
    return ["systemd-run", "--user", "--scope", "-q",
            "-p", f"MemoryMax={gb}G", "-p", "MemorySwapMax=0"] + argv


def pinned(cores, argv):
    return (["taskset", "-c", cores] if cores else []) + argv


def measured(cfg, stem, argv, timeout):
    """Run `argv` once under `measure.py` and return what it cost.

    The counters are opened per command, so each measured command needs its own
    wrapper; the memory scope wraps the tool alone, inside the measurement, so
    the cell's cgroup holds nothing whose death would lose the measurement.
    """
    stem = Path(stem)
    machine = stem.with_suffix(".measure.json")
    wrapper = [sys.executable, str(MEASURE), str(machine), "--timeout", str(timeout),
               "--stdout", str(stem.with_suffix(".stdout")),
               "--stderr", str(stem.with_suffix(".stderr")), "--"]
    subprocess.run(pinned(cfg.get("cores"), wrapper + capped(cfg, argv)), check=False)
    rec = json.loads(machine.read_text())
    rec["stdout_path"] = str(stem.with_suffix(".stdout"))
    return rec


def derna_front_argv(cfg, cell, arm, extra):
    """The shared shape of every `derna front` arm.

    `--forbid` is written out on every arm that has flanks even though it is
    now the default, because the default differs between builds: a command line
    that states the constraint cannot be misread, and the one in the row is the
    one that ran.

    A4 takes its own thread count, `a4_threads`: its search checks designs on
    every thread it is given and runs a few cells for hours, where the walks
    run many cells briefly and gain more from one core each.
    """
    out = cell_dir(cfg, cell) / f"{arm.lower()}.jsonl"
    threads = cfg["a4_threads"] if arm == "A4" else cfg["threads"]
    argv = [arm_cfg(cfg, arm)["bins"]["derna"], "front", "--protein", str(cell["fasta_path"]),
            "--threads", str(threads), "--out", str(out),
            "--report", str(cell_dir(cfg, cell) / f"{arm.lower()}.report.json"), "--force"]
    # Always spell the stop out, including `none`. Omitting the flag used to
    # mean "no stop codon", which is what `none` selects, so the two were
    # interchangeable; once the stop becomes a slot of the codon lattice an
    # omitted flag means "design it" instead, and a run that meant to add no
    # stop would quietly compute a different molecule under the same config.
    # `--stop none` is accepted by every build the grid uses and names the
    # molecule on the command line, where the header can be checked against it.
    argv += ["--stop", cfg["stop"]]
    return argv + extra, out


def walk_budget(cfg):
    """The walk's own limit, for the arms that walk a front.

    A walk stopped by `--max-secs` returns its partial front with the boxes it
    left open, and `--checkpoint` rewrites `--out` after each settled pair of
    hull points, so even the outer kill leaves the front found so far. Without
    both, a walk killed from outside writes no file at all, and a cell whose
    A2 wrote nothing has no A3 or A4 either.
    """
    return ["--max-secs", str(cfg["walk_budget_s"]), "--checkpoint"]


# ---- the arms ----------------------------------------------------------------------

def run_a1(cfg, cell, results):
    """The unconstrained arm: the UTRs inside the DP, no pair forbidden.

    The same molecule as A2 under the same energy model, with `--forbid none`,
    so the two fronts differ only by the constraint. `score` then applies the
    constraint to designs that were not made under it, which is what this arm
    is for. Its row carries `flanks_designed`, because an A1 that designed the
    coding sequence alone ran a different arm and does not count as this one.
    """
    u5, u3 = flank_files(cfg, cell)
    argv, out = derna_front_argv(cfg, cell, "A1",
                                 ["--utr5", u5, "--utr3", u3, "--forbid", "none",
                                  *walk_budget(cfg)])
    stem = cell_dir(cfg, cell) / "a1"
    rec = measured(cfg, stem, argv, cfg["walk_budget_s"] + cfg["emit_s"])
    results.append(row_for(cell, "A1", rec,
                           {"front": str(out), "constrained": False, "flanks_designed": True,
                            **completeness(out)}))
    return rec["exit"]


def run_a2(cfg, cell, results):
    """The constrained arm: the UTRs inside the DP.

    The walk runs under its own `--max-secs` and `--checkpoint`, so a walk
    that does not finish returns a partial front, marked incomplete, rather
    than nothing. A3 and A4 read whatever front this arm wrote; only a cell
    where it wrote none records them as missing.
    """
    u5, u3 = flank_files(cfg, cell)
    argv, out = derna_front_argv(cfg, cell, "A2",
                                 ["--utr5", u5, "--utr3", u3, "--forbid", cfg["forbid"],
                                  *walk_budget(cfg)])
    stem = cell_dir(cfg, cell) / "a2"
    rec = measured(cfg, stem, argv, cfg["walk_budget_s"] + cfg["emit_s"])
    results.append(row_for(cell, "A2", rec,
                           {"front": str(out), "constrained": True, **completeness(out)}))
    return rec["exit"]


def run_a4(cfg, cell, results):
    """The structurally optimal arm, starting from A2's front.

    `--front-from` reads A2's own output in place of walking it again, and is
    refused unless the file's header names the same molecule down to protein,
    flanks, stop codon, constraint and codon table. So A2's front is the only
    thing that can feed this arm, and a mismatched pairing is a refusal rather
    than a different question silently answered.

    The outer limit is the budget plus an emission allowance, because
    `--max-secs` bounds the search and not the writing of what the search
    found: a run killed between those two stages would lose a front it had
    already computed.

    `--checkpoint` appends each point checked, design found and band slice
    finished to `a4.jsonl.search.jsonl` as it settles, so a run in progress
    shows what it has found, and a killed one resumes with `--resume-search`
    rather than starting its hours over.
    """
    a2 = cell_dir(cfg, cell) / "a2.jsonl"
    if not a2.exists():
        return "missing-a2"
    u5, u3 = flank_files(cfg, cell)
    extra = ["--utr5", u5, "--utr3", u3, "--forbid", cfg["forbid"],
             "--structurally-optimal", "--front-from", str(a2),
             "--max-secs", str(cfg["budget_s"]), "--checkpoint"]
    if cfg.get("rss_cap_gb"):
        extra += ["--max-rss-gb", str(cfg["rss_cap_gb"])]
    argv, out = derna_front_argv(cfg, cell, "A4", extra)
    stem = cell_dir(cfg, cell) / "a4"
    rec = measured(cfg, stem, argv, cfg["budget_s"] + cfg["emit_s"])
    results.append(row_for(cell, "A4", rec,
                           {"front": str(out), "constrained": True, "front_from": str(a2),
                            **completeness(out), **prefix_check(out)}))
    return rec["exit"]


def row_for(cell, arm, rec, extra):
    keep = ("wall_s", "user_s", "sys_s", "maxrss_kb", "instructions", "cycles", "exit", "cpus")
    return {"cell": cell["cell"], "instance": cell["instance"], "pair": cell["pair"],
            "panel": cell["panel"], "band": cell["band"],
            "arm": arm, "index": None, "n_aa": cell["n_aa"],
            "utr5_len": len(cell["utr5"]), "utr3_len": len(cell["utr3"]),
            **{k: rec.get(k) for k in keep}, **extra}


def completeness(path):
    """A result's own verdict on itself, which is not its exit category.

    These are different facts and a cell needs both. A structurally optimal run
    that reaches `--max-secs` returns the part of the front it found and exits
    **0**, so `exit: ok` beside `complete: false` is the normal partial case,
    not a contradiction. Only the header can say which it was, and only the
    `open` regions say what the partial answer leaves unbounded.
    """
    for line in Path(path).open() if Path(path).exists() else []:
        if not line.strip():
            continue
        rec = json.loads(line)
        if rec.get("record") != "header":
            break
        return {"complete": rec.get("complete"), "stopped_by": rec.get("stopped_by"),
                "open_regions": len(rec.get("open") or []), "open": rec.get("open"),
                "header_points": rec.get("points"),
                "constrained_points": rec.get("constrained_points"),
                "peak_rss_kb_self": rec.get("peak_rss_kb"), "wall_s_self": rec.get("wall_s")}
    return {"complete": None, "stopped_by": None, "open_regions": None}


def prefix_check(path):
    """Whether a result file holds every point its own header counts.

    A structurally optimal run writes its points as it emits them, so a file
    an outer limit killed is a prefix of the answer rather than a corrupt file.
    The header's `points` is how many the front holds, so a file with fewer
    point lines than that is short by exactly what was never written. Checked
    on every file and not only on the ones believed to have been killed: a
    complete file passes the same check, and a count that disagrees unnoticed
    is the one way this arm could report a truncated front as a whole one.
    """
    if not Path(path).exists():
        return {"points_written": 0, "points_header": None, "whole_file": False}
    header, written = None, 0
    for line in Path(path).open():
        if not line.strip():
            continue
        rec = json.loads(line)
        if rec.get("record") == "header":
            header = rec.get("points")
        else:
            written += 1
    return {"points_written": written, "points_header": header,
            "whole_file": header is not None and written == header}


# ---- the comparator ----------------------------------------------------------------

def mrnafold_table(cfg, work):
    """The canonical table in mRNAfold's `Ala GCT 0.26` format, written once.

    Same fractions as the table every other arm is scored on, so the same
    weights: the comparator's CAI is not a different quantity measured the same
    way.
    """
    path = Path(work) / "mrnafold_codon_table.txt"
    if not path.exists():
        rows = read_codon_table(cfg["tables"]["canonical"])
        lines = [f"{THREE_LETTER[aa]} {codon.replace('U', 'T')} {f!r}"
                 for aa, codon, f in rows if aa in THREE_LETTER]
        path.write_text("\n".join(lines) + "\n")
    return path


def mrnafold_config(cfg, cell, lam_mf, path):
    """One mRNAfold run's config file.

    The stop codon goes at the head of `utr_3p`, because mRNAfold has no flag
    for one and the molecule has to be the same molecule the other arms
    design. For a tool with no region classes that placement is exactly
    equivalent -- it folds what it is given -- while for `derna` it would not
    be, which is why the other arms spell it with `--stop` instead.
    """
    stop = "" if cfg["stop"].lower() == "none" else cfg["stop"]
    body = [f"aa_seq {cell['protein']}", f"lambda {lam_mf!r}", "num_subopt_traces 1",
            f"codon_freq_path {mrnafold_table(cfg, cell_dir(cfg, cell))}",
            f"utr_5p {cell['utr5']}", f"utr_3p {stop}{cell['utr3']}"]
    Path(path).write_text("\n".join(body) + "\n")
    return path


def one_match(pattern, text, what, where):
    """The first group of `pattern` in `text`, or a refusal naming what is missing.

    A comparator that printed something other than what is expected should stop
    the cell with the reason, not hand back a row built from whatever did parse.
    """
    found = re.search(pattern, text, re.M)
    if not found:
        raise SystemExit(f"{where}: no {what} in the comparator's output")
    return found.group(1)


def parse_mrnafold(stdout, lam_mf, weights, cell, stop):
    """One run's design, with its energy recovered from the reported score.

    mRNAfold prints the molecule, its structure, `Score`, `EFE` and `CAI`.
    `Score` is the objective -- the structure's energy plus `lambda' * sum of
    -ln w` -- so the energy comes back by adding the weight term in again.
    `EFE` is an ensemble free energy from a partition function run after the
    search: it is not the MFE, it sits below it, and nothing here reads it.

    The flanks are checked rather than assumed: the coding sequence is taken by
    offset, so a tool that returned the molecule with its UTRs altered, or in a
    different order, would otherwise hand back a frame-shifted coding sequence
    that still parses.
    """
    text = Path(stdout).read_text()
    lines = [x for x in text.splitlines() if x.strip()]
    molecule = one_match(r"^([ACGU]+)$", lines[0], "design", stdout)
    score = float(one_match(r"^Score: ([-0-9.eE+]+)", text, "Score", stdout))
    reported_cai = float(one_match(r"^CAI: ([0-9.]+)", text, "CAI", stdout))
    tail = stop + cell["utr3"]
    if not molecule.startswith(cell["utr5"]) or not molecule.endswith(tail):
        raise SystemExit(f"{stdout}: the emitted molecule's flanks are not the given ones")
    cds = molecule[len(cell["utr5"]):len(molecule) - len(tail)]
    if len(cds) != 3 * cell["n_aa"]:
        raise SystemExit(f"{stdout}: coding sequence is {len(cds)} nt, not {3 * cell['n_aa']}")
    cai_log = sum(math.log(weights[cds[i:i + 3]]) for i in range(0, len(cds), 3))
    return {"sequence": molecule, "structure": lines[1], "score": score,
            "cai_log": cai_log, "cai_reported": reported_cai,
            "cai": math.exp(cai_log / cell["n_aa"]),
            "mfe_kcal_mol": round(score + lam_mf * cai_log, 6)}


def run_a3(cfg, cell, results, already):
    """The comparator at A2's hull weights, one run per weight.

    The weights come from A2's report, so a cell whose A2 stopped early asks
    the comparator only about the part of the front that exists. They run in
    bisection order and are capped, because a long instance's hull has more
    vertices than a comparator can be run at: truncating a coarse-to-fine order
    keeps the whole range at lower resolution, where truncating a sorted one
    would keep one end of it at full resolution and drop the other. A cap of
    `null` runs every hull weight; the order is the same either way, so a cell
    run under a cap continues from where the cap stopped it.
    """
    weights = codon_weights(read_codon_table(cfg["tables"]["canonical"]))
    lambdas = hull_lambdas(cfg, cell)
    if cfg["a3_max_lambdas"] is not None:
        lambdas = lambdas[:cfg["a3_max_lambdas"]]
    if not lambdas:
        return "missing-a2"
    deadline = time.monotonic() + cfg["a3_budget_s"]
    for i, lam in enumerate(lambdas):
        if (cell["cell"], "A3", i) in already:
            continue
        if time.monotonic() >= deadline:
            return "budget"
        results.append(a3_one(cfg, cell, i, lam, weights))
    return "ok"


def a3_one(cfg, cell, index, lam, weights):
    """One (weight, cell) run of the comparator, recorded whatever it did."""
    lam_mf = (1.0 - lam) / lam
    # No dot before the index: `with_suffix` would replace `.v{index}` and
    # every weight of a cell would share one config, stdout and measurement.
    stem = cell_dir(cfg, cell) / f"a3_v{index}"
    conf = mrnafold_config(cfg, cell, lam_mf, stem.with_suffix(".cfg"))
    rec = measured(cfg, stem, [cfg["bins"]["mrnafold"], str(conf)], cfg["mrnafold_budget_s"])
    extra = {"index": index, "lambda": lam, "lambda_mrnafold": lam_mf}
    if rec["exit"] == "ok":
        extra["design"] = parse_mrnafold(rec["stdout_path"], lam_mf, weights, cell,
                                         "" if cfg["stop"].lower() == "none" else cfg["stop"])
    row = row_for(cell, "A3", rec, extra)
    row["index"] = index
    return row


# ---- scoring, off every wall -------------------------------------------------------

def designs_input(cfg, cell, points, name):
    """The sequences of `points`, one per line, for `derna check --designs`.

    The designs are what is scored and nothing else is read from them: the energy and CAI
    each came with are the producing tool's claim about its own molecule, and the whole
    point of this measurement is to compute both again under the constraint.
    """
    path = cell_dir(cfg, cell) / name
    path.write_text("".join(p["sequence"] + "\n" for p in points))
    return path


def reported_mfe(point):
    """The energy a point came with, whichever arm produced it.

    A `derna` point calls it `mfe` and a comparator row carries the energy recovered from
    its score as `mfe_kcal_mol`. What reads this is each arm's own claim about its own
    design, against which a fold computed here is compared.
    """
    return point["mfe"] if "mfe" in point else point["mfe_kcal_mol"]


def read_points(path):
    """The point lines of a result file, header dropped."""
    if not Path(path).exists():
        return []
    out = []
    for line in Path(path).open():
        if not line.strip():
            continue
        rec = json.loads(line)
        if rec.get("record") != "header":
            out.append(rec)
    return out


def score_designs(cfg, cell, designs, out):
    """Each given molecule folded freely and under the constraint.

    This scores the two arms whose points are not `derna` results: the constraint is
    applied to a design that was not made under it, which is the measurement A1 and A3
    exist for. `derna check --designs` takes the molecule on the command line, because
    these sequences carry no header to read it from, and refuses any sequence that is not
    a molecule of that description.
    """
    u5, u3 = flank_files(cfg, cell)
    argv = [cfg["bins"]["derna"], "check", "--designs", str(designs),
            "--protein", str(cell["fasta_path"]),
            "--utr5", u5, "--utr3", u3, "--forbid", cfg["forbid"],
            "--stop", cfg["stop"], "--out", str(out),
            "--threads", str(cfg["score_threads"]), "--force"]
    subprocess.run(pinned(cfg.get("score_cores"), argv), check=True)
    return read_scored(out)


def read_scored(path):
    """A check's point rows, in this experiment's vocabulary.

    The result names the constrained energy `mfe`, the free fold `free_mfe` and the
    verdict `modular`. They are renamed on the way in so that nothing downstream of here,
    and nothing in the paper's tables, has two names for one quantity, and so that a row
    from this route reads the same as one from the deposit.
    """
    rename = {"mfe": "constrained_kcal_mol", "free_mfe": "free_kcal_mol",
              "gap": "gap_kcal_mol", "modular": "separable"}
    out = []
    for i, rec in enumerate(read_points(path)):
        row = {rename.get(k, k): v for k, v in rec.items()}
        row.setdefault("index", i)
        out.append(row)
    return out


def check(cfg, front, out):
    """`derna check` on a front that was produced under the constraint.

    The header carries the molecule and the constraint, so a check cannot be
    pointed at a different molecule than the result was produced over. Threads
    change nothing it reports: each point's check reads only its own sequence.
    """
    argv = [cfg["bins"]["derna"], "check", str(front), "--out", str(out),
            "--threads", str(cfg["score_threads"]), "--force"]
    subprocess.run(pinned(cfg.get("score_cores"), argv), check=True)
    return read_points(out)


# The scored files keep the name the deposit's rows were written under, so a re-run's
# output sits beside a deposited copy under the same name rather than looking like a
# different artefact.
def score_a1(cfg, cell):
    points = read_points(cell_dir(cfg, cell) / "a1.jsonl")
    if not points:
        return []
    designs = designs_input(cfg, cell, points, "a1.designs.txt")
    return score_designs(cfg, cell, designs, cell_dir(cfg, cell) / "a1.census.jsonl")


def score_a3(cfg, cell, results):
    """The comparator's designs, scored as the molecules it emitted.

    Its own bytes are used rather than a molecule rebuilt from the coding
    sequence and the flanks: `parse_mrnafold` has already checked the two agree,
    and scoring what the tool actually returned keeps that check meaningful.
    """
    designs = [r["design"] for r in results.rows()
               if r["cell"] == cell["cell"] and r["arm"] == "A3" and r.get("design")]
    if not designs:
        return []
    seqs = designs_input(cfg, cell, designs, "a3.designs.txt")
    return score_designs(cfg, cell, seqs, cell_dir(cfg, cell) / "a3.census.jsonl")


def score_cell(cfg, cell, results, sink):
    """Every arm's verdicts for one cell, written as one record.

    A4's points carry their check already, by construction. It is re-run
    anyway: a point that came back from the structurally optimal walk and does
    not certify is the one failure this experiment could not otherwise see, and
    the cost is one fold per point on cores no measurement is using.
    """
    scored = {"cell": cell["cell"], "instance": cell["instance"], "pair": cell["pair"],
              "n_aa": cell["n_aa"], "residues": residues(cfg, cell)}
    scored["A1"] = so_share(score_a1(cfg, cell))
    scored["A3"] = so_share(score_a3(cfg, cell, results))
    for arm in ("A2", "A4"):
        front = cell_dir(cfg, cell) / f"{arm.lower()}.jsonl"
        if read_points(front):
            out = cell_dir(cfg, cell) / f"{arm.lower()}.checked.jsonl"
            scored[arm] = so_share(check(cfg, front, out))
        else:
            scored[arm] = so_share([])
    sink.append(scored)
    return scored


def so_share(rows):
    """The share of points that are structurally optimal, and the gaps.

    One function over both scorers, because the share is the headline and it has to be the
    same quantity for every arm. A scored given design carries the two energies; a checked
    `derna` point carries the verdict and the gap. A point
    whose gap a beam only bounded counts as not structurally optimal and is
    counted separately, since a bound that refutes is a refutation.
    """
    gaps, optimal, bounded = [], 0, 0
    for r in rows:
        gap = r.get("gap_kcal_mol", r.get("gap"))
        if gap is None and "gap_at_least" in r:
            gap, bounded = r["gap_at_least"], bounded + 1
        if gap is None:
            continue
        gaps.append(gap)
        optimal += 1 if abs(gap) < 1e-9 else 0
    if not gaps:
        return {"points": 0, "structurally_optimal": 0, "share": None}
    return {"points": len(gaps), "structurally_optimal": optimal,
            "share": optimal / len(gaps), "bounded_gaps": bounded,
            "gap_mean": sum(gaps) / len(gaps), "gap_max": max(gaps),
            "gap_median": sorted(gaps)[len(gaps) // 2]}


# ---- the price A4 pays -------------------------------------------------------------

def price(a2_points, a4_points):
    """A2's front against A4's, at matched CAI, and what a partial A4 bounds.

    The distance is reported per A4 point, which is the direction that stays
    meaningful when A4 is partial: every A4 point has an A2 point at its CAI or
    above, because A4's designs are a subset of A2's. The reverse is not true,
    so a mean over A2's points would read a missing region as a cost of zero.
    """
    if not a2_points or not a4_points:
        return {"matched": 0, "mfe_mean": None, "mfe_max": None}
    floor = sorted(((p["cai_log"], p["mfe"]) for p in a2_points), reverse=True)
    deltas = []
    for p in a4_points:
        best = min((m for c, m in floor if c >= p["cai_log"]), default=None)
        if best is not None:
            deltas.append(p["mfe"] - best)
    if not deltas:
        return {"matched": 0, "mfe_mean": None, "mfe_max": None}
    return {"matched": len(deltas), "mfe_mean": sum(deltas) / len(deltas),
            "mfe_max": max(deltas)}


def coding(point, utr5, utr3):
    """The coding region of a molecule point, flanks removed by length."""
    seq = point["sequence"]
    return seq[len(utr5):len(seq) - len(utr3)]


def binds(cell, a1_points, a2_points):
    """Whether the constraint changed the design, or only the energy label.

    A cost measured against a constraint means nothing until the constraint is
    shown to change something, and the check is a diff of the emitted designs
    rather than of their reported energies. It is reported per cell because the
    answer is a property of the length and the pair, not of the method: a
    separable constraint leaves A2's designs byte-identical to A1's and moves
    every energy by one constant, in which case A2 carries no design A1 does
    not and the arm's own share is the only thing it contributes.

    `shift_min` and `shift_max` are the per-point energy differences. Equal
    ones mean the constraint was separable on this cell; a spread means the
    constrained fold used something separability forbids -- a 5'-3' pair, which
    the two forbidden classes leave allowed -- whether or not that changed any
    design.
    """
    if not a1_points or not a2_points:
        return {"designs_compared": 0, "identical": None}
    a1 = [coding(p, cell["utr5"], cell["utr3"]) for p in a1_points]
    a2 = [coding(p, cell["utr5"], cell["utr3"]) for p in a2_points]
    shifts = [round(b["mfe"] - a["mfe"], 6) for a, b in zip(a1_points, a2_points)]
    return {"designs_compared": min(len(a1), len(a2)),
            "identical": sum(1 for x, y in zip(a1, a2) if x == y),
            "same_set": set(a1) == set(a2),
            "shift_min": min(shifts), "shift_max": max(shifts),
            "separable": len(set(shifts)) == 1}


# ---- commands ----------------------------------------------------------------------

def cmd_plan(cfg, args, results):
    """What the grid owes, against what it has.

    A merged output is not evidence a sweep ran: the counts are printed against
    the grid's own shape so a short table is distinguishable from a small one.
    """
    have = done_keys(results)
    grid = selected(cfg, args)
    owed = 0
    per_panel = {}
    for cell in grid:
        missing = [a for a in ("A1", "A2", "A4") if (cell["cell"], a, None) not in have]
        a3 = sum(1 for k in have if k[0] == cell["cell"] and k[1] == "A3")
        owed += len(missing)
        tally = per_panel.setdefault(f"{cell['panel']}/{cell['band']}", [0, 0])
        tally[0] += 1
        tally[1] += len(missing)
        print(f"{cell['cell']:<28} n={cell['n_aa']:<5} {cell['panel']:<6} {cell['band']:<6} "
              f"missing={','.join(missing) or '-':<12} A3 weights done={a3}")
    print()
    for name in sorted(per_panel):
        held, due = per_panel[name]
        print(f"{name:<18} cells={held:<5} timed arms owed={due}")
    print(f"\ncells={len(grid)} arms={len(ARMS)} timed cells owed={owed}")


def selected(cfg, args):
    """The cells this invocation is about: the grid, or one named cell of it."""
    grid = cells(cfg)
    if not args.cell:
        return grid
    chosen = [c for c in grid if c["cell"] == args.cell]
    if not chosen:
        raise SystemExit(f"no cell {args.cell} in this grid")
    return chosen


def cmd_run(cfg, args, results):
    """The grid, one cell at a time, arms in dependency order.

    A1 and A2 are independent; A3 and A4 both read A2, so a cell whose A2 left
    nothing records them as missing rather than inventing a front for them. The
    arms are not run concurrently: each is timed, and two timed cells on one
    host measure each other.
    """
    have = done_keys(results)
    for cell in selected(cfg, args):
        cell["fasta_path"] = fasta_for(cfg, cell)
        for arm in [a for a in ARMS if a in args.arms.split(",")]:
            if (cell["cell"], arm, None) in have:
                continue
            print(f"{cell['cell']} {arm}: {run_arm(cfg, cell, arm, results, have)}", flush=True)


def run_arm(cfg, cell, arm, results, have):
    """One arm of one cell.

    A3 is the one arm whose unit of work is not the cell: it records a row per
    weight, so it is never "done" as a whole and re-entering it continues from
    the weights it has rather than repeating them.
    """
    if arm == "A3":
        return run_a3(cfg, cell, results, have)
    return {"A1": run_a1, "A2": run_a2, "A4": run_a4}[arm](cfg, cell, results)


def fasta_for(cfg, cell):
    """The protein as a one-record FASTA beside the cell.

    Written from the sequence the driver read, not referenced where it came
    from, so the molecule every arm of a cell designs is one file and the row
    can name it.
    """
    path = cell_dir(cfg, cell) / "protein.fa"
    path.write_text(f">{cell['instance']}\n{cell['protein']}\n")
    return path


def cmd_score(cfg, args, results):
    """Every cell's verdicts, on the scoring cores, after the timing."""
    sink = Results(Path(cfg["out_dir"]) / "scored.jsonl", results.stamp)
    done = {r["cell"] for r in sink.rows()}
    for cell in selected(cfg, args):
        if cell["cell"] in done and not args.again:
            continue
        scored = score_cell(cfg, cell, results, sink)
        shares = " ".join(f"{a}={fmt_share(scored[a])}" for a in ARMS)
        print(f"{cell['cell']:<28} {shares}", flush=True)


def fmt_share(s):
    if not s or s["share"] is None:
        return "-"
    return f"{s['structurally_optimal']}/{s['points']}"


def cmd_summarize(cfg, args, results):
    """The per-pair table and the price, with the row counts checked first."""
    sink = Results(Path(cfg["out_dir"]) / "scored.jsonl", results.stamp)
    scored = {r["cell"]: r for r in sink.rows()}
    grid = selected(cfg, args)
    print(f"cells: {len(scored)} scored of {len(grid)} in the grid")
    for cell in grid:
        s = scored.get(cell["cell"])
        if not s:
            continue
        a1 = read_points(cell_dir(cfg, cell) / "a1.jsonl")
        a2 = read_points(cell_dir(cfg, cell) / "a2.jsonl")
        a4 = read_points(cell_dir(cfg, cell) / "a4.jsonl")
        p, b = price(a2, a4), binds(cell, a1, a2)
        shares = " ".join(f"{a}={fmt_share(s[a])}" for a in ARMS)
        mean = "-" if p["mfe_mean"] is None else f"{p['mfe_mean']:+.2f}"
        print(f"{cell['cell']:<28} {shares:<44} price(mean)={mean} matched={p['matched']} "
              f"A1=A2 designs {b['identical']}/{b['designs_compared']} "
              f"shift {b.get('shift_min')}..{b.get('shift_max')}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("plan", "run", "score", "summarize"))
    ap.add_argument("--config", required=True)
    ap.add_argument("--cell", help="one cell only, as `INSTANCE__pair`")
    ap.add_argument("--again", action="store_true", help="re-score cells already scored")
    ap.add_argument("--arms", default=",".join(ARMS),
                    help="the arms `run` runs, comma-separated (default: all, in order)")
    args = ap.parse_args()
    cfg = load_cfg(args.config)
    Path(cfg["out_dir"]).mkdir(parents=True, exist_ok=True)
    arm_stamps = {arm: provenance(arm_cfg(cfg, arm)) for arm in cfg["arm_builds"]}
    results = Results(Path(cfg["out_dir"]) / "results.jsonl", provenance(cfg), arm_stamps)
    {"plan": cmd_plan, "run": cmd_run, "score": cmd_score,
     "summarize": cmd_summarize}[args.command](cfg, args, results)


main()
