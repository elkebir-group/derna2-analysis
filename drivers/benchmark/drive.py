#!/usr/bin/env python3
"""Core paper benchmark driver: correctness, single-λ runtime and memory, and
the approximate modes' speed against quality, in one grid.

Two passes, split by what host load does to each measurement:

  work   correctness, peak RSS and instruction counts. None of these move with
         host load, so this pass runs now, in parallel, on whatever cores the
         config gives it. It is most of the grid's cost.
  time   wall-clock and CPU time. These do move with load, so this pass runs
         only on a reserved core set, waits until those cores are idle, runs one
         cell at a time with every arm interleaved inside the cell, and records
         how busy the rest of the reserved set was during each run. A run that
         shared its cores is kept and flagged, not dropped.

The time pass runs only the cells the work pass finished, so it never spends
the quiet window rediscovering a timeout.

Subcommands:
  plan       print the grid and the number of rows each pass should write
  selftest   run every enabled arm on a short protein and check what it emits
  epscheck   check ε against its bound, and that ε/2 leaves every exact design unchanged
  work       run the work pass (resumable)
  time       run the time pass (resumable)
  summarize  check row counts, exit categories and cross-arm agreement

    python3 drive.py <subcommand> --config config.json [--arms a,b] [--corners mfe] [--max-n N]

Every run is wrapped by `measure.py`, which opens its own hardware counters, so
concurrent cells never count each other's instructions.
"""

import argparse
import fcntl
import json
import math
import os
import re
import struct
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
MEASURE = HERE / "measure.py"
CENSORED = ("timeout", "oom")
# Left out of the work pass on purpose (too long to be worth its correctness
# check there) but still timed: the time pass re-scores its own output.
SKIPPED = "skipped"
TIMEABLE = ("ok", SKIPPED)

# The corners are run at ε and 1 − ε rather than at 0 and 1, so each corner is
# the lexicographic optimum: the second objective breaks ties and never trades
# against the first. derna2 and LinearCDSfold's DN objective weight energy by λ;
# LinearDesign and mRNAfold weight CAI by λ' with energy at weight 1, and
# λ' = (1 − λ) / λ is the same point.
# How small ε must be is a bound in n: see `epsilon_bounds`.


def corner_lambda(kind, corner, eps):
    lam = 1.0 - eps if corner == "mfe" else eps
    if kind in ("lineardesign", "mrnafold"):
        return (1.0 - lam) / lam
    return lam


# Peak RSS per n² (kB), measured in the Scope 24 benchmark; the admission gate's
# estimate when an arm does not set its own.
RSS_K = {"lineardesign": 33.0, "derna": 11.0, "lcds": 7.0, "cdsfold": 33.0, "mrnafold": 4.0}

THREE_LETTER = dict(A="Ala", R="Arg", N="Asn", D="Asp", C="Cys", Q="Gln", E="Glu", G="Gly", H="His", I="Ile",
                    L="Leu", K="Lys", M="Met", F="Phe", P="Pro", S="Ser", T="Thr", W="Trp", Y="Tyr", V="Val")


def read_codon_table(path):
    """[(one-letter amino acid, RNA codon, fraction)] from a `amino_acid,codon,fraction` CSV."""
    rows = []
    for line in open(path):
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


# ---- instances -----------------------------------------------------------------

def read_fasta(path):
    seq = "".join(line.strip() for line in open(path) if not line.startswith(">"))
    return seq.rstrip("*")


def instances(cfg):
    """(label, protein) for every instance: protein directories scanned (never
    enumerated) and cut at a length, named single proteins, and a prefix ladder
    of one protein."""
    spec = cfg["instances"]
    out = []
    for d in spec.get("protein_dirs", []):
        for path in sorted(Path(d["dir"]).glob("*.fasta")):
            prot = read_fasta(path)
            if len(prot) <= d.get("max_aa", 10**9):
                out.append((path.stem, prot))
    for path in spec.get("protein_files", []):
        out.append((Path(path).stem, read_fasta(path)))
    ladder = spec.get("ladder")
    if ladder:
        full = read_fasta(ladder["fasta"])
        stem = Path(ladder["fasta"]).stem
        out += [(f"{stem}_p{n}", full[:n]) for n in ladder["lengths"] if n <= len(full)]
    return out


# ---- arms ----------------------------------------------------------------------

def variants(name, arm):
    """Expand one arm into its variants: one per beam width or thread count."""
    if "widths" in arm:
        return [dict(arm=name, variant=f"b{w}", width=w, ncores=1) for w in arm["widths"]]
    if "threads" in arm:
        return [dict(arm=name, variant=f"t{t}", threads=t, ncores=t) for t in arm["threads"]]
    return [dict(arm=name, variant="-", ncores=1)]


def grid(cfg, args):
    """Every (instance, corner, arm variant) cell the filters keep, smallest first."""
    want_arms = set(args.arms.split(",")) if args.arms else None
    want_corners = args.corners.split(",") if args.corners else cfg["corners"]
    cells = []
    for label, prot in instances(cfg):
        if args.max_n and len(prot) > args.max_n:
            continue
        for corner in want_corners:
            for name, arm in cfg["arms"].items():
                if not arm.get("enabled", True) or (want_arms and name not in want_arms):
                    continue
                if corner not in arm.get("corners", cfg["corners"]):
                    continue
                for v in variants(name, arm):
                    cells.append(dict(v, label=label, prot=prot, n=len(prot), corner=corner))
    return sorted(cells, key=lambda c: (c["n"], c["label"], c["corner"]))


def ext(stem, suffix):
    """`stem` plus a suffix. Cell ids contain dots, so `with_suffix` would
    replace the corner instead of appending."""
    return Path(f"{stem}{suffix}")


def cell_id(c):
    return f"{c['arm']}.{c['variant']}.{c['label']}.{c['corner']}"


def command(cfg, c, work):
    """(argv, cwd, stdin_path, output_parser) for one cell."""
    arm = cfg["arms"][c["arm"]]
    kind = arm["kind"]
    bins = cfg["bins"]
    stem = work / cell_id(c)
    fasta = ext(stem, ".fasta")
    fasta.write_text(f">{c['label']}\n{c['prot']}\n")
    if kind == "derna":
        lam = ext(stem, ".lambda.jsonl")
        lam.write_text(json.dumps({"lambda": corner_lambda("derna", c["corner"], cfg["epsilon"])}) + "\n")
        out = ext(stem, ".out.jsonl")
        argv = [bins["derna"], "--fasta", str(fasta), "--codon-usage", cfg["tables"]["derna"],
                "--lambdas", str(lam), "--out", str(out), "--machine", str(ext(stem, ".machine.json")),
                "--threads", str(c.get("threads", 1))]
        argv += [a.replace("{w}", str(c.get("width", ""))) for a in arm.get("args", [])]
        return argv, str(REPO), None, lambda _stdout: parse_derna(out, ext(stem, ".machine.json"))
    if kind == "lcds":
        out = ext(stem, ".out.txt")
        argv = [bins["lcds"], "-m", arm["mode"], "-O", "DN", "-l", repr(corner_lambda("lcds", c["corner"], cfg["epsilon"])),
                "-c", bins["lcds_table"], "-o", str(out), str(fasta)]
        if "width" in c:
            argv += ["-b", str(c["width"])]
        return argv, bins["lcds_dir"], None, lambda _stdout: parse_lcds(out.read_text())
    if kind == "lineardesign":
        argv = [bins["lineardesign"], repr(corner_lambda("lineardesign", c["corner"], cfg["epsilon"])), "0", bins["lineardesign_table"]]
        return argv, bins["lineardesign_dir"], str(ext(stem, ".stdin")), parse_lineardesign
    if kind == "cdsfold":
        return [bins["cdsfold"], str(fasta)], str(work), None, parse_cdsfold
    if kind == "mrnafold":
        lam = corner_lambda("mrnafold", c["corner"], cfg["epsilon"])
        threads = c.get("threads", 1)
        conf = ext(stem, ".cfg")
        conf.write_text(f"aa_seq {c['prot']}\nlambda {lam!r}\nnum_subopt_traces 1\n"
                        f"codon_freq_path {mrnafold_table(cfg, work)}\nparallel {'true' if threads > 1 else 'false'}\n")
        binary = bins["mrnafold_parallel"] if threads > 1 else bins["mrnafold"]
        weights = codon_weights(read_codon_table(cfg["tables"]["canonical"]))
        return [binary, str(conf)], str(work), None, lambda stdout: parse_mrnafold(stdout, lam, weights)
    raise ValueError(f"unknown arm kind {kind}")


def mrnafold_table(cfg, work):
    """The canonical table in mRNAfold's `Ala GCT 0.26` format, written once per
    work directory. Same fractions, so the same weights."""
    path = work / "mrnafold_codon_table.txt"
    if not path.exists():
        lines = [f"{THREE_LETTER[aa]} {codon.replace('U', 'T')} {f!r}"
                 for aa, codon, f in read_codon_table(cfg["tables"]["canonical"]) if aa in THREE_LETTER]
        tmp = ext(path, f".{os.getpid()}")
        tmp.write_text("\n".join(lines) + "\n")
        os.replace(tmp, path)
    return path


def parse_derna(out, machine):
    rec = json.loads(out.read_text().splitlines()[0])
    row = dict(seq=rec["rna_seq"], mfe_reported=rec["mfe_kcal_mol"], cai_reported=rec["cai"])
    if machine.exists():
        row["vmhwm_kb"] = json.loads(machine.read_text()).get("vmhwm_kb")
    return row


def parse_lcds(txt):
    seq = re.search(r"\n([ACGU]+)\n", txt)
    mfe = re.search(r"Folding free energy: ([-0-9.]+)", txt)
    cai = re.search(r"CAI: ([0-9.]+)", txt)
    return dict(seq=seq.group(1), mfe_reported=float(mfe.group(1)), cai_reported=float(cai.group(1)))


def parse_lineardesign(stdout):
    seq = re.search(r"mRNA sequence:\s*(\S+)", stdout)
    mfe = re.search(r"folding free energy:\s*([-0-9.]+)", stdout)
    cai = re.search(r"mRNA CAI:\s*([0-9.]+)", stdout)
    return dict(seq=seq.group(1), mfe_reported=float(mfe.group(1)), cai_reported=float(cai.group(1)))


def parse_cdsfold(stdout):
    """CDSfold's output format is not yet verified against a build: `selftest`
    must pass for this arm before its rows are read. It reports MFE only, so
    its CAI is a tie-break and never a loss."""
    seqs = re.findall(r"^([ACGU]{30,})$", stdout, re.M)
    mfe = re.search(r"(-?\d+\.\d+)\s*kcal", stdout) or re.search(r"MFE[^-\d]*(-?\d+\.?\d*)", stdout)
    return dict(seq=max(seqs, key=len), mfe_reported=float(mfe.group(1)), cai_reported=None)


def parse_mrnafold(stdout, lam, weights):
    """mRNAfold prints the design, its structure, `Score` (the objective: the
    structure's energy plus λ'·Σ −ln w) and `EFE`, an ensemble free energy from
    a partition function it runs afterwards. EFE is not the MFE, so the design's
    energy is recovered from Score by adding the CAI term back."""
    seq = re.search(r"^([ACGU]+)$", stdout, re.M).group(1)
    score = float(re.search(r"^Score: ([-0-9.eE+]+)", stdout, re.M).group(1))
    cai = float(re.search(r"^CAI: ([0-9.]+)", stdout, re.M).group(1))
    log_w = sum(math.log(weights[seq[i:i + 3]]) for i in range(0, len(seq) - 2, 3))
    return dict(seq=seq, mfe_reported=round(score + lam * log_w, 2), cai_reported=cai)


# ---- re-scoring ------------------------------------------------------------------

class Rescorer:
    """Back-translation, CAI on the canonical table and MFE on the common folder,
    so every arm is read on one energy model and one codon table."""

    def __init__(self, cfg):
        self.fold = cfg["bins"]["mfe_fold"]
        self.cores = cfg.get("rescore_cores")
        rows = read_codon_table(cfg["tables"]["canonical"])
        self.aa = {codon: aa for aa, codon, _ in rows}
        self.w = codon_weights(rows)

    def __call__(self, seq, prot):
        codons = [seq[i:i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
        translated = "".join(self.aa.get(c, "?") for c in codons).rstrip("*")
        logs = [math.log(self.w[c]) for c in codons if c in self.w]
        argv = (["taskset", "-c", self.cores] if self.cores else []) + [self.fold]
        r = subprocess.run(argv, input=seq + "\n", capture_output=True, text=True, cwd=str(REPO))
        mfe = int(r.stdout.split("\t")[1]) / 100.0 if "\t" in r.stdout else None
        return dict(translates=translated == prot, cai_rescored=math.exp(sum(logs) / len(logs)) if logs else None,
                    mfe_rescored=mfe)


# ---- host occupancy ----------------------------------------------------------------

def cpu_jiffies():
    """{cpu: (busy, total)} from /proc/stat."""
    out = {}
    for line in open("/proc/stat"):
        m = re.match(r"cpu(\d+) (.*)", line)
        if m:
            f = [int(x) for x in m.group(2).split()]
            idle = f[3] + f[4]
            out[int(m.group(1))] = (sum(f[:8]) - idle, sum(f[:8]))
    return out


def busy_between(a, b, cpus):
    """Per-cpu busy fraction between two /proc/stat snapshots."""
    return {c: (b[c][0] - a[c][0]) / max(1, b[c][1] - a[c][1]) for c in cpus}


class Sampler(threading.Thread):
    """Samples per-core busy fractions and the load average at a fixed interval,
    appending each sample to a timeline file, so any run can be labelled by what
    else was on the reserved cores while it ran."""

    def __init__(self, path, cpus, interval):
        super().__init__(daemon=True)
        self.path, self.cpus, self.interval = path, cpus, interval
        self.samples, self.lock = [], threading.Lock()

    def run(self):
        prev, t_prev = cpu_jiffies(), time.time()
        while True:
            time.sleep(self.interval)
            cur, t = cpu_jiffies(), time.time()
            row = dict(t0=t_prev, t=t, load=open("/proc/loadavg").read().split()[:3],
                       busy=busy_between(prev, cur, self.cpus))
            prev, t_prev = cur, t
            with self.lock:
                self.samples.append(row)
            with open(self.path, "a") as f:
                f.write(json.dumps(row) + "\n")

    def window(self, t0, t1, cores_at, burst=1):
        """Over the sample intervals wholly inside [t0, t1], each read on the
        cores `cores_at(a, b)` names for that interval: the mean of each
        interval's busiest core, the highest such mean over any `burst`
        consecutive intervals, and the fewest cores any interval had to read.
        None for each when no interval fits. The burst mean is taken over
        several intervals because a host daemon spiking one core for a few
        seconds every few minutes would otherwise flag nearly every long run,
        while a neighbour that stays for a minute or more is still caught.
        The cores are chosen per interval because, with several runs of ours on
        the socket, which cores are ours changes many times over a long run."""
        with self.lock:
            inside = [s for s in self.samples if s["t0"] >= t0 and s["t"] <= t1]
        worst, fewest = [], None
        for s in inside:
            cores = cores_at(s["t0"], s["t"])
            fewest = len(cores) if fewest is None else min(fewest, len(cores))
            if cores:
                worst.append(max(s["busy"][c] for c in cores))
        if not worst:
            return None, None, fewest
        k = max(1, min(burst, len(worst)))
        peak = max(sum(worst[i:i + k]) / k for i in range(len(worst) - k + 1))
        return sum(worst) / len(worst), peak, fewest


class Holds:
    """Which reserved cores this driver's own runs held, and over which
    interval. With several cells on the socket at once, a neighbour of ours is
    expected to be busy; only a core none of our runs held during an interval
    says something else was there."""

    def __init__(self):
        self.spans, self.lock = [], threading.Lock()

    def start(self, cores):
        span = dict(cores=list(cores), t0=time.time(), t1=None)
        with self.lock:
            self.spans.append(span)
        return span

    def end(self, span):
        span["t1"] = time.time()

    def held(self, a, b, exclude):
        with self.lock:
            return {c for s in self.spans if s is not exclude and s["t0"] <= b and (s["t1"] is None or s["t1"] >= a)
                    for c in s["cores"]}


def wait_quiet(cpus, gate, log):
    """Block until every core `cpus()` names sits below `busy_max` for `samples`
    consecutive intervals. The gate reads the reserved cores themselves, not the
    host load average: with the host's other work confined to the other socket
    the load average never drops, and it also counts this driver's own threads.
    `cpus` is re-read each interval, so cores this driver's own runs hold are
    left out as they come and go."""
    streak, prev = 0, cpu_jiffies()
    while streak < gate["samples"]:
        time.sleep(gate["interval_s"])
        cur = cpu_jiffies()
        watched = cpus()
        worst = max(busy_between(prev, cur, watched).values()) if watched else 0.0
        prev = cur
        streak = streak + 1 if worst < gate["busy_max"] else 0
        if gate["samples"] > 2 or worst >= gate["busy_max"]:
            log(f"gate: busiest free reserved core {worst:.3f}, streak {streak}/{gate['samples']}")


# ---- running cells -----------------------------------------------------------------

class Results:
    """Append-only JSONL with a lock, keyed for resume."""

    def __init__(self, path):
        self.path, self.lock = path, threading.Lock()

    def rows(self):
        if not self.path.exists():
            return []
        return [json.loads(line) for line in open(self.path) if line.strip()]

    def append(self, row):
        with self.lock, open(self.path, "a") as f:
            f.write(json.dumps(row) + "\n")


class CorePool:
    """Hands out sets of cores, so a threaded cell gets as many as it has threads."""

    def __init__(self, cores):
        self.free, self.cv = list(cores), threading.Condition()

    def take(self, k):
        with self.cv:
            self.cv.wait_for(lambda: len(self.free) >= k)
            got, self.free = self.free[:k], self.free[k:]
            return got

    def give(self, got):
        with self.cv:
            self.free += got
            self.cv.notify_all()


def mem_available_kb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable:"):
            return int(line.split()[1])
    return 0


class MemGate:
    """Admits a cell only when its estimated peak fits both the driver's own
    budget and what the host actually has free, less a headroom. The host is
    shared: an over-admitted cell would put the kernel's OOM killer among other
    people's jobs, not just this sweep's. A cell larger than the whole budget is
    refused rather than admitted alone."""

    def __init__(self, budget_kb, headroom_kb):
        self.budget, self.headroom = budget_kb, headroom_kb
        self.used, self.cv = 0, threading.Condition()

    def take(self, est):
        if est > self.budget:
            return False
        with self.cv:
            while not (self.used + est <= self.budget and mem_available_kb() - est >= self.headroom):
                self.cv.wait(timeout=30)
            self.used += est
        return True

    def give(self, est):
        with self.cv:
            self.used -= est
            self.cv.notify_all()


class RssModel:
    """Peak RSS per n² for each arm variant, taken from the longest finished
    cell of that variant (the largest k at that length), or from the configured
    prior until one has run. Memory follows a clean power law in n where time
    does not, so the measured k is the better estimate. It is read at the
    longest length rather than maximized over all of them: a short run's peak
    is mostly the binary's fixed footprint, and its k overstates a long run's
    threefold (derna2: 25.9 kB at n = 30 against 8.4 at n = 549)."""

    def __init__(self, cfg, rows, margin=1.25):
        self.cfg, self.margin, self.k, self.lock = cfg, margin, {}, threading.Lock()
        for r in rows:
            if r.get("exit") == "ok":
                self.observe(r["arm"], r["variant"], r["n_aa"], r["maxrss_kb"])

    def observe(self, arm, variant, n, rss_kb):
        with self.lock:
            key = (arm, variant)
            self.k[key] = max(self.k.get(key, (0, 0.0)), (n, rss_kb / n**2))

    def estimate(self, c):
        arm = self.cfg["arms"][c["arm"]]
        with self.lock:
            seen = self.k.get((c["arm"], c["variant"]))
        k = seen[1] if seen else arm.get("rss_k", RSS_K[arm["kind"]])
        return int(self.margin * k * c["n"] ** 2)


def run_once(cfg, c, cores, work, tag):
    """Run one cell once under measure.py on `cores`; returns the row fields."""
    argv, cwd, stdin, parse = command(cfg, c, work)
    stem = work / f"{cell_id(c)}.{tag}"
    if stdin:
        Path(stdin).write_text(c["prot"] + "\n")
    meas = ext(stem, ".measure.json")
    wrapped = ["taskset", "-c", ",".join(map(str, cores)), sys.executable, str(MEASURE), str(meas),
               "--timeout", str(cfg["timeout_s"]), "--stdout", str(ext(stem, ".stdout")),
               "--stderr", str(ext(stem, ".stderr"))]
    wrapped += (["--stdin", stdin] if stdin else []) + ["--"] + argv
    subprocess.run(wrapped, cwd=cwd)
    m = json.loads(meas.read_text())
    row = {k: m[k] for k in ("exit", "wall_s", "user_s", "sys_s", "maxrss_kb", "instructions", "cycles")}
    row["cores"] = cores
    row["cpu_share"] = (m["user_s"] + m["sys_s"]) / max(1e-9, m["wall_s"] * len(cores))
    if row["exit"] == "ok":
        try:
            row.update(parse(ext(stem, ".stdout").read_text()))
        except (AttributeError, IndexError, KeyError, ValueError, OSError) as err:
            row.update(exit="unparsed", note=repr(err))
    return row


def ident(c, cfg):
    return dict(arm=c["arm"], variant=c["variant"], kind=cfg["arms"][c["arm"]]["kind"], protein=c["label"],
                n_aa=c["n"], corner=c["corner"], ncores=c["ncores"])


def run_work(cfg, args, log):
    out = Path(cfg["out_dir"])
    results = Results(out / "work.jsonl")
    done = {r["cell"] for r in results.rows() if r["exit"] in ("ok", "projected", SKIPPED) + CENSORED}
    cells = [c for c in grid(cfg, args) if cell_id(c) not in done]
    log(f"work pass: {len(cells)} cells to run")
    rescore, pool = Rescorer(cfg), CorePool(cfg["work"]["cores"])
    mem = MemGate(cfg["mem_budget_gb"] * 1024 * 1024, cfg.get("mem_headroom_gb", 30) * 1024 * 1024)
    rss = RssModel(cfg, results.rows())
    failed, flock = {}, threading.Lock()

    def one(c):
        key = (c["arm"], c["variant"], c["corner"])
        with flock:
            at = failed.get(key)
        if at is not None and c["n"] >= at:
            results.append(dict(ident(c, cfg), cell=cell_id(c), exit="projected",
                                note=f"censored at n={at} on this arm variant and corner"))
            return
        est = rss.estimate(c)
        if not mem.take(est):
            results.append(dict(ident(c, cfg), cell=cell_id(c), exit="deferred-memory",
                                note=f"estimated {est / 2**20:.0f} GB over the {cfg['mem_budget_gb']} GB budget"))
            log(f"[deferred-memory] {cell_id(c)} est {est / 2**20:.0f} GB")
            return
        cores = pool.take(c["ncores"])
        try:
            row = run_once(cfg, c, cores, out / "work", "work")
        finally:
            pool.give(cores)
            mem.give(est)
        if row["exit"] == "ok":
            rss.observe(c["arm"], c["variant"], c["n"], row["maxrss_kb"])
        if row["exit"] in CENSORED:
            with flock:
                failed[key] = min(c["n"], failed.get(key, 10**9))
        if row.get("seq"):
            row.update(rescore(row["seq"], c["prot"]))
        if args.drop_seq:
            row.pop("seq", None)
        results.append(dict(ident(c, cfg), cell=cell_id(c), load=os.getloadavg()[0], **row))
        log(f"[{row['exit']}] {cell_id(c)} wall={row['wall_s']:.1f}s rss={row['maxrss_kb']}kB "
            f"instr={(row['instructions'] or 0) / 1e9:.1f}G")

    with ThreadPoolExecutor(max_workers=len(cfg["work"]["cores"]) + 4) as ex:
        list(ex.map(one, cells))
    log("work pass: done")


def time_rows_by_key(rows):
    """{(cell, rep): [rows]} in file order."""
    out = {}
    for r in rows:
        out.setdefault((r["cell"], r["rep"]), []).append(r)
    return out


def best_time_row(rows):
    """The row a (cell, rep) is read from: its latest clean run, else its latest run."""
    clean = [r for r in rows if r["exit"] == "ok" and not r.get("contaminated")]
    return (clean or rows)[-1]


def run_time(cfg, args, log):
    """Time every cell the work pass finished, up to `time.concurrency` at once
    on the reserved cores. Cells of one instance and corner are queued together,
    arms rotated per rep, and instances longest first so the last cell to start
    is a short one. A contaminated run is kept, and the cell is re-run after the
    free reserved cores have been quiet again, up to `time.retries` times."""
    out = Path(cfg["out_dir"])
    spec = cfg["time"]
    reserved = spec["cores"]
    concurrency = spec.get("concurrency", 1)
    retries = spec.get("retries", 2)
    # --redo-contaminated grants a cell that used up its retries another round of them.
    max_attempts = (retries + 1) * (2 if args.redo_contaminated else 1)
    work_rows = {r["cell"]: r for r in Results(out / "work.jsonl").rows() if r["exit"] in TIMEABLE}
    results = Results(out / "time.jsonl")
    attempts = {k: len(v) for k, v in time_rows_by_key(results.rows()).items()}
    settled = {k for k, v in time_rows_by_key(results.rows()).items()
               if any(r["exit"] != "ok" or not r.get("contaminated") for r in v)
               or len(v) >= max_attempts}
    groups = {}
    for c in grid(cfg, args):
        if cell_id(c) in work_rows:
            groups.setdefault((c["label"], c["corner"]), []).append(c)
    # A skipped cell has no wall of its own; its row carries an estimate instead,
    # so the longest instance-corners still go first.
    order = sorted(groups, key=lambda g: -sum(work_rows[cell_id(c)].get("wall_est_s") or work_rows[cell_id(c)]["wall_s"]
                                              for c in groups[g]))
    queue = [(rep, c, g) for g in order for rep in range(spec["reps"]) for c in rotate(groups[g], rep)
             if (cell_id(c), rep) not in settled]
    for _, c, _ in queue:
        if c["ncores"] > len(reserved):
            log(f"skip {cell_id(c)}: needs {c['ncores']} cores, {len(reserved)} reserved")
    queue = [q for q in queue if q[1]["ncores"] <= len(reserved)]
    log(f"time pass: {len(queue)} runs over {len(order)} instance-corners, up to {concurrency} at once")

    sampler = Sampler(out / "time-occupancy.jsonl", reserved, spec["sample_s"])
    sampler.start()
    burst_n = max(1, round(spec.get("contamination_peak_window_s", spec["sample_s"]) / spec["sample_s"]))
    holds, free, flock = Holds(), list(reserved), threading.Lock()
    mem = MemGate(cfg["mem_budget_gb"] * 1024 * 1024, cfg.get("mem_headroom_gb", 30) * 1024 * 1024)
    rss = RssModel(cfg, list(work_rows.values()))
    rescore = Rescorer(cfg)

    def free_now():
        with flock:
            return sorted(free)

    def one(rep, c, cores, span, est):
        before = cpu_jiffies()
        try:
            row = run_once(cfg, c, cores, out / "time", f"r{rep}")
        finally:
            after = cpu_jiffies()
            holds.end(span)
            # Released here, not by the scheduler, which may itself be blocked
            # in the memory gate waiting for exactly this.
            mem.give(est)
            with flock:
                free.extend(cores)

        def others(a, b):
            busy_ours = holds.held(a, b, span)
            return [x for x in reserved if x not in cores and x not in busy_ours]
        mean, peak, fewest = sampler.window(span["t0"], span["t1"], others, burst_n)
        if mean is None:
            # Shorter than a sampler interval: read the run's two ends instead,
            # on the cores no run of ours touched at any point during it.
            checked = others(span["t0"], span["t1"])
            mean = max(busy_between(before, after, checked).values()) if checked else None
            fewest = len(checked)
        row["others_busy_max"], row["others_busy_peak"], row["others_checked"] = mean, peak, fewest
        row["contaminated"] = contaminated(row, c, spec)
        if row.get("seq"):
            row.update(rescore(row["seq"], c["prot"]))
        if args.drop_seq:
            row.pop("seq", None)
        return row

    wait_quiet(lambda: reserved, spec["gate"], log)
    running, last_group, regate = {}, None, None
    with ThreadPoolExecutor(max_workers=concurrency) as ex:
        while queue or running:
            while queue and len(running) < concurrency and queue[0][1]["ncores"] <= len(free_now()):
                rep, c, g = queue.pop(0)
                if regate or g != last_group:
                    wait_quiet(free_now, regate or spec["gate_between"], log)
                    regate, last_group = None, g
                est = rss.estimate(c)
                if not mem.take(est):
                    log(f"skip {cell_id(c)}: estimated {est / 2**20:.0f} GB over the budget")
                    continue
                with flock:
                    cores, free[:] = free[:c["ncores"]], free[c["ncores"]:]
                span = holds.start(cores)
                running[ex.submit(one, rep, c, cores, span, est)] = (rep, c, g)
            if not running:
                continue
            done = next(as_completed(running))
            rep, c, g = running.pop(done)
            row = done.result()
            key = (cell_id(c), rep)
            attempts[key] = attempts.get(key, 0) + 1
            results.append(dict(ident(c, cfg), cell=cell_id(c), rep=rep, attempt=attempts[key],
                                load=os.getloadavg()[0], **row))
            flag = ""
            if row["exit"] == "ok" and row["contaminated"]:
                flag = " CONTAMINATED"
                if attempts[key] < max_attempts:
                    queue.insert(0, (rep, c, g))
                    regate = spec.get("gate_after_contamination", spec["gate"])
                    flag += ", re-queued"
            log(f"[{row['exit']}] {cell_id(c)} rep {rep} wall={row['wall_s']:.2f}s cpu_share={row['cpu_share']:.3f} "
                f"others avg={row['others_busy_max']} peak={row['others_busy_peak']} "
                f"over {row['others_checked']} cores{flag}")
    log("time pass: done")


def contaminated(row, c, spec):
    """A run shared the reserved socket if a reserved core none of our own runs
    held was busy while it ran: on average over the sampler intervals inside
    the run (or, for a run shorter than one, between /proc/stat reads at its two
    ends), or over any `contamination_peak_window_s` stretch of it, so a burst is not
    diluted by a long run. A single-core
    run that got less than its whole core is flagged too, but only past
    `cpu_share_after_s`: below that, process start-up alone keeps the share
    under 1."""
    if (row["others_busy_max"] or 0) > spec["contamination_busy_max"]:
        return True
    if (row.get("others_busy_peak") or 0) > spec.get("contamination_peak_max", 1.0):
        return True
    return (c["ncores"] == 1 and row["wall_s"] >= spec["cpu_share_after_s"]
            and row["cpu_share"] < spec["min_cpu_share"])


def rotate(xs, k):
    """Each rep starts one arm later, so no arm always runs first in a cell."""
    k %= max(1, len(xs))
    return xs[k:] + xs[:k]


# ---- ε -------------------------------------------------------------------------------

def epsilon_bounds(table, n, kcal_per_nt=2.0):
    """The largest ε at which each corner is still lexicographic, for a protein
    of n residues, in derna2's objective λ·E[dcal] − (1 − λ)·100·Σ ln w.

    MFE corner: the CAI term must never pay for one dcal, and Σ ln w spans at
    most n·L, L the worst codon's log penalty. CAI corner: no energy gain may pay
    for leaving a CAI-optimal codon, which costs at least δ (the smallest
    non-zero log penalty); |MFE| is bounded by `kcal_per_nt` per nucleotide."""
    by_aa = {}
    for line in open(table):
        if line.startswith("amino_acid") or not line.strip():
            continue
        aa, _codon, frac = line.strip().split(",")
        by_aa.setdefault(aa, []).append(float(frac))
    penalties = [-math.log(f / max(fs)) for fs in by_aa.values() for f in fs]
    worst = max(penalties)
    delta = min(p for p in penalties if p > 1e-12)
    energy_span = kcal_per_nt * 100 * 3 * n
    return dict(mfe=1 / (100 * n * worst + 1), cai=100 * delta / (energy_span + 100 * delta))


def single_exact(x):
    return struct.unpack("f", struct.pack("f", x))[0] == x


def epscheck(cfg, args, log):
    """The configured ε against its bound at the longest instance, then every
    enabled exact arm at ε and ε/2 on the longest instance under --max-n: the
    re-scored (MFE, CAI) must not move. A tool that rounds λ internally fails
    here even when the bound holds.

    LinearCDSfold parses λ in single precision, which steps by 2⁻²⁴ just below
    1, so 1 − ε must be exactly representable there or the MFE corner silently
    runs at a different ε, or at ε = 0 with no tie-break at all. A power of two
    at or above 2⁻²³ satisfies that, and halving it stays representable."""
    eps = cfg["epsilon"]
    n_max = max(len(p) for _, p in instances(cfg))
    bound = epsilon_bounds(cfg["tables"]["canonical"], n_max)
    log(f"ε = {eps:g}; bound at n = {n_max}: MFE corner {bound['mfe']:.2e}, CAI corner {bound['cai']:.2e}")
    ok = eps < min(bound.values())
    for x in (eps, eps / 2):
        if not (single_exact(1.0 - x) and single_exact(x)):
            log(f"FAIL ε = {x!r}: 1 − ε is not exact in single precision")
            ok = False
    work = Path(cfg["out_dir"]) / "epscheck"
    work.mkdir(parents=True, exist_ok=True)
    rescore = Rescorer(cfg)
    cap = args.max_n or 300
    label, prot = max((i for i in instances(cfg) if len(i[1]) <= cap), key=lambda i: len(i[1]))
    for name, arm in cfg["arms"].items():
        if not arm.get("enabled", True) or not arm.get("exact") or "threads" in arm:
            continue
        for corner in arm.get("corners", cfg["corners"]):
            got = []
            for eps in (cfg["epsilon"], cfg["epsilon"] / 2):
                c = dict(variants(name, arm)[0], label=label, prot=prot, n=len(prot), corner=corner)
                row = run_once(dict(cfg, epsilon=eps), c, cfg["work"]["cores"][:1], work, f"eps{eps:g}")
                got.append((row["exit"], rescore(row["seq"], prot) if row.get("seq") else None))
            same = got[0][0] == "ok" and got[0] == got[1]
            ok &= same
            log(f"{'PASS' if same else 'FAIL'} {name} {corner} n={len(prot)}: {got[0][1]} vs {got[1][1]}")
    sys.exit(0 if ok else 1)


# ---- selftest, plan, summarize -------------------------------------------------------

def selftest(cfg, args, log):
    """Every enabled arm on one short protein: it must exit, parse and translate
    back, and its re-scored MFE must match what it reported (exact arms) or not
    exceed it (approximate arms)."""
    work = Path(cfg["out_dir"]) / "selftest"
    work.mkdir(parents=True, exist_ok=True)
    rescore = Rescorer(cfg)
    n = args.max_n or 30
    label, prot = min((i for i in instances(cfg) if len(i[1]) >= n), key=lambda i: len(i[1]))
    prot = prot[:n]
    bad = 0
    for name, arm in cfg["arms"].items():
        if not arm.get("enabled", True):
            continue
        v = variants(name, arm)[0]
        corner = arm.get("corners", cfg["corners"])[0]
        c = dict(v, label=f"{label}_selftest", prot=prot, n=n, corner=corner, ncores=1, threads=1)
        row = run_once(cfg, c, cfg["work"]["cores"][:1], work, "selftest")
        if row.get("seq"):
            row.update(rescore(row["seq"], prot))
        good = row["exit"] == "ok" and row.get("translates") and row.get("mfe_rescored") is not None
        if good:
            # An exact arm's reported MFE is the molecule's. A beam arm reports
            # the fold its beam found, which the molecule's true MFE can only
            # undercut, so only the re-scored energy compares across arms.
            gap = row["mfe_rescored"] - row["mfe_reported"]
            good = abs(gap) < 0.051 if arm.get("exact") else gap < 0.051
        bad += not good
        log(f"{'PASS' if good else 'FAIL'} {name}: exit={row['exit']} translates={row.get('translates')} "
            f"mfe reported={row.get('mfe_reported')} rescored={row.get('mfe_rescored')} {row.get('note', '')}")
    sys.exit(1 if bad else 0)


def plan(cfg, args, log):
    cells = grid(cfg, args)
    by_arm = {}
    for c in cells:
        by_arm[(c["arm"], c["variant"])] = by_arm.get((c["arm"], c["variant"]), 0) + 1
    for (a, v), k in sorted(by_arm.items()):
        log(f"{a:24s} {v:6s} {k:5d} cells")
    log(f"work pass rows: {len(cells)}; time pass rows: at most {len(cells) * cfg['time']['reps']}")


def summarize(cfg, args, log):
    """Row counts against the grid, exit categories, contamination, and the
    re-scored optimum of every exact arm against the reference arm."""
    out = Path(cfg["out_dir"])
    cells = {cell_id(c) for c in grid(cfg, args)}
    # A deferred cell that later ran has two rows; the latest is the cell's state.
    latest = {r["cell"]: r for r in Results(out / "work.jsonl").rows() if r["cell"] in cells}
    work = list(latest.values())
    by_key = time_rows_by_key(r for r in Results(out / "time.jsonl").rows() if r["cell"] in cells)
    times = [best_time_row(v) for v in by_key.values()]
    log(f"work rows {len(work)} / {len(cells)} planned")
    for status in sorted({r["exit"] for r in work}):
        log(f"  {status}: {sum(r['exit'] == status for r in work)}")
    ok_cells = sum(r["exit"] in TIMEABLE for r in work)
    log(f"time rows {len(times)} / {ok_cells * cfg['time']['reps']} expected, "
        f"still contaminated {sum(r.get('contaminated', False) for r in times)}, "
        f"runs repeated after contamination {sum(len(v) > 1 for v in by_key.values())}")
    for r in slow_per_instruction(times, cfg["time"]):
        log(f"  slow per instruction: {r['cell']} rep {r['rep']} {r['ratio']:.2f}x its neighbours' median")
    ref_arm = cfg.get("reference_arm")
    ref = {(r["protein"], r["corner"]): r for r in work if r["arm"] == ref_arm and r["exit"] == "ok"}
    for r in work:
        if r["exit"] != "ok" or r["arm"] == ref_arm or not cfg["arms"].get(r["arm"], {}).get("exact"):
            continue
        base = ref.get((r["protein"], r["corner"]))
        if base and r.get("mfe_rescored") != base.get("mfe_rescored") and r["corner"] == "mfe":
            log(f"  MFE differs from {ref_arm}: {r['cell']} {r.get('mfe_rescored')} vs {base.get('mfe_rescored')}")
        if not r.get("translates"):
            log(f"  does not translate back: {r['cell']}")
        # Every amino acid has a codon of weight 1, so the CAI corner's optimum
        # is CAI = 1 exactly; an exact arm short of it traded CAI for energy.
        if r["corner"] == "cai" and abs((r.get("cai_rescored") or 0) - 1.0) > 1e-9:
            log(f"  CAI corner below 1: {r['cell']} CAI {r.get('cai_rescored')}")


def slow_per_instruction(times, spec):
    """Timed runs that took longer per instruction than the same arm variant at
    nearby lengths. Instructions retired do not move with load, so time per
    instruction should change smoothly with n; a run well above its neighbours
    was slowed by something the core-busy checks cannot see, such as a
    neighbour's memory traffic on the same socket. Runs under `cpu_share_after_s`
    are left out: start-up dominates them."""
    ok = [r for r in times if r["exit"] == "ok" and r.get("instructions") and r["wall_s"] >= spec["cpu_share_after_s"]]
    by_variant = {}
    for r in ok:
        by_variant.setdefault((r["arm"], r["variant"]), []).append(r)
    span, limit = spec.get("tpi_window", 1.25), spec.get("tpi_max", 1.2)
    flagged = []
    for rows in by_variant.values():
        for r in rows:
            near = sorted(x["wall_s"] / x["instructions"] for x in rows
                          if x is not r and r["n_aa"] / span <= x["n_aa"] <= r["n_aa"] * span)
            if len(near) >= 3:
                ratio = (r["wall_s"] / r["instructions"]) / near[len(near) // 2]
                if ratio > limit:
                    flagged.append(dict(r, ratio=ratio))
    return sorted(flagged, key=lambda r: -r["ratio"])


def resolve_paths(cfg):
    """Tables and instance paths may be given relative to the repository root.
    Binaries' own working files (LinearDesign's table) stay as written."""
    def at_repo(p):
        return str(p if Path(p).is_absolute() else REPO / p)
    cfg["tables"] = {k: at_repo(v) for k, v in cfg["tables"].items()}
    spec = cfg["instances"]
    for d in spec.get("protein_dirs", []):
        d["dir"] = at_repo(d["dir"])
    spec["protein_files"] = [at_repo(p) for p in spec.get("protein_files", [])]
    if spec.get("ladder"):
        spec["ladder"]["fasta"] = at_repo(spec["ladder"]["fasta"])
    return cfg


# ---- entry -----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cmd", choices=["plan", "selftest", "epscheck", "work", "time", "summarize"])
    ap.add_argument("--config", required=True)
    ap.add_argument("--arms")
    ap.add_argument("--corners")
    ap.add_argument("--max-n", type=int)
    ap.add_argument("--drop-seq", action="store_true", help="do not keep designed sequences in the results")
    ap.add_argument("--redo-contaminated", action="store_true", help="time pass: re-run flagged reps")
    args = ap.parse_args()
    cfg = resolve_paths(json.loads(Path(args.config).read_text()))
    out = Path(cfg["out_dir"])
    for sub in ("work", "time"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    def log(msg):
        print(f"{time.strftime('%F %T')} {msg}", flush=True)

    if args.cmd in ("work", "time"):
        # One driver per pass per output directory: two runs of one grid write
        # the same work files and the same results rows.
        lock = open(out / f".{args.cmd}.lock", "w")
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            sys.exit(f"another {args.cmd} pass holds {out}/.{args.cmd}.lock")
        (out / f"{args.cmd}.config.json").write_text(json.dumps(cfg, indent=1))
    {"plan": plan, "selftest": selftest, "epscheck": epscheck, "work": run_work, "time": run_time, "summarize": summarize}[args.cmd](cfg, args, log)


if __name__ == "__main__":
    main()
