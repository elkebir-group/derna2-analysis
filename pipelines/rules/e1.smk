# E1: single-λ correctness, runtime and memory. The grid is drivers/benchmark/e1.config.json;
# drivers/benchmark/BENCHMARK.md explains the two passes and the socket isolation.

E1 = config["e1"]
E1_RUNS = RUNS / "e1"
E1_DRIVE = DRIVERS / "benchmark" / "drive.py"


rule e1_figures:
    input:
        DERIVED / "e1" / "correctness.tsv",
        FIGURES / "e1_scaling.pdf",
        FIGURES / "e1_speed_quality.pdf",


rule e1_runs:
    input:
        E1_RUNS / "summary.txt",


# ---- figures ------------------------------------------------------------------------

rule e1_tables:
    input:
        ok=source_ok("e1"),
        work=source("e1") / "work.jsonl",
        time=source("e1") / "time.jsonl",
    output:
        cells=DERIVED / "e1" / "cells.tsv",
        correctness=DERIVED / "e1" / "correctness.tsv",
    params:
        rows=E1["rows"]["work"],
    shell:
        "{PY} {SCRIPTS}/e1/tables.py {input.work} {input.time} {output.cells} "
        "{output.correctness} {params.rows}"


rule e1_plot_scaling:
    input:
        DERIVED / "e1" / "cells.tsv",
    output:
        FIGURES / "e1_scaling.pdf",
        FIGURES / "e1_scaling.png",
    shell:
        "{PY} {SCRIPTS}/e1/plot_scaling.py {input} {output[0]}"


rule e1_plot_speed_quality:
    input:
        DERIVED / "e1" / "cells.tsv",
    output:
        FIGURES / "e1_speed_quality.pdf",
        FIGURES / "e1_speed_quality.png",
    shell:
        "{PY} {SCRIPTS}/e1/plot_speed_quality.py {input} {output[0]}"


# ---- runs -----------------------------------------------------------------------------

rule e1_config:
    input:
        template=DRIVERS / "benchmark" / "e1.config.json",
        bins=STAMPS / "bins_e1.ok",
    output:
        E1_RUNS / "config.json",
    params:
        bin=Path(E1["bin_dir"]),
    shell:
        "{PY} {SCRIPTS}/render_driver_config.py {input.template} {output} "
        "bin={params.bin} inputs={INPUTS} out={E1_RUNS}"


# Every enabled arm on a 30-residue protein: exits, parses, translates back, re-scores.
rule e1_selftest:
    input:
        E1_RUNS / "config.json",
    output:
        E1_RUNS / "selftest.log",
    shell:
        "{DRIVER_PY} {E1_DRIVE} selftest --config {input} > {output} 2>&1"


# ε against its bound at the longest instance, and every exact arm unchanged at ε/2.
rule e1_epscheck:
    input:
        config=E1_RUNS / "config.json",
        selftest=E1_RUNS / "selftest.log",
    output:
        E1_RUNS / "epscheck.log",
    shell:
        "{DRIVER_PY} {E1_DRIVE} epscheck --config {input.config} --max-n 300 > {output} 2>&1"


# Correctness, peak RSS, instructions: load-insensitive, runs on any cores. Resumable.
rule e1_work:
    input:
        config=E1_RUNS / "config.json",
        epscheck=E1_RUNS / "epscheck.log",
    output:
        E1_RUNS / "work.jsonl",
    log:
        E1_RUNS / "work.log",
    shell:
        "{DRIVER_PY} {E1_DRIVE} work --config {input.config} > {log} 2>&1"


# Wall-clock on the reserved socket, gated on it being idle. The driver runs on the other
# socket so its own threads never read as a neighbour. Resumable.
rule e1_time:
    input:
        config=E1_RUNS / "config.json",
        work=E1_RUNS / "work.jsonl",
    output:
        E1_RUNS / "time.jsonl",
    log:
        E1_RUNS / "time.log",
    params:
        cpus=E1["driver_cpus"],
    shell:
        "taskset -c {params.cpus} {DRIVER_PY} {E1_DRIVE} time --config {input.config} > {log} 2>&1"


# Row counts against the plan, exit categories, contamination, exact arms against the
# reference. Read the counts before any result.
rule e1_summarize:
    input:
        config=E1_RUNS / "config.json",
        time=E1_RUNS / "time.jsonl",
    output:
        E1_RUNS / "summary.txt",
    shell:
        "{DRIVER_PY} {E1_DRIVE} summarize --config {input.config} > {output} 2>&1"
