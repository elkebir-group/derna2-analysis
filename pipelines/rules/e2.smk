# E2: what the whole Pareto front costs each method. The grid is
# drivers/benchmark/e2.config.json; the design is drivers/benchmark/E2-PLAN.md.

E2 = config["e2"]
E2_RUNS = RUNS / "e2"
E2_DRIVE = DRIVERS / "benchmark" / "drive_e2.py"


rule e2_figures:
    input:
        FIGURES / "e2_fronts.pdf",


rule e2_runs:
    input:
        E2_RUNS / "summary.txt",


# ---- figures ------------------------------------------------------------------------

# The E2 grid is still running, so its deposit has no checksum manifest yet; the figure
# reads only rows stamped with this pin.
rule e2_plot_fronts:
    input:
        source("e2") / "sweep.jsonl",
    output:
        pdf=FIGURES / "e2_fronts.pdf",
        png=FIGURES / "e2_fronts.png",
        cells=DERIVED / "e2" / "cells.tsv",
    params:
        pin=E2["pin"],
    shell:
        "{PY} {SCRIPTS}/e2/plot_fronts.py {input} {params.pin} {output.pdf} {output.cells}"


# ---- runs -----------------------------------------------------------------------------

rule e2_config:
    input:
        template=DRIVERS / "benchmark" / "e2.config.json",
        bins=STAMPS / "bins_e2.ok",
        comparators=STAMPS / "bins_e1.ok",
    output:
        E2_RUNS / "config.json",
    params:
        bin=Path(E2["bin_dir"]),
        bin_e1=Path(E2["comparator_bin_dir"]),
        derna1_params=Path(E2["derna1_params"]),
    shell:
        "{PY} {SCRIPTS}/render_driver_config.py {input.template} {output} "
        "bin={params.bin} bin_e1={params.bin_e1} derna1_params={params.derna1_params} "
        "inputs={INPUTS} out={E2_RUNS} root={ROOT}"


# derna2's chord sweep per instance: hull vertices and the replay λs every arm runs.
rule e2_hull:
    input:
        E2_RUNS / "config.json",
    output:
        E2_RUNS / "hull.done",
    log:
        E2_RUNS / "hull.out",
    shell:
        "{DRIVER_PY} {E2_DRIVE} hull --config {input} > {log} 2>&1 && touch {output}"


# Each comparator's λ map checked against derna2 at interior λs, before the grid.
rule e2_verify:
    input:
        config=E2_RUNS / "config.json",
        hull=E2_RUNS / "hull.done",
    output:
        E2_RUNS / "verify.done",
    log:
        E2_RUNS / "verify.out",
    shell:
        "{DRIVER_PY} {E2_DRIVE} verify --config {input.config} > {log} 2>&1 && touch {output}"


rule e2_sweep:
    input:
        config=E2_RUNS / "config.json",
        verify=E2_RUNS / "verify.done",
    output:
        E2_RUNS / "sweep.jsonl",
    log:
        E2_RUNS / "sweep.out",
    params:
        cpus=E2["driver_cpus"],
    shell:
        "taskset -c {params.cpus} {DRIVER_PY} {E2_DRIVE} sweep --config {input.config} "
        "> {log} 2>&1"


rule e2_summarize:
    input:
        config=E2_RUNS / "config.json",
        sweep=E2_RUNS / "sweep.jsonl",
    output:
        E2_RUNS / "summary.txt",
    params:
        pin=E2["pin"],
    shell:
        "{DRIVER_PY} {E2_DRIVE} summarize --config {input.config} --pin {params.pin} "
        "> {output} 2>&1"
