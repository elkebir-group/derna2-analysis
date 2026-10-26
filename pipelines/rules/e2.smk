# E2: what the whole Pareto front costs each method. The grid is drivers/e2/e2.config.json;
# the design is drivers/e2/E2-PLAN.md.

E2 = config["e2"]
E2_RUNS = RUNS / "e2"
E2_DRIVE = DRIVERS / "e2" / "drive_e2.py"


rule e2_figures:
    input:
        FIGURES / "e2_fronts.pdf",
        DERIVED / "e2" / "fig2.json",


rule e2_runs:
    input:
        E2_RUNS / "summary.txt",


# ---- figures ------------------------------------------------------------------------

# The E2 grid is still running, so its deposit has no checksum manifest yet; the figure
# reads only rows stamped with this pin. Drawn from this checkout's runs, it waits for the
# run's summary.
rule e2_plot_fronts:
    input:
        sweep=source("e2") / "sweep.jsonl",
        ok=lambda w: [] if config["figures_from"] == "deposit" else E2_RUNS / "summary.txt",
    output:
        pdf=FIGURES / "e2_fronts.pdf",
        png=FIGURES / "e2_fronts.png",
        cells=DERIVED / "e2" / "cells.tsv",
    params:
        pin=E2["pin"],
    shell:
        "{PY} {SCRIPTS}/e2/plot_fronts.py {input.sweep} {params.pin} {output.pdf} {output.cells}"


# Figure 2's data: per arm on the cells every arm completed, the cost quartiles; per front,
# the share of its points no weight reaches, from hull vertices computed on the points.
rule e2_fig2_data:
    input:
        sweep=source("e2") / "sweep.jsonl",
        ok=lambda w: [] if config["figures_from"] == "deposit" else E2_RUNS / "summary.txt",
    output:
        DERIVED / "e2" / "fig2.json",
    params:
        pin=E2["pin"],
    shell:
        "{PY} {SCRIPTS}/e2/fig2_data.py {input.sweep} {params.pin} {output}"


# ---- runs -----------------------------------------------------------------------------

rule e2_config:
    input:
        template=DRIVERS / "e2" / "e2.config.json",
        bins=STAMPS / "bins_e2.ok",
    output:
        E2_RUNS / "config.json",
    params:
        bin=Path(E2["bin_dir"]),
        lineardesign_dir=Path(E2["lineardesign_dir"]),
        lcds_dir=Path(E2["lcds_dir"]),
        pin=E2["pin"],
    shell:
        "{PY} {SCRIPTS}/render_driver_config.py {input.template} {output} "
        "bin={params.bin} lineardesign_dir={params.lineardesign_dir} "
        "lcds_dir={params.lcds_dir} pin={params.pin} "
        "inputs={INPUTS} out={E2_RUNS} root={ROOT}"


# Each comparator's λ map checked against derna2 at interior λs, before the grid.
rule e2_verify:
    input:
        E2_RUNS / "config.json",
    output:
        E2_RUNS / "verify.done",
    log:
        E2_RUNS / "verify.out",
    shell:
        "{DRIVER_PY} {E2_DRIVE} verify --config {input} > {log} 2>&1 && touch {output}"


# The grid's four passes in the order they depend on each other (single-core walk,
# single-core replays, eight-core walk, eight-core replays), refused while the timing cores
# or their hyperthread siblings are busy. Resumable: a re-run skips the cells this pin
# already recorded. sweep.jsonl is not a declared output, so a failed pass never makes
# Snakemake delete the rows the passes before it paid for.
rule e2_grid:
    input:
        config=E2_RUNS / "config.json",
        verify=E2_RUNS / "verify.done",
    output:
        E2_RUNS / "grid.done",
    log:
        E2_RUNS / "grid.out",
    shell:
        "bash {DRIVERS}/e2/run_e2_grid.sh {input.config} > {log} 2>&1 && touch {output}"


rule e2_summarize:
    input:
        config=E2_RUNS / "config.json",
        grid=E2_RUNS / "grid.done",
    output:
        E2_RUNS / "summary.txt",
    params:
        pin=E2["pin"],
    shell:
        "{DRIVER_PY} {E2_DRIVE} summarize --config {input.config} --pin {params.pin} "
        "> {output} 2>&1"
