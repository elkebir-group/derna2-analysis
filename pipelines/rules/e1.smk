# E1: what the whole Pareto front costs each method. The grid is drivers/e1/e1.config.json.

E1 = config["e1"]
E1_RUNS = RUNS / "e1"
E1_DRIVE = DRIVERS / "e1" / "drive_e1.py"


rule e1_figures:
    input:
        FIGURES / "e1_fronts.pdf",
        FIGURES / "e1_overview.pdf",
        DERIVED / "e1" / "fig2.json",


rule e1_runs:
    input:
        E1_RUNS / "summary.txt",


# ---- figures ------------------------------------------------------------------------

# The E1 grid is still running, so its deposit has no checksum manifest yet; the figure
# reads only rows stamped with this pin. Drawn from this checkout's runs, it waits for the
# run's summary.
rule e1_plot_fronts:
    input:
        sweep=source("e1") / "sweep.jsonl",
        ok=lambda w: [] if config["figures_from"] == "deposit" else E1_RUNS / "summary.txt",
    output:
        pdf=FIGURES / "e1_fronts.pdf",
        png=FIGURES / "e1_fronts.png",
        cells=DERIVED / "e1" / "cells.tsv",
    params:
        pin=E1["pin"],
    shell:
        "{PY} {SCRIPTS}/e1/plot_fronts.py {input.sweep} {params.pin} {output.pdf} {output.cells}"


# The overview figure: the hull of each front against the whole front, and what the whole
# front costs. The pruned hull search's wall time is read from the walk reports in hull/.
rule e1_plot_overview:
    input:
        sweep=source("e1") / "sweep.jsonl",
        walks=source("e1") / "hull",
        ok=lambda w: [] if config["figures_from"] == "deposit" else E1_RUNS / "summary.txt",
    output:
        pdf=FIGURES / "e1_overview.pdf",
        png=FIGURES / "e1_overview.png",
    params:
        pin=E1["pin"],
    shell:
        "{PY} {SCRIPTS}/e1/plot_overview.py {input.sweep} {input.walks} {params.pin} {output.pdf}"


# Figure 2's data: per arm on the cells every arm completed, the cost quartiles; per front,
# the share of its points no weight reaches, from hull vertices computed on the points.
rule e1_fig2_data:
    input:
        sweep=source("e1") / "sweep.jsonl",
        ok=lambda w: [] if config["figures_from"] == "deposit" else E1_RUNS / "summary.txt",
    output:
        DERIVED / "e1" / "fig2.json",
    params:
        pin=E1["pin"],
    shell:
        "{PY} {SCRIPTS}/e1/fig2_data.py {input.sweep} {params.pin} {output}"


# ---- runs -----------------------------------------------------------------------------

rule e1_config:
    input:
        template=DRIVERS / "e1" / "e1.config.json",
        derna=DERNA,
        fold=MFE_FOLD,
    output:
        E1_RUNS / "config.json",
    params:
        bin=Path(E1["bin_dir"]),
        derna2_bin=lambda w: Path(DERNA).parent,
        lineardesign_dir=Path(E1["lineardesign_dir"]),
        lcds_dir=Path(E1["lcds_dir"]),
        pin=E1["pin"],
    shell:
        "{PY} {SCRIPTS}/render_driver_config.py {input.template} {output} "
        "bin={params.bin} derna2_bin={params.derna2_bin} "
        "lineardesign_dir={params.lineardesign_dir} "
        "lcds_dir={params.lcds_dir} pin={params.pin} "
        "inputs={INPUTS} out={E1_RUNS} root={ROOT}"


# Each comparator's λ map checked against derna2 at interior λs, before the grid.
rule e1_verify:
    input:
        E1_RUNS / "config.json",
    output:
        E1_RUNS / "verify.done",
    log:
        E1_RUNS / "verify.out",
    shell:
        "{DRIVER_PY} {E1_DRIVE} verify --config {input} > {log} 2>&1 && touch {output}"


# The grid's four passes in the order they depend on each other (single-core walk,
# single-core replays, eight-core walk, eight-core replays), refused while the timing cores
# or their hyperthread siblings are busy. Resumable: a re-run skips the cells this pin
# already recorded. sweep.jsonl is not a declared output, so a failed pass never makes
# Snakemake delete the rows the passes before it paid for.
rule e1_grid:
    input:
        config=E1_RUNS / "config.json",
        verify=E1_RUNS / "verify.done",
    output:
        E1_RUNS / "grid.done",
    log:
        E1_RUNS / "grid.out",
    shell:
        "bash {DRIVERS}/e1/run_e1_grid.sh {input.config} > {log} 2>&1 && touch {output}"


rule e1_summarize:
    input:
        config=E1_RUNS / "config.json",
        grid=E1_RUNS / "grid.done",
    output:
        E1_RUNS / "summary.txt",
    params:
        pin=E1["pin"],
    shell:
        "{DRIVER_PY} {E1_DRIVE} summarize --config {input.config} --pin {params.pin} "
        "> {output} 2>&1"


# ---- Figure 4 ----------------------------------------------------------------------------
# Figure 4 reads its own frozen deposit, `results/e1_fig4`: a snapshot of the E1 sweep taken
# before rows were added to it on 2026-10-30, and the walk's hull/ outputs it was drawn from
# (pins/e1_fig4-PROVENANCE.md). Every table below is rebuilt from those files, and
# e1_fig4_verify checks the tables and the figure against pins/e1_fig4.outputs.

F4 = config["e1_fig4"]
F4_DIR = DERIVED / "e1" / "fig4"


rule e1_fig4_figures:
    input:
        FIGURES / "fig4_e1_overview.pdf",
        STAMPS / "e1_fig4.verified",


# One cell row per (arm, protein) of the pin: the best attempt, not the last line.
rule e1_fig4_cells:
    input:
        sweep=RESULTS / "e1_fig4" / "sweep.jsonl",
        ok=STAMPS / "results_e1_fig4.ok",
    output:
        F4_DIR / "cells.csv",
    params:
        pin=F4["pin"],
    shell:
        "{PY} {SCRIPTS}/e1/fig4_cells.py {input.sweep} {params.pin} {output}"


# Each arm's exit and wall time, and the dichotomic hull search with and without pruning.
rule e1_fig4_per_protein:
    input:
        sweep=RESULTS / "e1_fig4" / "sweep.jsonl",
        ok=STAMPS / "results_e1_fig4.ok",
    output:
        F4_DIR / "per_protein.csv",
    params:
        hull=RESULTS / "e1_fig4" / "hull",
    shell:
        "{PY} {SCRIPTS}/e1/fig4_per_protein.py {input.sweep} {params.hull} {output}"


# Every front and its sweep hull, and the steps between neighbouring designs of each.
rule e1_fig4_steps:
    input:
        per=F4_DIR / "per_protein.csv",
        ok=STAMPS / "results_e1_fig4.ok",
    output:
        F4_DIR / "steps.npz",
    params:
        hull=RESULTS / "e1_fig4" / "hull",
    shell:
        "{PY} {SCRIPTS}/e1/fig4_steps.py {params.hull} {input.per} {output}"


# The figure. Fails, keeping the files, if any text covers another text, a mark or another
# panel. `focus` picks panel (a)'s protein.
rule e1_fig4:
    input:
        steps=F4_DIR / "steps.npz",
        per=F4_DIR / "per_protein.csv",
        cells=F4_DIR / "cells.csv",
    output:
        pdf=FIGURES / "fig4_e1_overview.pdf",
        png=FIGURES / "fig4_e1_overview.png",
    params:
        focus=F4["focus"],
    shell:
        "{PY} {SCRIPTS}/e1/fig4_overview.py {input.steps} {input.per} {input.cells} "
        "{output.pdf} {params.focus}"


# The tables and the figure against the digests they were published with. The tables match
# on any host with requirements.lock; the figure's bytes also need the TeX the manifest names.
rule e1_fig4_verify:
    input:
        manifest=ROOT / "pins" / "e1_fig4.outputs",
        built=[F4_DIR / "cells.csv", F4_DIR / "per_protein.csv", F4_DIR / "steps.npz",
               FIGURES / "fig4_e1_overview.pdf", FIGURES / "fig4_e1_overview.png"],
    output:
        STAMPS / "e1_fig4.verified",
    shell:
        "{PY} {SCRIPTS}/check_manifest.py {input.manifest} {DERIVED} {output}"
