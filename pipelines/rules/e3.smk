# E3: the price of motif-freedom, as five nested exact fronts per instance. The design is
# drivers/e3/E3-PLAN.md. The grid needs `derna front --motif`, which is not in a pinned
# build yet, so only what has run is here: the forced-motif census and the hull probe.

E3 = config["e3"]
E3_DRIVERS = DRIVERS / "e3"


E3_BANDS = ["grid", "mid", "long"]


rule e3_tables:
    input:
        DERIVED / "e3" / "forced_motifs.jsonl",
        DERIVED / "e3" / "probe.tsv",


rule e3_figures:
    input:
        DERIVED / "e3" / "fig3_panels.json",
        DERIVED / "e3" / "fig3_planned.json",
        DERIVED / "e3" / "fig3_strip.json",


# Which motif entries the codon table forces, on which peptides: recomputed by the vendored
# script from the committed sets and codon table, and checked against the plan's record.
rule e3_forced_motifs:
    input:
        driver=E3_DRIVERS / "forced_motifs.py",
        sets=expand(str(INPUTS / "motifs" / "{name}.txt"), name=E3["motif_sets"]),
        table=INPUTS / "codon_usage_human.csv",
        recorded=source("e3") / "e3-probe" / "forced.jsonl",
    output:
        DERIVED / "e3" / "forced_motifs.jsonl",
    params:
        motifs=INPUTS / "motifs",
    shell:
        "{PY} {SCRIPTS}/e3/forced_motifs.py {input.driver} {params.motifs} {input.table} "
        "{input.recorded} {output}"


# The hull probe (plan § 3a): one row per rung, the peak-RSS fit and its extrapolation to
# the long arm, and the wall's consecutive-rung exponents. The hull lower-bounds the front.
rule e3_probe_table:
    input:
        [source("e3") / d / "summary.csv" for d in E3["probe_dirs"]],
    output:
        tsv=DERIVED / "e3" / "probe.tsv",
        fit=DERIVED / "e3" / "probe_fit.txt",
    params:
        dirs=[source("e3") / d for d in E3["probe_dirs"]],
    shell:
        "{PY} {SCRIPTS}/e3/probe_table.py {output.tsv} {output.fit} {params.dirs}"


# ---- Figure 3 ---------------------------------------------------------------------------

# Every cell of the three bands (`grid` is the short band), from their summaries, with the
# rows the driver wrote with an unquoted core mask rejoined.
rule e3_fig3_cells:
    input:
        [source("e3") / b / "summary.csv" for b in E3_BANDS],
    output:
        DERIVED / "e3" / "fig3_cells.json",
    params:
        root=source("e3"),
    shell:
        "{PY} {SCRIPTS}/e3/fig3_cells.py {params.root} {output}"


# Coverage per protein and nest, peak RSS fitted per build, points found relative to F1.
rule e3_fig3_panels:
    input:
        DERIVED / "e3" / "fig3_cells.json",
    output:
        DERIVED / "e3" / "fig3_panels.json",
    shell:
        "{PY} {SCRIPTS}/e3/fig3_panels.py {input} {output}"


# Per short-band cell: ΔMFE at matched CAI against the protein's unconstrained front,
# whether its sets bind, and its cost against that front.
rule e3_fig3_planned:
    input:
        source("e3") / "screen" / "binding.jsonl",
        source("e3") / "grid" / "summary.csv",
    output:
        DERIVED / "e3" / "fig3_planned.json",
    params:
        root=source("e3"),
    shell:
        "{PY} {SCRIPTS}/e3/fig3_planned.py {params.root} {output}"


# The constraint strip: per blacklist entry, proteins admitting it against proteins whose
# front spells it; one instance's nested fronts and their prices; every cell's coverage.
rule e3_fig3_strip:
    input:
        [source("e3") / "screen" / "binding.jsonl"]
        + [source("e3") / b / "summary.csv" for b in E3_BANDS],
    output:
        DERIVED / "e3" / "fig3_strip.json",
    params:
        root=source("e3"),
    shell:
        "{PY} {SCRIPTS}/e3/fig3_strip.py {params.root} {output}"
