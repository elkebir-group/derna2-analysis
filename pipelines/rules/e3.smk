# E3: the price of motif-freedom, as five nested exact fronts per instance. The design is
# drivers/e3/E3-PLAN.md. The grid needs `derna front --motif`, which is not in a pinned
# build yet, so only what has run is here: the forced-motif census and the hull probe.

E3 = config["e3"]
E3_DRIVERS = DRIVERS / "e3"


rule e3_tables:
    input:
        DERIVED / "e3" / "forced_motifs.jsonl",
        DERIVED / "e3" / "probe.tsv",


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
