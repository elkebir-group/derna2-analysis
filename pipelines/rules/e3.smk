# E3: what a UTR-aware design practice delivers, per practical UTR pair (30 proteins on five
# UTR pairs, 150 cells). The deposit is the figure bundle's paper-e4/, against
# pins/results-e3.sha256. Two figures are rules: the main figure (e3_fig_story) and a supplement
# figure on where the CDS-UTR contacts fall (e3_fig_supp_contacts).
#
# Two kinds of target, as everywhere: e3_figures reads the deposit and the refolds below, and
# e3_refolds regenerates the refolds with `derna check`, which needs the binary of config
# the built `derna`. The refolds carry the positions of each design's
# CDS-UTR pairs, which the census rows do not: they hold the gap and nothing about where the
# pairs fall. They are deposited for that reason, so a host that has the deposit can draw both
# figures without the binary; e3_refold_deposit takes them from it, and `--force-run e3_refold`
# recomputes one instead.

E3 = config["e3"]
E3_CELLS = [line.split()[0] for line in open(INPUTS / "e3" / "cells.txt") if line.strip()]
E3_REFOLD = DERIVED / "e3" / "hullcheck"
E3_CAND = DERIVED / "e3" / "cand"
E3_STORY = E3["story"]
E3_SCRIPTS = SCRIPTS / "e3"
# The prefixes the deposited rows recorded, which the readers rewrite to this host's trees.
E3_ENV = " ".join(
    f'E3_RECORDED_{name.upper()}="{prefix}"'
    for name, prefix in {"e2": "e2_fronts", "deposit": "deposit"}.items()
    for prefix in [E3["recorded_prefixes"][prefix]]
)


rule e3_figures:
    input:
        FIGURES / "fig_e3_story.pdf",
        FIGURES / "fig_supp_contacts.pdf",
        DERIVED / "e3" / "story_stats.txt",


rule e3_refolds:
    input:
        expand(str(E3_REFOLD / "{cell}.out.jsonl"), cell=E3_CELLS),


# The deposited refold of one cell, which is what the figures read on a host that has the
# deposit. e3_refold recomputes the same file from the binary; the order below prefers the
# deposit, so regenerating is an explicit `--force-run e3_refold`.
ruleorder: e3_refold_deposit > e3_refold


rule e3_refold_deposit:
    input:
        ok=source_ok("e3"),
        out=lambda w: source("e3") / "hullcheck" / f"{w.cell}.out.jsonl",
    output:
        E3_REFOLD / "{cell}.out.jsonl",
    wildcard_constraints:
        cell="[A-Za-z0-9]+__[a-z0-9]+",
    shell:
        "cp {input.out} {output}"


# One cell's concatenation hull vertices as a derna result file.
rule e3_hull_input:
    input:
        ok=source_ok("e3"),
        seed=lambda w: next(q for q in (source("e3") / "seeds" / f"{w.cell}.check.{n}.jsonl"
                                        for n in ("freestop.ec083208", "pooled-fixed")) if q.exists()),
        a0=lambda w: next(q for q in (source("e3") / t / "cells" / w.cell / "a0.jsonl"
                                      for t in ("a0", "a0-all")) if q.exists()),
    output:
        E3_REFOLD / "{cell}.in.jsonl",
    wildcard_constraints:
        cell="[A-Za-z0-9]+__[a-z0-9]+",
    params:
        bundle=source("e3"),
        e2=source("e2"),
    shell:
        "{E3_ENV} {PY} {E3_SCRIPTS}/hull_input.py {params.bundle} {params.e2} {wildcards.cell} {output}"


# `derna check` of those vertices: the structure of each, with the positions of its CDS-UTR
# pairs. Written under a temporary name so that a stopped refold leaves nothing behind.
rule e3_refold:
    input:
        E3_REFOLD / "{cell}.in.jsonl",
    output:
        E3_REFOLD / "{cell}.out.jsonl",
    wildcard_constraints:
        cell="[A-Za-z0-9]+__[a-z0-9]+",
    log:
        E3_REFOLD / "{cell}.log",
    threads: 2
    params:
        derna=DERNA,
    shell:
        "{params.derna} check {input} --threads {threads} --out {output}.tmp 2> {log} && mv {output}.tmp {output}"


# ---- The main figure -----------------------------------------------------------------------
# Besides the refolds of every cell's hull vertices, panel (e) draws the best design at the CAI
# floor of the instance for concatenation and for mRNAfold, refolded on the whole mRNA.

rule e3_cand_input:
    input:
        ok=source_ok("e3"),
    output:
        conc=E3_CAND / (E3_STORY["cell"] + ".conc.in.jsonl"),
        mrf=E3_CAND / (E3_STORY["cell"] + ".mrf.in.jsonl"),
    params:
        bundle=source("e3"),
        e2=source("e2"),
        cell=E3_STORY["cell"],
        floor=E3_STORY["cai_floor"],
        prefix=E3_CAND / E3_STORY["cell"],
    shell:
        "{E3_ENV} {PY} {E3_SCRIPTS}/cand_input.py {params.bundle} {params.e2} {params.cell} {params.floor} {params.prefix}"


rule e3_cand_refold:
    input:
        E3_CAND / "{cell}.{kind}.in.jsonl",
    output:
        E3_CAND / "{cell}.{kind}.out.jsonl",
    wildcard_constraints:
        kind="conc|mrf",
    log:
        E3_CAND / "{cell}.{kind}.log",
    threads: 2
    params:
        derna=DERNA,
    shell:
        "{params.derna} check {input} --threads {threads} --out {output}.tmp 2> {log} && mv {output}.tmp {output}"


rule e3_fig_story:
    input:
        ok=source_ok("e3"),
        refolds=expand(str(E3_REFOLD / "{cell}.out.jsonl"), cell=E3_CELLS),
        cand=expand(str(E3_CAND / (E3_STORY["cell"] + ".{kind}.out.jsonl")), kind=["conc", "mrf"]),
    output:
        FIGURES / "fig_e3_story.pdf",
    params:
        bundle=source("e3"),
        e2=source("e2"),
        refolds=E3_REFOLD,
        cand=E3_CAND,
        cell=E3_STORY["cell"],
        floor=E3_STORY["cai_floor"],
    shell:
        "{E3_ENV} {RNA_PY} {E3_SCRIPTS}/fig_story.py {params.bundle} {params.e2} {params.refolds} {params.cand} "
        "{params.cell} {params.floor} {output}"


# The numbers the paper's UTRs subsection quotes, which the figure's panels draw or sit beside.
rule e3_story_stats:
    input:
        ok=source_ok("e3"),
        refolds=expand(str(E3_REFOLD / "{cell}.out.jsonl"), cell=E3_CELLS),
        cand=expand(str(E3_CAND / (E3_STORY["cell"] + ".{kind}.out.jsonl")), kind=["conc", "mrf"]),
    output:
        DERIVED / "e3" / "story_stats.txt",
    params:
        bundle=source("e3"),
        e2=source("e2"),
        refolds=E3_REFOLD,
        cand=E3_CAND,
        cell=E3_STORY["cell"],
        floor=E3_STORY["cai_floor"],
    shell:
        "{E3_ENV} {PY} {E3_SCRIPTS}/story_stats.py {params.bundle} {params.e2} {params.refolds} {params.cand} "
        "{params.cell} {params.floor} {output}"


# The numbers against the digest they were published with.
rule e3_story_stats_verify:
    input:
        manifest=ROOT / "pins" / "e3_story_stats.outputs",
        built=[DERIVED / "e3" / "story_stats.txt"],
    output:
        STAMPS / "e3_story_stats.verified",
    shell:
        "{PY} {SCRIPTS}/check_manifest.py {input.manifest} {DERIVED} {output}"


# The main figure against the digest it was published with.
rule e3_fig_story_verify:
    input:
        manifest=ROOT / "pins" / "e3_story.outputs",
        built=[FIGURES / "fig_e3_story.pdf"],
    output:
        STAMPS / "e3_story.verified",
    shell:
        "{PY} {SCRIPTS}/check_manifest.py {input.manifest} {DERIVED} {output}"


# ---- The supplement figure ------------------------------------------------------------------

# The histograms the figure is drawn from, so a layout is redrawn in a second.
rule e3_contacts_data:
    input:
        cells=INPUTS / "e3" / "cells.txt",
        ok=source_ok("e3"),
        refolds=expand(str(E3_REFOLD / "{cell}.out.jsonl"), cell=E3_CELLS),
    output:
        DERIVED / "e3" / "supp_contacts.json",
    params:
        bundle=source("e3"),
        e2=source("e2"),
        refolds=E3_REFOLD,
    shell:
        "{E3_ENV} {PY} {E3_SCRIPTS}/contacts_data.py {params.bundle} {params.e2} {params.refolds} {input.cells} {output}"


rule e3_fig_supp_contacts:
    input:
        DERIVED / "e3" / "supp_contacts.json",
    output:
        FIGURES / "fig_supp_contacts.pdf",
    shell:
        "{PY} {E3_SCRIPTS}/fig_supp_contacts.py {input} {output}"


# The histograms against the digest they were published with.
rule e3_supp_contacts_verify:
    input:
        manifest=ROOT / "pins" / "e3_supp_contacts.outputs",
        built=[DERIVED / "e3" / "supp_contacts.json"],
    output:
        STAMPS / "e3_supp_contacts.verified",
    shell:
        "{PY} {SCRIPTS}/check_manifest.py {input.manifest} {DERIVED} {output}"
