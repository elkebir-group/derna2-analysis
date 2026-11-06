# E2: the price of motif-freedom, as five nested exact fronts per instance. The grid needs `derna front --motif`, which is not in a pinned
# build yet, so only what has run is here: the forced-motif census and the hull probe.

E2 = config["e2"]
E2_DRIVERS = DRIVERS / "e2"


E2_BANDS = ["grid", "mid", "long"]


rule e2_tables:
    input:
        DERIVED / "e2" / "forced_motifs.jsonl",
        DERIVED / "e2" / "probe.tsv",


rule e2_figures:
    input:
        DERIVED / "e2" / "fig3_panels.json",
        DERIVED / "e2" / "fig3_planned.json",
        DERIVED / "e2" / "fig3_strip.json",


# Which motif entries the codon table forces, on which peptides: recomputed by the vendored
# script from the committed sets and codon table, and checked against the plan's record.
rule e2_forced_motifs:
    input:
        driver=E2_DRIVERS / "forced_motifs.py",
        sets=expand(str(INPUTS / "motifs" / "{name}.txt"), name=E2["motif_sets"]),
        table=INPUTS / "codon_usage_human.csv",
        recorded=source("e2") / "e3-probe" / "forced.jsonl",
        ok=source_ok("e2"),
    output:
        DERIVED / "e2" / "forced_motifs.jsonl",
    params:
        motifs=INPUTS / "motifs",
    shell:
        "{PY} {SCRIPTS}/e2/forced_motifs.py {input.driver} {params.motifs} {input.table} "
        "{input.recorded} {output}"


# The hull probe (plan § 3a): one row per rung, the peak-RSS fit and its extrapolation to
# the long arm, and the wall's consecutive-rung exponents. The hull lower-bounds the front.
rule e2_probe_table:
    input:
        [source("e2") / d / "summary.csv" for d in E2["probe_dirs"]],
        ok=source_ok("e2"),
    output:
        tsv=DERIVED / "e2" / "probe.tsv",
        fit=DERIVED / "e2" / "probe_fit.txt",
    params:
        dirs=[source("e2") / d for d in E2["probe_dirs"]],
    shell:
        "{PY} {SCRIPTS}/e2/probe_table.py {output.tsv} {output.fit} {params.dirs}"


# ---- Figure 3 ---------------------------------------------------------------------------

# Every cell of the three bands (`grid` is the short band), from their summaries, with the
# rows the driver wrote with an unquoted core mask rejoined.
rule e2_fig3_cells:
    input:
        [source("e2") / b / "summary.csv" for b in E2_BANDS],
        ok=source_ok("e2"),
    output:
        DERIVED / "e2" / "fig3_cells.json",
    params:
        root=source("e2"),
    shell:
        "{PY} {SCRIPTS}/e2/fig3_cells.py {params.root} {output}"


# Coverage per protein and nest, peak RSS fitted per build, points found relative to F1.
rule e2_fig3_panels:
    input:
        DERIVED / "e2" / "fig3_cells.json",
    output:
        DERIVED / "e2" / "fig3_panels.json",
    shell:
        "{PY} {SCRIPTS}/e2/fig3_panels.py {input} {output}"


# Per short-band cell: ΔMFE at matched CAI against the protein's unconstrained front,
# whether its sets bind, and its cost against that front.
rule e2_fig3_planned:
    input:
        source("e2") / "screen" / "binding.jsonl",
        source("e2") / "grid" / "summary.csv",
        ok=source_ok("e2"),
    output:
        DERIVED / "e2" / "fig3_planned.json",
    params:
        root=source("e2"),
    shell:
        "{PY} {SCRIPTS}/e2/fig3_planned.py {params.root} {output}"


# The constraint strip: per blacklist entry, proteins admitting it against proteins whose
# front spells it; one instance's nested fronts and their prices; every cell's coverage.
rule e2_fig3_strip:
    input:
        [source("e2") / "screen" / "binding.jsonl"]
        + [source("e2") / b / "summary.csv" for b in E2_BANDS],
        ok=source_ok("e2"),
    output:
        DERIVED / "e2" / "fig3_strip.json",
    params:
        root=source("e2"),
    shell:
        "{PY} {SCRIPTS}/e2/fig3_strip.py {params.root} {output}"


# ---- The composite figure: the price of motif-freedom, in seven panels -------------------

# The cached measurements the plate is drawn from. Four steps write them and the drawing
# step reads them back, so a layout can be redrawn in seconds without touching the deposit.
E2_CACHE = DERIVED / "e2" / "composite"
# The example the fronts and the four structures are taken from: protein, CAI floor, and
# which alternative design (`none` = the rung optima themselves).
E2_EXAMPLE = E2["example"]
E2_SLUG = f"{E2_EXAMPLE['cell']}-c{round(E2_EXAMPLE['floor'] * 100):02d}-{E2_EXAMPLE['alt']}"
# Where each step reads and writes. Paths go in the environment rather than on the command
# line because two of these steps exec a prefix of another, which argv cannot reach; the
# files themselves are declared as inputs and outputs above the shell line, as usual.
#
# E2_E1_SWEEP and E2_E1_EXTRA are panel (a)'s only reach outside E2's deposit: the screen
# is scored over the proteins E1 ran, and its motif census over an auxiliary tree of E1
# instances screened for this figure. Both are deposited and checksummed since 2026-11-05
# (pins/results-e1.sha256, pins/results-e1_extra.sha256), so they resolve through `source`
# like every other deposit and need no host override. The rules reading them declare both
# stamps.
E2_ENV = (
    f'E2_DEPOSIT="{source("e2")}" E2_CACHE="{E2_CACHE}" '
    f'E2_E1_SWEEP="{source("e1") / "sweep.jsonl"}" '
    f'E2_E1_EXTRA="{source("e1_extra")}" '
    f'E2_MOTIF_NAMES="{INPUTS / "motifs" / "names.json"}" '
    f'E2_CORPUS="{INPUTS / "e2" / "corpus.json"}"'
)


# The unit-square AUC per protein and rung over the finished set, which is also the list of
# proteins every later step is computed over.
rule e2_composite_auc:
    input:
        finished=source("e2") / "FINISHED.tsv",
        ok=source_ok("e2"),
    output:
        E2_CACHE / "auc01_fin.json",
    shell:
        "{E2_ENV} {PY} {SCRIPTS}/e2/auc01_fin.py"


# Panels (a) and (b): the motif screen over the catalog, and what each rung costs against
# the unconstrained design at four CAI floors, per protein.
rule e2_composite_data:
    input:
        # Asserted against, not read from: the population this step derives from E1's sweep
        # is the one the supplement's panel (d) counts over, and the two must not drift.
        corpus=INPUTS / "e2" / "corpus.json",
        auc=E2_CACHE / "auc01_fin.json",
        sets=expand(str(INPUTS / "motifs" / "{name}.txt"), name=E2["motif_sets"]),
        names=INPUTS / "motifs" / "names.json",
        # This step reads the deposit through E2_ENV rather than through a declared path,
        # so the stamp is what makes its checksums a precondition of the figure.
        ok=source_ok("e2"),
        ok_e1=source_ok("e1"),
        ok_e1_extra=source_ok("e1_extra"),
    output:
        data=E2_CACHE / "composite-data.json",
        # The host's census at build time: not pinned, and nothing reads it. See the file's
        # own note in composite_data.py.
        census=E2_CACHE / "composite-census.json",
    shell:
        "{E2_ENV} {PY} {SCRIPTS}/e2/composite_data.py"


# The other direction of (b): the CAI given up at four floors on the MFE axis.
rule e2_composite_dcai:
    input:
        auc=E2_CACHE / "auc01_fin.json",
    output:
        E2_CACHE / "dcai-at-mfe.json",
    shell:
        "{E2_ENV} {PY} {SCRIPTS}/e2/dcai_at_mfe.py"


# Panels (c)-(g): one protein's four fronts, and the four designs at the floor with their
# folds laid out and aligned. This step needs ViennaRNA.
rule e2_composite_example:
    input:
        # `example_coords.py` reads the deposit through E2_ENV, and picks its structure
        # files with a glob rather than by name, so no declared path would cover them.
        ok=source_ok("e2"),
    output:
        E2_CACHE / f"example-{E2_SLUG}.json",
    params:
        cell=E2_EXAMPLE["cell"],
        floor=E2_EXAMPLE["floor"],
        alt=E2_EXAMPLE["alt"],
    shell:
        "{E2_ENV} {RNA_PY} {SCRIPTS}/e2/example_coords.py {params.cell} {params.floor} "
        "{params.alt}"


# The plate. Fails, writing nothing, if any text covers another text, a mark or a panel it
# does not belong to. The geometry is solved in inches at the size given here, so the two
# numbers are part of the figure and not a convenience.
rule e2_fig_composite:
    input:
        data=E2_CACHE / "composite-data.json",
        dcai=E2_CACHE / "dcai-at-mfe.json",
        example=E2_CACHE / f"example-{E2_SLUG}.json",
        names=INPUTS / "motifs" / "names.json",
    output:
        pdf=FIGURES / "fig3_e2_composite.pdf",
        png=FIGURES / "fig3_e2_composite.png",
    params:
        slug=E2_SLUG,
        width=E2["figure"]["width_in"],
        height=E2["figure"]["height_in"],
        subset=E2["figure"]["subset"],
    shell:
        "{E2_ENV} {PY} {SCRIPTS}/e2/fig_composite.py 4 {params.slug} {params.subset} "
        "mfe {params.height} {params.width} {E2_CACHE} {output.pdf}"


# The figure and its cache against the digests they were published with.
rule e2_fig_composite_verify:
    input:
        manifest=ROOT / "pins" / "e2_composite.outputs",
        built=[E2_CACHE / "composite-data.json", E2_CACHE / "dcai-at-mfe.json",
               E2_CACHE / f"example-{E2_SLUG}.json",
               FIGURES / "fig3_e2_composite.pdf", FIGURES / "fig3_e2_composite.png"],
    output:
        STAMPS / "e2_composite.verified",
    shell:
        "{PY} {SCRIPTS}/check_manifest.py {input.manifest} {DERIVED} {output}"


# ---- The motif supplement plate -----------------------------------------------------------
# What a motif constraint can do to a protein, and what reaching each level costs. Three
# cached measurements and one drawing step, the same shape as the composite above, so a
# layout can be redrawn in a second without rescanning the deposit.
E2_SUPP = DERIVED / "e2" / "supp"
E2_SUPP_ENV = (
    f'E2_DEPOSIT="{source("e2")}" E2_SUPP="{E2_SUPP}" '
    f'E2_MOTIF_DIR="{INPUTS / "motifs"}" '
    f'E2_MOTIF_NAMES="{INPUTS / "motifs" / "names.json"}" '
    f'E2_MOTIF_EFFECTS="{INPUTS / "motifs" / "effects.json"}" '
    f'E2_F0_TIMINGS="{INPUTS / "e2" / "f0_timings_imported.json"}" '
    f'E2_CORPUS="{INPUTS / "e2" / "corpus.json"}" '
    f'E2_E1_EXTRA="{source("e1_extra")}" '
    # The grid ran on two hosts, and its early rounds recorded the host only through the
    # home directory in each line's environment.
    f'E2_FIRST_HOST="{E2["hosts"]["first"]}" '
    f'E2_SECOND_HOST="{E2["hosts"]["second"]}" '
    f'E2_SECOND_HOST_MARKER="{E2["hosts"]["second_marker"]}" '
    # The set files and their order come from the config, so the Roman numeral a panel
    # prints and the file a rule declares cannot drift apart.
    f'E2_MOTIF_SETS="{",".join(E2["motif_sets"][:3])}" '
    f'PYTHONPATH="{SCRIPTS / "e2"}"'
)


# Which motifs each protein can spell at all, over the codon lattice. Exact: every window,
# every synonymous choice.
rule e2_supp_admits:
    input:
        finished=source("e2") / "FINISHED.tsv",
        ok=source_ok("e2"),
        ok_e1_extra=source_ok("e1_extra"),
        sets=expand(str(INPUTS / "motifs" / "{name}.txt"), name=E2["motif_sets"][:3]),
        # The population panel (d) counts over, which the composite plate derives from E1's
        # sweep and asserts against.
        corpus=INPUTS / "e2" / "corpus.json",
    output:
        E2_SUPP / "admits.json",
    shell:
        "{E2_SUPP_ENV} {PY} {SCRIPTS}/e2/supp_admits.py"


# How much of each motif set survives on the front as levels are added.
rule e2_supp_presence:
    input:
        finished=source("e2") / "FINISHED.tsv",
        ok=source_ok("e2"),
        ok_e1_extra=source_ok("e1_extra"),
        sets=expand(str(INPUTS / "motifs" / "{name}.txt"), name=E2["motif_sets"][:3]),
    output:
        E2_SUPP / "presence.json",
    shell:
        "{E2_SUPP_ENV} {PY} {SCRIPTS}/e2/supp_presence.py"


# What reaching each constrained front costs against the unconstrained run.
rule e2_supp_cost:
    input:
        finished=source("e2") / "FINISHED.tsv",
        ok=source_ok("e2"),
        ok_e1_extra=source_ok("e1_extra"),
        # The denominator for an unconstrained front E2 imported rather than ran, which
        # therefore carries points and no timing.
        timings=INPUTS / "e2" / "f0_timings_imported.json",
    output:
        E2_SUPP / "cost.json",
    shell:
        "{E2_SUPP_ENV} {PY} {SCRIPTS}/e2/supp_cost.py"


rule e2_fig_supp_motif:
    input:
        admits=E2_SUPP / "admits.json",
        presence=E2_SUPP / "presence.json",
        cost=E2_SUPP / "cost.json",
        # Panel (d)'s 24 tick labels are read from here, so a rename there is a change to
        # the figure and has to rebuild it.
        names=INPUTS / "motifs" / "names.json",
        ok_e1_extra=source_ok("e1_extra"),
    output:
        pdf=FIGURES / "fig_supp_motif.pdf",
        png=FIGURES / "fig_supp_motif.png",
        prov=FIGURES / "fig_supp_motif.prov.json",
    params:
        width=E2["supp_figure"]["width_in"],
        height=E2["supp_figure"]["height_in"],
    shell:
        "{E2_SUPP_ENV} {PY} {SCRIPTS}/e2/fig_supp_motif.py "
        "{params.width} {params.height} {output.pdf}"


# The motif catalog as a table. The plate labels panel (d) by name because a
# name-over-sequence tick does not fit there, so the sequences are published here.
rule e2_supp_table:
    input:
        admits=E2_SUPP / "admits.json",
        names=INPUTS / "motifs" / "names.json",
        # The table's two editorial columns, so neither is written into the script.
        effects=INPUTS / "motifs" / "effects.json",
        ok_e1_extra=source_ok("e1_extra"),
    output:
        FIGURES / "tab_supp_motifs.tex",
    shell:
        "{E2_SUPP_ENV} {PY} {SCRIPTS}/e2/supp_table.py {output}"


# The supplement plate and its caches against the digests they were published with.
rule e2_fig_supp_motif_verify:
    input:
        manifest=ROOT / "pins" / "e2_supp_motif.outputs",
        built=[E2_SUPP / "admits.json", E2_SUPP / "presence.json", E2_SUPP / "cost.json",
               FIGURES / "fig_supp_motif.pdf", FIGURES / "fig_supp_motif.png",
               FIGURES / "tab_supp_motifs.tex"],
    output:
        STAMPS / "e2_supp_motif.verified",
    shell:
        "{PY} {SCRIPTS}/check_manifest.py {input.manifest} {DERIVED} {output}"
