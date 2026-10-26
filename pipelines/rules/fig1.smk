# The overview figure: GLPE_HUMAN (78 aa) with the hbb UTRs. drivers/fig1/FIG1.md describes
# the instance, the guard and how the structurally optimal front was found.

FIG1 = config["fig1"]
FIG1_RUNS = RUNS / "fig1"
GLPE = INPUTS / "glpe_hbb"


rule fig1_figures:
    input:
        FIGURES / "fig1_glpe_hbb.pdf",
        DERIVED / "fig1" / "so_front.jsonl",


# Re-runs of the two fronts the deposit holds, each checked against it.
rule fig1_runs:
    input:
        FIG1_RUNS / "cdsonly" / "check.txt",
        FIG1_RUNS / "symmetric" / "check.txt",


# ---- figures ------------------------------------------------------------------------

rule fig1_plot:
    input:
        ok=source_ok("fig1"),
        so=DERIVED / "fig1" / "so_front.jsonl",
    output:
        FIGURES / "fig1_glpe_hbb.pdf",
        FIGURES / "fig1_glpe_hbb.png",
    params:
        results=source("fig1"),
    shell:
        "{PY} {SCRIPTS}/fig1/plot.py {params.results} {output[0]}"


# The structurally optimal front reassembled from its window reports; fails unless it
# equals the deposited front point for point.
rule fig1_so_front:
    input:
        ok=source_ok("fig1"),
        front=source("fig1") / "d_so_front.jsonl",
    output:
        DERIVED / "fig1" / "so_front.jsonl",
    params:
        windows=source("fig1") / "so_windows",
    shell:
        "{PY} {SCRIPTS}/fig1/so_front.py {params.windows} {input.front} {output}"


# The paper's overview figure: the modular front against UTR masking and concatenation, and
# the structures of three designs. Its designs need the pinned fold_seq (the modular front's
# window walks report energies only), so it is not part of fig1_figures or `all`.
rule fig1_overview_designs:
    input:
        ok=source_ok("fig1"),
        ok_fronts=source_ok("fig1_fronts"),
        bins=STAMPS / "bins_fig1.ok",
    output:
        DERIVED / "fig1" / "overview_designs.json",
    params:
        results=source("fig1"),
        mrnafold=source("fig1_fronts") / "mrnafold" / "compare.jsonl",
        motif=FIG1["overview_motif"]["name"],
        site=FIG1["overview_motif"]["site"],
        fold=lambda w: Path(FIG1["bin_dir"]) / "fold_seq",
        work=DERIVED / "fig1" / "overview_fold",
    shell:
        "{PY} {SCRIPTS}/fig1/overview_designs.py {params.results} {params.mrnafold} {params.motif} "
        "{params.site} {params.fold} {params.work} {output}"


rule fig1_overview:
    input:
        ok=source_ok("fig1"),
        ok_fronts=source_ok("fig1_fronts"),
        designs=DERIVED / "fig1" / "overview_designs.json",
        fasta=GLPE / "GLPE_HUMAN_78.fasta",
        utr5=GLPE / "hbb5.txt",
        utr3=GLPE / "hbb3.txt",
    output:
        FIGURES / "fig1_overview.pdf",
        FIGURES / "fig1_overview.png",
    params:
        results=source("fig1"),
        mrnafold=source("fig1_fronts") / "mrnafold" / "compare.jsonl",
    shell:
        "{PY} {SCRIPTS}/fig1/overview.py {params.results} {params.mrnafold} {input.designs} "
        "{input.fasta} "
        "{input.utr5} {input.utr3} {output[0]}"


# ---- runs -----------------------------------------------------------------------------

def fig1_flags(w):
    utr5 = (GLPE / "hbb5.txt").read_text().strip()
    utr3 = (GLPE / "hbb3.txt").read_text().strip()
    return FIG1["fronts"][w.arm].format(utr5=utr5, utr3=utr3)


rule fig1_front:
    input:
        bins=STAMPS / "bins_fig1.ok",
        fasta=GLPE / "GLPE_HUMAN_78.fasta",
    output:
        front=FIG1_RUNS / "{arm}" / "front.jsonl",
        report=FIG1_RUNS / "{arm}" / "report.json",
    log:
        FIG1_RUNS / "{arm}" / "stderr.log",
    params:
        walk=lambda w: Path(FIG1["bin_dir"]) / "scope35_walk_front",
        flags=fig1_flags,
        n=FIG1["n_aa"],
    shell:
        "/usr/bin/time -v {params.walk} --fasta {input.fasta} --n {params.n} "
        "--out {output.report} --front {output.front} {params.flags} 2> {log}"


FIG1_DEPOSITED = {
    "cdsonly": "a_cdsonly_front.jsonl",
    "symmetric": "cd_symmetric_guarded_front.scored.jsonl",
}


rule fig1_check_front:
    input:
        run=FIG1_RUNS / "{arm}" / "front.jsonl",
        ok=STAMPS / "results_fig1.ok",
        deposited=lambda w: RESULTS / "fig1" / FIG1_DEPOSITED[w.arm],
    output:
        FIG1_RUNS / "{arm}" / "check.txt",
    shell:
        "{PY} {SCRIPTS}/fig1/check_front.py {input.run} {input.deposited} {output}"
