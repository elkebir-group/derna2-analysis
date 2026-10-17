# Peak memory of each method by protein length (a supplementary figure of the paper). It reads
# Figure 4's frozen deposit and per-protein table, so it shows the same proteins, length
# groups and finished runs as Figure 4d.


rule e1_memory_figure:
    input:
        sweep=RESULTS / "e1_fig4" / "sweep.jsonl",
        per=F4_DIR / "per_protein.csv",
        ok=STAMPS / "results_e1_fig4.ok",
    output:
        pdf=FIGURES / "fig_e1_memory.pdf",
        png=FIGURES / "fig_e1_memory.png",
    params:
        pin=F4["pin"],
    shell:
        "{PY} {SCRIPTS}/e1/fig_memory.py {input.sweep} {params.pin} {input.per} {output.pdf}"


rule e1_memory_figures:
    input:
        FIGURES / "fig_e1_memory.pdf",
