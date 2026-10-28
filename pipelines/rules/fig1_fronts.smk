# Figure 1's three fronts on GLPE_HUMAN with the hbb UTRs: the plain front, the front under
# each of three CDS–UTR constraints, and a structurally optimal design per constrained hull
# vertex. drivers/fig1-fronts/README.md describes the runs and reads the tables below.


rule fig1_fronts_tables:
    input:
        ok=source_ok("fig1_fronts"),
    output:
        DERIVED / "fig1_fronts" / "tables.md",
    params:
        results=source("fig1_fronts"),
    shell:
        "{DRIVER_PY} {DRIVERS}/fig1-fronts/make_tables.py {params.results} "
        "{params.results}/mrnafold/compare.jsonl > {output}"
