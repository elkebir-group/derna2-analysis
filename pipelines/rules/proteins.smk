# The supplement's table of every protein the paper designs, with the experiments run on each:
# the Pareto-front comparison (all of them, each also without forbidden motifs), the motif
# levels (E2's finished set) and the UTR experiment (E3's cells).

PROTEIN_DIRS = [INPUTS / d for d in ("uniprot", "e1-natives", "prt", "e3")] + [INPUTS]
# The accession each sequence is, where its run label names another entry; see supp_proteins.py.
RENAMED = {"P01707": "P01706"}
PROTEINS = [RENAMED.get(p, p) for p in json.load(open(INPUTS / "e2" / "corpus.json"))["members"]]


def protein_fasta(protein):
    return next(d / f"{protein}.fasta" for d in PROTEIN_DIRS if (d / f"{protein}.fasta").exists())


rule supp_proteins:
    input:
        corpus=INPUTS / "e2" / "corpus.json",
        motif_levels=DERIVED / "e2" / "composite" / "auc01_fin.json",
        utr_cells=INPUTS / "e3" / "cells.txt",
        fasta=[protein_fasta(p) for p in PROTEINS],
    output:
        FIGURES / "tab_supp_proteins.tex",
    params:
        dirs=" ".join(str(d) for d in PROTEIN_DIRS),
    shell:
        "{PY} {SCRIPTS}/supp_proteins.py {input.corpus} {input.motif_levels} {input.utr_cells}"
        " {output} {params.dirs}"


# The table against the digest it was published with.
rule supp_proteins_verify:
    input:
        manifest=ROOT / "pins" / "supp_proteins.outputs",
        built=[FIGURES / "tab_supp_proteins.tex"],
    output:
        STAMPS / "supp_proteins.verified",
    shell:
        "{PY} {SCRIPTS}/check_manifest.py {input.manifest} {DERIVED} {output}"
