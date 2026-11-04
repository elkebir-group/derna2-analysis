# Protein structures — where each model comes from

`AF-P15421-F1-model_v6.pdb` — `sha256 b7ba0f760f4fb31dad93008c63eac69fa36c8d83efbaa4c1ba32867ae26a846f`,
retrieved 2026-10-27 from the AlphaFold Protein Structure Database
(`https://alphafold.ebi.ac.uk/files/AF-P15421-F1-model_v6.pdb`, the `pdbUrl` the database's
prediction API gives for P15421; model version 6). AlphaFold DB data are released under CC BY 4.0.

It is the AlphaFold prediction for glycophorin-E (UniProt P15421, GLPE_HUMAN), residues 1–78,
the whole protein. UniProt lists no experimental structure for P15421 (no PDB cross-reference),
so this predicted model is the only structure available. Its 78 C-alpha atoms spell the sequence
of `../glpe_hbb/GLPE_HUMAN_78.fasta` exactly, and consecutive C-alpha atoms are 3.3–3.9 Å apart.
Per-residue confidence (pLDDT, in the B-factor column) is 37–88: the two helices are 70–88, the
linker between them below 50.

The overview figure draws its C-alpha trace as a ribbon (`pipelines/scripts/fig1/overview2.py`).
It is an illustration of the protein, not a structural claim.
