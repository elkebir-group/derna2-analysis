# Protein records added or corrected after the first corpus snapshot

P22301 (interleukin-10, 178 residues) and P41159 (leptin, 167 residues) were added to E1's
corpus after its 2026-10-27 sweep. These files are the FASTA records E1 was given for them;
the residues equal the translation of E1's designs for both (`e2-out-1031/hull/*.front.jsonl`).

`P01706.fasta` (immunoglobulin lambda variable 2-11) was run under the label `P01707`, which
the DERNA dataset it comes from uses for this record. The runs are filed under that label and
the frozen deposits still carry it in their file names, so the label remains the key that joins
a result to its protein; everything a reader sees names the entry the sequence is. The rename
lives in `pipelines/scripts/supp_proteins.py` and the rule beside it.
