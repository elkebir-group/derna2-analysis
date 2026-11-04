# PRT-05: two sequences, and which results each belongs to

`PRT-05.fasta` here is the sequence E1 designed: its residues are the translation of every
design in E1's deposited front (`hull/PRT-05.front.jsonl`). It replaced a copy that differed at
8 of its 1,082 residues, which is still what E2's screen holds.

The divergence is real and not yet resolved:

- **E1** (the Pareto-front comparison) ran this sequence.
- **E2** (the motif work) ran the other one. Its front header, which is where the motif scripts
  read residues from, carries it: `screen/f0-xl/PRT-05.jsonl`, matching `screen/corpus-xl/`.

Measured consequence, over the three motif sets the supplement plate reports: **none**. Both
sequences admit all 24 motifs, so PRT-05 contributes "admissible" to every motif either way and
no count in the motif table moves. PRT-05 is also not among the 45 proteins whose fronts are
finished at every motif level, so it enters no cost figure.

What does depend on the sequence is E2's own per-cell run record for PRT-05 (wall, memory,
points, hull vertices). Those describe a run on the other protein.

Correcting E2's side means recomputing its PRT-05 front, which is a change to a frozen deposit
and belongs to whoever owns that lane. The file here is left as E1's sequence because that is
the one the paper's front comparison reports.
