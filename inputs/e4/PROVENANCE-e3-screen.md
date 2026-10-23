# The six therapeutic natives

The screen's 55 UniProt proteins are the LinearDesign set, capped at 500 aa and
reproducible from the `_cpp` mirror. These six are the realistic mRNA targets
added beside them, and they are **not** in that mirror, so their provenance is
recorded here rather than being implicit in a path.

| input | residues | source |
|---|---|---|
| `P01588.fasta` | 193 | UniProt `P01588` (EPO_HUMAN), full entry |
| `P09919.fasta` | 207 | UniProt `P09919` (CSF3_HUMAN, G-CSF), full entry |
| `P01241.fasta` | 217 | UniProt `P01241` (SOMA_HUMAN, hGH), full entry |
| `P01563.fasta` | 188 | UniProt `P01563` (IFNA2_HUMAN, IFN-α2), full entry |
| `eGFP.fasta` | 239 | see below |
| `RBD.fasta` | 223 | `P0DTC2` residues 319–541, UniProt numbering |

Fetched from `https://rest.uniprot.org/uniprotkb/<accession>.fasta` on
2026-10-19.

**Full entries, not mature chains.** Four of these are secreted proteins whose
therapeutic form is the mature chain, without the signal peptide — EPO's is 166
of its 193 residues. The screen uses the **full UniProt entry** anyway, because
an accession is reproducible and a hand-trimmed chain is not, and because what
the screen measures is which motifs a composition can spell: a signal peptide
contributes residues like any other segment. Any price quoted for these proteins
is therefore a price for the entry, and a row that claims to be about a
therapeutic product would have to be recomputed on the mature chain. That is a
statement about scope, not an approximation to be corrected later.

**eGFP** is a laboratory construct with no UniProt accession of its own. Its
sequence here is the `aa_seq` of `e2/mRNAfold/configs/egfp_example.config`,
truncated at the eGFP terminus `…MDELYK`, which yields the standard 239-residue
construct. The config's own sequence continues for 45 further residues
(`RSRDISHGFPPAVAAQ…`), a fusion extension belonging to that paper's benchmark;
it is dropped here so the row is about eGFP rather than about a fusion. Both
forms are available from the same config if the comparison is ever wanted.

**RBD** has no accession either — it is a domain of the spike. The window 319–541
is the receptor-binding domain as conventionally delimited, and the extracted
sequence begins `RVQPTESIVRFP` and ends `NKCVNF`, which is the check that the
numbering was applied 1-indexed and inclusive.
