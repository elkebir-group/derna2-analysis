# The eight middle-band proteins

E3's population is 30 proteins: the 22 of the short band and eight of E2's middle band
(500-1000 aa), each on the five UTR pairs, 150 cells.

| input | residues |
|---|---|
| `Q8TAA9.fasta` | 524 |
| `Q8N1G4.fasta` | 583 |
| `P31645.fasta` | 630 |
| `Q9UNN5.fasta` | 650 |
| `Q06187.fasta` | 659 |
| `P05160.fasta` | 661 |
| `Q9BXB4.fasta` | 747 |
| `Q9HCI6.fasta` | 787 |

Source: `data/uniprotSeq/<accession>.fasta` in the `uniprotSeq` corpus of DERNA
(https://github.com/elkebir-group/derna, commit `adecda82`, as for the short band), by way
of E2's middle-band selection of 2026-10-24, with CRLF line ends stripped; residues unchanged.
Q9BXB4 and Q9HCI6 reached that selection with bare `>accession` headers, so their files here
are taken from that corpus directly, which keeps the UniProt header that names the protein.
Each sequence equals the `protein.sequence` of the cell headers in the E3 deposit.

E3 also started three other middle-band proteins, Q96NI6, Q9H3R1 and Q9Y5G9. They were dropped
before their arms completed and are not in the population.
