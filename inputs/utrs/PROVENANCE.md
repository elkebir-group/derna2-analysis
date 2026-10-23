# UTR pin — where each sequence comes from

`utrs.fasta` — `sha256 dc8b4f7eb8077f5910839b4f09cb95a82a119d6ee447e201299f6969876f4058`, pinned
2026-10-19. Every record was extracted programmatically from the source document named
below; none was transcribed by hand. Each header carries the record's own length and
the first 16 hex of its SHA-256, so a record can be checked without the whole file.

Sequences are written as DNA (`ACGT`). `derna` reads `T` as `U` on `--utr5` / `--utr3`;
mRNAfold's config wants RNA, so its driver converts. BNT162b2's `U` is in fact
N1-methylpseudouridine and mRNA-1273's likewise; base identity is unaffected and no
energy model here distinguishes them, which is a limitation to state rather than model.

## Sources

| record | document | retrieved |
|---|---|---|
| `hbb_5p`, `hbb_3p`, `hbb132_3p` | RefSeq **NM_000518.5** (HBB), `KEYWORDS RefSeq; MANE Select.`, 628 bp, dated 06-SEP-2026, via NCBI efetch | 2026-10-19 |
| `hba_5p`, `hba_3p` | RefSeq **NM_000558.5** (HBA1), `KEYWORDS RefSeq; MANE Select.`, 577 bp, dated 06-SEP-2026, via NCBI efetch (`sha256 1c3361e9…`) | 2026-10-19 |
| `bnt_5p`, `bnt_3p` | **FDA GSRS** record for tozinameran, UNII `5085ZFP6SJ`, `substanceClass nucleicAcid`, `nucleicAcidType "mRNA VACCINE"`, `sequenceType COMPLETE`, subunit 1 = 4283 nt (`sha256 db6cc790…`). Corroborated by BioNTech CTD Module 3.2.S.1.2 released under UK FOI 22-1116 (`sha256 91a950e1…`) | 2026-10-19 |
| `mrna1273_5p`, `mrna1273_3p` | **FDA GSRS** record for elasomeran, UNII `EPK39PL4R4`, same classes, subunit 1 = 4100 nt (`sha256 6b2b2243…`) | 2026-10-19 |
| `curevac_5p`, `curevac_3p` | **US 11,241,493 B2** (CureVac, *Coronavirus vaccine*) sequence listing, SEQ ID NO: 148, the complete 4094-nt mRNA. Element identity corroborated by **US 2022/0233568 A1** (the US national phase of WO 2019/077001, CureVac's UTR patent) SEQ ID NO: 1 / 23, and by RefSeq **NM_000414.4** (HSD17B4, MANE Select) and **NM_002795.3** (PSMB3) | 2026-10-19 |

## How each boundary is established

The boundary, not the sequence, is where a pin like this goes wrong, so each record
says whether its source states the span or whether arithmetic produced it.

- **`hbb`, `hba` — from the record's own CDS feature.** `NM_000518.5` annotates
  `CDS 51..494`, so the 5′ UTR is 1–50 and the 3′ UTR 495–628; `NM_000558.5` annotates
  `CDS 38..466`, giving 1–37 and 467–577. Both coding sequences check out: lengths 444
  and 429 nt, divisible by three, beginning `ATGGTG` and ending in a stop. HBA1's
  record also annotates `regulatory 556..561 /note="hexamer: AATAAA"` and
  `polyA_site 577`, both inside the 3′ UTR span as expected.

- **`bnt` — stated by the primary source.** The GSRS record annotates the sequence it
  prints: a `NUCLEIC ACID FEATURE` of type `Site Range`, value `1_1-1_53`, named *5′-UTR
  sequence of the human alpha-globin mRNA with an optimized "Kozak sequence"*, and a
  second, `1_54-1_3878`, named *SARS-COV2 SPIKE CODING SEQUENCE*. Both strings were
  confirmed present in the retrieved record. So the 5′ UTR is 1–53 and the 3′ UTR
  begins at 3879. Its **end** is arithmetic, pinned on both sides: the A30 tract begins
  at 4174 in the printed sequence (`A×30 + GCATATGACT + A×70`), so the element is
  3879–4173, 295 nt, beginning `CTCGAG` and ending `GCTAGC`. The BioNTech Module 3
  document states the same elements with coordinates offset by exactly 15 nt, that copy
  being the BA.4/BA.5 drug substance whose spike is 15 nt shorter; every boundary
  reconciles under that offset, including total length 4268 + 15 = 4283 and the stated
  `FI element (nucleotides 3864 to 4158)` = 295 nt. The element as the manufacturer
  delimits it **excludes both stop codons**.
  One secondary source (Mamaghani et al., *Cells* 2024;13:1242) gives the 5′ UTR as 52
  nt; its 295 for the 3′ UTR agrees. 53 is pinned here, on the primary annotation.

- **`mrna1273` — no manufacturer annotation exists.** The elasomeran GSRS record
  carries **zero** `Site Range` features, confirmed by inspection, so this is the one
  pair whose boundary rests on secondary sources plus arithmetic. The 5′ UTR is
  unambiguous: the first `ATG` is at 58, so the UTR is 1–57. The 3′ UTR starts at 3886,
  after three tandem stops `TGA TAA TAG` at 3877–3885 — Xia (*Vaccines* 2021;9:734)
  places the element "between the last stop codon and a poly(A) tail", and Mamaghani's
  Table S1 prints it starting `GCTGGAGCC…`.
  Its **end** is a real disagreement of one adenosine at the poly(A) junction: 111 nt
  (Mamaghani, which prints the sequence) against 110 (Xia). **111 is pinned**, on
  evidence independent of both: at 111 nt the element aligns 1:1 to the native HBA1
  3′ UTR pinned here as `hba_3p` and differs at **exactly 2 of 111 positions** — 19–20,
  `AT` → `TA`, turning the native in-UTR `ATG` at 19–21 into `TAG` and so removing a
  spurious start codon. That is a designed substitution, and it only lines up under the
  111-nt framing; at 110 the two cannot be aligned 1:1 at all. For the 110-nt framing,
  drop the trailing `A`.

- **`curevac` — spans read off the construct, element identity stated by the patent.**
  Two framings exist and the choice changes what sits at the boundary the constraint
  acts on, so both are recorded and one is pinned.
  The patent defines the *elements* as discrete SEQ IDs: a 62-nt 5′ UTR element
  "derived from a HSD17B4 gene" (SEQ ID NO: 1 in US 2022/0233568 A1 = SEQ ID NO: 231 in
  US 11,241,493 B2) and a 57-nt 3′ UTR element "derived from a PSMB3 gene" (SEQ ID
  NO: 23 = 253). Both check out against the native genes: the 62-mer is `NM_000414.4`
  positions **18–79**, the native 5′ UTR with its first 17 nt — the pyrimidine-rich
  5′TOP motif the patent says to exclude — removed, ending exactly where the native CDS
  begins; the 57-mer is `NM_002795.3` positions **733–789**, starting immediately after
  the native stop. Neither contains an internal `ATG`; the 3′ element carries `AATAAA`
  at 38.
  **What is pinned is the as-built framing**: everything between the coding sequence and
  the poly(A) tail, which is the same rule the other four pairs already follow. Read off
  SEQ ID NO: 148, whose coding sequence runs 78–3899 (3822 nt, 1274 codons, stop `TGA`)
  and whose final 100 nt are pure `A`: the 5′ UTR is 1–77 and the 3′ UTR 3900–3994.
  The as-built 5′ UTR adds a 6-nt leader, a HindIII site and the **Kozak `ACC`** that
  the element framing omits — and every other pinned 5′ UTR ends in that `ACC`, so the
  element framing would have given this pair alone a start codon with no Kozak context,
  at exactly the junction where the `5-cds` class acts. The as-built 3′ UTR likewise
  carries a SpeI site, a BglII site and the 24-nt histone stem-loop that the patent
  enumerates as an element separate from the 3′ UTR; it is kept because it is
  physically present between the stop codon and the tail, and because `bnt_3p` keeps
  its own cloning scar (`CTCGAG`, XhoI) on the same reasoning.
  The element framings, if ever wanted: 62 nt `sha256 834b2dbc3c995a9b`, 57 nt
  `sha256 0829ff91cfdccf86`.
  **Attribution:** this pair belongs to **CV2CoV**, CureVac's second-generation
  candidate, not to CVnCoV. Checked rather than assumed: the FDA GSRS record for
  zorecimeran (CVnCoV, UNII `5TP24STD1S`, 4039 nt) contains **neither** element, and
  uses a 15-nt leader and a mutated α-globin 3′ element instead. No retrieved document
  states CV7202's UTRs. The chain for this pin is therefore: the literature states
  CV2CoV uses an HSD17B4 5′ UTR and a PSMB3 3′ UTR, and CureVac's own patents define
  those named elements as these sequences — one link weaker than `bnt`, where the
  regulator's record annotates the shipped molecule.
  **Overlap check**, since the point of this pair is independent material: the longest
  common substring against every other record in this file is **7 nt**, which is
  chance-level at these lengths. The unrelated-lineage premise holds.

## `hbb` carries two 3′ UTRs, and the difference is measured

`hbb132_3p` is the sequence pinned in the repo since Scope 19
(`experiments/scope-19-utr-support/m5/hbb3.txt`), which had no recorded source document.
Checked against `NM_000518.5`: the 5′ UTR is **identical**, and the 3′ UTR is the MANE
3′ UTR **truncated by 2 nt**, the terminal `AA` absent — positions 495–626 rather than
495–628. The truncation sits at the 3′ terminus, where a poly(A) tail would begin.

**E4 uses `hbb_3p`, the 134-nt sequence the source document delimits** (decided
2026-10-19). `hbb132_3p` stays pinned because every earlier result for `hbb` — Scope
37's shift of −31.40 kcal/mol, its census verdicts, and this experiment's own
three-class control — was computed on it, so any comparison back to those numbers must
name which sequence it used. A cell records the sequence it ran. The two differ only in
the terminal `AA`, so the three-class control is worth re-running on the 134-nt pair to
get this experiment's own shift figure rather than inheriting Scope 37's.

## Not pinned

**`min`** — no pair was sourced. The published minimal construct (Mamaghani et al.,
*Cells* 2024;13:1242, "minRNA") states a 4-nt 5′ UTR (`CAAG`, the TISU element) or a
Kozak variant (`ACC`), and **deliberately has no 3′ UTR at all**: "UAA stop codon
directly followed by a poly-A tail". Its Table S1's 3′ UTR section lists only long
comparators. A search of four further candidate papers found no construct printing both
a short 5′ and a short 3′ UTR. No minimal 3′ UTR is synthesized to fill the slot; the
grid runs with the pairs above until a sourced one exists.
