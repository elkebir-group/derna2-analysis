# Paper benchmark E3 — the price of motif-freedom

Written 2026-10-19. Pre-registered: every selection rule, metric and decision
rule below is fixed before any cell of the headline grid runs.

E1 measures one λ per cell; E2 measures what the whole front costs. **E3 measures
what a constraint costs**: the exact Pareto front of (MFE, CAI) is computed over
the unconstrained design set, then over the designs that avoid a set of forbidden
motifs, then over those that also avoid a second set, a third and a fourth. Four
nested feasible sets, five nested fronts, one curve per instance.

The question is **not** whether derna2 can exclude motifs — Scope 17 built that and
the paper states it as MF-BRD (`prob:motif`). It is what exclusion *costs*: how far
the front moves, how many Pareto points disappear, and what the exact computation
costs against the unconstrained front.

**E3 has no comparator arms, and that is a result rather than a gap.** mRNAfold has
only `banned_codon` / `forced_codon` (position × codon) plus a post-hoc `FindMotif`
check on its output; LinearDesign, LinearCDSfold and DERNA v1 have nothing. A
constraint whose occurrences straddle codon boundaries is not expressible in any of
them at all, so there is nothing to race. E3 measures one method against its own
unconstrained baseline.

## The nesting

| front | feasible set | label |
|---|---|---|
| `F0` | all of `V(w)` | baseline |
| `F1` | `F0` minus designs containing a set-I motif | + manufacturing |
| `F2` | `F1` minus set II | + processing |
| `F3` | `F2` minus set III | + repeats |
| `F4` | `F3` minus set IV | + initiation |

The sets are **cumulative**: `F4`'s run forbids `I ∪ II ∪ III ∪ IV`, passed as four
repeated `--motif` files rather than as a concatenated fifth file, so each entry
keeps its set label in the output. The order is fixed and is not a free parameter of
the report — it runs from the constraint a protocol imposes on you, outward to the
ones biology argues for, with binding density rising along the way so that
consecutive fronts separate progressively.

**Nesting is a free soundness gate.** The feasible sets shrink monotonically, so at
every λ the scalarized optimum is monotone non-decreasing:

    obj_λ(F0) ≤ obj_λ(F1) ≤ obj_λ(F2) ≤ obj_λ(F3) ≤ obj_λ(F4)

A violation is a bug, never a finding. It is checked at every λ of every cell.

## The four motif sets

Files are in `experiments/paper-benchmark/motifs/`, one pattern per line with `#`
comments and a source note per entry. `T` is read as `U`; IUPAC codes are expanded
at trie build time (Scope 17 M1, D3). Each entry carries its own `removed` count in
the output — a motif that forbade nothing must be visible as such.

| set | file | content | shape |
|---|---|---|---|
| I | `set1-manufacturing.txt` | 8 Type IIS sites and reverse complements (the `p8` set every prior speed measurement used) + 9 common six-cutters | GC-rich 6–8mers, the sparsest set |
| II | `set2-processing.txt` | cryptic poly(A) signals `AATAAA`, `ATTAAA`; 5′ splice donor consensus `GGTRAGT` | AU-rich 6–7mers |
| III | `set3-repeats.txt` | homopolymer runs `AAAAAA`, `TTTTTT`, `CCCCCC`; `GGGG` for G-quadruplex propensity | the densest set |
| IV | `set4-initiation.txt` | internal Kozak-like initiation context `RCCATGG` | one 7-mer |

Two things the files record, because both are easy to get wrong:

- **A fixed 6-mer is exactly equivalent to the run-length rule `A{6,}`** on the
  exclusion axis, since every longer run contains the 6-mer. The homopolymer family
  looks like it needs a different instrument and does not.
- **Of two nested patterns, the shorter is the constraint.** Set IV lists `RCCATGG`
  and *not* the fuller `GCCRCCATGG`, because every occurrence of the 10-mer contains
  the 7-mer: listing both would leave the 10-mer unable to be the entry that removed
  anything. The same check was run within and across all four sets.

### What the substring instrument cannot state

Set IV's single entry is the whole of the translation category that hard substring
exclusion can express. The rest is reported as a limitation **with evidence**, which
is the stronger contribution:

- **Out-of-frame AUG needs frame.** Excluding `ATG` as a substring forbids every Met
  codon, so every protein is infeasible. The product lattice does carry the codon
  position, so a frame-aware axis is buildable — named as future work, not claimed.
- **Rare-codon clusters are a density rule** ("≥ k rare codons in a window of w"),
  not a fixed string.
- **A G-run threshold of 4 is infeasible, and the genetic code sets the floor at 5.**
  Trp has only `UGG` and Gly is `GGN`, so a Trp immediately followed by a Gly spells
  `UGG|GGN` and carries a 4-G run at every synonymous choice. `GGGG` is therefore
  infeasible for **any protein containing WG** — measured on P31645 at residue 82,
  where the walk reports infeasibility at lattice construction in 8 ms. Only Trp and
  Met are single-codon, so 4 is the longest run the code can force, and `GGGGG` is
  always satisfiable; the instance confirms it. So a G-quadruplex constraint at the
  biologically argued threshold of 4 is **not expressible as hard exclusion at all**
  and needs a window or a soft penalty, while set III carries `GGGGG`. The threshold
  is set by the code, not by the biology, which is the point worth reporting.
- **Poly-proline stalls are expressible and self-defeating.** All four Pro codons are
  `CCN`, so forbidding `CCNCCN` makes **every protein containing adjacent PP
  infeasible**, with no design to find. Measured once, as its own row, on a
  PP-bearing and a PP-free protein from the screen — a one-line result, not a front.
- **The 3′ splice acceptor** is a composition rule over a window (its polypyrimidine
  tract), not a k-mer. Omitted from set II by design, not by oversight.
- **CpG density and uridine content** (innate immunity) are composition targets for
  the same reason, and are out of scope for the exclusion axis entirely.
- **Shine–Dalgarno `GGAGG`** is bacterial and irrelevant to an mRNA therapeutic in
  human cells.

### How motif matching meets the flanks

Settled from the code rather than left to a test: **every occurrence overlapping the
coding sequence is forbidden, including one straddling a CDS/UTR boundary**, since
the codons on the coding side are being chosen (ADR 0043). Every `derna` coordinate
is a CDS position; the UTRs are fixed sequence. Two consequences bear on E3:

- An occurrence lying **wholly inside a fixed UTR** cannot be removed by any codon
  assignment. The run proceeds, warns, and records it in the header's
  `motifs_in_fixed_flanks` — a field the screen reads rather than re-deriving.
- An occurrence straddling the **5′** boundary reaches into the start codon, which
  the protein fixes. **If the flank spells the rest of it, the run is refused.**

That second clause is a live risk for exactly one cell of this plan: set IV's
`RCCATGG` straddles the junction if the 5′ UTR ends in `RCC`. So the flank
cross-check below carries a **pre-check**: before it runs, test the `hbb` 5′ UTR's
last 6 nt against all four sets. A refusal there is a property of that UTR and the
cross-check moves to a flank pair that does not trip it, recording which and why.
The headline grid runs without flanks, so nothing there can hit this.

## Instances

Chosen by measurement, not by assumption. Two repo facts set up the tension that
selection has to resolve:

- **A motif that cannot be spelled forbids nothing.** Scope 27 found **4 of 6**
  standard six-cutters have **no spellable placement at all** in a 78-residue
  protein. A Δobjective of zero from a constraint that forbade nothing is not a
  price (AGENTS.md, *check the baseline binds*); it reads identically to "exclusion
  is cheap" and is a different claim.
- **The opposite pull:** exclusion is the expensive direction, and the exact
  motif-excluded *front* has never been computed at any length (see *Dependency*).

Short instances do not carry the motifs; long ones may not finish. The reasoning here
was that binding is driven by **composition** rather than length — a GC-rich protein
spells Type IIS sites, an AU-rich one spells poly(A) signals and U-runs — so the
screen would read composition while the ladder handled length.

**The screen refutes that, and the ladder is dropped, so both halves of the sentence
have changed.** Length is the driver: of 59 proteins with a complete `F0`, those
whose front carries all four sets number 0 of 20 at ≤ 211 aa, 2 of 20 at 212–343 and
3 of 19 above 343. Composition does not select into that pool — its GC propensity
median is 0.6547 against the corpus's 0.6680. So the screen reads both axes and
reports length as the one that separates them, and there is no ladder to hand length
off to.

### 1. The binding screen (cheap, runs first, no new code)

Over all 55 `uniprotSeq` proteins at native length plus six therapeutic natives —
EPO (P01588), G-CSF (P09919), hGH (P01241), IFN-α2 (P01563), eGFP, SARS-CoV-2 RBD —
compute the **baseline front `F0` only** and record per protein per set:

- `occurs` — motif occurrences in the designs **on the baseline front**, the only
  occurrences that can move the front;
- `spellable` — placements the lattice admits at all (the Scope 27 quantity), which
  separates *this protein's optimum happens to avoid it* from *no design could
  contain it*;
- both counted **per entry**, so a set is never credited for a motif that never
  appeared;
- `pp_adjacent` — whether the protein contains adjacent PP, which selects the
  instances for the poly-proline row: every grid member carrying it, however many
  that turns out to be.

### Which entries the code forces, computed exhaustively

`forced_motifs.py` answers this from the codon table alone, with no runs. An entry is
**forced** by a peptide when every synonymous encoding of that peptide contains it, so
any protein carrying the peptide is infeasible under a set holding that entry. A motif
of length `m` spans at most `ceil((m + 4) / 3)` residues, so the peptide search is
complete at that length rather than truncated.

Result (`e3-probe/forced.jsonl`): **all 25 entries of the four sets are satisfiable on
every protein** — no forcing peptide exists for any of them. The two rejected
thresholds are each forced by exactly one peptide: `GGGG` by **WG**, `CCNCCN` by **PP**.
Those two peptides are the whole of the obstruction, which is why raising the G-run
threshold to 5 removes it.

The same script, run over the union of all 25 patterns, finds **no peptide of length
≤ 4 that forces the union**, so the cumulative set has no *local* obstruction either.

That is weaker than global feasibility and is stated as such: a local window forcing a
violation is sufficient for infeasibility but not necessary, since avoiding a motif in
one window could in principle force one in the next, and no single window would show
it. The lattice check is what settles an instance, and it is why:

**Step 0 of the screen is a feasibility pre-check**, per instance per set, and it is
cheap: infeasibility is detected at lattice construction, measured at 8 ms, so the
check costs nothing beside the front it precedes. It is mandatory because two of the
four sets have already been found to carry an infeasible entry on sight — `CCNCCN`
by construction, and `GGGG` on any WG-bearing protein — and an infeasible cell
reported as a front of zero points would read as a price.

`occurs` is computed by an independent substring scan over the emitted `sequence`
field of each point line, sharing no code with the DP. `spellable` needs a lattice
query, which is what `motif_outside_census` already does.

### 2. Selection rule — ten instances in each of three length bands

**Restated 2026-10-20 against the experiment's goal.** What the paper needs is a
constrained front for **at least ten instances in each of three length bands**:
under 500 residues, 500–1000, and 1000 and above. The earlier rule — six instances
stratified by length tertile and GC propensity, with a five-instance short block
beside them — answered a narrower question and is superseded in full. What survives
from it is the discipline: every band's membership is fixed by a rule stated before
its cells run, and no instance is substituted after its numbers are seen.

The bands are not equally reachable, and the rule differs per band for reasons that
are measured rather than assumed.

**The short band: the 18 proteins that bind sets I–III**, 163–478 aa, which is the
whole of that pool rather than a sample of it. Six of the eighteen bind all four sets
and run under all four nests; the other twelve run under the first three, with set IV
reading *not binding* beside the screen's `occurs = 0`. That is **60 cells** — 24 plus
36 — against a floor of 22.0 h single-threaded, summing the proteins' own measured
baseline walls. Set IV therefore remains a six-instance claim while the first three
nests carry eighteen, and the band's length coverage is a property of the pool rather
than of a stratification: no protein at or below 211 aa binds all four sets, which is
why no stratified rule could have reached ten instances here.

**The middle band: ten of the 32 mirror proteins at 500–1000 aa**, selected once the
501–700 screen now running reports which of its nineteen bind. The selection rule
within the band is the same as the short band's — take every I–III binder, ascending
in length until ten — and it cannot be written more tightly until the screen reports.

The 13 proteins at 701–1000 were excluded on cost, each being 3–12 h of baseline
front. **Five of them arrive free instead** (2026-10-20): E2's pass 1 re-run walks
twelve natives of this band — every third by length of the 32 — and five sit above
E3's 700 aa ceiling: Q96NI6 (719), Q9BXB4 (747), Q9HCI6 (787), Q9H3R1 (872) and
Q9Y5G9 (962), about 30 h of front E3 does not have to spend. E2's walks emit the
same design sets E3's screen does, so the occurrence counts can be read straight off
them.

**The build difference is checked on the quantity E3 consumes, not on a proxy.** The
re-run is not E3's pinned `8603d37d` but that commit merged with a walk fix, and E2's
own A/B for it reports identical *objective multisets* on three proteins near 470 aa.
That is the wrong relation for this purpose: two designs can share an (MFE, CAI) pair
and spell different nucleotides, and `occurs` is a substring scan over nucleotides.
Checked directly on E2's A/B artifacts, two of whose three proteins are short-band
members: on Q9NWZ3, Q9UBM8 and Q9Y512 the two builds and E3's own `8603d37d` fronts
agree as **sets of design sequences**, and the per-set occurrence totals are equal
across all three arms — 7708/501/14192/1957 on Q9NWZ3, 7732/374/8397/3671 on Q9UBM8,
14513/0/12540/3316 on Q9Y512. So the inherited fronts are admissible for selection.

They are admissible for **selection only**. No published price comes from them: every
priced front is computed here on the pinned commit, and the inherited five contribute
binding counts that decide membership. The band's candidates therefore
reach **962 aa rather than 694**, which matters because its upper half is where a
complete constrained front stops being reachable.

Seven of E2's twelve overlap the nineteen already screening here, and E3 does not drop
them. A pre-registered screen with six holes in it has a completeness that depends on
another experiment's schedule, and E2's timing run needs its own walls regardless, so
its half of the duplication is unavoidable and only E3's could be skipped — at a
saving of 12 h, against a dependency the experiment carries across a host move.

**The long band: ten instances at 1000 aa and above, as partial fronts.** Decided
2026-10-20: a complete constrained front is out of reach there and the band publishes
**anytime curves** instead, exactly as E2 does for a censored cell — each cell runs
under a fixed wall and memory budget and reports the front it reached, which is a
lower bound on the cost of completeness and a complete record of what the budget
bought. The supply is 15 mirror proteins plus the four Cas orthologs at 1053–1629 aa,
Cas9 at 984 aa falling in the middle band rather than this one.

**Selection in the long band cannot use a baseline front**, one costing 12 h at
1000 aa and 65 h at 1629 on the measured `n^3.53`, which would spend longer screening
the band than running it. It uses the **lattice census** instead: `spellable` is a
lattice query, cheap at any length and load-immune, so eligibility is *the lattice
admits every entry of the nest*, and `occurs` is read off the partial front the cell
produces anyway.

That is a genuinely weaker instrument and the asymmetry is the point: on a partial
front `occurs > 0` still proves the entry appears among optimal designs and so can
move the front, while `occurs = 0` proves nothing, the vertex that would have carried
it possibly lying in the part of the front the budget never reached. So **a
non-binding cell in the long band is reported as unresolved, not as vacuous**, and no
vacuity claim rests on a partial front.

**And the census is not accurate enough to carry selection on its own, measured
against the screen that does: 33% for `processing`.** On the 61 screened proteins the
census admits 57 as able to spell a processing motif while only 19 have optima that
contain one, so two selections in three would land on a cell whose constraint
forbids nothing. The same instrument is 93–100% predictive for `manufacturing` and
`repeats`, which is what makes the processing figure a property of that nest rather
than of the census.

**Decided 2026-10-20: ten long candidates are screened properly, and all three nests
keep their claim in all three bands.** The alternative was to restrict the long
band to `manufacturing` and `repeats`, where the census stands unaided, and spend no
compute on selection at all. It is rejected because a nest the paper reports in two
bands and withdraws in the third reads as a result about long proteins when it is a
result about the selection instrument, and the price of motif-freedom under processing
motifs at length is the row the table exists to carry.

The screen costs roughly 150–200 h and precedes every long-band cell. **That estimate
is derived from tamazirt's E-core walls and does not transfer**: it is re-derived on
the new host by the three cost-and-memory runs named in the handoff's *Moving to
another host* section, together with every other budget here, and the figure above
stands only as the order of magnitude that made the decision.

### 2a. What the three bands cost, and why the host had to change

Floors, single-threaded, assuming a constrained front costs what a free one does —
which it does not, exclusion being the expensive direction, so each is a lower bound:

| band | instances | cells | floor | note |
|---|---|---|---|---|
| < 500 aa | 18 | 60 | 22.0 h | complete fronts throughout |
| 500–1000 aa | 10 | 40 | 170 h | complete below n ≈ 554, budgeted above |
| ≥ 1000 aa | 10 | 40 | 240–960 h | budget × 40; never complete |

The long band's figure is not a floor but a *choice*: 40 cells times whatever budget
each is given, 10 days at 6 h per cell and 40 at 24 h. And it cannot be bought back
with cores, because memory serialises it — on the plan's `n²` model a cell needs
about 17 GB at 1000 aa and 45 GB at 1629, against a 61 GB host on which E2 reserves
up to 53 GB.

**So the experiment moves to another host** (decided 2026-10-20). Every budget in this
plan was derived on tamazirt's E-cores under contention and none of them transfers:
the per-protein caps, the 12 h front budget, the concurrency proposal and the memory
guards are all re-derived on the new host before any cell of the middle or long band
runs. The three decisions above are about *what is measured* and stand independently
of where it runs.

**Concurrency is set from the calibration cell, not fixed in advance** (decided
2026-10-20). The earlier proposal of three concurrent cells outside E2's windows is
withdrawn rather than carried to the new host: it was sized against tamazirt's single
core and against a memory model whose `n²` fit was never checked on a constrained
run. The order is therefore fixed — **`P54852` under set I runs first**, the shortest
instance of the short band at 163 aa and the cheapest nest, and its measured wall and
`peak_rss_kb` set the concurrency for every cell after it.

That cell is also what resolves the width of Phase B's own bracket. Phase B sits
between a 14 h floor and a 468 h ceiling, a factor of 33, and the quantity that
separates them — the ratio of a constrained front's cost to a free one's on the same
instance — has never been measured on any instance at any length. The cell reports it
in under an hour against a baseline that already exists, so no concurrency number,
budget or band order is worth deriving before it has run.

### 3. The ladder — dropped

A truncation ladder on one protein was specified here: rungs at roughly 0.5×, 0.7×
and 1.0× native plus one above, under all four nested sets, fitted against `n`. **It
is dropped.** Its two jobs are reassigned:

- **the scaling fit** is already on record from the hull probe below, which ran those
  rungs under the exact bounded hull rather than the front. That probe is closed and
  its exponents stand; it was always the instrument that produced the numbers the
  plan quotes, the front ladder being the version that had not run.
- **the frontier length** moves to Phase B's own completeness column, where it is a
  property of the instances the paper actually reports rather than of prefixes.

Dropping it loses one thing and gains another, and the balance is why it goes. What
is lost is a length series at *fixed composition*: real proteins confound length with
composition, and prefixes do not. What is gained is that the frontier stops resting
on an instrument the plan itself flagged as weak — rungs are prefixes of one protein,
so a frontier read off them is a frontier for that protein's composition.

The consequence is that **nothing between the grid's ceiling and the shortest Cas
ortholog is covered by a real instance.** With the corpus capped at 500 aa that gap
runs 494–984 aa, and the ladder was what nominally spanned it. Extending the screen
above 500 aa is therefore no longer an enhancement but the only way to locate the
frontier on real proteins; the mirror holds 19 proteins at 501–700 aa and 13 more at
701–1000.

**Decided 2026-10-20: the screen extends over 501–700 aa.** The 19 mirror proteins in that band run as a second pass of
the same screen — baseline fronts, counts only, core 24 — at roughly 2–3 h each, so
40–55 h sequentially. The band above 700 aa is not run: each cell there costs several
hours more and needs the memory gate, and 501–700 already reaches into the gap.

The extension was first taken on terms that forbade it to change the instance
selection, so that no headline instance was conditional on a run that had not
happened. **The three-band restatement later the same day overrode those terms**, and
the middle band is now drawn from exactly these nineteen. The two are consistent in
what they protect: the short band was fixed before the screen reported and stays
fixed, while the middle band has no members until it reports. What the extension also
feeds is the frontier column and the choice of long-band candidates.

### 3a. The hull probe — pre-registered 2026-10-19, before it runs

This probe was registered as a substitute for the ladder, at a time when
`front --motif` did not exist. It is now the whole of the length series: the ladder
is dropped, and `front --motif` landing on `scope/38-motif-front` does not revive it.
The exact motif **hull** is the cheaper half of the front — the front's work is the
hull vertices **plus** the off-hull gap solves — so a hull measurement is a **strict
lower bound** on the front's cost at the same length, and the probe is therefore a
**one-sided** test of the long arm:

- if the fit puts 984 aa **outside** the budget, the front is certainly outside, and
  the Cas arm is settled as non-exact without the scope-38 build;
- if it puts 984 aa **inside**, the verdict stays open and only Phase B can close it.

**Instrument.** `motif_hull --fasta F --motifs FILE --arms walk --workers 1
--epsilon 0`, the frozen binary `motif-work/bin/tip` at commit `5d12df5f`, no UTRs —
the headline setting. Its motif file is read with `split_whitespace()` and supports
**no comments**, so the probe passes the stripped `motifs/probe-all-sets.bare.txt`
(25 patterns, all four sets, deduplicated) rather than the annotated set files.

**Rungs.** P31645 at `n` = 100, 200, 300, 400, then 600 if 400 lands inside the
budget. One rung at a time, `--workers 1`, pinned to the E-cores, under a memory cap
and a 1 h per-rung wall budget. A rung that hits either limit is recorded with the
limit it hit, and the fit uses the rungs that completed.

**Primary readout is peak RSS, not wall.** Memory is what binds on this host (52 GB
cap; `n²`, measured 5.7 GB at n = 300 and 24.1 GB at n = 600 at 6 workers), RSS is
load-immune, and memory extrapolates more reliably than time. The walls here are
E-core walls taken while another session builds and tests on those same cores, so
they are recorded as **load-qualified lower bounds** with the load at launch, and no
published frontier rests on them. Note also that the n = 600 figure on record came
from a **20-minute capped run**: it is a memory measurement, not a completed hull,
and completed walls exist only to about n = 300.

**Decision rule — amended 2026-10-19, after the first rungs.** The original rule
tested predicted peak RSS against the 52 GB cap. The rungs show that memory is not
the binding axis at `--workers 1` (≈15 GB projected at 984 aa, ≈40 GB at 1629 aa,
both inside the cap) while the wall is: the n = 300 rung exceeded an hour for the
hull alone, which lower-bounds the front.

**No instance is dropped for being expensive, and no verdict gates the long arm.**
Every instance in the grid runs, including all five Cas orthologs. A cell that cannot
complete within its budget **reports the partial front it reached and is labelled
incomplete**, carrying the budget it hit and the point count it had. Completeness is
a reported column, never a selection criterion.

The **frontier length** — the largest `n` whose cell completed — is therefore reported
as a measured quantity rather than used as a gate, with whichever axis bound first
named beside it. With the ladder dropped it comes from Phase B's completeness column
over real instances, and the hull probe's rungs below bound it from one side only.
The fitted exponents are reported for both axes, and the fit
is reported as what it is: the rungs so far are not a clean power law (100→150
implies 3.2, 150→200 implies 1.4), consistent with the walk's cost being set by its
hardest single probe rather than by a smooth function of length.

**The load explanation was tested and refuted.** The n = 300 overrun launched at load
11.9 against 3.9 for n = 200, so n = 200 was re-run as a control under comparable
load: 354.9 s at load 6.03 against 344.4 s at 3.92, **+3%**. The overrun is the
instance, not the host.

### 4. The long band's Cas orthologs (984–1629 aa)

Addressed by ID (`PRT-03`…`PRT-07`), as in E2, so no ortholog name reaches a row.
These five are part of the long band rather than the whole of it — the band takes ten
instances and the mirror supplies fifteen more above 1000 aa — and `PRT-03` at 984 aa
belongs to the middle band by length. **Their baseline front for `PRT-03` is already
on disk** in `e2/out/hull`, so the long band inherits one baseline it does not have to
recompute, subject to a check that the build that wrote it is compatible.
**The route is decided by the hull probe's fit, not assumed** — the front ladder that
was to refine it is dropped, so the probe's one-sided bound is what there is. What is
on record:

- the exact bounded motif **hull** walk has run at **n = 600** (P31645) at 24.1 GB
  peak RSS with 6 anchor workers; memory goes as `n²` and `--workers 1` restores
  serial memory, so 1629 aa extrapolates to roughly 45 GB serial — inside
  tamazirt's 52 GB cap, as an extrapolation;
- **no wall time is on record:** the n = 600 figure came from a 20-minute capped
  run, and at n = 300 one anchor scan is already ~66 s;
- the hull is the cheaper half. The whole front adds the off-hull gap solves, and
  **nothing has been measured for those under motif exclusion at any length.**

So: if the fit puts 984 aa inside the frontier, the Cas arm runs **exact**. If it
does not, it runs under a wall and memory budget and records the partial front as an
anytime curve, exactly as E2 does for a censored cell — a lower bound on cost and a
complete record of what was produced by the limit. Either way the published reason
is the fit, not a judgement call. Memory extrapolates more reliably than time, so
the memory prediction may be quoted and the time prediction must be checked.

### What the rungs measured: the fold, not the search

The hull probe's rungs (`e3-probe/`, `e3-probe2/`) separate cleanly into a smooth axis
and an erratic one, and the erratic one is localized:

| n | wall | peak RSS | vertices | probes | folds/probe | s/fold |
|---|---|---|---|---|---|---|
| 100 | 62.3 s | 0.173 GB | 26 | 51 | 2.51 | 0.35 |
| 150 | 231.9 s | 0.374 GB | 40 | 79 | 1.92 | 1.21 |
| 200 | 344.4 / 354.9 s | 0.668 GB | 49 | 97 | 1.52 | 1.82 |
| 250 | 2281.9 s | 1.020 GB | 62 | 123 | 1.50 | 11.47 |
| 300 | 5097.2 s | 1.509 GB | 68 | 135 | 1.58 | 22.85 |

Fitted over all five rungs: **RSS exponent 1.97 at r² = 0.9997**, projecting 15.4 GB at
984 aa and 41.4 GB at 1629 aa, both inside the 52 GB cap. **Wall exponent 3.97 at
r² = 0.939**, with consecutive-rung exponents 3.24, 1.48, 8.34, 4.41 — a poor power
law with a discontinuity between n = 200 and n = 250 rather than a smooth curve. The
wall projects to ≈ 117 h at 984 aa and ≈ 871 h at 1629 aa.

**The frontier, measured.** The exact motif hull completes within 1 h up to n ≈ 297
(calibrated: n = 300 took 85 min) and within a 12 h cell up to **n ≈ 485–554**,
the range spanning the global exponent and the top-of-range one.

**What that settles.** The hull is the cheaper half of the front, so the long arm's
cells at 984–1629 aa **cannot complete** — they will report partial fronts under the
completeness policy, and the paper can say so in advance rather than discovering it
per cell. The 6 selected instances are all ≤ 500 aa and so sit at or below the
frontier, which is where the price itself is measured.

- **Peak RSS is a clean `n²`** across every rung (fitted exponents 1.9–2.0).
- **The wall is not a power law**: consecutive rungs imply exponents 3.24, 1.47, 8.33.
- **The search is not degrading.** Folds per probe *falls* (2.51 → 1.50), probe count
  grows about linearly (≈ n/2) and vertices about n/4. The pruning counters rise.
- **The cost of one fold is what explodes**: 0.33 → 1.15 → 1.71 → 11.47 → 22.85 s for
  the ordinary folds, where `n³` predicts 1.95× for the 200 → 250 step and the
  measurement shows 6.7×. The step after it is 2.0×, so the jump is a discontinuity
  between 200 and 250 rather than a steeper exponent throughout.
- **It is not the known λ = 0.0001 critical path.** Stratified by λ, that single fold
  per rung costs 3.36 / 9.54 / 14.78 / 34.79 s — smooth, about `n^2.4`. The
  **ordinary** folds at λ ≥ 0.01 are the ones that jump: 0.33 / 1.15 / 1.71 / 11.47 /
  22.85 s against 3.36 / 9.54 / 16.92 / 34.79 / 61.99 s for the λ = 0.0001 fold.
- **Load was excluded as the cause** by the control described above.

This is a result about the fold internals, which the motif speed-up round owns; E3
records it and does not chase it. Two consequences for E3 itself: a frontier fitted
from rungs at or below n = 200 would have badly under-predicted the cost, and the
runtime price of a *longer blacklist* is a quantity this experiment should report,
since the discriminator between length and AC-axis density is one rung under the
8-pattern `p8` set against the same rung under these 25.

### Blacklist size costs a step change, and the frontier is per blacklist

The hull rungs were run again under the 8-pattern `p8` set — same binary, same cores, same
instance, 8 patterns instead of 25 (`e3-probe3/`, `e3-probe4/`). Every figure below is
regenerated by `python3 -I fit_probe.py`.

| n | `p8` wall | `p8` s/fold | 25-pattern wall | 25-pattern s/fold | ratio |
|---|---|---|---|---|---|
| 200 | 142.6 s | 0.45 | 354.9 s | 1.71 | 2.5× |
| 250 | 293.1 s | 0.82 | 2281.9 s | 11.47 | 7.8× |
| 300 | 722.8 s | 2.92 | 5097.2 s | 22.85 | 7.1× |
| 400 | 1577.7 s | 5.00 | — | — | — |

| | wall exponent (r²) | RSS exponent (r²) | 12 h frontier | wall at 984 aa |
|---|---|---|---|---|
| `p8` | 3.55 (0.987) | 2.04 (0.9993) | n ≈ 995 | ≈ 11.5 h |
| four sets | 3.97 (0.939) | 1.97 (0.9997) | n ≈ 554 | ≈ 117 h |

**Both blacklists scale as roughly `n^3.5`–`n^4`.** The difference between them is not
the exponent: `p8`'s power law fits at r² = 0.987, while the 25-pattern fit is poor at
0.939 *because of a step change between n = 200 and n = 250* — its per-fold cost jumps
6.7× over that step where `p8`'s rises 1.8×. Above the step the two curves run roughly
parallel, 7–8× apart. The earlier reading of these rungs as "`p8` smooth, dense set
super-cubic" was wrong and is corrected here: it is a step, then parallel curves.

**Peak RSS is blind to the blacklist**: 1.0206 GB against 1.0211 GB at n = 250, and the
RSS exponent is ≈ 2 under both. Memory does not see how much is forbidden; time does.

Three consequences for E3:

- **The runtime price of exclusion is a function of how much you forbid**, and the
  nested grid produces that curve for free, `F1…F4` being four blacklists of increasing
  size on the same instances. It belongs in the paper beside the objective price.
- **The frontier is per blacklist, so every reported frontier names its set.** Under all
  four sets the hull completes a 12 h cell up to n ≈ 554; under `p8`, up to n ≈ 995.
- **The shortest Cas ortholog is a near miss, not a rout.** At 984 aa the `p8` *hull*
  projects to ≈ 11.5 h against a 12 h cell, and the front adds the gap solves on top, so
  a complete front at that length is out of reach even under the lightest blacklist —
  but by a factor of a few rather than of a hundred.

**Load was excluded twice, not assumed.** The n = 200 control put load's cost at +3%
(354.9 s at load 6.0 against 344.4 s at 3.9), and the n = 400 `p8` rung launched at load
22.5 — three times any other rung — and still landed below its predecessor's implied
step. Walls here are E-core walls and remain qualified as such, but load does not
explain any of the structure above.

### What reporting a partial front requires of the build

Measured on the censored rung: **`motif_hull` is not anytime.** It buffers its output
and prints after the walk completes, so the n = 300 rung killed at 3600 s produced
**zero** lines — no partial hull, no vertices, nothing to report. The budget also has
to be enforced *inside* the measurement wrapper, or the censored cell loses even its
own wall and RSS, which is what happened on that rung.

That made the completeness policy a requirement on `front --motif` before Phase B
could honour it: stream or checkpoint the points found so far, stop on its own budget
rather than being killed from outside, and say which gaps were left unsolved — or an
incomplete cell is indistinguishable from a failed one.

**Checked against `origin/scope/38-motif-front` on 2026-10-19 (after `eda42cf0`,
`6222df86`): all three are built, so Phase B is unblocked.** `front` carries:

| flag | what it gives the policy |
|---|---|
| `--max-secs S`, `--max-rss-gb G` | the walk stops itself and returns the partial front |
| `--checkpoint` | rewrites `--out` with the front so far after each settled gap |
| `--resume FRONT.jsonl` | walks only the gaps a budgeted front left open, returning its points with theirs |

and the result record carries `complete` (true when no gap is left open) and `gaps`,
with the unfinished gaps listed per box. So a censored cell reports its own partial
front and its own incompleteness, and `--resume` means a cell stopped by an E2 window
is continued rather than restarted. E3 therefore runs every front cell with
`--max-secs`, `--max-rss-gb` and `--checkpoint` set, and no cell is ever killed from
outside the process.

The probe's own measurement stands: it used `motif_hull`, which has none of this, and
its censored rung is the evidence that the flags are load-bearing rather than
decorative.

### What the branch owner has in flight, and what it costs us (2026-10-20)

Four changes are local and unpushed on `scope/38-motif-front`. `--resume`'s format is
unchanged in all of them, and `complete` keeps its meaning. Two are reporting fixes we
need and one makes the pin stricter than "any commit with the flags".

| change | consequence for E3 |
|---|---|
| `stopped_by` is null whenever nothing is left open (before, an idle worker could set `"memory"` on a front that had completed) | **read `complete`, never `stopped_by`.** A pre-fix commit can label a complete cell as memory-stopped, and that is exactly the false cost claim the screen's `censored()` bug would have made |
| a checkpoint can no longer list an already-settled box as open | a pre-fix commit **overstates `gaps`**, so the incompleteness we report would be worse than the truth. No front point changes either way |
| `--checkpoint` with a non-regular `--out` is refused at parse | no effect: every cell writes a regular file |
| a separate branch settles the highest-tie-weight box first | **the one that binds the pin.** Same front when complete, *different* partial fronts under a budget — so two cells of the grid built from different commits are not comparable as partials, and neither is a `--resume` of a cell against its own earlier run |

So the pin is not "a commit that has the flags" but **one commit for the entire grid**,
including every resume and every re-run of a cell, carried in each cell's provenance.
And we want a commit *after* the first two fixes, because both move numbers that go in
the paper's incompleteness column rather than its front.

**The pin is `60ddd0f9`, the `xg` tip as of 2026-10-20** (decided that day, once the
stack merged). The whole of it is now in `xg` and reachable from that one commit:
#429 motif front (`b3adb049`), #431 `--pin` (`60ddd0f9`), #432 the plain walk's
budget and resume (`f500f44a`), #433 (`473b6516`), and the two perf merges #434
jemalloc (`e523c001`) and #436 capfit (`f584a090`).

The earlier form of this pin was the gated `scope/38-motif-front` head, on the
grounds that #431 and #432 were being retargeted to `xg` behind it, so an `xg`
commit's contents depended on when it was read while a branch head could not acquire
them. **The merges retire that reasoning rather than override it**: a named `xg`
commit is immutable once the branches merging into it have merged, and it is now the
branch *head* that moves. The pin is a SHA either way, which is what the grid needs.

Two consequences of pinning the tip rather than #429 alone:

- **#429 alone could not run the long band.** The budget and checkpoint surface
  arrived with #432, not with the motif front — `front.rs` carries 34 references to
  it at `b3adb049` against 55 at `f500f44a` — and every long-band cell is a budgeted
  partial.
- **#436 changes the cap fit, so it changes what a budgeted partial reaches.** That
  is a reason to pin it, not to avoid it: the comparability the grid needs is across
  its own cells, which one SHA delivers, and the alternative gave up #434's 1.41×
  across a band budgeted at 240–960 h.

One unit question is settled and does not apply. #431 makes `Dfa::codon_offset`
non-zero, and the walk had been reporting vertices and open boxes in lattice rather
than front units — but that conversion is reached only under `--pin`. Unpinned,
`codon_offset` is exactly zero and the walk adds that zero to points, vertices and
boxes alike; #431's unpinned front on `P0DTC2[:100]` is point-for-point identical to
the front branch's, and pre-#431 code refuses `--pin` with `--motif` outright. **No
E3 phase uses `--pin`** — the flank work fixes flank *sequence*, which is a different
mechanism — so E3's `gaps` carry no unit caveat, before or after that merge.

### `--max-rss-gb` is a stop, not a ceiling

It is read only between folds, and the two end folds run before the first check. At
n = 500 they peak at 4.7 GB side by side; at n = 1000 a 6 GB setting was passed inside
two minutes, reaching 10.7 GB. Our instances are all ≤ 500 aa, so the unchecked peak is
the 4.7 GB one — but it means the flag cannot be the only guard.

Every front cell therefore runs under a hard cgroup cap as well: `systemd-run --user
--scope -p MemoryMax=12G -p MemorySwapMax=0` with `--max-rss-gb 8`. The two numbers are
set rather than described, so that a later reader can check them: 8 GB leaves the
graceful stop clear of the 4.7 GB unchecked peak at n = 500, and 12 GB leaves the hard
cap clear of the stop. Both are re-derived, not reused, if an instance over 500 aa ever
enters the grid — the peak grows with length, and at n = 1000 it is already 10.7 GB. The two guards fail differently and that is the point:
the flag stops the walk and returns a partial front with its `gaps`, while the cap
kills it outright. `--checkpoint` is what makes the second survivable — the last
settled gap is already on disk, so even an OOM kill leaves a partial front to report,
and the cell resumes from it. The rewrite is atomic (written beside `--out`, synced,
renamed, one checkpoint at a time), so a kill during a rewrite leaves the previous
complete checkpoint rather than a truncated line.

A cell killed by the cap is reported as killed, not as budget-stopped, and the two are
told apart by the header's `checkpoint` field rather than by whether a file exists:
`false` means the run ended itself, whether complete or budget-stopped, and `true`
means it did not. Only a kill before the first settled task leaves no file at all. The
reading has one precondition — `true` also describes a cell that is *still running*, so
a cell is classified only once its process is gone, never from the file alone.

## Flanks

**Decided 2026-10-19: the headline grid runs without UTRs**, to isolate the motif
axis. The paper's other contribution already measures what fixed flanks cost, and a
cell carrying both constraints cannot attribute a front shift to either.

One screened instance is then repeated with the `hbb` UTR pair and the region
constraint, as a cross-check that the price of motif-freedom does not depend on the
flank setting. If it does, that is a finding and gets its own row. The cross-check
runs the junction pre-check described above first.

**Which instance, fixed 2026-10-20 before any price is visible: `P47211`, 349 aa.**
The rule is the **median-length instance of the short band's six all-four binders**,
taking the lower of the two middle instances ascending in length — the six being
P26436 (265), Q9HC38 (313), P47211 (349), P15086 (417), Q9NWZ3 (460) and Q9UBM8
(478), so the rule selects the third. Stating the rule now is the point of it: the
cross-check is the single row that licenses quoting every headline price without a
flank qualifier, and an instance chosen once the prices are known cannot carry that
weight, whatever it then shows.

Two properties of the rule, each measured rather than assumed:

- **It runs over the all-four binders, not over all eighteen.** A cross-check on an
  instance that reads `occurs = 0` in `initiation` yields flank-vs-no-flank rows at
  `F1`–`F3` only, and `F4` is the one nest the flank setting can interact with at
  all, `RCCATGG` being the sole entry able to straddle the 5′ junction. The rule over
  all eighteen would have selected `Q6PF06` at 316 aa, whose `initiation` occurrence
  count is zero, so the junction pre-check would stand in front of a cell that
  carries no price.
- **Little representativeness is given up for it.** The eighteen-instance band's own
  median is 326.5 aa — an exact tie, 316 and 337 lying 10.5 aa either side, which is
  a second reason the wider pool needed a tiebreak of its own. At 349 aa `P47211` sits
  23 aa above that median and inside the band's middle third.

`P47211` also carries adjacent PP, so it is already a member of the poly-proline row,
and the two cells share an instance without sharing a claim.

## Metrics, per cell per front

Field names are the ones the result JSONL already carries, so nothing is re-derived:
a header line records `motifs`, `motif_window`, `motifs_in_fixed_flanks`, `wall_s`
and `peak_rss_kb` (the process high-water mark `VmHWM`, a true peak), and each point
line carries `sequence`, `structure`, `mfe` (the energy **under the constraint** —
the objective), `cai`, `cai_log` and `on_hull`.

**Front geometry** (the price):

- `points`, `vertices` — front points and hull vertices (`on_hull`), and how many of
  `F_{k-1}`'s points are gone from `F_k`;
- `dMFE_at_CAI` — the vertical gap between consecutive fronts: at each CAI attained
  on `F_{k-1}`, the best MFE still reachable on `F_k`. Min, median, max over the
  front, in kcal/mol;
- `dCAI_at_MFE` — the same read horizontally;
- `dobj_lambda` — the scalarized gap at each λ of `F_{k-1}`'s hull, the quantity the
  monotonicity gate checks;
- `area` — the region between consecutive fronts, as one summary number.

**Binding** (what makes the price meaningful): `occurs` and `removed` per set and per
entry, carried on every row. A cell with `removed = 0` is labelled **non-binding**
and its `dMFE` of zero is reported as such, not averaged into a price. Per AGENTS.md,
a non-binding cell is worthless to a *price* and is a genuine result for a *runtime*,
so the runtime table keeps it and the price table excludes it, on the same column.

**Cost**: `wall_s`, CPU user+system, `peak_rss_kb`, instructions retired, exit
category, and the ratio of each to the baseline front on the same instance.

**Completeness**: every row carries `complete`, and an incomplete row carries the
budget it hit and the points it had reached. A price computed against an incomplete
front has a **direction**, and it is reported with it: the true front has at least
the points found and each missing point is weakly closer to `F_{k-1}`, so when
`F_{k-1}` is complete and `F_k` is not, the measured shift is an **upper bound** on
the price. When both are incomplete, no direction is claimed and the row carries the
gap counts instead of a bound.

**Feasibility**: an instance made infeasible by a set is recorded as `infeasible`
with the set that did it. That is an outcome, not a failed run.

## Soundness gates

1. **Monotonicity** — `obj_λ` non-decreasing along `F0…F4` at every λ, every cell.
2. **Independent motif-freedom certification** — every emitted design is scanned for
   every motif of its own set by a checker sharing no code with the DP. Zero hits
   required. A design on `F_k` is also scanned for sets `1…k−1`.
3. **Empty-set pass-through** — `front` with an empty motif set must reproduce the
   baseline front byte-identically (the Scope 17 invariant, re-checked here).
4. **Non-vacuity by mutation** — the oracle test is proved non-vacuous by mutating
   the exclusion code and watching it fail, **and** by confirming the mutation
   applied: a `sed` once silently matched nothing after `cargo fmt` moved the
   indentation, and the test "passed".
5. **Re-score** — every front point's MFE and CAI are recomputed by `derna check` /
   `mfe_fold`, independently of the walk that produced it.
6. **Baseline binds** — gate 2 is correctness; this one is meaning. A set whose
   `removed` is zero across every cell is reported as vacuous on this corpus, and no
   price is quoted against it.
7. **Row count before reading** — the merged output's row count is checked against
   `instances × fronts` before any number is read from it, and no `_*` glob is run
   inside the experiment directory while a sweep is live (AGENTS.md).

## Reporting

**Decided 2026-10-19: the vacuous-entry table is a paper row**, not supporting
material. It publishes, per entry, which standard blacklist entries have no spellable
placement and which are spellable but never appear on the baseline front. It costs
nothing beyond the screen, it is the evidence that justifies every price quoted
beside it, and it is a claim no comparator paper can make, since none of those tools
can express the constraint.

**Decided 2026-10-20: vacuity is asserted from the screen, not demonstrated.** The
confirmation runs that would have shown an identical front under and without a
non-occurring entry are not run. `occurs = 0` is the stronger statement already: an
entry absent from every design on the baseline front cannot move that front, so the
table rests on the screen's integers rather than on a pair of fronts agreeing. The
two mechanisms it separates stay separate in the reporting, since they are different
claims — `processing` is spellable on 55 of 59 and occurs on 18, the lattice
admitting what the optimum avoids, while `initiation` occurs on 22 with only 7
spellable-but-absent, the lattice forbidding it outright.

Everywhere a price is reported, the sentence that qualifies it travels with it: the
flank setting, the corpus, and the fact that a non-binding cell is excluded from the
price and kept in the runtime.

The twelve short-band instances that bind only I–III read *not binding* in set IV's
column, carrying the screen's
`occurs = 0` beside it, which is the stronger statement: the entry cannot move that
protein's front, so the absence of a price is a property of the instance rather than
a measurement that came back small. A 0.00 in that cell would be unreadable against
a genuine zero price elsewhere in the same table.

## Phases, and the dependency

**`front --motif` is on `xg`** (PR #429, merged), with `--motif-window`, the budget
flags, `--resume` and `--checkpoint` described above, and the pin resolved to
`60ddd0f9` as set out under *The pin* — so Phase B builds from `xg` at that commit
and records it in every cell's provenance.

**This was Phase B's one blocking dependency and it is discharged.** The phase
previously waited on another session's unmerged branch for a gated head SHA, and
nothing of it could start without that handover. What remains before the first cell
is not a dependency on anyone: the host move, the three cost-and-memory runs that
re-derive every budget there, and the `P54852` calibration cell that sets the
concurrency.

- **Phase A′ — the binding screen.** Baseline fronts only, so it runs today, and the
  selection rule resolves before the dependency lands. Steps: write the screen
  harness; run `F0` over 61 proteins; the independent substring scan; the
  `spellable` census; publish the table including non-binders.

  **Decided 2026-10-19: A′ runs under contention, and publishes counts only.** The
  E-cores it is confined to carry another session's build, so its walls are not
  reportable — but `occurs`, `spellable` and `pp_adjacent` are exact integers that
  load cannot perturb, and they are the whole of what the selection rule consumes. No
  runtime number from A′ enters any table; E3's runtime claims come from the probe
  and from Phase B, which are budgeted and load-qualified. **All 61 proteins run at
  native length**, with no length cap: a cap would drop exactly the proteins most
  likely to spell Type IIS sites, and a prefix would make the screen a measurement of
  the prefix.

  That decision was taken on the premise that binding is driven by composition rather
  than length. **The screen refutes the premise and strengthens the decision.** Of 59
  proteins with a complete `F0`, those whose front carries all four sets number 0 of
  20 at ≤ 211 aa, 2 of 20 at 212–343, and 3 of 19 above 343 — mean sets binding 2.30,
  2.75, 2.79. Length is the driver. Composition is not: the all-four pool's GC
  propensity median is **0.6547**, *below* the 0.6680 of the proteins it was selected
  out of. So a length cap would have removed the only band that binds the full nest at
  all, which is a stronger argument against capping than the one originally given —
  and the stratification the selection rule was built on does not survive it. See the
  selection rule below.

  **Re-capped 2026-10-20, from 3600 s to 10800 s per protein, on the measured tail.**
  The original cap was set against a fit over all completed cells, which gives
  `wall ≈ 3.2e-4 · n^2.58` and projects every remaining protein inside an hour. Fit
  instead on the 17 cells at or above 300 aa — the only ones that say anything about
  the tail — and the exponent is `3.88`: 494 aa projects to 70 min centrally, not 48,
  and eight of the eleven remaining proteins put their +1 sd band over 3600 s. Short
  proteins dominate a whole-corpus regression and flatten exactly the region the cap
  acts on, which is the under-scoping direction of *measure before capping*.

  Raising it costs the screen nothing it reports: no wall from A′ enters a table, so
  the cap is only a guard against a pathological run, and at 10800 s it still sits at
  2.2× the worst remaining +1 sd projection. What it buys is the difference between a
  protein that contributes `occurs` counts and one that contributes nothing: plain
  `front` has no `--checkpoint`, so a protein cut by `timeout` leaves no partial front
  to report, only an exit code. The grid would have lost up to eight of its longest,
  GC-richest proteins — the ones most likely to spell Type IIS sites, and so the ones
  the set IV bottleneck most depends on.
- **Phase A — the per-λ price. Derived from Phase B, not run.**

  **Decided 2026-10-20: Phase A stops being a measurement phase.** The per-λ price is
  a table computed from the fronts Phase B and Phase A′ already produce, at no compute
  cost, and the 2,108 `design --motif` calls it was specified as are not run. What the
  phase reports does not change: Δobjective, ΔMFE and ΔCAI at each hull λ of the
  baseline front, per nested set.

  The reason is that **a per-λ optimum is a lookup on a front rather than a solve.**
  Weighted-sum scalarization can only return a vertex of the convex hull of the
  feasible set, so the optimum at λ is the point of the front minimizing
  `λ·mfe − (1−λ)·cai_log`, and the price at λ is the difference between that argmin
  over the constrained front and over the baseline front. Verified on the screen's own
  output before the decision was taken: over P26436's 796-point baseline front the
  argmin landed on an `on_hull` vertex at **99 of 99** λs tested across (0, 1). Both
  inputs are in hand or will be — A′ emitted `F0` for all 61 proteins, and Phase B
  emits the constrained front per cell.

  So the earlier reasoning that the two phases measure different quantities was right
  about the quantities and wrong about what it implied. A per-λ price does hold λ fixed
  at the baseline's vertices while a constrained walk returns its own, and those are
  different numbers — but the walk's output *contains* the per-λ answer for every λ,
  so computing one from the other needs arithmetic rather than a second experiment.

  **Two things the lookup cannot give, and what happens to each:**

  - **An incomplete front breaks it.** The lookup is exact only where Phase B's walk
    finished: a budget-stopped front may not hold the vertex optimal at some λ, and
    the argmin would then return a suboptimal point *silently*. So the derived table
    is emitted **only for cells whose checkpoint header reports `complete`**, which is
    every grid cell on the measured frontier (n ≈ 554 against a longest instance of
    478 aa) and no dense Cas cell. Where the front is partial the row is withheld
    rather than approximated; running Phase A there instead is not an option either,
    at `n^1.97` a single call at 1629 aa being about 470 s.
  - **The per-λ wall and RSS ratio is not recoverable** — it is the cost of one
    constrained solve against one unconstrained solve, where Phase B only ever reports
    the cost of a whole walk. That column is dropped from the phase and, if the claim
    is wanted, comes from a handful of `design` calls on two or three instances, which
    is a runtime probe rather than a phase.

  What the decision gives up is the one thing Phase A did that Phase B cannot: it runs
  without the unmerged branch, so it would have put a prior on Phase B's numbers before
  Phase B could start. That is project risk rather than a result, and the measured cost
  makes the trade a poor one — 19.2 h over the grid to preview cells budgeted at 12 h
  each. The costs, for the record: the call count is **`4 nests × hull λs` summed over
  proteins**, which was 2,108 over the six-instance grid this rule replaced and 16,508
  over the whole corpus; a single `design` call takes
  14.2 s at 265 aa and 45.5 s at 478 aa, an exponent of 1.97, putting the grid at
  19.2 h and the corpus at 112 h; and `--lambda L,L,…` does not amortize, eight λs
  costing 7.8× one, so there was no saving to be had from batching. An earlier
  estimate of 9–13 h was self-inconsistent at the 0.66 s per front point it quoted —
  0.66 s over 15,496 calls is 2.8 h — and no corpus-wide average per call can hold
  when the cost goes as `n^1.97`.

  **Consequence for breadth, which is why the open decision about it closes too.** The
  wide version over all 61 proteins was attractive as a corpus-scale result robust to
  the grid being wrong, and at 112 h it was not affordable. Derived, it costs nothing,
  so the table is emitted for **every protein whose baseline and constrained fronts are
  both complete** — which the short band's cells are and the long band's never are, the
  middle band's depending on where its instances fall against the n ≈ 554 frontier. The breadth question is therefore no longer a
  trade-off against compute.
- **Phase B — the nested fronts.** `front --motif` × {I, I∪II, I∪II∪III, I∪…∪IV} ×
  the three bands of the selection rule: 60 cells in the short band, 40 in the middle
  and 40 in the long, the long band's being budgeted and never complete. The ladder
  that sat between them is dropped.

  **Decided 2026-10-19: 12 h per front cell**, as `--max-secs 43200`, with
  `--max-rss-gb` set from the headroom at launch and `--checkpoint` always on. The
  figure is the measured frontier, not a guess: 12 h reaches n ≈ 554 under all four
  sets and n ≈ 995 under `p8`, so the six screened instances complete under every
  set, `p8` very nearly completes at 984 aa, and only the dense Cas cells report
  partial fronts. Doubling to 24 h was rejected on the measurement — at `n^3.5–4` it
  moves the frontier about 20% and collides with two E2 windows. Those frontiers were
  themselves measured on E-core walls at load 4–22, so they already price contention
  rather than assuming a quiet host.
- **Phase C — the limitation rows.** PP-infeasibility on every short-band instance
  carrying adjacent PP, the vacuous-entry table, and the flank cross-check with its
  pre-check. Seven of the eighteen carry adjacent PP — P30281, P47211, P01860,
  Q92734, Q9NWZ3, Q6B0B8 and Q9UBM8, three of them among the six that bind all four
  sets and none shorter than 292 aa. **No
  vacuity row is drawn from the long band**, where `occurs = 0` on a partial front is
  unresolved rather than vacuous.

A′ and C's first row do not depend on B and are ordered first. Phase A now does
depend on B, being derived from it, and costs no compute of its own. Nothing in the
plan waits on a run it does not consume.

## Host

tamazirt, shared with the E2 grid and the motif walk speed work. E2 holds the
P-cores (0–14 even, siblings idle); **E3 screens and builds on the E-cores (16–27)**
and takes P-cores only when E2 is not timing. Memory: `systemd-run --user --scope -p
MemoryMax=… -p MemorySwapMax=0`, and `--workers 1` for any motif walk running beside
another walk — each anchor in flight holds a whole scan's tables (~3.7× serial RSS).

**Three windows are reserved for the E2 timing run.** Its longest walk is projected at
53.4 GB on a 61 GB host, so anything of E3's over ~1 GB resident would push it into the
host OOM killer rather than its own cap. Hold off during, local time, as earliest start
→ latest end:

| instance | window |
|---|---|
| PRT-06, 1368 aa | 06:00 Oct 20 → 15:00 Oct 20 |
| PRT-07, 1629 aa | 13:15 Oct 20 → 03:00 Oct 21 (the ~53 GB one) |
| P0DTC2_p1273 | 07:45 Oct 21 → 16:45 Oct 21 |

Nothing large runs before 06:00 Oct 20. These are the third set of times given for
the same three walks — they have moved both later and earlier as the grid was
rebuilt — which is the standing argument for gating on the check below rather than
on any of them.

E2's rescoring folds have moved off the walk path and now wait for its first pass
to end (~Wed afternoon), so **CPUs 16–27 carry no E2 load until then**. That makes a
quiet-host wall available on the E-cores for the first time this week, which matters
for Phase A, where walls are a reported quantity. It does not change Phase A′: its
first pass ran under contention, and a table mixing loaded and quiet walls would be
worse than one uniformly loaded, so A′ stays counts-only as decided.

The projections carry about ±20%, so the gate is the check rather than the clock:
`pgrep -af 'scope35_walk_front.*--n (1273|1368|1629)($| )'` — if it prints anything,
wait for it to exit before launching anything large. The alternation needs the
end-or-space anchor; a trailing literal space misses the match when the length is the
last argument.

Wall times from a loaded host qualify nothing. Instruction counts survive load, but
`perf_event_paranoid=4` here, so `icount.py` falls back to wall/CPU/RSS and the
primary cost measure is a quiet-host wall with E2's contamination checks, or
instructions once counters are available. Before launching any grid, check that
nothing else is already running it (`ps -eo args | grep …`).

## Running rule

Keep going until at least the time named at launch. Before ending any turn, run
`date`; if it is earlier, start the next step or extend the plan. Do not wait on a
run while any independent item remains — A′, the oracle work and the derived tables'
harness are independent of each other. Stop early only for a decision that belongs to the user,
for anything destructive or outward-facing, or when no meaningful task remains.
After that time: finish the step, commit, record results here, leave a handoff brief.

Extension steps, in order: the screen above 500 aa, which is what locates the
frontier now that the ladder is dropped; the remaining therapeutic natives from the
screen; the flank cross-check on a second instance; the frame-aware exclusion axis,
which would turn out-of-frame AUG from a limitation into a set.
