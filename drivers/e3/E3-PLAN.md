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

Short instances do not carry the motifs; long ones may not finish. Binding is driven
by **composition**, not length — a GC-rich protein spells Type IIS sites, an AU-rich
one spells poly(A) signals and U-runs — so the screen reads composition, and length
is handled separately by the ladder.

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
- `pp_adjacent` — whether the protein contains adjacent PP, which selects the two
  instances for the poly-proline row.

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

### 2. Selection rule (fixed now)

From the screen, keep proteins where **all four sets bind on the baseline front**
(`occurs ≥ 1` for each of I–IV). Stratify the survivors by length tertile and by GC
propensity (the fraction of residues whose codon menu admits a G or C at the wobble
position), and take **6**: two per length tertile, within a tertile preferring the
extreme GC-propensity cells. Ties broken by UniProt accession, ascending.

If fewer than 6 survive, the requirement drops to sets I–III binding and set IV is
reported on whatever subset binds — never by substituting a different protein after
seeing its numbers.

### 3. The ladder (3–4 rungs, one protein)

Truncation rungs on **one** screened protein — the median-length survivor — at
roughly 0.5×, 0.7× and 1.0× native plus one rung above it, under all four nested
sets. This is the scaling probe: wall, peak RSS and instructions fitted against `n`.

Its job is to say **where the exact nested-front route stops finishing**, which is
the number a reviewer will ask for and the only honest basis for the long arm below.
Rungs are prefixes, so they share composition — acceptable here because the ladder
measures *cost*, and the Cas arm is what separates length from identity.

### 3a. The hull probe — pre-registered 2026-10-19, before it runs

The ladder as specified needs `front --motif`, which does not exist yet. The exact
motif **hull** does exist, and the front's work is the hull vertices **plus** the
off-hull gap solves, so a hull measurement is a **strict lower bound** on the front's
cost at the same length. That makes a hull-only ladder a **one-sided** test of the
long arm, available today:

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

The ladder therefore reports a **frontier length** — the largest `n` whose cell
completed — as a measured quantity rather than as a gate, with whichever axis bound
first named beside it. The fitted exponents are reported for both axes, and the fit
is reported as what it is: the rungs so far are not a clean power law (100→150
implies 3.2, 150→200 implies 1.4), consistent with the walk's cost being set by its
hardest single probe rather than by a smooth function of length.

**The load explanation was tested and refuted.** The n = 300 overrun launched at load
11.9 against 3.9 for n = 200, so n = 200 was re-run as a control under comparable
load: 354.9 s at load 6.03 against 344.4 s at 3.92, **+3%**. The overrun is the
instance, not the host.

### 4. The long arm (the 5 Cas orthologs, 984–1629 aa)

Addressed by ID (`PRT-03`…`PRT-07`), as in E2, so no ortholog name reaches a row.
**The route is decided by the ladder's fit, not assumed.** What is on record:

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
| 250 | 2281.9 s | 1.020 GB | 62 | 123 | 1.50 | 11.60 |
| 300 | > 3600 s | — | — | — | — | — |

- **Peak RSS is a clean `n²`** across every rung (fitted exponents 1.9–2.0).
- **The wall is not a power law**: consecutive rungs imply exponents 3.24, 1.47, 8.33.
- **The search is not degrading.** Folds per probe *falls* (2.51 → 1.50), probe count
  grows about linearly (≈ n/2) and vertices about n/4. The pruning counters rise.
- **The cost of one fold is what explodes**: 0.35 → 1.21 → 1.82 → 11.60 s, where `n³`
  predicts 1.95× for the last step and the measurement shows 6.4×.
- **It is not the known λ = 0.0001 critical path.** Stratified by λ, that single fold
  per rung costs 3.36 / 9.54 / 14.78 / 34.79 s — smooth, about `n^2.4`. The
  **ordinary** folds at λ ≥ 0.01 are the ones that jump: 0.33 / 1.15 / 1.73 / 11.47 s.
- **Load was excluded as the cause** by the control described above.

This is a result about the fold internals, which the motif speed-up round owns; E3
records it and does not chase it. Two consequences for E3 itself: a frontier fitted
from rungs at or below n = 200 would have badly under-predicted the cost, and the
runtime price of a *longer blacklist* is a quantity this experiment should report,
since the discriminator between length and AC-axis density is one rung under the
8-pattern `p8` set against the same rung under these 25.

### What reporting a partial front requires of the build

Measured on the censored rung: **`motif_hull` is not anytime.** It buffers its output
and prints after the walk completes, so the n = 300 rung killed at 3600 s produced
**zero** lines — no partial hull, no vertices, nothing to report. The budget also has
to be enforced *inside* the measurement wrapper, or the censored cell loses even its
own wall and RSS, which is what happened on that rung.

So the completeness policy is a requirement on `front --motif` before Phase B can
honour it:

- stream points as they are found, as the plain walk's `--stream` already does;
- accept `--max-secs` and `--max-rss-gb`, as `front --structurally-optimal` already
  does, so the command stops itself rather than being killed from outside;
- report, on stopping, which gaps were left unsolved, so an incomplete front is
  distinguishable from a complete one with few points.

Without those, an incomplete cell is indistinguishable from a failed one. This is the
one thing E3 needs from the scope-38 build beyond `--motif` itself.

## Flanks

**Decided 2026-10-19: the headline grid runs without UTRs**, to isolate the motif
axis. The paper's other contribution already measures what fixed flanks cost, and a
cell carrying both constraints cannot attribute a front shift to either.

One screened instance is then repeated with the `hbb` UTR pair and the region
constraint, as a cross-check that the price of motif-freedom does not depend on the
flank setting. If it does, that is a finding and gets its own row. The cross-check
runs the junction pre-check described above first.

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

Everywhere a price is reported, the sentence that qualifies it travels with it: the
flank setting, the corpus, and the fact that a non-binding cell is excluded from the
price and kept in the runtime.

## Phases, and the dependency

**`derna front` has no `--motif` today.** `front --help` lists only
`--structurally-optimal`; `--motif` / `--motif-window` are `design`-only, and the
guide says why: "the front walk runs over the codon lattice, and the motif lattice
has only its hull so far." The motif-excluded exact front is being built on
`scope/38-motif-front` (worktree `motif-front-wt`), per
`experiments/scope-38-derna-cli/MOTIF-FRONT-NEXT.md`. E3's headline grid consumes it.

- **Phase A′ — the binding screen.** Baseline fronts only, so it runs today, and the
  selection rule resolves before the dependency lands. Steps: write the screen
  harness; run `F0` over 61 proteins; the independent substring scan; the
  `spellable` census; publish the table including non-binders.
- **Phase A — the per-λ price. Needs no new code; runs today.** `design --motif` at
  the hull λs of each baseline front, over the screen's survivors: Δobjective, ΔMFE,
  ΔCAI, wall and RSS ratio per λ, per nested set. A real row of the paper on its own
  (the price at a *given* trade-off), and it de-risks Phase B: the λ-wise gaps bound
  what the front-wise gaps can be, so Phase B's numbers have a prior before they are
  computed.
- **Phase B — the nested fronts.** `front --motif` × {I, I∪II, I∪II∪III, I∪…∪IV} ×
  the 6 selected instances; then the ladder; then the long arm on the fit's verdict.
- **Phase C — the limitation rows.** PP-infeasibility on the two screen-selected
  proteins, the vacuous-entry table, and the flank cross-check with its pre-check.

Phases A′, A and C's first row do not depend on B and are ordered first. Nothing in
the plan waits on a run it does not consume.

## Host

tamazirt, shared with the E2 grid and the motif walk speed work. E2 holds the
P-cores (0–14 even, siblings idle); **E3 screens and builds on the E-cores (16–27)**
and takes P-cores only when E2 is not timing. Memory: `systemd-run --user --scope -p
MemoryMax=… -p MemorySwapMax=0`, and `--workers 1` for any motif walk running beside
another walk — each anchor in flight holds a whole scan's tables (~3.7× serial RSS).

Wall times from a loaded host qualify nothing. Instruction counts survive load, but
`perf_event_paranoid=4` here, so `icount.py` falls back to wall/CPU/RSS and the
primary cost measure is a quiet-host wall with E2's contamination checks, or
instructions once counters are available. Before launching any grid, check that
nothing else is already running it (`ps -eo args | grep …`).

## Running rule

Keep going until at least the time named at launch. Before ending any turn, run
`date`; if it is earlier, start the next step or extend the plan. Do not wait on a
run while any independent item remains — Phase A, A′ and the oracle work are
independent of each other. Stop early only for a decision that belongs to the user,
for anything destructive or outward-facing, or when no meaningful task remains.
After that time: finish the step, commit, record results here, leave a handoff brief.

Extension steps, in order: the deferred ladder rungs; the remaining therapeutic
natives from the screen; the flank cross-check on a second instance; the frame-aware
exclusion axis, which would turn out-of-frame AUG from a limitation into a set.
