# The paper's figures

One plan for every figure, written before any of them is drawn. Each entry names the
panel, the quantity, the file the panel reads and the state of that file today, so a
figure that cannot be built yet is visibly waiting on data rather than on a decision.

Three conventions are fixed here and are not re-argued per figure.

- **Main text is front-led.** Figure 1 and E2–E5 carry main-text figures; E1 contributes
  a correctness-and-cost table to the main text and its scaling panels to the supplement.
  The headline is the exact whole front and exact constrained design, and E1's single-λ
  speed against LinearDesign is 0.75–1.02×, which is not a figure's worth of claim.
- **Every front panel has the same axes:** x = CAI on 0–1, y = MFE in kcal/mol, optimum
  to the lower right, hull vertices as filled marks and off-hull points as small open
  marks. Where convexity is the argument being made (E2's recovery panel), a log-CAI
  inset is added rather than switching the main axes, since the hull is convex in
  (cai_log, mfe) and not in (cai, mfe).
- **Every plotted row carries its build and its full configuration**, not its table.
  Two rows that differ in a flank length, a stop codon, a default or a commit are not
  comparable, and a cross-build pair is only visible at the point it is written down if
  the stamp is per row. This is what a provenance block owes a figure.
- **Figures build in this repo**, as Snakemake rules under `pipelines/rules/`, from a
  shared renderer in `derna2_analysis/`. `figures.py` already carries the style and the
  method registry; a front renderer joins it rather than being written per experiment.

## The figure set

| # | experiment | what the figure asserts |
|---|---|---|
| 1 | Fig 1 data (GLPE_HUMAN + hbb) | attaching UTRs breaks a CDS-only design, and the structurally optimal front is computable |
| 2 | E2 | a λ sweep is not a route to the front: it reaches hull vertices only, and not many of those at length |
| 3 | E3 | what forbidding motifs costs, against a baseline shown to bind |
| 4 | E4 | what a UTR-aware practice delivers, per practical UTR pair |
| 5 | E5 | the front is what a designer filters; a handful of structural classes summarize it |
| S1–S4 | E1 and the overflow panels of 2–5 | correctness, single-λ cost, memory, host qualification |

## Figure 1 — the overview, on one instance

Data is on disk (`data/results/fig1_fronts`, tables built by `fig1_fronts_tables`).

| panel | content | source |
|---|---|---|
| a | CDS-only front, 145 points, 21 hull vertices — the object the field computes today | `a_cdsonly_front.jsonl` |
| b | the isolation gap of each of those designs once the UTRs are attached: **145 of 145 above zero**, 1.0–11.0 kcal/mol, median 3.3 | `b_cdsonly_with_utrs.scored.jsonl` |
| c | the constrained front (117 points) and the structurally optimal front (114) on shared axes, 0–10.3 kcal/mol apart, median 3.0 | `cd_symmetric_guarded_front.scored.jsonl`, `d_so_front.jsonl` |

**Panel b plots the gap, not the raw drop** (decided 2026-10-24). Attaching the hbb UTRs
lowers each design's free MFE by 33.0–42.5 kcal/mol, of which **31.4–33.9 is the UTRs
folding on their own** — a constant that costs the designer nothing. The crossing
penalty, which is what "the molecule did not keep the structure you designed" means, is
`symmetric_gap`: 1.0–11.0, median 3.3. The large number is a true measurement of the
wrong quantity, and its qualifier is a dimension rather than a bound, so it would travel
without it. The categorical claim is what the panel asserts and it is airtight: the gap
is positive on every one of the 145.

**The panel names its constraint table on its face.** Structural optimality is
table-dependent and the three tables in the deposit disagree on this instance: 0 of 145
under `symmetric`, 0 of 145 under `unpaired`, **14 of 145** under `no_clamp`. The paper's
constraint is `symmetric` (`5-cds,cds-3`, ADR 0044), so without the label "none of them"
and "14 of them" are both quotable from one file.

Two mechanical notes for the renderer. The SO front records carry `cai_log` only and no
`cai`, so the conversion at n = 78 is asserted against panel a's own `cai`/`cai_log`
pairs rather than assumed — a length off by one silently shifts a whole curve. And one
caveat travels with the figure: the most adapted structurally optimal design gives
up 0.58 log-CAI against the maximum, so panel c's SO curve stops short of the top of the
CAI axis. That is a property of the instance and is labelled, not cropped away.

## Figure 2 — E2, the cost of the whole front

| panel | content | source | state |
|---|---|---|---|
| a | Pareto points recovered against *n*, per arm: the walk returns the whole front, every replay arm returns hull vertices only (11 against 31 at n = 32) | `e2` grid rows | grid running |
| b | anytime curve — hypervolume against wall, one long instance, every arm, censoring marked where a cell hit 12 h | per-λ rows with completion walls | grid running |
| c | wall against *n* per arm with the fitted exponent, and the **finishing length** per arm as a rug | `e2` cell rows | grid running |

Panel a compares **arms within a cell**, never point counts across cells: E2 runs
without flanks and every arm of one cell sees one configuration, which is what makes the
11-against-31 comparison a statement about the methods. A yield compared across cells
would be measuring the configuration instead.

Three things the panels must carry as text, because each changes what the numbers mean:

- the replay arms are **given** their λs for free, which understates their cost; the
  gap between `derna1-own` and `derna1-hull` is the size of that gift;
- walls on a contended host are upper bounds, and instruction counts are the axis that
  survives load — if the grid runs under the raised gate, panel c is drawn on
  instructions and the wall version moves to the supplement;
- `derna1-own` spends 1.71–2.33 λs per front point and still misses vertices (9 of 11
  at n = 32); that ratio is a small inset, not a fourth panel.

**As drawn for the paper** (its Figure 4, rule `e2_fig4`): (a) one protein's front with a
magnified hull step, the ΔMFE and ΔCAI between neighbouring designs of its hull and of its
front, and each method's wall time on it; (b) the share of each front a weight reaches, by
length group; (c) those steps pooled per length group; (d) wall time per run, with the runs
each method did not complete; (e) Carrillo-Lipman pruning's speed-up of the dichotomic hull
search, and that search's time over the whole front's. Source: `data/results/e2_fig4`.

## Figure 3 — E3, the price of motif-freedom

| panel | content | source | state |
|---|---|---|---|
| a | the five nested fronts `F0…F4` on one short-band instance, chosen by the rule the plan already fixed | Phase B cells | short band running |
| b | `dMFE_at_CAI` (median, with min–max) per nest across the short band's binding cells, **non-binding cells excluded and counted in the caption** | Phase B metrics | partial |
| c | runtime ratio against the baseline front as the blacklist grows, `F1…F4` being four blacklists on one instance | Phase B cost fields | partial |

The vacuous-entry table is a paper **table**, not a panel: per entry, spellable against
occurring on the baseline front. It is what licenses every price in panels b and c, and
it is the claim no comparator can make.

Panel b's exclusion rule is the one that is easy to get wrong in the plotting code: a
cell with `removed = 0` is dropped from the price panel and **kept** in the runtime
panel, on the same column. Two different filters over one table.

## Figure 4 — E4, UTR-aware design

| panel | content | source | state |
|---|---|---|---|
| a | `A1`–`A4` on one instance at one pair, four curves on the shared axes | A1–A4 cells | partial; some A1 fronts known damaged |
| b | the headline: share of each arm's points that are structurally optimal, per UTR pair — a grouped bar per pair with `A1`, `A2`, `A3` | `free_check` per point | partial |
| c | the isolation-gap distribution of `A1`, per pair, in kcal/mol | `free_check` | partial |
| d | the price of structural optimality: `A2` → `A4` ΔMFE at matched CAI against length, each point labelled complete or bounded | A2/A4 fronts | partial |

Panel b is the figure. Panels c and d are supporting and can move to the supplement if
space binds.

Two labelling obligations the renderer enforces rather than the caption: an `A4` point
is drawn differently when `provisional` rather than `certified`, and a price taken from
a partial `A4` is drawn as an interval, since the open regions have a floor and so bound
the price from both sides.

The pair axis is narrower than five pairs suggests — three of the five are α-globin on
at least one side — so panel b's pairs are ordered by independence of material, with the
shared-material pairs adjacent and marked, rather than alphabetically.

## Figure 5 — E5, the use case

| panel | content | state |
|---|---|---|
| a | `R0…R3` on shared axes, points coloured by MASS cluster, the filter-surviving region shaded | plan only |
| b | `selected_k` against τ, with the knee marked by the pre-registered rule | plan only |
| c | the filter table as a bar: share of the front surviving each practical filter, and the MFE given up at matched CAI | plan only |

Panel a is the only figure in the paper where colour carries a categorical variable that
is not a method, so the MASS clusters take their own ramp and the method palette is not
reused for them.

## Supplement

- **S1 — E1 cost.** Wall and peak RSS against *n* per arm, log-log, fitted exponents.
  Already built as `e1_scaling`; it moves from main text to supplement unchanged.
- **S2 — E1 speed against quality.** The beam arms' instructions against ΔMFE from the
  exact optimum. Already built as `e1_speed_quality`.
- **S3 — host qualification.** What load costs a wall on the shared host, which is the
  evidence behind every "upper bound" label in Figures 2–4.
- **S4 — E3's probe.** The hull probe's rungs, the step change between n = 200 and 250,
  and the per-blacklist frontier.

E1's correctness result is a **table**: every exact arm's re-scored MFE against the
reference, the loop-cap-20 cells, and the 17 cells where LinearCDSfold scores below
Turner 2004.

## Caption hazards, from the certificate lane (2026-10-24)

Four claims were retracted on the day this plan was written, and each is the shape that
reaches a caption intact. They are recorded here so the figures are drawn against them
rather than corrected after review. Source: the certificate lane's section P0,
`handoff/prompts/scope-37-s28-back-on-ada.md` on `chore/s28-p0-margin-is-the-gap`
at `8e0acae4`. The earlier `8b304755` is **not** the citation: its ladder paragraph
carries a figure corrected at `8e0acae4`, so a reference to it would reproduce the
error this section exists to prevent.

- **Seedlessness is a property of a cell and a stop codon, never of a UTR pair.**
  `P54852__hbb` at a pinned `UGA` returns 303 seeds and 314 certified points, so "`hbb`
  is seedless" is false. The surviving statement names the cell and the stop:
  `P01563__hbb`, zero at `UAA`, `UAG`, `UGA` and `none`.
- **A zero-point row is not an empty front.** Every A4 row on record is `stopped_by:
  time`, `complete: false`, so a plotted zero is budget-bound on the recovery axis. The
  budget-free quantity is the **seed** count, a property of the constrained front alone:
  `points ≤ seeds + window + band`, and only the first addend needs no search. A panel
  that shows zeros says which of the two kinds it is.
- **An open region's width is not a front's span**, and they differ 3–4× — 153.48 and
  171.03 against 56.88 and 44.26 on the same cells, with 63% and 74% of each region
  holding no front point. Figure 4d's intervals are drawn from the front's span.
- **A stop choice does not remove a base pair.** `BasePair` includes `GU` and `UG`, and
  0 of 9 anchors is made illegal by any stop, so no caption frames a stop convention as
  closing an anchor.
- **Recovered-point counts are not comparable across cells unless the flanks match to
  the nucleotide — and even then a factor of two carries no signal.** A same-build pair
  differing in one thing, a 3′ UTR of 132 against 134 nt, same stop, same protein, the
  same 117-point constrained front and identical best emitted adaptation, returns **16
  points against 8**. Two nucleotides of flank halve the walk's yield. So any panel
  comparing yields across cells, proteins or UTR variants is measuring the
  configuration, not the method. The stable cross-cell quantity is the **seed count**, a
  property of the constrained front alone that needs no search. Cited to `6bebee78`.
- **Withdrawn: the 78-residue matched-flank ladder comparison.** The claim that the walk
  recovers 8 structurally optimal points on GLPE and none on `P54852`'s first 78, both
  with empty seed sets, is retracted in both halves: the two rows ran on builds **156
  commits apart**, a range containing PR #477 — the outside-pass fix on record as having
  damaged E4 fronts, biting on both stop lanes and so on the `--stop none` setting both
  arms use — plus two default changes to seeding and gap folding, with no same-build
  control. The seed half fails separately and is cross-build too. No figure or caption
  makes a length claim from that pair.

## S5 — withdrawn (the shortfall panel is not drawn)

**Dropped 2026-10-24, on its author's recommendation, before anything was drawn.** The
five-row table's arithmetic reproduces and its rows are sound; the defect is one level
up, in what the number is a number *about*, which no caption fixes.

The figures measure the **decoupled** condition `sup_{d∈R} E_d(S) < inf_R A`, while a
witness certificate requires only the **pointwise** one, `∀d ∈ R: E_d(S) < A(d)`.
Decoupled implies pointwise strictly: it compares the *worst* design's witness energy
against the *best* design's constrained energy — two different designs — so it is false
by inspection on any region holding one mediocre design, and needs no measurement. The
`inf_R A` in the `S`-independent form is where the relaxation enters, which is why the
11.3×–21.5× range quantifies the relaxation rather than the obstruction. It therefore
cannot support "a witness certificate cannot fire", "the regional emptiness bound is
closed" or "certification is impossible".

Two further defects in the same table, independent of that one: the numerator is an
all-designs quantity while the 12.70 ceiling is measured on the 1,483-point **front**
census, an asymmetry of populations; and the 78-aa row borrows the full-length protein's
ceiling, so it cannot separate "the shortfall grows with length" from "both quantities
are extensive and the ratio is scale-free".

A re-captioned version would have to say it quantifies a relaxation strictly stronger
than the certificate's requirement and independently false in one line — at which point
it is redundant rather than wrong, while a skimming reader still leaves with the
impossibility claim. So it is not drawn. Retraction cited to the certificate lane;
the superseded specification is in this file's history.

**What the lane's state actually is, and the ceiling on what any figure may say.** The
closure claim is withdrawn and the question is reopened. Measured on a 482-point census
row with the evaluator checked against our own folds 482/482: a witness from one design's
free fold is a legal structure of **5.0%** of the other designs; where legal, the
pointwise test fires on **25.2%**; 48 witnesses together refute **33.4%** of the front.
Coverage grows sublinearly, so this is not a route to emptiness either. No figure may say
anything stronger than *reopened, with coverage measured*.

**And no figure presents an ε value as final.** The ε-front's soundness is untouched — a
lower bound over all completions also bounds the structurally optimal completions — but a
sound, strictly stronger instrument is in view (a Lagrangian floor pricing the
structural-optimality predicate in, recovering the present floor at zero weights). ε is
captioned as an **upper bound that will shrink**, never as a result. The
*incomplete, not empty* convention stands and is strengthened by this.

## Build order

1. **The front renderer** — one function drawing a front on the fixed axes, with hull
   marks, censoring, provisional points and the optional log-CAI inset. Every front panel
   in Figures 1, 3, 4 and 5 calls it. Buildable today, against Fig 1's data.
2. **Figure 1** — the only main-text figure whose data is complete and frozen.
3. **S1, S2** — re-emit E1's two existing figures at supplement size.
4. **Figures 2–4 as skeletons wired to the real schema**, each with a rule that fails
   loudly on missing rows rather than drawing a partial panel silently. The row-count
   check runs before any panel is drawn, as the experiment plans require of their own
   summaries.
5. **Figure 5** when E5 runs. S5 is withdrawn and is not in the build order.

## What each figure is waiting on

| figure | waiting on |
|---|---|
| 1 | nothing |
| 2 | the E2 grid; the quiet-host decision, which decides whether panel c is walls or instructions |
| 3 | Phase B's short band, then the middle and long bands' selection |
| 4 | the A1 re-run on the cells it was wrong on, and E3's selection for the upper bands |
| 5 | E5 has not run |
| S1, S2 | nothing |
| S3 | the host qualification numbers the grids record anyway |
| S4 | nothing — the probe is closed |

## State of the data, 2026-10-25

Verified this date against the deposits and the branch tips. Each entry names what it
changes in the specs above.

### Figure 1 — two reserved panels are now run, and one was framed backwards

The specs above reserve the motif-excluded CDS-only front and the motif-free modular
front as unrun. `experiment/glpe-hbb-arms` at `97428137` runs both.

- **The motif-excluded CDS-only front exists**: 153 points and 17 gaps, closed in 28.2 s
  at 0.11 GB, against A0's 145 on the same coding sequence with nothing forbidden.
- **Forbidding motifs raises the point count.** 153 against 145 without UTRs, 131 against
  A2's 117 with them. Fewer feasible designs, more lambda intervals to describe what
  remains. Any panel drawing the motif-excluded series as a *subset* of the 145 is wrong;
  it is a larger, separate front and needs its own curve. A reader expects a constraint to
  shrink a solution set, so the direction has to be stated on the panel rather than left
  to be inferred.
- **A5 exists for hba and bnt, not for hbb**, the pair Figure 1 is specified on:

  | arm | hba | bnt |
  |---|---|---|
  | A0 CDS only, no UTRs | 145, complete | 145, complete (copied) |
  | A1 UTRs in the DP, nothing forbidden | 101, complete | 178, complete |
  | A2 UTRs in the DP, `5-cds` and `cds-3` forbidden | 117, complete | 124, complete |
  | A3 competitors at A2's hull vertices | 44 + 44 scored | 44 + 44 scored |
  | A4 modular under A2 | 88, 1 open region, 1.9239 open | 95, 5 regions, 54.8067 |
  | A5 the same over motif-free designs | 80, 1 open region, 4.3444 | 85, 8 regions, 2.9273 |

  `hba` is the pre-registered hold-out and carries no analysis, so bnt is the only
  retarget candidate. A0 is pair-independent and byte-identical across pairs, so the
  124-of-145 claim survives a pair change untouched.
- **hbb's own A4 is complete and exact**, and the arms ladder does not change that. Every
  window over `cai_log` [−77.64, −0.5816] closed, the stretch above −0.5816 (156,698
  designs) was settled exhaustively with each design free-folded and guard-folded, and all
  114 points verify free MFE = guarded MFE = stated energy. **No A4 or A5 in the ladder
  above is exact**: every one has open regions, and the open *measure* rather than the
  region count is the statistic — bnt A4 is 5 regions but 54.8067 log-units open, against
  bnt A5's 8 regions and 2.9273, so bnt's A4 is the most uncertain front in the ladder and
  not the safer of the two. Retargeting Figure 1 from hbb to bnt therefore trades an exact
  front for a badly incomplete one on the arm the figure most needs. The branch's reading is that
  three independent routes fail to find a modular design in the well-adapted corner, which
  is evidence the corner holds nothing and not proof, since none is exhaustive. The
  *incomplete, not empty* convention therefore binds any A5 panel.
- The run outputs are **not on ada**. The branch reads them through
  `glpe-arms-runs/so_state.py <rundir> --front <front.jsonl>`, and no such directory
  exists under any readable path here, so these counts are the markdown's and are not
  re-derived.
- Figure 1's existing inputs were re-checked against `pins/results-fig1.sha256` and
  `pins/results-fig1_fronts.sha256`: all six match, `unguarded/front.scored.jsonl`
  included.

### Figure 2 — the grid's rows are on the second host and unreadable from ada

The three panels above are specified against a grid that, in the deposit reachable here,
produced **one cell with all seven arms**. Distinct rescored points per protein per arm,
from the 437-row sweep of 2026-10-16:

| protein | n_aa | derna2-walk | derna2-hull | derna1-own | derna1-hull | LinearDesign | LCDS | mRNAfold |
|---|---|---|---|---|---|---|---|---|
| P0DTC2_p32 | 32 | 31 | 11 | 9 | 11 | 11 | 11 | 11 |
| P0DTC2_p41 | 41 | . | . | 11 | . | . | . | . |
| P0DTC2_p51 | 51 | . | . | 10 | . | . | . | . |
| P0DTC2_p64 | 64 | . | . | 14 | . | . | . | . |
| P15421 | 78 | . | . | 13 | . | . | . | . |
| P0DTC2_p81 | 81 | . | . | 13 | . | . | . | . |
| P13501 | 91 | . | . | 14 | . | . | . | . |
| Q99525 | 98 | . | . | 11 | . | . | . | . |
| P0DTC2_p102 | 102 | . | . | 14 | . | . | . | . |

133 of 437 rows carry `wall_s`, and `derna2-walk` carries 2 of them. Panel a needs several
arms at several lengths and has them at one length; panel c needs walls against length per
arm and has two points for the reference arm.

What that one cell does establish is sharper than the specs claim. Rescoring every arm's
designs and testing them against the walk front's own vertex set: **all five replay arms
return the identical 11-point set, and that set is exactly the hull vertex set, with zero
off-hull points returned by any lambda arm.** `derna1-own`, choosing its own weights,
returns 9 — a strict subset of those 11, so it misses 2 vertices as well as all 20 off-hull
points. `derna2-walk` returns 31, being the 11 vertices plus 20 points no weight reaches.
No competitor design is dominated under the rescoring model, so the panel carries a claim
about reachability and not a confound about energy models.

A second thing the specs drop and should not: `derna2-walk` completes at **one of nine**
lengths and carries `wall_s` on 2 rows. A completion strip — arm against length, filled
where the arm returned — discloses that in the main figure instead of leaving the exact
route's own censoring to the supplement, and it costs nothing, being the dot pattern of the
table above.

`experiment/e2-tamazirt` has since advanced to `1eaf4892`, over `42be727e` (pass-1 re-run
on the xg build) and `4a304a05` (concurrent hull pass with a backfilling walk scheduler,
walk flags by length band, `--rewalk`). Those commits carry drivers, configs, plans and 56
A/B records across `occ-ab`, `occxg-ab`, `ownloose-ab` and `control` — **not grid rows**.
`config.e2.tamazirt.json` writes to `/home/melkebir/projects/derna2/e2/out`, which on ada
is a local filesystem with mode 700 and no access rule, and tamazirt does not resolve from
ada. **Figure 2 is blocked on a deposit, not on a run**: the rows have to reach a readable
path before any panel is drawn, and `data/results/e2/` already symlinks into scratch for
exactly this.

### Off-hull share against length, and the convention it needs

Computed with `fronts.hull_vertices`, which is the lower hull in (cai_log, MFE) and so the
lambda-reachable set:

| instance | front points | lambda-reachable | off-hull |
|---|---|---|---|
| n = 32, P0DTC2_p32 | 31 | 11 | 64.5% |
| n = 78, P15421 | 145 | 21 | 85.5% |
| n = 100 | 258 | 33 | 87.2% |
| n = 200 | 763 | 59 | 92.3% |

The instrument was checked before the numbers were used: it reproduces the n = 78 file's
own `hull_vertex` flag exactly (21) and the n = 100 run's own `vertices` list exactly (33).

One qualifier travels with this table: the four lengths are **four different proteins**, so
it is a trend across instances and not a length series on one — the rule above about
comparing arms within a cell applies to yields across cells as well.

Two questions the table raised were settled on 2026-10-25 by refolding.

**The 2026-10-07 energies are current.** The n = 100 and n = 200 runs predate the
folding-scan fix, and known-stale siblings exist at n = 300 and n = 500. Nine designs
spread across both fronts were refolded with `fold_seq` from the Figure 1 deposit
(`1b3739b`, `xg` `7bb14ec` + step 0) and **all nine reproduce the stored energy exactly**,
so the scan fix does not move this instance's energies. This tests **pricing, not
coverage**: a scan defect that caused the walk to *miss* front points would not show up in
a refold of the points it did return, so the completeness of these two fronts remains
untested and a point count from them is the weaker half of the claim.

**At n = 200 the vertex count is not resolvable from the stored energies.** The run records
58 where `fronts.hull_vertices` computes 59. The disputed point (`cai_log` −1.579181,
−227.8 kcal/mol, index 717 of 763, an interior point and not an endpoint) refolds to its
stored energy, and its turn is not shallow in normalized terms (0.597 against its two edge
lengths) — but it lies only **0.066 kcal/mol below the chord joining its neighbours**, and
energies are stored to 0.1 kcal/mol. Its vertex status is therefore below the resolution of
the recorded data, and a tolerance in the run's own routine is the sufficient explanation.
Quote the share as ~92%, not 92.3%, and read any run's stored `vertices` as what its search
found at its own tolerance rather than as what the front contains.

There is **no cross-deposit replicate** at any length: the 2026-10-16 sweep carries n = 32
and the 2026-10-07 deposit carries n = 30, 50, 100 and 200, so no instance appears in both
and the two cannot validate each other.

### Figure 3 — the bands are stopped on purpose

E3's middle band (Q9UNN5, Q06187, 650-659 aa, four depths each) is stopped until the
gap/slice stack carrying gap-rest lands on `xg`: five of eight cells settled zero gaps in
2.5 h while two cells of the same length settled 26-27, so at that length a gap can exceed
a 4 h round and an aborted round banks nothing. The long band (>= 850 aa) has never
started. The grid ran on pin `06605630` against `xg` at `ab620b96`, and the decision on
file is to re-pin at the restart. `mid/` holds what was banked; Q06187 F1/F2 hold 26-27
settled gaps as a resume point. Figure 3's waiting row above is unchanged, but the reason
is now a pin boundary rather than queue depth.
