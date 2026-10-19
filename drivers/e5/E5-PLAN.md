# E5 — one design, start to finish: the combined ladder and what the front says

> Status: **plan, not run.** Written 2026-10-19. Nothing below has been measured
> except where a line names an earlier scope's result or a `--help` output and
> says so.

E1 measures one λ. E2 measures what a whole front costs. E3 measures what a motif
constraint costs. E4 measures what attaching the UTRs costs. **E5 is the use
case**: one protein, one UTR pair, every practice applied at once, and then the
question a designer actually has — *given this front, which molecule do I order?*

It contributes the cell neither E3 nor E4 contains — **UTRs and motifs
constraining the same design** — and then summarizes the result with
[MASS](https://github.com/elkebir-group/MASS) (Maximum Agreement Secondary
Structures), so the front is reported as a handful of structural design classes
rather than as a few hundred points.

## The claim E5 exists to make

A single-λ design cannot be filtered. A front can.

The practical constraints a designer brings — GC window, uridine content, CpG
count, a Kozak preference, a manufacturer's composition spec — are mostly **not
expressible as hard substring exclusions**, and E3 says so entry by entry
(CpG density, U content and rare-codon clusters are composition rules, not
k-mers). They can therefore only be applied *after* optimization, as filters on a
set of optimal designs. That is an argument for computing the exact front rather
than one scalarized point, and it is the argument E5 makes with numbers:

> for each practical filter, the share of the front that survives it, and the
> MFE/CAI the surviving designs give up against the unfiltered front.

Everything else in E5 serves that sentence.

## The ladder

Four nested feasible sets, four fronts, in this order (decided 2026-10-19):

| rung | feasible set | practice added | flags added |
|---|---|---|---|
| `R0` | all of `V(w)` | none — the baseline | — |
| `R1` | `R0` minus designs pairing with a UTR | the molecule is designed with its UTRs attached | `--utr5 --utr3 --stop UAA --forbid 5-cds,cds-3` |
| `R2` | `R1` minus designs containing a blacklisted motif | the manufacturing / processing / repeat / initiation blacklist | `--motif` × E3 sets I–IV |
| `R3` | `R2` minus designs that pair in the start window | the cap-proximal start codon is kept accessible | `--access UP,DOWN` |

Cumulative, as in E3: `R3` forbids all three practices at once. The order runs
from the outermost design context inward, and it puts the paper's UTR
contribution first; the alternative ordering (motifs first, making E3's `F4` this
ladder's rung 1) was considered and rejected, and the cost of the choice is that
**attribution is order-dependent** — a design removed by both the UTR constraint
and the blacklist is charged to whichever comes first. The plan states the order
rather than hiding it, and the per-rung `removed` counts are reported per
practice *and* per design, so a reader can see the overlap.

`--forbid 5-cds,cds-3` is explicit on every rung from `R1` on, for E4's reason:
under the three-class constraint `R1` is `R0` shifted by a constant and carries
no design `R0` does not, so the two-class setting is the experiment and the
three-class setting is a control. **E4's step-2 exactness gate for the two-class
constraint is a prerequisite of E5 and is not re-derived here** — E5 consumes it.

`--access UP,DOWN` needs its window fixed before any run. **Pre-registered:
`--access 4,30`** — four nucleotides upstream of the A of AUG (the Kozak −4…−1
positions, which lie in the fixed 5′ UTR) through 30 nucleotides downstream (the
start codon plus nine codons), the window the ribosome-scanning literature treats
as the start-accessibility footprint. It is one number chosen once, written here,
and not tuned after seeing a front.

**Monotonicity is a free soundness gate**, exactly as in E3: at every λ,
`obj_λ(R0) ≤ obj_λ(R1) ≤ obj_λ(R2) ≤ obj_λ(R3)`. A violation is a bug, never a
finding, and it is checked at every λ of every cell.

## Instances — two, reported as two

Decided 2026-10-19: one short and one long, never pooled.

| role | instance | length | UTR pair | why |
|---|---|---:|---|---|
| short | SARS-CoV-2 RBD | ~220 aa | `bnt` | the antigen with the UTRs of an approved construct; inside the lengths Scope 37 ran, so complete exact fronts are reachable and the full summary is computable |
| long | one Cas ortholog, `PRT-03` | 984 aa | `hbb` | a realistic editor-sized transcript with the standard IVT globin pair; the shortest of E3's long arm, so it is the long instance most likely to produce something |

Each instance gets the UTR pair a real construct of its kind uses, so each reads
as a design rather than as a grid cell. Sequences come from E4's `utrs.fasta`
pin with its `PROVENANCE.md` SHA-256s; E5 adds no new sequence source, and a pair
E4 could not source publicly is not reconstructed here either.

**The long instance will be partial, and that is planned for, not discovered.**
E3's hull probe projects the *hull alone* outside an hour at n = 300 and ≈15 GB
at 984 aa, and the front is the hull plus its off-hull gap solves. `PRT-03`'s
rungs are therefore expected to report partial fronts, and the MASS summary then
runs over a partial front: **every design answer taken from it is a bound, signed
and labelled as one.** The short instance carries the complete-front claims.

Addressed by ID as in E2, E3 and E4, so no ortholog name reaches a row.

## Phases — Phase A runs today

`derna design --help` (frozen bin `e2/bin/derna`, `df835619`, verified
2026-10-19) accepts `--utr5 --utr3 --stop --forbid --access --motif
--motif-window --pin` **together**. `derna front --help` accepts `--utr5 --utr3
--stop --forbid --access` but has **no `--motif`**, and its `--max-secs` /
`--max-rss-gb` are listed only under the structurally-optimal walk, with no
`--stream` at all.

So the plan splits where the code splits:

- **Phase A — the ladder at fixed λs. Needs no new code; runs today.** Run the
  four rungs with `design` at the hull λs of `R0`'s front (from E2's grid, or a
  fresh `front` run on the short instance), on both instances. Yields Δobjective,
  ΔMFE, ΔCAI, wall and RSS per λ per rung, and — the part that matters for E5 —
  **a set of designs large enough to run MASS on**. A ladder of per-λ designs is
  not a front, is never called one, and the summary computed from it is
  explicitly labelled as the λ-sampled summary.
- **Phase B — the exact combined fronts.** `front --motif` from
  `scope/38-motif-front`, with its UTR and access paths exercised at the same
  time; Phase A's per-λ gaps bound what Phase B's front-wise gaps can be, so
  Phase B's numbers arrive with a prior. Phase B is where `R2` and `R3` become
  fronts rather than λ samples.
- **Phase C — the summary.** MASS over the pooled designs, the filter table, and
  the cluster annotations. Phase C runs on Phase A's output and runs again on
  Phase B's; its code is written and validated against Phase A, so nothing waits.
- **Phase D — the side cells.** The Kozak pin, the three-class control, and the
  no-UTR / no-motif cross-checks below.

Nothing in E5 waits on a run it does not consume.

### What E5 needs from the scope-38 build

The same three things E3 needs, for the same reason, plus one E3 does not:

1. `front --motif`, at all;
2. `--stream`, `--max-secs` and `--max-rss-gb` on `front`, so a stopped long cell
   leaves its hull and its finished gaps on disk instead of nothing. E3 measured
   the failure mode: `motif_hull` buffers and prints at the end, so a rung killed
   at 3600 s produced **zero** lines;
3. **`--motif` and `--utr5/--utr3/--forbid` exercised together on the front walk.**
   No run in this project has combined them. Phase B's first step is the
   combination check below, and it is a gate.

## Run order, with the gate for each step

1. **Pin.** One commit, one frozen bin directory, `SHA256SUMS` + `PIN`, as E2 and
   E4 do. MASS is pinned too: its commit, `--algorithm`, `--tau` range and
   `requirements.txt` versions go in the same `PIN`.
2. **Combination check** (Phase B's gate, and runnable in part today with
   `design`): on a 30-residue protein with the `bnt` pair and all four motif
   sets, the emitted design must (a) translate back to the protein, (b) leave the
   UTRs byte-identical, (c) contain no motif of any set under an independent
   scan, (d) form no `5-cds` or `cds-3` pair under an independent classifier, and
   (e) keep the `--access 4,30` window unpaired. **Gate: all five, on the same
   design.** A constraint that silently stops applying when another is present is
   the specific failure this catches.
3. **Feasibility pre-check, per instance per rung.** Infeasibility is detected at
   lattice construction (E3 measured 8 ms), so it costs nothing. `R3` is the live
   risk: the access window and the blacklist both constrain the first ten codons,
   and `GGGGG` plus a pinned unpaired window could in principle have no joint
   solution. An infeasible rung is recorded as `infeasible` with the practice that
   did it — an outcome, not a failed run.
4. **The three-class control**, one rung, short instance: `R1` under
   `--forbid 5-cds,cds-3,5-3` must equal `R0` shifted by `MFE(u5) + MFE(u3)`,
   point for point, byte-identical designs. **Gate: exact equality.** It verifies
   the UTR plumbing against a closed-form answer; if it fails, no other cell is
   readable.
5. **Phase A, short instance**, all four rungs.
6. **Phase C against Phase A**, so the summary code is finished and validated
   before the expensive fronts exist.
7. **Phase A, long instance**, all four rungs.
8. **Phase B** when scope-38 lands: short instance first, longest-first within an
   instance, one cell at a time.
9. **Summarize** — row counts against the plan first, then exit categories and
   completeness, then numbers. No number is read before the row counts match.

## MASS — how the front gets summarized

MASS decomposes each dot-bracket structure into **structure units** (hairpin,
interior/bulge, multiloop, stem run, as sets of base pairs plus unpaired
nucleotides) and selects a minimum-size set of units whose presence pattern
partitions the collection into at most τ clusters. It reports `selected_units`,
`cluster_map`, `selected_k` and `objective_value`. `mstp` is the default and
needs no license; `ilp` is exact and needs Gurobi; `mstp_beam` sits between.

Two properties of this front make it an unusually clean MASS input, and both are
worth stating in the paper:

- **Every design on one instance's front has identical length** — same protein,
  same UTRs, same stop. MASS's gapping/alignment path is therefore never
  exercised, and the unit coordinates are directly comparable across designs.
  This is why the pooling in E5 stops at the instance boundary.
- **The collection is a Pareto front, not a Boltzmann ensemble.** MASS's usual
  input is one molecule's structural ensemble; here it is many molecules' single
  optimal structures. The clusters therefore answer *how many structurally
  distinct ways there are to be optimal*, which is a different question from
  ensemble diversity, and the plan does not blur them.

### What MASS is run on

**Pooled across the ladder, labelled by rung** (decided 2026-10-19). One MASS run
per instance over the union of `R0…R3`'s designs, each structure ID carrying its
rung and its λ (or its front index). One feature basis across all four rungs
makes cluster identities comparable, which is the only way to say *which
structural class each practice removes*. Per-rung runs are kept as a secondary
check, and the two are reported side by side. Pooling across the two instances
was considered and rejected: different lengths would force the gapping path for
no gain.

**Both folds, as paired runs** (decided 2026-10-19). Each design has two
structures, and MASS runs on each collection separately:

- the **constrained fold** the walk emits — what the method promised;
- the **free re-fold** of the whole molecule by `mfe_fold` — what the molecule
  actually does.

The paired cluster maps are then compared. Where a design's constrained cluster
and free cluster disagree, the molecule did not keep the structure it was
designed for; that is E4's isolation gap told in structural terms rather than in
kcal/mol, and it is the first structural account of it this project has. The
agreement between the two cluster maps is reported as a confusion table per rung,
not as a single index.

### τ, and the curve that replaces it

τ is not a parameter to pick; it is an axis. **Pre-registered: `--tau-range 2 10`**
on every collection, with the reported quantity the **`selected_k` vs. τ curve**
— how few structural features suffice to resolve the front into τ classes.
Its shape is the summary:

- a curve that stays flat and low means a front with a few structural themes and
  a designer needs to look at τ molecules, not hundreds;
- a curve that climbs linearly means every design is its own structural class and
  the front genuinely cannot be summarized — a negative result E5 is willing to
  report, and the reason the headline is a *filter* table rather than a cluster
  count.

The headline τ for the figures is chosen by the **knee of that curve, computed
by a rule fixed here** (largest τ whose `selected_k` increase over τ−1 is below
the mean increase over the whole range; ties to the smaller τ), not by eye.

### MASS soundness gates

1. **Exactness where it fits.** On the short instance, `mstp` and `mstp_beam`
   (beam 5 and 10) are run against `ilp` at every τ. `mstp`'s `objective_value`
   must be ≥ `ilp`'s, and where it is equal the heuristic is exact on this input.
   If no Gurobi license is available the comparison is dropped **with that stated**,
   and no claim of minimality is made — `selected_k` is then reported as an upper
   bound at every τ, which it is.
2. **Non-vacuity by permutation.** A control collection — the same structures with
   their unit columns permuted independently per structure — must need
   substantially more units at the same τ. If a shuffled front summarizes as well
   as the real one, the clusters carry no information and the summary is not
   reported. This is E3's mutation-oracle discipline applied to the summary.
3. **Determinism.** Every MASS run is repeated once from the same input; the
   `selected_units` and `cluster_map` must match. `mstp_beam`'s beam order is a
   known place for this to fail.
4. **Input provenance.** The structures MASS reads come from the emitted
   `structure` field (constrained) and from `mfe_fold`'s output (free), with no
   re-folding in between, and their lengths are asserted equal to
   `|u5| + 3·(|P|+1) + |u3|` before the run.
5. **Row count before reading**, as everywhere: structure count in equals
   structure count in `cluster_map`.

## The practical design questions

Four groups, all four selected 2026-10-19. Every metric is computed post hoc from
the emitted `sequence` and `structure` by code sharing nothing with the DP, so
none of them costs a run. Each is reported **per front point**, **per MASS
cluster**, and **as a filter** (what share of the front survives, and what the
survivors give up).

**Composition.** GC fraction and GC3 (wobble-position GC) over the CDS, globally
and in a sliding 50-nt window with the window extremes reported; uridine
fraction, which is the immunogenicity axis and is not expressible as an
exclusion; CpG count and density, likewise; longest homopolymer run; and the
per-entry occurrence counts of E3's four motif sets, so `R0`'s and `R1`'s designs
can be read for what the blacklist would have removed. The filter form is a
manufacturer's spec: **GC in [40 %, 60 %]** and **GC3 in [40 %, 70 %]**, fixed
here, reported as surviving share and as the MFE/CAI cost at the filter's edge.

**Initiation.** The Kozak context's −3 and −6…−1 positions lie in the fixed 5′
UTR and are not design variables; **the +4 position is** — it is the first
nucleotide of codon 2, so `+4 = G` is attainable only if residue 2's codon menu
admits a G there. Reported per instance: whether the front contains any design
with `+4 = G` at all, the share that do, and what the preference costs in MFE at
matched CAI. Because `design --pin` exists today, the cost is also *measured*
rather than only read off: **Phase D pins `4:G`** and runs the ladder again on the
short instance. Alongside it: the unpaired fraction of the `--access 4,30` window
under the free fold (which `R0`–`R2` do not constrain and `R3` does), the local
MFE of the first 50 nt of the molecule, and the MFE of the 25-nt window centred
on the 5′ UTR–CDS junction.

**Structure.** The count of CDS–UTR pairs by class (`5-cds`, `cds-3`, `5-3`)
under the free fold; the longest stem; the hairpin count; MFE per nucleotide; the
number of structure units; and the isolation gap (free MFE minus constrained
energy, positive when the molecule prefers a fold the constraint forbids), which
is the quantity E4 defines and E5 reuses unchanged.

**Decision support.** This is where MASS does the work. The ≤ τ clusters at the
knee are the **shortlist**: for each cluster, its representative (the point
closest to the cluster's centroid in unit space, ties to the better MFE), the
share of front points it holds, the MFE and CAI range it spans, and its
composition and initiation profile. Then: how much of the front's MFE/CAI range
the τ representatives cover; the pairwise diversity of front designs (Hamming
distance between sequences, base-pair distance between structures) within and
between clusters; and the distance from each representative to three reference
answers — **the native coding sequence**, **LinearDesign's single design**, and
**DERNA v1's design at the matched λ** — all three of which exist as binaries in
`e2/bin`. "How far is the shortlist from what you would have shipped" is the
question a reader brings, and it is answerable with no new runs.

## What each cell records

Per point: `mfe` (re-scored by `mfe_fold` on the emitted molecule), `cai`,
`cai_log`, `on_hull`, `sequence`, `structure` (constrained), the free fold's
structure and MFE, the isolation gap, the per-rung constraint verdicts, and every
metric of the four groups above.

Per cell: rung, instance ID, UTR key, constraint classes, motif sets, access
window, points, hull vertices, gaps solved and gaps left open, wall, CPU, peak
RSS (`VmHWM`), exit category, completeness verdict, and the open regions of a
partial answer.

Per MASS run: instance, fold (constrained or free), pooled or per-rung,
algorithm, τ, `selected_k`, `objective_value`, `num_clusters`, `selected_units`,
`cluster_map`, runtime, status, and the determinism and non-vacuity verdicts.

## The numbers the paper takes from it

1. **The filter table** — per practical filter (GC window, GC3 window, U
   fraction, CpG count, `+4 = G`), the share of the front that survives and the
   MFE at matched CAI that the survivors give up. The argument for the front.
2. **The ladder's prices** — ΔMFE at matched CAI and ΔCAI at matched MFE between
   consecutive rungs, plus the points removed and which practice removed them.
   The combined cell neither E3 nor E4 contains.
3. **Whether the practices interact** — the combined price of `R3` against the
   sum of the practices applied alone. A super-additive cost is the result that
   justifies solving them jointly instead of in sequence, and it is the one
   number that needs both E3's and E4's machinery in the same run.
4. **How few molecules a designer has to look at** — `selected_k` vs. τ, the knee,
   and the share of the front's objective range the τ representatives cover.
5. **Designed structure against actual structure** — the paired cluster maps'
   confusion table per rung, with the isolation-gap distribution beside it.
6. **The distance to what you would have shipped** — representative vs. native,
   LinearDesign and DERNA v1.

One figure per instance: the four rungs' fronts on shared axes, points coloured
by MASS cluster, with the filter-surviving region shaded. The figure's instances
were chosen above, before any front was read.

## Risks, and what each one would cost

- **No rung of the long instance completes.** Expected for `R2`/`R3`. The long
  instance then reports Phase A's λ-sampled ladder, completeness and cost, and
  every summary number from it is labelled a bound. The short instance carries
  claims 1, 4, 5 and 6 unchanged. E5 still answers its own question.
- **`R3` is infeasible.** The access window and the blacklist both bite on the
  first ten codons. Step 3 catches it before any front runs. The ladder then
  stops at `R2`, and the infeasibility is reported with the pair of practices that
  caused it — a sharper result than a price, and the same kind of finding as E3's
  `GGGG`-on-WG and `CCNCCN`-on-PP rows.
- **`--motif` and the UTR constraint do not compose on the front walk.** Step 2
  catches it. Phase B then reports whichever composition does hold, and the
  failure to compose is stated as a limitation of the build with its evidence.
- **The front does not summarize** — `selected_k` climbs linearly in τ. Then
  claim 4 is reported as the negative result it is, the clusters are dropped from
  the figure, and the filter table carries E5 alone. Gate 2 exists so this is
  distinguishable from a summary that works.
- **A filter is non-binding** — every design on the front already satisfies the GC
  window. Reported as non-binding, which is a real property of that protein's
  codon menu, not a null result. Per E3's rule, a non-binding filter is worthless
  to a *price* and genuine for a *characterization*, and the two tables treat it
  accordingly.
- **The λ-sampled ladder is read as a front.** Phase A's output is labelled
  `lambda_sampled` on every row, and no figure plots it as a curve.

## Host

tamazirt, after E2's `par` phase, sharing E3's and E4's rules: `--threads 1` for
anything timed, `systemd-run --user --scope -p MemoryMax=48G -p
MemorySwapMax=0`, E-cores (16–27) for builds, screens, MASS and all post-hoc
metrics, P-cores (0–14 even, siblings idle) only when nothing else is timing.
MASS and the metric code are cheap and can run beside anything.

Before launching any grid, check nothing else is already running it
(`ps -eo args | grep …`), and run no `_*` glob inside the experiment directory
while a sweep is live.

## Running rule

Keep going until at least the time named at launch. Before ending any turn, run
`date`; if it is earlier, start the next step or extend the plan. Do not wait on a
run while any independent item remains — Phase A, the summary code of Phase C,
and the Phase D side cells are independent of each other and of scope-38. Stop
early only for a decision that belongs to the user, for anything destructive or
outward-facing, or when no meaningful task remains.

Extension steps, in order: the second Cas ortholog; the `min` low-structure UTR
pair as the non-binding control on both instances; the reversed ladder order
(`R0 → +motif → +UTR → +access`) on the short instance, which measures how much
of the attribution the ordering carries; and MASS over the structurally-optimal
front of E4's `A4`, which would summarize designs whose own free fold respects the
constraint and so needs only one fold rather than two.
