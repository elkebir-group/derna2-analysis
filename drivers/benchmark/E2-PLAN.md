# Paper benchmark E2 — enumerating the Pareto front

The core grid (E1) measures one λ per cell. This experiment measures what it
costs to produce the **whole front**: runtime and peak memory for every method
that can only reach the front by sweeping its scalarization parameter, against
derna2's front enumeration.

The question is not which method is faster at one λ — E1 answers that, and the
answer is within a factor of about two. It is whether sweeping λ is a viable
route to the front at all at therapeutic protein lengths, and what the exact
front costs by comparison.

## What each method is given

Four of the five comparators optimize a scalarization of the same two
objectives, so each has a parameter that plays λ's role. None of them
enumerates a front: a sweep is the only route they have.

| arm | binary | λ it is given |
|---|---|---|
| `derna2-hull` | `scope35_walk_front` | its own chord sweep |
| `derna2-walk` | `scope35_walk_front` | its own chord sweep plus gap solves |
| `derna1-own` | DERNA v1 `-s 3` | its own adaptive sweep |
| `derna1-hull` | DERNA v1 `-s 2 -l λ` | derna2's hull λs |
| `lineardesign-hull` | LinearDesign | derna2's hull λs, mapped |
| `lcds-exact-hull` | LinearCDSfold `-m exact` | derna2's hull λs, mapped |
| `mrnafold-hull` | mRNAfold | derna2's hull λs, mapped |

DERNA v1 appears twice on purpose. `derna1-own` is the only arm that both
chooses its own λs and reports a front, so it is the honest sweep baseline;
`derna1-hull` isolates the cost of the folds from the cost of the λ search by
giving v1 the same λ list every other comparator gets.

`derna2-walk` is the one arm whose unit of work is not a λ: it enumerates the
whole front in a single invocation, hull vertices and the off-hull points
between them alike. It is therefore timed in the `hull` pass rather than as a
cell of `sweep` — the same measurement wrapper, the same row schema, one cell
row and one row per front point, with `lambdas_done` carrying the point count
because a per-λ cost is a category error for it. That pass has to run first
anyway, since every replay arm takes its λs from its output, and timing it
there is also what keeps the most expensive arm from being paid for twice.

### The λ map is per tool and is read from its source

Three of the comparators state the trade-off as a ratio rather than a convex
weight. For a derna2 weight λ ∈ (0,1), the equivalent parameter is

    λ' = (1 − λ) / λ

for LinearDesign and mRNAfold (`corner_lambda` in `drive.py` already maps both
this way, cross-validated at 168/168 hull vertices for mRNAfold).

**DERNA v1 and LinearCDSfold take the convex weight on a different scale,
measured, not assumed.** Both score `λ·E_dcal − (1 − λ)·Σ ln w` where derna2
scores `λ·E_dcal − (1 − λ)·100·Σ ln w`. The energy is in dekacal in all three;
those two leave the CAI term unscaled, so **at the same λ they weight CAI 100×
less**, and the map is

    λ_tool = 1 / (1 + 100·(1 − λ)/λ)

so derna2's 0.9, 0.7, 0.5, 0.3, 0.1 become 0.08256881, 0.02280130, 0.00990099,
0.00426743, 0.00110988.

Handed derna2's λ unchanged, v1 at λ = 0.5 returns CAI 0.066 where derna2
returns 0.881, and **LinearCDSfold returns the MFE-corner design at every
interior λ tried** (0.5, 0.2 and 0.05 all gave −56.20 kcal/mol at CAI 0.733 on
a 32-residue instance). Under the map above both reproduce derna2 exactly at
all five λs: −56.20, −55.90, −53.60, −50.60, −41.90 kcal/mol at CAI 0.733,
0.754, 0.831, 0.902, 0.972. v1's arithmetic checks too — it printed
`O: -1319.3` against `0.5·(−2720) − 0.5·(−81.4029) = −1319.29855`, and its
`cai_s` is the geometric mean, `exp(−81.4029/30) = 0.066309` against the
printed `0.0663089`.

**The core grid could not have caught this, and neither can any corner check.**
At a corner the CAI term's scale cannot change the answer, because the second
objective only breaks ties — which is what the corners are chosen for. E1 ran
corners only, so `corner_lambda` passing LinearCDSfold the weight unchanged is
correct there and wrong for every interior λ. This is why the verification
below is specified at interior λs, and it has now failed twice at exactly that
step.

**v1 also needs a codon-table conversion.** It takes `-c`, but in its own
format — two lines per amino acid, a codon header row then a count row padded
to six fields — while the canonical human table used in E1 is `CODON,AA,freq`
triples. The table v1 ships (`data/InputFiles/sample_codon_usage.csv`) is not
the human one either: its Ala counts are 10/65/4/0 against the canonical
0.26/0.40/0.23/0.11. Until the canonical table is converted into v1's format,
v1 optimizes a different objective than every other arm, and the λ map above
cannot be checked end to end on designs.

**Every map is verified, not assumed.** For each comparator:

1. An ε-vs-ε/2 toggle at both corners, as E1 does — a tool that ignores the
   change is rounding λ (LinearCDSfold parses with `stof`, which is why
   ε = 2⁻²³).
2. At a sample of interior λs, the arm must return the hull vertex derna2
   reports as optimal there. A map off by a transform that happens to agree at
   the corners fails here and nowhere else.

A map that cannot be verified disqualifies the arm from the replay comparison;
it does not get run with an unverified λ.

### Which λs the replay arms run

One λ per hull vertex, at the **midpoint of the λ interval over which that
vertex is optimal**, plus the two ε corners. At a breakpoint λ two vertices tie
and each tool may return either one, so a breakpoint λ makes the recovered
vertex a function of tie-breaking rather than of the method; an interior λ has
one answer. V vertices therefore cost V + 2 runs.

The λs run **corners first, then in bisection order** over the vertex index, so
a cell stopped at its limit has covered the whole λ range coarsely rather than
one end finely. This is what makes a truncated cell's partial front comparable
to a complete one.

**The replay arms are given their λs for free.** derna2's chord sweep paid to
find them; the comparators do not. This understates their cost and must be
stated wherever the numbers are reported. `derna1-own` is the arm that pays its
own λ search, and the gap between `derna1-own` and `derna1-hull` is the size of
what the others are being given.

## Instances

77, in three blocks:

- **55 UniProt proteins ≤ 500 aa** — the LinearDesign set, as in E1
  (`uniprotSeq`, restored by `just sync-cpp`; discovered by scanning the
  directory, never enumerated here).
- **Every 3rd rung of the P0DTC2 ladder**, counted down from the full protein:
  32, 41, 51, 64, 81, 102, 128, 161, 203, 255, 321, 404, 508, 640, 804, 1012,
  1273 aa. A ~27% step locates the length at which each method stops finishing
  to within about 27%, at a third of the long-rung cost of all 50.
- **The 5 Cas orthologs**, `PRT-03` … `PRT-07` at 984, 1053, 1082, 1368 and
  1629 aa. These extend the length range past the spike and are the first
  instances in this benchmark that are not one protein's prefixes, so they
  separate length from sequence identity: a ladder rung and its predecessor
  share a prefix, two orthologs do not.

  They are addressed by ID, not by ortholog name, and the config reads them
  from `/scratch/xinyugu4/paper-e2-inputs/prt/PRT-0N.fasta` rather than from
  the Scope 37 sequence directory, whose filenames are the original names. That
  keeps the name out of the config, out of the instance labels, and so out of
  every row, table and figure downstream. Each copy's length was checked
  against the ID table on write.

**Deferred:** the other 33 ladder rungs (30, 35, 38, 44, 47, 55, 60, 70, 75,
88, 94, 110, 119, 139, 150, 174, 188, 219, 237, 276, 298, 347, 374, 436, 471,
549, 592, 690, 745, 868, 937, 1092, 1179). They are a pure addition — same
arms, same λ rule, same limits — so they can be run later into the same output
without redesigning anything. Adding them roughly doubles the ladder's cost and
is only worth it if the finishing-length crossover needs to be resolved more
finely than 27%.

## Limits, and what a cell that hits one records

- **12 h per (instance, arm).** A cell is one method's whole sweep over one
  instance, not one λ.
- **200 GB per cell.** Socket-local memory is about 250 GB, so a cell above the
  cap would reach across the interconnect and change its own walls. A cell over
  the cap records `oom`.
- **reps = 1.**

A cell that times out or OOMs **keeps every λ it finished**, each stamped with
the wall time at which it completed. That record is the point, not a
consolation: it gives an **anytime curve** — front quality against elapsed time
— for every method on every instance, which is the comparison that survives
when most comparators do not finish at length. A censored cell is a lower bound
on cost and a complete measurement of what that method had produced by the
limit.

## Parallel arms

At 8 cores, on a subset of about 6–8 instances spanning the length range and
including the 1273 rung and PRT-06 — enough for a speedup-against-length curve
without paying for the grid twice.

- **Threaded:** mRNAfold's TBB build, and `derna2-hull` / `derna2-walk` with
  `--threads 8`.
- **Parallel-λ:** each replay arm running its λs side by side. This is the
  comparison the sweep methods should win most clearly — independent λs are
  embarrassingly parallel, while derna2's walk shares work between gaps — so
  leaving it out would flatter derna2. The worker count is limited by the
  200 GB cap rather than by the core count for the memory-heavy arms, and the
  count each cell actually ran is recorded.

## Measurement

**One timed pass**, not E1's work-then-time split. At 12 h per cell a second
pass would double the grid's cost for correctness fields that the timed run
records anyway. The consequence is that a contaminated cell costs a full rerun,
which is what the retry budget is for.

Per cell and per λ: wall, CPU user + system, peak RSS, instructions retired,
exit category, and the designed sequence. E1's gate, contamination checks
(two-end average and worst 60 s stretch), slow-per-instruction check and retry
logic apply unchanged.

### The reserved socket moves to the even CPUs

E1 timed on the odd CPUs (socket 1). Another account now runs single-threaded
jobs there that this account cannot move, so **E2 reserves socket 0, the even
CPUs**, and this account's own work moves to the odd CPUs.

**What the reservation can and cannot deliver.** The premise that the other
account's jobs are confined to the odd socket no longer holds: two samples at
the end of this preparation put 4–5 of the 14 even cores under its pegged
single-threaded jobs, which set no affinity and so migrate across both
sockets. The reservation therefore removes only *this* account's neighbours. A
6 h waiter requiring all 14 even cores under 5% busy on three consecutive 30 s
samples never fired — 720 samples, mean 63–67% throughout.

`isolate.sh` was dry-run in the even-socket orientation and reports exactly
this: this account's jobs listed with their current masks, the other account's
listed as unmovable. One consequence to time deliberately — `--apply` moves
*every* process this account owns, so it would move a held-out walk currently
running on the even socket at ~589% CPU, changing that run's own timings from
that point on. Either let it finish first or accept the perturbation knowingly.

**A socket cannot be reserved against unpinned foreign jobs: they follow free
capacity.** This fell out as a natural experiment during the preparation rather
than being set up. Earlier, this account's own 8-thread walk held the even
socket and the other account's 8 jobs sat 4–5 even and 4 odd. That walk then
finished, and within the hour the same 8 jobs had migrated to **8 even and 0
odd**, stable across five samples 20 s apart, while this account's remaining
work sat 5–6 odd and 0 even.

They did not choose the even socket; they filled what was freed. That explains
why the quiet-socket waiter never fired, and it predicts what flipping sockets
will and will not buy: pinning E2 to a socket does not evict anyone, and 8
unpinned jobs facing 14 pinned cores on one side and 14 on the other will
spread across both. **Reserving a socket dilutes a foreign neighbour; it does
not exclude one.** Of the three ways out above, only the third — having those
jobs pinned off, which needs their owner or root — actually delivers a quiet
socket, and the choice between the other two is a choice about how to qualify
the walls, not about whether to.

It also means the even/odd orientation is not the lever it looked like. At the
time of writing the even socket holds 8 of their cores and the odd socket holds
none, so the orientation the plan specifies is the worse one *today* — but by
the same mechanism that is a consequence of where this account's own jobs are
pinned, and it will move again when they do. The orientation should be chosen
at launch from a fresh sample, not from this paragraph.

So the quiet gate as specified (`busy_max` 0.05) will not open while those jobs
run, and that is left as it is on purpose: it fails loudly rather than
admitting contaminated cells. Three ways out, and which one to take is the open
decision before launch — it changes what every wall in the grid means:

1. Wait for the other account's jobs to drain and keep the gate.
2. Raise the gate to the floor the host actually offers and carry the
   neighbour's share as a column on every row, treating walls as upper bounds.
3. Have the other account's jobs pinned off the even socket, which needs
   either their cooperation or root.

**What the load is worth, measured on this experiment's own arm.** The walk at
32 residues, same binary, same instance, same single core: **6.1 s at 10:01 and
12.3 s at 10:48**. Nothing about the work changed — a 2.0× spread in wall from
the host alone, which is larger than most of the differences between arms this
grid exists to measure.

Instruction counts are the fallback that survives this: they reproduce under
load where walls do not, and every row already carries them. A grid run under
option 2 can still answer the ordering questions on instructions while its
walls stay qualified.

The re-pin was also checked for cost at the short end, since its levers were
measured at n = 1000 and the grid starts at 32. Paired and alternating on one
core, the new defaults and `--no-fold-by-width --no-slice-width` run within 1%
of each other at 32 residues (12.25/12.38/12.40 s against 12.25/12.14/12.33 s),
and both return the same 11 vertices and 31 points. The levers are free where
they do not help, so there is no short-instance reason to turn them off.

The flip happens **at E2's launch**, not before: `isolate.sh --apply
0,2,4,6,8,10,12,14,16,18,20,22,24,26 1,3,5,7,9,11,13,15,17,19,21,23,25,27`,
with the config's `cores` and `rescore_cores` set to even CPUs. Other runs are
then pinned with `taskset -c 1-27:2`. Until launch the existing rule stands, so
nothing in flight is disturbed.

The reservation is a whole socket because the last-level cache and memory
bandwidth are per socket: a neighbour on the same socket slows a run without
touching its cores.

## Scoring the fronts

Every design from every arm is **re-scored with one folder and one codon
table**, then filtered for dominance. Self-reported numbers are kept beside the
re-score but are not compared across tools: the comparators disagree with each
other and with Turner 2004 (LinearDesign and LinearCDSfold cap loops at 20
against derna2's 30; LinearCDSfold reports energies below Turner for 17 E1
cells, documented and deliberately not debugged).

Re-scoring is also what makes a recovered point meaningful: a tool's own front
can contain points that are dominated once every design is scored on the same
model.

**Reference front** = `derna2-walk` where it finishes, which is the exact whole
front. Where it does not, the reference is the union of every arm's re-scored
points, Pareto-filtered — a lower bound on the true front, so recovery
fractions against it are upper bounds and are labelled as such.

## Metrics

1. **Hull vertices recovered** out of the exact hull's — the supported points,
   which is all a λ sweep can reach in principle.
2. **Pareto points recovered** out of the reference's — includes the
   non-supported points, which no λ sweep reaches at any resolution. This is
   the structural gap, not a speed difference.
3. **Hypervolume** against the reference, from a fixed reference point per
   instance (recorded, since hypervolume is meaningless without it).
4. **Hypervolume against wall time** — the anytime curve, per instance and arm.
5. **Wall and peak RSS against n**, per arm, with the fitted exponent.
6. **Finishing length** per arm: the longest instance it completed inside both
   limits.

## Before the grid runs

Done:

- **DERNA v1 builds.** Not through its own CMake — that wants cmake ≥ 3.1
  against this host's 2.8.12 and C++17 against g++ 4.8.5. `build_derna1.sh`
  compiles CMakeLists' source list directly with the conda g++ E1 used for
  mRNAfold. Its `-s` dispatch is crossed against its own help text: the mode
  the help calls the λ sweep (3) is the adaptive search, which is what the
  `derna1-own` arm runs, and the undocumented 4 is a uniform grid that aborts
  without `-l` and `-a`.
- **The codon-table converter.** `derna1_table.py` writes the canonical human
  table in v1's row format. v1 normalizes per amino acid
  (`codon_cai = log(usage / max_usage)`), so the absolute values are free and
  frequencies serve as well as counts.
- **The λ maps, verified at interior λs** for all six comparator arms against
  derna2's own objective: every arm matches exactly at 0.9, 0.7, 0.5, 0.3, 0.1.
  This is where both 100× scale errors were found.
- **The cost table** below, from E1's timing rows.
- **The driver**, `drive_e2.py`, reusing the core grid's measurement,
  isolation, core pool and re-scoring. `verify` is the map check above;
  `hull`, `sweep` and `summarize` are the grid, and `plan` states the cells the
  grid owes before it runs and checks them off after.

  `plan` exists because a short output file is indistinguishable from a small
  grid: it names all 539 cells (7 arms over 77 instances), says how many are
  present, how many are missing, how many of those are blocked on a hull that
  has not run, and whether anything is present that this config does not plan
  — which is what another grid writing into the same output root would look
  like. The counts go in the log before any result is read.
- **The binaries, frozen from one pinned commit.** derna2 at `acfaae60`
  (origin/xg tip) built in the detached worktree
  `/scratch/xinyugu4/paper-e2-src`; the comparators are E1's freeze, their
  hashes re-verified against that directory's own record and all six matching.
  `PROVENANCE.md` in `/scratch/xinyugu4/paper-e2-bin` carries every hash. E1
  pinned derna2 at `17ee9d8b`, so the two experiments' derna2 numbers are not
  interchangeable.

  The pin moved once already, from `d3aa667a` to `acfaae60`, when
  `experiment/walk-at-length` merged into xg (#414). That merge carries
  `--fold-by-width` and `--slice-width 250` — 0.69× the instructions at
  n = 1000 — which the earlier pin could not include, and they are the walk's
  **defaults** there, so the arm now runs the walk as it ships rather than a
  slower stand-in. `mfe_fold` rebuilt byte-identical across the two pins, so
  re-scored values are comparable either way; the other two binaries changed.

  One lever is still outside the pin: the hull-walk anchor schedule
  (`ad963f76`, unmerged). It only makes the walk cheaper, so the arm remains an
  upper bound by that one lever.

  **Every row records the build that produced it** — the pin from the config
  and the first 12 hex of each derna2 binary's hash from disk, so a config
  naming a pin its binaries were not built from is visible in the rows. A
  re-pin is the normal case here, not an exception, and a row that does not say
  which build produced it cannot be separated from one that does afterwards.
  `sweep` treats a cell from another pin as missing rather than done, `plan`
  counts them separately, and `summarize` refuses to read a mixed file as one
  experiment.

- **The whole pipeline re-run on the frozen binaries**, hull through
  summarize, reproducing the table below.

Still to do:

1. **A DERNA v1 cost pilot** — running at the time of writing over every
   instance up to 102 residues. v1's cost at these lengths is measured nowhere,
   and its adaptive sweep chooses its own λs, so the λ count per instance is
   unknown too; the count is deterministic and survives a loaded host even
   where the walls do not. Its cells are not affected by the re-pin: v1's own
   binary and codon table are unchanged and `mfe_fold` is byte-identical across
   both pins, so nothing those rows depend on moved. They are stamped with the
   pin in force when they ran all the same, which for the cells already
   finished is no pin at all.

2. **Re-running the 32-residue cells under the new pin.** The eight cells
   measured before the re-pin used the earlier walk and chord binaries, so they
   are measurements of a different derna2 and `plan` already counts them as
   missing. The table below is from those cells and is superseded as soon as
   the re-run lands.
3. **A concurrency pilot**, alone against the chosen concurrency on a few
   cells, as E1 did (E1 measured a 1.09× median slowdown at 4 and chose 2).
   This one needs a quiet host and has not run.
4. **The quiet-socket decision** above. Both sockets were loaded throughout
   this preparation, by this account's jobs and another's; nothing here was
   timed on a quiet host.

## The pipeline, run end to end

On the shortest instance (32 residues), every command and all seven arms:

| arm | λs | wall | Pareto points | hull vertices |
|---|---|---|---|---|
| mrnafold-hull | 11 | 0.08 min | 11 | 11/11 |
| derna2-hull | 11 | 0.09 min | 11 | 11/11 |
| **derna2-walk** | **one run** | **0.10 min** | **31** | **11/11** |
| lineardesign-hull | 11 | 0.18 min | 11 | 11/11 |
| lcds-exact-hull | 11 | 0.45 min | 11 | 11/11 |
| derna1-hull | 11 | 2.02 min | 11 | 11/11 |
| **derna1-own** | **21** | **3.76 min** | **9** | **9/11** |

The walk ran in 6.3 s on one core at 31 MB peak, giving 11 vertices, 10
breakpoints and a front of 31 points, and every point re-scores to the energy
and CAI it reported.

Two rows carry the experiment's two questions, at the only length so far.
**The walk returns 31 Pareto points where every λ-replay arm returns 11, for
the same wall** — a λ sweep can reach a hull vertex and nothing else, so the
20 off-hull points are not slow for the comparators, they are unreachable. And
**v1's own λ search spends 21 λs and 3.76 min to recover 9 of the 11 vertices,
where the same tool handed derna2's 11 λs recovers all 11 in 2.02 min.**

Paying for the λs costs ~1.9× the wall and loses two vertices. Its sweep
bisects from hardcoded corners at `EPSILON = 1e-5`
(`params/constants.h`) and stops on its `-t`/`-p` thresholds, so neither its
corners nor its stopping rule are the ones the replay arms use — at n beyond
about 560 its MFE-corner λ of 1 − 1e-5 is over the ε bound this plan derives,
so its corner is not the lexicographic optimum there. One instance is not a
trend; the shape is what the grid will test.

## Three ways DERNA v1 fails silently

All three were found by running it, and all three would have corrupted results
rather than stopped the run.

- **It dumps its DP tables to fixed names in its working directory**
  (`E1.bin`, `M1.bin`, `TM.bin`, … via `save_all_vectors("./")`). Two v1 runs
  sharing a cwd overwrite each other with no diagnostic, which the parallel-λ
  arm and any concurrency above 1 would do by default. Every v1 run gets its
  own directory.
- **Its `-s` dispatch is crossed against its own help text.** `-h` lists modes
  1–3 and calls 3 "lambda sweep"; `main.cpp` wires 3 to `lambda_sweep_2`, the
  adaptive search, and the undocumented 4 to `lambda_sweep`, a uniform grid
  that needs `-l` and `-a`. The `derna1-own` arm runs 3.
- **Its sweep CSV is flushed only at close, while its stdout is flushed per
  λ.** A run killed at the limit can leave an empty CSV however many λs it
  finished — and this is the arm most likely to be killed, so the anytime
  curve it exists to provide would be the one thing lost. The driver falls
  back to the λs in stdout, which carry no design and so cannot be re-scored:
  those rows are marked `reported_only` and are kept for the anytime curve on
  reported values alone. Its two modes also disagree on case for the same
  fields (`mfe`/`cai_s` against `MFE`/`sCAI`), which drops every row of one
  mode if the parser is written against the other.

## Cost, from E1's measured per-λ walls

Fitted on E1's single-core `ok` rows, 99 points per arm over n = 30–1273, the
slower corner kept per length. Per-λ wall at the lengths E2 adds, and the
number of λs that fit inside one 12 h cell:

| arm | wall exp | per λ at 1273 | at 1629 | λs in 12 h at 1273 | at 1629 |
|---|---|---|---|---|---|
| derna-exact-30 | 2.05 | 7 min | 12 min | 99 | 60 |
| lineardesign | 1.99 | 20 min | 32 min | 37 | 23 |
| mrnafold | 2.26 | 17 min | 30 min | 43 | 24 |
| lcds-exact | 2.16 | 1.5 h | 2.6 h | 8 | 5 |
| cdsfold | 2.43 | 47 min | 1.4 h | 15 | 8 |

Peak RSS over the same fit, per concurrent process at 1629: LinearDesign 79 GB,
derna2 exact 19 GB, LinearCDSfold exact 14 GB, mRNAfold 6 GB. **LinearDesign is
what the 200 GB cap binds on** — two concurrent cells of it at Cas length is
158 GB, so the parallel-λ arm cannot run many workers there whatever the core
count says.

This settles the shape of the result and is why the anytime curve is the
measurement rather than a fallback: a front with tens to low hundreds of
vertices is **out of reach inside 12 h for every comparator at Cas length** —
LinearCDSfold exact gets 5 λs, LinearDesign 23, mRNAfold 24. The experiment
will mostly report partial fronts there, which is a finding about sweeping as a
route to the front, not a missing measurement. The fit residuals are large in
log (0.27–0.89 across arms, worst for `lcds-exact`), so treat every figure as
an order of magnitude for admission and ordering, not a prediction.

The derna2 arms are not in the table: E1 never ran a front walk, so their cost
comes from the walk's own measurements (n = 1000 in 2 h 21 min at 83.5 GB,
growing about n^2.85), which is the one arm whose cost per *front* rather than
per λ is already known.

### DERNA v1: the λ count cannot be costed, the per-λ wall can

The pilot ran v1's own sweep over every instance up to 102 residues. Two
quantities, and they behave oppositely.

**The λ count does not track length.** It ranges 21–32 over n = 32–102 with no
order in n at all:

| n | 32 | 41 | 51 | 64 | 78 | 81 | 91 | 98 | 102 |
|---|---|---|---|---|---|---|---|---|---|
| λs | 21 | 22 | 21 | 28 | 25 | 27 | 24 | 21 | **32** |

64 residues spends 28 and 98 spends 21. The count tracks where v1's own
thresholds stop bisecting on that instance's front, so there is nothing to fit
and the 12 h limit is *enforced* for this arm rather than planned around. Two
probes outside the pilot agree: n = 20 stops at 15, and n = 50 had already
passed 23.

**The per-λ wall fits cleanly.** Dividing the varying count out leaves a
tighter fit than most arms in the table above — exponent **2.01**, worst
residual 0.24 in log over 8 points (the `P13501` cell is dropped: its wall is
the one this experiment's own probe contaminated).

| arm | per-λ exponent | per λ at n = 102 | at 1273 | at 1629 | λs in 12 h at 1629 |
|---|---|---|---|---|---|
| derna2 exact | 2.05 | 2.5 s | 0.1 h | 0.2 h | 60 |
| mRNAfold | 2.26 | 3.3 s | 0.3 h | 0.5 h | 24 |
| LinearDesign | 1.99 | 7.6 s | 0.3 h | 0.5 h | 23 |
| LinearCDSfold exact | 2.16 | 23.2 s | 1.5 h | 2.6 h | 5 |
| **DERNA v1** | **2.01** | **124.1 s** | **5.5 h** | **9.1 h** | **1.3** |

So v1 is not asymptotically worse than anything here — its exponent sits
between LinearDesign's and derna2's. It is **50× slower per λ** at the same
length, and that constant is what decides the experiment: against a demand of
21–32 λs, **v1 affords between one and two of them in 12 h at Cas length.**
It is the arm that times out hardest, worse than LinearCDSfold exact.

That conclusion survives the obvious objection. The pilot's walls were measured
on a loaded host, which inflates them, and the extrapolation reaches 16× past
the fitted range. Halving every figure for load still leaves 2–3 λs against a
demand of 21–32. The anytime curve is v1's entire output at long length, which
is a finding about adaptive sweeping and not a missing measurement.

**And about half of what v1 spends buys nothing.** Its λs convert to distinct
front points at a ratio that barely moves across the pilot:

| n | 32 | 41 | 51 | 64 | 78 | 81 | 91 | 98 | 102 |
|---|---|---|---|---|---|---|---|---|---|
| λs | 21 | 22 | 21 | 28 | 25 | 27 | 24 | 21 | 32 |
| points | 9 | 11 | 10 | 14 | 13 | 13 | 14 | 11 | 14 |
| λs per point | 2.33 | 2.00 | 2.10 | 2.00 | 1.92 | 2.08 | 1.71 | 1.91 | 2.29 |

1.71–2.33 over nine instances is a mechanism rather than noise: a bisection has
to evaluate both sides of a breakpoint to establish it, so roughly every second
λ returns a design the sweep already holds. None of its points are dominated —
distinct designs and Pareto points coincide — so this is pure duplication, not
wasted exploration.

That is the arithmetic the two v1 arms were built to separate. Handed derna2's
λs, v1 gets **one hull vertex per λ**. Left to its own search it gets **one
point per two λs** and still misses vertices, 9 of 11 at 32 residues. The λ
search costs about 2× in folds and loses coverage at the same time.

## Three ways DERNA v1 fails silently

All three were found by running it, and all three would have corrupted results
rather than stopped the run.

- **It dumps its DP tables to fixed names in its working directory**
  (`E1.bin`, `M1.bin`, `TM.bin`, … via `save_all_vectors("./")`). Two v1 runs
  sharing a cwd overwrite each other with no diagnostic, which the parallel-λ
  arm and any concurrency above 1 would do by default. Every v1 run gets its
  own directory.
- **Its `-s` dispatch is crossed against its own help text.** `-h` lists modes
  1–3 and calls 3 "lambda sweep"; `main.cpp` wires 3 to `lambda_sweep_2`, the
  adaptive search, and the undocumented 4 to `lambda_sweep`, a uniform grid
  that needs `-l` and `-a`. The `derna1-own` arm runs 3.
- **Its sweep CSV is flushed only at close, while its stdout is flushed per
  λ.** A run killed at the limit can leave an empty CSV however many λs it
  finished — and this is the arm most likely to be killed, so the anytime
  curve it exists to provide would be the one thing lost. The driver falls
  back to the λs in stdout, which carry no design and so cannot be re-scored:
  those rows are marked `reported_only` and are kept for the anytime curve on
  reported values alone. Its two modes also disagree on case for the same
  fields (`mfe`/`cai_s` against `MFE`/`sCAI`), which drops every row of one
  mode if the parser is written against the other.

## Cost, from E1's measured per-λ walls

Fitted on E1's single-core `ok` rows, 99 points per arm over n = 30–1273, the
slower corner kept per length. Per-λ wall at the lengths E2 adds, and the
number of λs that fit inside one 12 h cell:

| arm | wall exp | per λ at 1273 | at 1629 | λs in 12 h at 1273 | at 1629 |
|---|---|---|---|---|---|
| derna-exact-30 | 2.05 | 7 min | 12 min | 99 | 60 |
| lineardesign | 1.99 | 20 min | 32 min | 37 | 23 |
| mrnafold | 2.26 | 17 min | 30 min | 43 | 24 |
| lcds-exact | 2.16 | 1.5 h | 2.6 h | 8 | 5 |
| cdsfold | 2.43 | 47 min | 1.4 h | 15 | 8 |

Peak RSS over the same fit, per concurrent process at 1629: LinearDesign 79 GB,
derna2 exact 19 GB, LinearCDSfold exact 14 GB, mRNAfold 6 GB. **LinearDesign is
what the 200 GB cap binds on** — two concurrent cells of it at Cas length is
158 GB, so the parallel-λ arm cannot run many workers there whatever the core
count says.

This settles the shape of the result and is why the anytime curve is the
measurement rather than a fallback: a front with tens to low hundreds of
vertices is **out of reach inside 12 h for every comparator at Cas length** —
LinearCDSfold exact gets 5 λs, LinearDesign 23, mRNAfold 24. The experiment
will mostly report partial fronts there, which is a finding about sweeping as a
route to the front, not a missing measurement. The fit residuals are large in
log (0.27–0.89 across arms, worst for `lcds-exact`), so treat every figure as
an order of magnitude for admission and ordering, not a prediction.

The derna2 arms are not in the table: E1 never ran a front walk, so their cost
comes from the walk's own measurements (n = 1000 in 2 h 21 min at 83.5 GB,
growing about n^2.85), which is the one arm whose cost per *front* rather than
per λ is already known.

### DERNA v1 cannot be costed from length at all

v1 is not in the table either, and the pilot says it cannot be put there. Its
adaptive sweep picks its own λs, and the count is **instance**-dependent rather
than length-dependent:

| n | instance | λs | wall | per λ |
|---|---|---|---|---|
| 20 | ladder | 15 | — | — |
| 32 | ladder | 21 | 3.8 min | 10.7 s |
| 50 | ladder | **≥ 23** | — | — |
| 98 | UniProt | 21 | 33.5 min | 95.8 s |
| 102 | ladder | **32** | 61.0 min | 114.3 s |

The count is not monotone in n and does not track it: 50 residues had already
passed 23 λs where both 32 and 98 stop at 21. It tracks where the tool's own
thresholds stop bisecting on that instance's front. A length-fitted model of
v1's cost would therefore be fitting the wrong variable, and the 12 h limit has
to be *enforced* per cell rather than planned around for this arm.

Two of those rows are bounds rather than values, and the reason is worth
recording: the n = 50 probe was pinned to a core the pilot's own pool was
handing out, so it spent 25 minutes splitting one core with the `P13501` cell
at 57%/50%. The probe was stopped at 23 λs — the finding needs only that 23
exceeds 21 — and `P13501`'s **wall** is contaminated by this experiment's own
probe and must be re-measured before it is quoted. Its λ count is unaffected,
being deterministic. The mistake is the ordinary one: a side probe pinned
inside the pool the main run draws from is a neighbour, not a bystander.

It also means the λ count and the vertex count are unrelated quantities. At 32
residues v1 spends 21 λs — nearly twice the 11 vertices — and still recovers
only 9 of them, because its λs are bisected on its thresholds rather than
placed one per hull interval. Spending more λs than there are vertices while
recovering fewer is not a paradox; it is what an adaptive sweep with no access
to the hull does.
