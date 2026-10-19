# Paper benchmark — the core grid

Correctness, single-λ runtime and memory, and the approximate modes' speed
against quality, over one grid of instances × corners × arms. The arm set is
`e1.config.json`; what each cell records is the row format `drive.py` writes.

## Two passes, split by what load does to a measurement

| pass | measures | survives host load | where it runs |
|---|---|---|---|
| `work` | correctness (re-scored MFE and CAI, back-translation), peak RSS, `VmHWM`, instructions, cycles | yes, except cycles (±4%) | any cores, in parallel, now |
| `time` | wall-clock, CPU user + system | no | a reserved socket, one cell at a time, gated |

Instructions retired repeat to ~0.03% at load 40–50 and peak RSS is set by
what is resident, so the `work` pass needs no quiet host and carries most of
the grid's cost. Wall-clock time drifts 9–11% on a busy host even when phase
shares do not, so the `time` pass runs only once the reserved cores are idle
and runs only the cells the `work` pass finished, plus any it records as
`skipped`: a cell too long to be worth its work-pass run is still timed, and
the `time` pass re-scores its output. Its row carries `wall_est_s` so the
longest-first order still holds.

## Files

- `measure.py` wraps one run: hardware counters opened by the wrapper itself
  (so concurrent cells never count each other), `wait4` rusage, a process-group
  timeout, and exit categories `ok` / `timeout` / `oom` / `crash`.
- `drive.py` runs the grid: `plan`, `selftest`, `work`, `time`, `summarize`.
- `isolate.sh` moves this account's processes off the reserved socket.
- `config.example.json` is a full configuration; copy it and fill in the
  binaries.

## Corners

The corners run at λ = ε and 1 − ε, not 0 and 1, so each is the
lexicographic optimum: the second objective breaks ties and never trades
against the first. LinearDesign's λ is mapped as λ_LD = (1 − λ)/λ.

- **ε has an upper bound in n.** At the MFE corner the CAI term must not pay
  for one dcal; at the CAI corner no energy gain may pay for leaving a
  CAI-optimal codon. On the human table that is ε < 2.0e-6 at n = 2828
  (`epsilon_bounds`). The 1e-5 used before is over the MFE-corner bound past
  n ≈ 560.
- **ε has a floor set by the tools.** LinearCDSfold parses λ with `stof`, so
  1 − ε is rounded to single precision, which steps by 2⁻²⁴ below 1: at
  ε = 1e-8 it reads λ = 1 and drops the CAI tie-break (CAI 0.694 against the
  optimum's 0.723 at n = 102, same MFE).
- **ε = 2⁻²³ ≈ 1.19e-7** sits inside both: exact in single and double
  precision, ε/2 still exact, 17× under the bound at n = 2828.
- At the CAI corner every exact arm must re-score to CAI = 1 exactly, since
  every amino acid has a codon of weight 1; `summarize` reports any that do not.

## mRNAfold

mRNAfold reads a config file rather than flags, so the driver writes one per
cell, and writes the canonical table once in its `Ala GCT 0.26` format.
Its λ maps like LinearDesign's, λ' = (1 − λ)/λ, and it parses λ in double
precision, so ε has no tool floor there.

- It prints `Score`, the objective (energy plus λ'·Σ −ln w), and `EFE`, an
  ensemble free energy from a partition function run on the design after the
  search. EFE is not the MFE: the reported MFE is recovered as
  Score + λ'·Σ ln w. The partition function is part of its wall time.
- Only the TBB build runs `parallel true` in parallel; `std::execution` falls
  back to serial in libstdc++ without TBB. The threaded arm uses that build,
  and TBB sizes its pool from the affinity mask the driver sets.

## Isolation without root

The cpuset cgroup is root-owned, so the only lever is `taskset` on processes
this account owns. The host has two sockets of 14 cores, no SMT, one NUMA node
per socket. Reserve a **whole socket** (odd CPUs, node 1) rather than a few
cores: the last-level cache and memory bandwidth are per socket, so a neighbour
on the same socket slows a run without touching its cores.

`isolate.sh --apply` moves **every** process this account owns to socket 0,
drivers included, because a sweep driver idling between cells would otherwise
launch its next binary onto the reserved socket. It moves only processes that exist
when it runs: a job this account starts later is placed on all cores again,
and one landing on the reserved socket stalls the time pass at its gate until
it ends. Another account's process can land there too; the `time` pass catches that by reading the reserved cores
themselves:

- **Gate:** before the first run, every reserved core below 5% busy on 20
  consecutive 30 s samples; a 2 × 5 s re-check of the reserved cores not held
  by a run of ours before each new instance and corner. The gate does not read
  the load average, which stays high while the other socket is busy and also
  counts the driver's own threads.
- **Contamination:** a run is flagged if a reserved core that none of this
  driver's own runs held was busy while it ran. Two readings: the average over
  the run, from `/proc/stat` at its two ends (over 5%), and the busiest 60 s
  stretch of 5 s sampler intervals inside it (over 25%). The second catches a
  burst that a long run would average away: a 5-minute neighbour is 8% of a
  1-hour run but 2% of a 4-hour one. It is a 60 s mean rather than a single
  interval because a host daemon that cannot be moved (root's configuration
  management) spikes one core to 20–70% for 5–15 s every few minutes, which
  flagged nearly every run of more than a few minutes. A single-core run longer than 5 s that got under 95% of its
  core is flagged too.
- **Retry:** a flagged run is kept, the free reserved cores must be quiet for
  6 × 10 s, and the cell runs again, up to `time.retries` times.
  `--redo-contaminated` gives cells that ran out of retries another round.
- **Slow per instruction:** instructions retired do not move with load, so a
  method's time per instruction should change smoothly with n. `summarize`
  lists runs over 1.2× the median of the same method within ±25% of their
  length: a slowdown the core-busy checks cannot see, such as a neighbour's
  memory traffic on the same socket.
- **Order:** cells of one instance and corner run together, arms rotated per
  rep, and instances longest first, so the run that starts last is short.

`time.concurrency` runs that many cells at once on the reserved socket. They
share its last-level cache and memory bandwidth, so every wall reads somewhat
above a cell run alone; comparisons between methods hold only if the slowdown
is similar across methods. Measure it on a few cells, alone and at the chosen
concurrency, before committing to a level.

## Run order

1. Build every binary from one pinned commit into a frozen directory on
   `/scratch` (the home volume is full) and point the config at it.
2. `python3 drive.py selftest --config C` — every enabled arm on a 30-residue
   protein must exit, parse, translate back, and re-score to its reported MFE
   (exact arms) or below it (beam arms, whose report is the fold the beam found).
   CDSfold's parser is unverified until this passes.
3. `python3 drive.py epscheck --config C --max-n 300` — ε against its bound
   at the longest instance, and every exact arm unchanged at ε/2.
4. `python3 drive.py plan --config C` — the rows each pass must write.
5. `python3 drive.py work --config C` — now, under load.
6. `bash isolate.sh` (dry run), then `bash isolate.sh --apply`.
7. `taskset -c 0,2,4,…,26 python3 drive.py time --config C` — waits for the
   gate, then runs. The driver itself goes on the other socket so its own
   threads never read as a neighbour; each run is pinned to the reserved cores
   explicitly.
8. `python3 drive.py summarize --config C` — row counts against the plan, exit
   categories, contamination, and exact arms' re-scored MFE against the
   reference arm. Check the row counts before reading any result.

## Memory on a shared host

The `work` pass admits a cell only when its estimated peak fits both
`mem_budget_gb` and the host's `MemAvailable` less `mem_headroom_gb`: an
over-admitted cell would bring the kernel's OOM killer down on other people's
jobs. The estimate is `k·n²` per arm variant, starting from a prior and raised
to the largest `k` any finished cell of that variant has shown. A cell whose
estimate exceeds the whole budget is recorded `deferred-memory` and is not
run alone; it is not censored, and a later run with a larger budget picks it
up.

Both passes are resumable and hold a lock per output directory, so a second
driver on the same grid refuses to start instead of overwriting rows.

## Known tool discrepancies

Every arm's design is re-folded by one folder (`mfe_fold`, Turner 2004,
internal loops up to 30 unpaired), and correctness is judged on that re-scored
MFE, not on the energy a tool reports.

- **Loop cap 20 (LinearDesign, LinearCDSfold).** Both tools cap internal loops
  at 20 unpaired nucleotides, as `derna-exact-20` does. A capped tool's
  re-scored MFE can only be equal to or lower than its reported one, and on
  some instances the capped optimum is worse than the cap-30 one: Q9HAE3 at the
  MFE corner, where LinearDesign, LinearCDSfold exact and `derna-exact-20` all
  return the same design at −369.4 kcal/mol and the cap-30 arms reach −369.9.
- **LinearCDSfold scores some structures below the Turner 2004 model.** In 17
  cells of the grid a LinearCDSfold design re-scores *higher* than LinearCDSfold
  reported. A larger cap searches a superset of structures, so the cap cannot
  cause this; the two energy models disagree on some loop. No other tool shows
  it. The cells:
  - exact and beam width 2000, CAI corner, P0DTC2 prefixes n = 508, 549, 640,
    690, 745, 804, 868: +0.2 kcal/mol each;
  - Q5SWX8, MFE corner: exact −858.5 reported, −857.1 re-scored (+1.4); beam
    width 500 −843.1 → −842.1 (+1.0);
  - P0DTC2 n = 1012, MFE corner, beam width 100: −1869.0 → −1867.7 (+1.3).

  The loop responsible is not identified. Each structure is in the cell's
  `work/<cell>.out.txt`. On Q5SWX8 the effect makes LinearCDSfold exact miss the
  optimum that `derna-exact-20` and LinearDesign reach at the same cap
  (−858.5), so `summarize` lists it as an MFE mismatch.

## The front-enumeration experiment

This grid measures one λ per cell. A second experiment measures what the whole
front costs: every comparator reaches the front only by sweeping its
scalarization parameter, while derna2 enumerates it. Its design — arms, the λ
map per tool and how it is verified, the instance set (the 55 UniProt proteins,
every 3rd ladder rung, and five Cas orthologs to 1629 aa), the 12 h / 200 GB
limits with partial fronts kept as anytime curves, and the move of the reserved
socket to the even CPUs — is in [`E2-PLAN.md`](E2-PLAN.md). It reuses this
driver's measurement, isolation and re-scoring machinery.
