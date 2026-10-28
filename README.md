# derna2-analysis

The analyses and figures of the derna2 paper: exact Pareto-optimal mRNA codon design, its cost
against the other design tools, and constrained design with UTRs. The method itself lives in
[derna2](https://github.com/elkebir-group/derna2); this repository holds what the paper
measured with it and how.

- **`derna2_analysis/`**: the shared library: result rows (`io.py`), Pareto fronts, hull
  vertices and hypervolume (`fronts.py`), checksums against `pins/` (`provenance.py`), one module
  per experiment (`e1.py`, `fig1.py`) and the figure style (`figures.py`).
- **`pipelines/`**: a Snakemake pipeline from the inputs, or from the deposited results, to the
  figures. See "Running the pipeline".
- **`drivers/`**: the benchmark drivers that produced the results, one directory per
  experiment, copied from the derna2 repository (`drivers/README.md` names the commits and the
  changes).
- **`pins/`**: what every result depends on, as checksums: the frozen binaries of each
  experiment, their provenance, and the deposited results.
- **`inputs/`**: the protein sequences, UTRs, codon table and E3's motif sets, committed.

## Install

Python 3.11–3.13. On a host whose glibc is older than 2.28 (CentOS 7), the numpy cap in
`pyproject.toml` keeps pip on wheels it can install. ViennaRNA has no wheel for such a host and
its source build fails there; the overview figures, the only targets that import it, then need a
Python whose ViennaRNA 2.7.2 bindings come from elsewhere (e.g. conda's `viennarna`), with the
rest of `requirements.lock` installed beside it.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.lock
.venv/bin/pip install -e '.[dev,pipeline]' --no-deps
```

`requirements.lock` pins every package version the published figures were made with. Install
from it, not from `pyproject.toml`'s ranges, to reproduce them.

## Data

The deposited results go under `data/results/` (or set `results_root` in
`pipelines/config/host.yaml`). The figure targets check every file against
`pins/results-<experiment>.sha256` before drawing anything.

```
data/results/
├── e1/
│   ├── config.json          # the grid as run (arms, instances, ε, limits, host paths)
│   ├── work.jsonl           # correctness, peak RSS, instructions: 3,465 rows, one per cell
│   ├── time.jsonl           # wall-clock and CPU time on a reserved socket: 3,918 rows (retries)
│   └── structures.jsonl     # every design: sequence, tool structure, re-folded structure
├── e2/                      # in progress: no checksum manifest yet
│   ├── sweep.jsonl          # one row per λ, per front point and per cell, stamped with its pin
│   └── hull/                # derna2's hull vertices and replay λs per instance (partial if the
│                            # walk was stopped)
├── e3/                      # in progress: no checksum manifest yet
│   ├── e3-probe/            # hull probe round 1: per rung the walk's rows, its measurement
│   │                        # record and summary.csv; forced.jsonl, the forced-motif census
│   ├── e3-probe2/           # round 2: the n = 200 load control, n = 150, 250 and 300
│   ├── grid/, mid/, long/   # the short, middle and long bands: summary.csv and each cell's
│   │                        # front by round (`<protein>_<nest>*.jsonl`)
│   └── screen/              # binding.jsonl, the binding screen; f0/, unconstrained fronts
├── fig1/
│   ├── a_cdsonly_front.jsonl                    # exact CDS-only front, hull vertices marked
│   ├── b_cdsonly_with_utrs.scored.jsonl         # those designs folded with the UTRs
│   ├── cd_symmetric_guarded_front.scored.jsonl  # exact guarded front
│   ├── d_so_front.jsonl                         # exact structurally optimal front
│   ├── so_windows/          # the window walks that found it, by round
│   └── top_exhaustive/      # the stretch above the most adapted seed, by enumeration
├── fig1_hba/                # not deposited yet: the overview figure for the hba UTR pair
│   ├── run3-hba/            # the GLPE arms' results for that pair (a0 CDS-only front and its
│   │                        # refolds, a2 masking front, a3 mRNAfold)
│   └── xhoi/so.jsonl        # the structurally optimal front with XhoI forbidden
└── fig1_fronts/
    ├── <arm>/               # unguarded, symmetric, no_clamp, unpaired: report.json, the front
    │                        # and its hull vertices re-scored; front3.jsonl and its designs
    ├── mrnafold/compare.jsonl  # mRNAfold at each hull vertex's weight
    ├── plot/                # plot_data.json and fronts.png
    └── raw/                 # every run output untrimmed, and the binaries' commit and SHA-256
```

Nothing under `data/` is committed.

## Running the pipeline

`pipelines/run.sh` runs `pipelines/Snakefile` with `pipelines/config/config.yaml`, then
`pipelines/config/host.yaml` if it exists (copy `host.example.yaml`). A target is a rule name or
an output path relative to the repository root.

```bash
pipelines/run.sh --dry-run                                     # what would run
pipelines/run.sh                                               # every figure, from the deposit
pipelines/run.sh -- fig1_figures                               # one experiment's figures
pipelines/run.sh -- data/derived/figures/e1_scaling.pdf        # one output
pipelines/run.sh -- e1_runs                                    # regenerate E1 from binaries
```

`run.sh` passes `--rerun-triggers mtime`, so a code or config change never reruns an existing
output by itself: delete it, or add `--forcerun <rule>`. Outputs go to `data/derived/`.

There are two kinds of target. **Figure targets** read results and need no binaries; `all` is
these. Set `figures_from: runs` to draw from this checkout's own runs instead of the deposit.
**Run targets** regenerate the results from frozen binaries. They take hours to days, and the
timed passes need a reserved socket, so they run only when named.

### Provenance

| Rule | Output | Content |
|---|---|---|
| `check_results` | `stamps/results_{exp}.ok` | every deposited file of an experiment against `pins/results-{exp}.sha256` |
| `check_bins` | `stamps/bins_{exp}.ok` | the frozen binaries in the experiment's `bin_dir` against `pins/{exp}.sha256` |
| `build_derna2` | `bin/derna2-{pin}/` | derna2's benchmark binaries rebuilt from a fresh clone at a commit |
| `check_derna2_build` | `bin/derna2-{exp}.compare.txt` | whether that rebuild is byte-identical to the frozen binaries of the experiment |
| `package_deposit` | `deposit/{exp}.tar` | the deposit of one experiment: exactly the files its manifest names, checked, without the run logs and working scripts beside them |

### E1: single-λ correctness, runtime and memory

Every method at both corners of the objective (λ = ε and 1 − ε, ε = 2⁻²³) on 55 UniProt
proteins up to 500 aa and a 50-rung prefix ladder of P0DTC2 up to 1,273 aa. Arms: derna2 exact
(internal-loop cap 30 and 20), its beam modes, 8 threads; LinearDesign; LinearCDSfold exact and
beam; CDSfold (MFE corner only); mRNAfold serial and threaded. Every design is re-folded by one
folder, and correctness is judged on that re-scored MFE. `drivers/e1/BENCHMARK.md`
describes the two passes, the socket isolation and the known tool discrepancies. derna2 is
pinned at `17ee9d8b`.

| Rule | Output | Content |
|---|---|---|
| `e1_tables` | `e1/cells.tsv`, `e1/correctness.tsv` | one row per cell; per method and corner: cells, re-scored MFE equal to / above / below derna2 exact, the largest gap, CAI = 1 at the CAI corner, cells without a clean timed run |
| `e1_plot_scaling` | `figures/e1_scaling.pdf` | wall time and peak RSS against length, exact serial methods, both corners |
| `e1_plot_speed_quality` | `figures/e1_speed_quality.pdf` | the beam modes: MFE gap to the exact optimum and wall time relative to derna2 exact, per width |
| `e1_figures` | | the three above |
| `e1_config` | `runs/e1/config.json` | `drivers/e1/e1.config.json` with this host's paths |
| `e1_selftest` | `runs/e1/selftest.log` | every arm on a 30-residue protein: exits, parses, translates back, re-scores |
| `e1_epscheck` | `runs/e1/epscheck.log` | ε against its bound, and every exact arm unchanged at ε/2 |
| `e1_work` | `runs/e1/work.jsonl` | the work pass: correctness, memory, instructions; load-insensitive |
| `e1_time` | `runs/e1/time.jsonl` | the time pass on the reserved socket, gated on it being idle |
| `e1_summarize` | `runs/e1/summary.txt` | row counts against the plan, exits, contamination, exact arms against the reference |
| `e1_runs` | | the run chain above |

### E2: the whole front

What the whole Pareto front costs each method. derna2 enumerates it (`derna2-walk`, the walk
`derna front` ships: `--pass-grid auto --seed-sweep --cut-probes --lend-idle`); every other
method sweeps its scalarization weight at derna2's hull λs (one per hull vertex plus the
corners), and DERNA v1, through its Rust port, also by its own adaptive sweep. derna2's walk
and mRNAfold run once on one core, against every other single-core arm, and once on eight
(`-par`). Instances: the 55 UniProt proteins up to 500 aa, twelve native proteins of
511–962 aa (`inputs/e2-natives/`), P0DTC2 whole (1,273 aa) and five Cas orthologs (PRT-03 to
PRT-07, 1,000–1,629 aa). Limits: 2 h per cell below 1,000 aa (soft for the walk, which runs on
to 12 h), 6 h above (12 h for `derna2-hull` and the walk); a cell is admitted only while its
estimated peak memory fits the host's budget. A walk the limit stops keeps its partial front,
and the replay arms then run at its partial hull. Design: `drivers/e2/E2-PLAN.md`. The grid
runs on a second host; **the pin is `pending`** and the binaries are not frozen yet (the walk
and the replay arms currently run different derna2 builds), so `check_bins` refuses
`pins/e2.sha256` and the run targets stop there.

| Rule | Output | Content |
|---|---|---|
| `e2_plot_fronts` | `figures/e2_fronts.pdf`, `e2/cells.tsv` | wall time and Pareto points recovered per cell against length, eight-core arms hollow and dashed; rows of the configured pin only |
| `e2_fig2_data` | `e2/fig2.json` | Figure 2's data: per arm, wall, peak RSS, designs and seconds per design on the cells every arm completed, with completions, timeouts and the longest protein finished; per walk front, the share of its points no weight reaches, from hull vertices computed on the points, and the fronts whose stored `vertices` disagrees; rows of the configured pin only |
| `e2_figures` | | the two above |
| `e2_config` | `runs/e2/config.json` | `drivers/e2/e2.config.json` with this host's paths and the pin |
| `e2_verify` | `runs/e2/verify.done` | each comparator's λ map against derna2 at interior λs |
| `e2_grid` | `runs/e2/grid.done` | `run_e2_grid.sh`: refuses to start while a timing core or its hyperthread sibling is busy, then the four passes in order (single-core walk, single-core replays, eight-core walk, eight-core replays) into `runs/e2/sweep.jsonl`; resumable |
| `e2_summarize` | `runs/e2/summary.txt` | fronts, recovery against the reference, anytime curves, for the configured pin |
| `e2_runs` | | the run chain above |

### E3: the price of motif-freedom

What excluding motifs costs: the exact Pareto front over every design, then over the designs
that avoid each of four cumulative motif sets (`inputs/motifs/`: manufacturing sites,
processing signals, repeats, internal initiation), five nested fronts per instance, no UTRs.
No other tool can express the constraint, so there are no comparator arms. Design and
pre-registered rules: `drivers/e3/E3-PLAN.md`. The grid runs `derna front --motif` at derna2
`06605630` (`drivers/e3/E3-PIN.md`), in 4 h rounds that stop at a gap boundary and resume.
Pass 1 (six short-band proteins × F1–F4) is closed with two complete fronts and 22 anytime ones
(`drivers/e3/E3-GRID-PASS1.md`); pass 2 and the middle band's selection are pending, and the
long band was chosen from the lattice census (`drivers/e3/calib/`). The three bands' fronts
are read in place, while they run, by Figure 3's rules (`grid/`, `mid/`, `long/` and `screen/`
under `data/results/e3/`, no checksum manifest yet). What is here as tables:

- **the forced-motif census**: no entry of the four sets is forced by any peptide, so every
  set is satisfiable on every protein; the rejected thresholds `GGGG` and `CCNCCN` are forced
  by `WG` and `PP`;
- **the hull probe** (plan § 3a): the exact motif hull on P31645 prefixes, a lower bound on the
  front's cost. Peak RSS grows as n^1.95 (≈15 GB at 984 aa, ≈39 GB at 1,629 aa, extrapolated);
  the wall is not a power law, and n = 300 did not finish inside an hour.

| Rule | Output | Content |
|---|---|---|
| `e3_forced_motifs` | `e3/forced_motifs.jsonl` | the census recomputed from the motif sets and the codon table, and checked entry by entry against the recorded one |
| `e3_probe_table` | `e3/probe.tsv`, `e3/probe_fit.txt` | one row per rung: exit, wall, peak RSS, vertices, probes, folds per probe and the cost of one fold (all probes, the λ = 10⁻⁴ probe, the ordinary ones); the peak-RSS fit and the wall's consecutive-rung exponents |
| `e3_tables` | | the two above |
| `e3_fig3_cells` | `e3/fig3_cells.json` | every cell of the three bands' summaries (`grid/`, `mid/`, `long/`): state, wall, peak RSS, points, hull, completeness; rows split by an unquoted core mask are rejoined, and a stopped cell keeps its wall and RSS |
| `e3_fig3_panels` | `e3/fig3_panels.json` | coverage per protein and nest; peak RSS fitted per build where a build spans enough length (the pooled fit mixes builds sitting in different bands); points found relative to F1 |
| `e3_fig3_planned` | `e3/fig3_planned.json` | per short-band cell: ΔMFE at matched CAI against the protein's unconstrained front, whether its sets bind (removed from the price if not, kept in the cost), and wall against that front |
| `e3_fig3_strip` | `e3/fig3_strip.json` | per blacklist entry, proteins admitting it against proteins whose front spells it; one instance's nested fronts with computed hull vertices and the price of each step; every cell's coverage against the retirement gate |
| `e3_figures` | | Figure 3's data, the three above |

### E4: UTR-aware design (running)

Four fronts per protein and UTR pair, all re-scored on the whole molecule by one folder: the
CDS designed alone with the UTRs pasted on afterwards, the UTRs inside the DP with CDS–UTR pairs
forbidden (`--forbid 5-cds,cds-3`; 5′–3′ pairs allowed), mRNAfold with the UTRs in its energy
model, and the structurally optimal front. Instances: the short band's proteins
(`inputs/e4/`, with their source checksums), each with up to five UTR pairs used in practice
(`inputs/utrs/utrs.fasta` and its provenance); the middle and long bands wait on E3's
selections. The stop codon is a recorded dimension of every cell. Design:
`drivers/e4/E4-PLAN.md`; driver `drivers/e4/drive_e4.py`. The grid runs on derna2 builds that
follow `xg` arm by arm, so E4 has no single pin, and its results stay on the host until it
does; it has no rule here.

### E5: one design, start to finish (planned)

A ladder of four cumulative practices on one short instance (SARS-CoV-2 RBD with the `bnt`
pair) and one long one (`PRT-03` with `hbb`): no constraint, the UTRs attached, E3's motif sets,
the start window kept accessible (`--access 4,30`). The pooled designs are summarized with MASS
into structural classes, and practical filters (GC and GC3 windows, uridine, CpG, `+4 = G`) are
applied to the front. Design: `drivers/e5/E5-PLAN.md`. **Nothing has run**; E5 consumes E3's
motif sets and E4's UTR pin and gate, and has no driver or rule.

### Figure 1: GLPE_HUMAN with the hbb UTRs

One real instance (P15421, 78 aa; hbb UTRs, 50 nt 5′ and 132 nt 3′) carrying both
contributions: the exact CDS-only front, what attaching the UTRs does to its designs, a guarded
front with no CDS–UTR pairs, and the exact structurally optimal front under that guard (free
MFE = guarded MFE). `drivers/fig1/FIG1.md` describes the method.

| Rule | Output | Content |
|---|---|---|
| `fig1_so_front` | `fig1/so_front.jsonl` | the structurally optimal front reassembled from its window reports; fails unless the windows settle the whole log-CAI range and the front equals the deposited one |
| `fig1_plot` | `figures/fig1_glpe_hbb.pdf` | panels a–d |
| `fig1_figures` | | the two above |
| `fig1_overview_designs` | `fig1/overview_designs.json` | the overview figure's three designs at CAI ≈ 0.87: mRNAfold's design (CDS–UTR pairs discouraged, from the `fig1_fronts` deposit) with the forbidden XhoI site, a hull vertex of the concatenation front (what a weighted sum returns) and the nearest modular design, refolded with the pinned `fold_seq`; fails unless the refold matches the deposited energy with no CDS–UTR pair and no modular design contains XhoI |
| `fig1_overview` | `figures/fig1_overview.pdf` | the paper's overview figure; needs a LaTeX install (Times) for its text |
| `fig1_front` | `runs/fig1/{arm}/front.jsonl` | the CDS-only (`cdsonly`) or symmetric-guarded (`symmetric`) front, re-walked |
| `fig1_check_front` | `runs/fig1/{arm}/check.txt` | that front against the deposited one, energy at every log-CAI level |
| `fig1_runs` | | both fronts, re-walked and checked |
| `fig1_overview2_prep_hbb` | `fig1/overview2_hbb/` | the revised overview's inputs for the hbb pair from the fig1 and fig1_fronts deposits, hull vertices computed, and its three designs picked as `fig1_overview_designs` does; the structurally optimal design refolded with the pinned `fold_seq` and checked against its recorded energy, no CDS–UTR pair and no XhoI site |
| `fig1_overview2_prep_hba` | `fig1/overview2_hba/` | the same for the hba pair, from `fig1_hba/` |
| `fig1_overview2_plot` | `figures/fig1_overview2_{pair}.pdf` | the revised overview figure: (a) the problem, with the AlphaFold model of the protein (`inputs/structures/`); (b) the fronts with a table of what each front guarantees (UTRs, motifs, complete, modular) and a one-row marker legend (square modular, circle not; filled on the hull, open off it); (c)–(e) the three designs' RNApuzzler structures; layout settings fixed per pair in `config.yaml` (`fig1.overview2`) |
| `fig1_overview2` | | the revised overview figure for both pairs; needs ViennaRNA and a LaTeX install with Times and the `pifont` ZapfDingbats font for the table's marks (on EL7, `texlive-times` and `texlive-zapfding`); the PDF carries no creation date, so a rebuild on the same TeX fonts is byte for byte the same |

### Figure 1's three fronts

The same instance with three fronts drawn together: the exact plain front (114 points, 23 hull
vertices); the exact front under each of three CDS–UTR constraints (`5-cds,cds-3`: 117 points,
none structurally optimal; `5-cds,5-3`: 119 points, 42 structurally optimal; `5-5,5-cds,5-3`:
119 points, none, since the hbb 5′ UTR pairs with itself); and, per constrained hull vertex,
the first structurally optimal design a lazy k-best walk finds (21 of 22 vertices under
`5-cds,cds-3`, 22 of 22 under `5-cds,5-3`). mRNAfold is run at each vertex's weight. The deposit
is derna2's `experiment/glpe-hbb-fronts-data` at `f67a04a6` (merged as derna2 PR #398), built
from `1b3739b`; `drivers/fig1-fronts/README.md` describes the runs. The scripts there carry the
paths of the host they ran on and are not pipeline steps.

| Rule | Output | Content |
|---|---|---|
| `fig1_fronts_tables` | `fig1_fronts/tables.md` | per front: points, hull vertices, structurally optimal points and wall; per front-3 walk: λ, the vertex's energy, rank found, the design's energy and the walk's cost; mRNAfold against each vertex |

## Tests

```bash
.venv/bin/pytest
```

The pipeline tests dry-run every target the tables above name, and run `run.sh` with no
arguments.

## Open items

- **TODO: run only clean `derna` executables.** The pipeline should call only the `derna` CLI,
  built from the integration branch `xg` at a recorded commit, never an experiment binary.
  Figure 1 does not yet: its fronts come from `scope35_walk_front` and `scope37_window`, its
  scoring from `fold_seq` and `guard_fold`, and the overview figure refolds its modular design
  with `fold_seq` because the window walks record no structures. For each of these, either
  `derna` already does the job (then switch the rule to it) or the CLI lacks it, which is a CLI
  gap to fix in derna2, not a reason to keep the custom binary. Then rerun Figure 1 with
  `derna` and drop `pins/fig1.sha256`.
- **Figure 1's binaries have no recorded commit.** `pins/fig1.sha256` pins the binaries
  themselves, so the figure is reproducible from the freeze, but they cannot be rebuilt from
  source until the commit is identified. The test is mechanical: the window reports stamp
  `binary_fnv64` = `1ff256ead61c92ba`, which is `provenance.fnv1a64` of the pinned
  `scope37_window` and of no other pinned binary, so a candidate commit is confirmed by
  building it and hashing. See `pins/fig1-PROVENANCE.md`.
- **Figure 1 panel b's scoring and the window search are not pipeline rules.** The scored files
  (free and guarded folds of each design) and the window rounds were produced by the scripts in
  `drivers/fig1/`, which carry absolute paths to a derna2 checkout. The pipeline checks their
  output (`fig1_so_front`) but does not rerun them.
- **Figure 1's three fronts** have no run target and no frozen binaries here. Their raw outputs
  record the host paths of the checkout, the cargo target directory and the run directory they
  ran in — 97 files, every one of them named by `pins/results-fig1_fronts.sha256` and so every
  one published. Three substitutions cover all of them:

  ```sh
  .venv/bin/python pipelines/scripts/scrub_host_paths.py /scratch/xinyugu4/fig1-fronts-deposit \
      '/home/xinyugu4/derna2/.claude/worktrees/snug-cooking-turtle/=derna2/' \
      '/scratch/xinyugu4/derna2-target-glpe-hbb/=target/' \
      '/scratch/xinyugu4/glpe-hbb-fronts/=run/'
  ```

  The script reports what it would rewrite and exits non-zero on any host path the
  substitutions miss; `--apply` rewrites the files in place. **Rewriting the deposit
  invalidates its manifest**, so the manifest has to be regenerated afterwards and the
  substitutions recorded beside it. Whether to rewrite the deposit or to scrub only a
  published copy is open — but all 191 files are byte-identical to derna2 `f67a04a6`
  (checked file by file; `pins/fig1_fronts-PROVENANCE.md`), so a rewrite loses nothing
  that cannot be recovered from the commit.
- **E2 is in progress.** Its deposit has no checksum manifest until the grid finishes.
- **E2's freeze.** When the second host chooses its commit, set `e2.pin`, fill
  `pins/e2.sha256` with every binary in `e2.bin_dir` (derna2's and the comparators') and add
  `pins/e2-PROVENANCE.md`. Rows of the first host's pilot (`acfaae60`) are not part of the grid.
- **E3's grid** has no rule: its binary (`06605630`) is not frozen here and its fronts are not
  deposited. The probe's binary (`motif_hull` at `5d12df5f`) is not frozen either, so the probe
  is a record, not a run target. The union rows of the forced-motif census came from code the
  plan does not include, and are not checked.
- **E4** has a driver, inputs and UTRs here but no pin or rule; **E5** needs its driver and a
  pin, and also pins MASS.
- **The comparator builds** (LinearDesign, LinearCDSfold, CDSfold, mRNAfold) are described in
  `pins/e1-PROVENANCE.md`, not automated.
- **A rebuild of derna2 at a pin has not been compared** with the frozen binaries
  (`check_derna2_build`).
- **License**: none chosen yet.
