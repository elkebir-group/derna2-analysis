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
- **`inputs/`**: the protein sequences, UTRs, codon table and E2's motif sets, committed.

## Install

Python 3.11–3.13. On a host whose glibc is older than 2.28 (CentOS 7), the numpy cap in
`pyproject.toml` keeps pip on wheels it can install. ViennaRNA has no wheel for such a host and
its source build fails there. Two rules fold a sequence and so import it — `fig1_overview` and
`e2_composite_example` — and on such a host they need a Python whose ViennaRNA 2.7.2 bindings
come from elsewhere (e.g. conda's `viennarna`), with the rest of `requirements.lock` installed
beside it. Name that interpreter as `rna_python` in `host.yaml`; the other rules keep running in
the environment snakemake itself runs in, and nothing else in the pipeline needs it.

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
├── e1/                      # in progress: no checksum manifest yet
│   ├── sweep.jsonl          # one row per λ, per front point and per cell, stamped with its pin
│   └── hull/                # derna2's hull vertices and replay λs per instance (partial if the
│                            # walk was stopped)
├── e1_fig4/                 # Figure 4's frozen snapshot of e1 (pins/e1_fig4-PROVENANCE.md)
│   ├── sweep.jsonl          # the sweep after the grid's last pass, 2026-10-31: 258,083 rows
│   └── hull/                # <P>.front.jsonl and <P>.walk.json for the 85 proteins
├── e2/                      # cut 2026-11-05: pins/results-e2.sha256, pins/e2-PROVENANCE.md
│   ├── e3-probe/            # hull probe round 1: per rung the walk's rows, its measurement
│   │                        # record and summary.csv; forced.jsonl, the forced-motif census
│   ├── e3-probe2/           # round 2: the n = 200 load control, n = 150, 250 and 300
│   ├── grid/, mid/, long/   # the short, middle and long bands: summary.csv and each cell's
│   │                        # front by round (`<protein>_<nest>*.jsonl`)
│   └── screen/              # binding.jsonl, the binding screen; f0/, unconstrained fronts
├── fig1/
│   ├── a_cdsonly_front.jsonl                    # exact CDS-only front, hull vertices marked
│   ├── b_cdsonly_with_utrs.scored.jsonl         # those designs folded with the UTRs
│   ├── cd_symmetric_guarded_front.scored.jsonl  # exact constrained front
│   ├── d_so_front.jsonl                         # exact modular front
│   ├── so_windows/          # the window walks that found it, by round
│   └── top_exhaustive/      # the stretch above the most adapted seed, by enumeration
├── fig1_hba/                # not deposited yet: the overview figure for the hba UTR pair
│   ├── run3-hba/            # the GLPE arms' results for that pair (a0 CDS-only front and its
│   │                        # refolds, a2 masking front, a3 mRNAfold)
│   └── xhoi/so.jsonl        # the modular front with XhoI forbidden
└── fig1_fronts/
    ├── <arm>/               # unguarded, symmetric, no_clamp, unpaired: report.json, the front
    │                        # and its hull vertices re-scored; front3.jsonl and its designs
    ├── mrnafold/compare.jsonl  # mRNAfold at each hull vertex's weight
    ├── plot/                # plot_data.json and fronts.png
    └── raw/                 # every run output untrimmed, and the binaries' commit and SHA-256
```

Nothing under `data/` is committed.

## The figures

| # | from | what the figure asserts | target |
|---|---|---|---|
| 1 | Figure 1's instance (GLPE_HUMAN + hbb) | attaching UTRs breaks a design made without them, and the modular front is computable | `fig1_figures`, `fig1_overview2` |
| 2 | E1 | a λ sweep is not a route to the front: it reaches hull vertices only, and not many of those at length | `e1_figures` |
| 3 | E2 | what forbidding motifs costs, against a baseline shown to bind | `e2_figures`, `e2_fig_composite` |
| 4 | E3 | what a UTR-aware practice delivers, per practical UTR pair | `e3_figures` |
| S1–S4 | the overflow panels of 2–4 | memory, host qualification, the motif supplement, the protein tables | `e1_memory_figures`, `e2_fig_supp_motif`, `e3_fig_supp_contacts`, `supp_proteins` |

Three conventions hold across them and are not re-argued per figure.

- **Every front panel has the same axes:** x = CAI on 0–1, y = MFE in kcal/mol, the optimum
  to the lower right, hull vertices as filled marks and off-hull points as small open marks.
  Where convexity is the argument being made, a log-CAI inset is added rather than switching
  the main axes, since the hull is convex in (`cai_log`, `mfe`) and not in (`cai`, `mfe`).
- **Every plotted row carries its build and its full configuration.** Two rows that differ in
  a flank length, a stop codon, a default or a commit are not comparable, and a cross-build
  pair is visible only where the stamp is per row.
- **The figures are built here**, as Snakemake rules under `pipelines/rules/`, from the shared
  style and method registry in `derna2_analysis/figures.py`.

## Running the pipeline

`pipelines/run.sh` runs `pipelines/Snakefile` with `pipelines/config/config.yaml`, then
`pipelines/config/host.yaml` if it exists (copy `host.example.yaml`). A target is a rule name or
an output path relative to the repository root.

```bash
pipelines/run.sh --dry-run                                     # what would run
pipelines/run.sh                                               # every figure, from the deposit
pipelines/run.sh -- fig1_figures                               # one experiment's figures
pipelines/run.sh -- data/derived/figures/e1_fronts.pdf         # one output
pipelines/run.sh -- e1_runs                                    # regenerate E1 from binaries
```

`run.sh` passes `--rerun-triggers mtime`, so a code or config change never reruns an existing
output by itself: delete it, or add `--forcerun <rule>`. Outputs go to `data/derived/`.

There are two kinds of target. **Figure targets** read results and need no binaries; `all` is
these. Set `figures_from: runs` to draw from this checkout's own runs instead of the deposit.
**Run targets** regenerate the results. They take hours to days, and the timed passes need a
reserved socket, so they run only when named. Every one of them runs the `derna` command line
and `mfe_fold`, built by `build_derna2` at the commit `derna2.pin` names, or taken from a
directory given as `derna2.bin_dir`; the comparators come from `bin_dir` of the experiment
that runs them (`pins/comparators-PROVENANCE.md`).

### Provenance

| Rule | Output | Content |
|---|---|---|
| `check_results` | `stamps/results_{exp}.ok` | every deposited file of an experiment against `pins/results-{exp}.sha256` |
| `build_derna2` | `bin/derna2-{pin}/` | `derna` and `mfe_fold` built from a fresh clone at `derna2.pin`, with their checksums and the commit beside them |
| `package_deposit` | `deposit/{exp}.tar` | the deposit of one experiment: exactly the files its manifest names, checked, without the run logs and working scripts beside them |

### E1: the whole front

What the whole Pareto front costs each method. derna2 enumerates it (`derna2-walk`, the walk
`derna front` ships: `--pass-grid auto --seed-sweep --cut-probes --lend-idle`); every other
method sweeps its scalarization weight at derna2's hull λs (one per hull vertex plus the
corners), and DERNA v1, through its Rust port, also by its own adaptive sweep. derna2's walk
and mRNAfold run once on one core, against every other single-core arm, and once on eight
(`-par`). Instances: the 55 UniProt proteins up to 500 aa, twelve native proteins of
511–962 aa (`inputs/e1-natives/`), P0DTC2 whole (1,273 aa) and five Cas orthologs (PRT-03 to
PRT-07, 1,000–1,629 aa). Limits: 2 h per cell below 1,000 aa (soft for the walk, which runs on
to 12 h), 6 h above (12 h for `derna2-hull` and the walk); a cell is admitted only while its
estimated peak memory fits the host's budget. A walk the limit stops keeps its partial front,
and the replay arms then run at its partial hull. The grid runs on a second host and
**`e1.pin` reads `pending`**: the commit its rows will be stamped with is chosen when it
runs. The comparators it is measured against are in `pins/comparators-PROVENANCE.md`.

| Rule | Output | Content |
|---|---|---|
| `e1_plot_fronts` | `figures/e1_fronts.pdf`, `e1/cells.tsv` | wall time and Pareto points recovered per cell against length, eight-core arms hollow and dashed; rows of the configured pin only |
| `e1_fig2_data` | `e1/fig2.json` | Figure 2's data: per arm, wall, peak RSS, designs and seconds per design on the cells every arm completed, with completions, timeouts and the longest protein finished; per walk front, the share of its points no weight reaches, from hull vertices computed on the points, and the fronts whose stored `vertices` disagrees; rows of the configured pin only |
| `e1_figures` | | the two above |
| `e1_config` | `runs/e1/config.json` | `drivers/e1/e1.config.json` with this host's paths and the pin |
| `e1_verify` | `runs/e1/verify.done` | each comparator's λ map against derna2 at interior λs |
| `e1_grid` | `runs/e1/grid.done` | `run_e1_grid.sh`: refuses to start while a timing core or its hyperthread sibling is busy, then the four passes in order (single-core walk, single-core replays, eight-core walk, eight-core replays) into `runs/e1/sweep.jsonl`; resumable |
| `e1_summarize` | `runs/e1/summary.txt` | fronts, recovery against the reference, anytime curves, for the configured pin |

#### The paper's Figure 4

Drawn from its own frozen deposit, `data/results/e1_fig4`: a snapshot of the E1 sweep and the
walk's `hull/` outputs, pinned by `pins/results-e1_fig4.sha256` and described in
`pins/e1_fig4-PROVENANCE.md`. Every table is rebuilt from those files, and the tables and the
figure are checked against the digests in `pins/e1_fig4.outputs`. The protein of panel (a) is
`e1_fig4.focus` in `config.yaml` (P47211); any protein of the corpus can be drawn there.

| Rule | Output | Content |
|---|---|---|
| `e1_fig4_cells` | `e1/fig4/cells.csv` | one cell row per (arm, protein) of the replay pin `8603d37d`: the delivered and then faster attempt |
| `e1_fig4_per_protein` | `e1/fig4/per_protein.csv` | per protein, each arm's exit and wall time, and the dichotomic hull search with Carrillo-Lipman pruning (measured) and without (projected) |
| `e1_fig4_steps` | `e1/fig4/steps.npz` | every front Pareto-filtered on raw (energy, log CAI), its sweep hull vertices, and the ΔMFE and ΔCAI between neighbouring designs of each |
| `e1_fig4` | `figures/fig4_e1_overview.pdf` | the figure; fails if any text covers another text, a mark or another panel |
| `e1_fig4_verify` | `stamps/e1_fig4.verified` | the three tables and the figure against `pins/e1_fig4.outputs` (the figure's bytes also need the TeX that file names) |
| `e1_fig4_figures` | | the two above |
| `e1_runs` | | the run chain above |

### E2: the price of motif-freedom

What excluding motifs costs: the exact Pareto front over every design, then over the designs
that avoid each of four cumulative motif sets (`inputs/motifs/`: manufacturing sites,
processing signals, repeats, internal initiation), five nested fronts per instance, no UTRs.
No other tool can express the constraint, so there are no comparator arms. The grid's rules
were pre-registered before it ran, in a plan document that is not published here. It runs `derna front --motif` at derna2
`06605630` (`drivers/e2/E2-PIN.md`), in 4 h rounds that stop at a gap boundary and resume.
Pass 1 (six short-band proteins × F1–F4) is closed with two complete fronts and 22 anytime ones
(`drivers/e2/E2-GRID-PASS1.md`); pass 2 and the middle band's selection are pending, and the
long band was chosen from the lattice census (`drivers/e2/calib/`). The three bands' fronts
are read by Figure 3's rules from the cut deposit (`grid/`, `mid/`, `long/` and `screen/`
under `data/results/e2/`, checksummed by `pins/results-e2.sha256`), not in place while they
run, since the cut of 2026-11-05. What is here as tables:

- **the forced-motif census**: no entry of the four sets is forced by any peptide, so every
  set is satisfiable on every protein; the rejected thresholds `GGGG` and `CCNCCN` are forced
  by `WG` and `PP`;
- **the hull probe** (plan § 3a): the exact motif hull on P31645 prefixes, a lower bound on the
  front's cost. Peak RSS grows as n^1.95 (≈15 GB at 984 aa, ≈39 GB at 1,629 aa, extrapolated);
  the wall is not a power law, and n = 300 did not finish inside an hour.

| Rule | Output | Content |
|---|---|---|
| `e2_forced_motifs` | `e2/forced_motifs.jsonl` | the census recomputed from the motif sets and the codon table, and checked entry by entry against the recorded one |
| `e2_probe_table` | `e2/probe.tsv`, `e2/probe_fit.txt` | one row per rung: exit, wall, peak RSS, vertices, probes, folds per probe and the cost of one fold (all probes, the λ = 10⁻⁴ probe, the ordinary ones); the peak-RSS fit and the wall's consecutive-rung exponents |
| `e2_tables` | | the two above |
| `e2_fig3_cells` | `e2/fig3_cells.json` | every cell of the three bands' summaries (`grid/`, `mid/`, `long/`): state, wall, peak RSS, points, hull, completeness; rows split by an unquoted core mask are rejoined, and a stopped cell keeps its wall and RSS |
| `e2_fig3_panels` | `e2/fig3_panels.json` | coverage per protein and nest; peak RSS fitted per build where a build spans enough length (the pooled fit mixes builds sitting in different bands); points found relative to F1 |
| `e2_fig3_planned` | `e2/fig3_planned.json` | per short-band cell: ΔMFE at matched CAI against the protein's unconstrained front, whether its sets bind (removed from the price if not, kept in the cost), and wall against that front |
| `e2_fig3_strip` | `e2/fig3_strip.json` | per blacklist entry, proteins admitting it against proteins whose front spells it; one instance's nested fronts with computed hull vertices and the price of each step; every cell's coverage against the retirement gate |
| `e2_figures` | | Figure 3's data, the three above |

#### Figure 3, the composite plate (`e2_fig_composite`)

Seven panels at the paper's text width, drawn by `pipelines/scripts/e2/fig_composite.py` from a
cache that four rules write. The geometry is solved in inches at the size `e2.figure` names, so
the width and height are part of the figure's definition rather than a rendering convenience —
a different size is a different layout and has to be re-measured.

| Rule | Output | Content |
|---|---|---|
| `e2_composite_auc` | `e2/composite/auc01_fin.json` | the finished set and each cell's AUC per rung, which is the population every later step is computed over |
| `e2_composite_data` | `e2/composite/composite-data.json`, `composite-census.json` | panel (a)'s screen catalog and panel (b)'s per-rung cost at four CAI floors, per protein, with the provenance of every front it read; the census file is the host's own state at build time (below) |
| `e2_composite_dcai` | `e2/composite/dcai-at-mfe.json` | panel (b)'s other direction: the CAI given up at four floors on the MFE axis |
| `e2_composite_example` | `e2/composite/example-<cell>-c<floor>-<alt>.json` | panels (c)–(g): one protein's four nested fronts, and the four designs at the floor folded and laid out on one shared scale. **Imports ViennaRNA** — see `rna_python` under Install |
| `e2_fig_composite` | `figures/fig3_e2_composite.pdf`, `.png` | the plate |
| `e2_fig_composite_verify` | `stamps/e2_composite.verified` | the five published files against `pins/e2_composite.outputs` |

```bash
pipelines/run.sh -- e2_fig_composite          # the plate, from a cold cache: about 35 s
pipelines/run.sh -- e2_fig_composite_verify   # and the same bytes as the published figure
```

**The plate fails rather than ships a collision.** Before saving, the script rasterizes itself
at 200 dpi and tests every text against every other text, against every panel it is not part of,
and against the canvas edge, asking 4 device pixels of clearance; one failure aborts the save,
so a layout regression cannot reach the paper as ink. The caption widths are checked the same
way and raise, so `python -O` cannot drop the check.

**The bytes are stable, and the manifest passes while the bands are still running.** The strip's
jitter is seeded, the PDF's `CreationDate` and the PNG's `Software` stamp are suppressed, and the
cache carries no build time, no repository commit and no checksum of the finished table — each of
those moved between two builds whose every plotted number was identical, and a self-referential
commit field can never agree with a rebuild in any case. What identifies the code is the digest
of every script the cache was computed by, which `composite-data.json` records. The one thing
that does move with the host rather than with the figure — the finished table's write time and
checksum, and which cells have a round file on disk and no row in that table — is written beside
the cache as `composite-census.json`, which nothing reads and `pins/e2_composite.outputs` does
not name.

The chain is **not closed on an arbitrary host**: panel (a) reads E1's sweep and an auxiliary
tree of screened E1 instances, neither of which is in E2's deposit, and the short2/mid bands'
fronts are read at their host paths. `e2.e1_sweep` and `e2.e1_extra` name the first two in
`host.yaml` until they are deposited.

#### The motif supplement, and what a constraint can do to a protein (`e2_fig_supp_motif`)

Nine panels in two rows at the text width, drawn by `pipelines/scripts/e2/fig_supp_motif.py`
from a cache three rules write, with the motif catalog published beside it as a table. A
protein fixes its residues and a design chooses only among synonymous codons, so each pattern
stands in one of three relations to a given protein — **inadmissible**, **admissible**,
**unavoidable** — and which one holds is a property of the protein. (a)–(c) draw one residue
window as two columns of the codon lattice, one edge per synonymous pair, marked where the
pair spells the pattern; (d) is the same question over every window and every synonymous
choice of every protein; (e), (f) are how much of each set survives on the front as levels
are added; (g)–(i) are what reaching each front costs against the same protein's F0 run. The
geometry is solved in inches at the size `e2.supp_figure` names, so the width and height are
part of the figure's definition rather than a rendering convenience — a different size is a
different layout and has to be re-measured.

| Rule | Output | Content |
|---|---|---|
| `e2_supp_admits` | `e2/supp/admits.json` | per protein and pattern, whether any synonymous coding sequence spells it, and the first residue window that does; exact, over every window and every synonymous choice |
| `e2_supp_presence` | `e2/supp/presence.json` | per protein, level and set: the share of that level's designs spelling the set, the mean and median sites per design, and the sites in the most stable design at CAI ≥ 0.80 |
| `e2_supp_cost` | `e2/supp/cost.json` | per protein and level: wall summed over the level's rounds and peak RSS maxed over them, each carrying the chain below it where the run narrowed from the level below, and the front's size — all three as a ratio to the same protein's F0 run |
| `e2_fig_supp_motif` | `figures/fig_supp_motif.pdf`, `.png`, `.prov.json` | the plate, and a sidecar carrying the median and the n behind every box it draws |
| `e2_supp_table` | `figures/tab_supp_motifs.tex` | the 24 patterns by set: name, sequence, and the proteins each is inadmissible for. The plate labels (d) by name, because a name-over-sequence tick does not fit at that width, so the sequences and the counts are published here at full size |
| `e2_fig_supp_motif_verify` | `stamps/e2_supp_motif.verified` | the six published files against `pins/e2_supp_motif.outputs` |

```bash
pipelines/run.sh -- e2_fig_supp_motif e2_supp_table   # from a cold cache: about 32 s
pipelines/run.sh -- e2_fig_supp_motif_verify          # and the same bytes as the published figure
```

**The plate fails rather than ships a collision**, as Figure 3 does. Before saving, the gate
in `pipelines/scripts/e2/supp_gate.py` draws the figure and tests every label against every
other label, against every panel it is not part of, and against the canvas edge, asking 2.9 pt
of clearance; one failure aborts the save, and the check is an `assert` so `python -O` cannot
drop it. The clearance is given in points and converted at the figure's own dpi, because a
threshold written in device pixels means a different distance on the page at every dpi. Two
labels are compared even when they share a panel — nothing but the drawing code places those —
unless they declare themselves one block by sharing a `gid`, which is how the two lines under
each lattice say their leading is deliberate. The gate is a module rather than inline, so
`tests/test_supp_motif.py` can drive it: it is pinned to fire on touching text, to pass text
set apart, to catch a letter on its own panel's title and a label over a neighbour, and not to
change its answer when the dpi doubles.

**The bytes are stable, and the pin holds while the grid is still running.** Every jitter is
seeded per panel, the PDF's `CreationDate` and the PNG's `Software` stamp are suppressed, and
the caches carry no build time and no repository commit. The layout is measured off the drawn
figure, so the script fixes `figure.dpi` and refuses a tight save box: both decide the measured
geometry and neither is the host's to choose. The finished table is digested parsed rather than
as a file, over the complete and retired rows only and by each round file's basename — its
header carries a write time, it is regenerated hourly, an open cell's round advances every
round the grid runs, and the round file is spelled as an absolute path, so a digest over the
file, over the open rows, or over the paths moves while every plotted number stands still. All
three were seen to move exactly that way. `admits.json` carries no digest of the finished set
at all, its population being the 103 proteins with an F0 front rather than the 47 closed cells.
Each cache names its members, since a count cannot report a one-in-one-out swap.

The chain is **not closed on an arbitrary host**: `presence.json` and `cost.json` read each
cell's fronts and round logs at the absolute paths `FINISHED.tsv` records, the same reach as
the composite's short2 and mid bands. `admits.json` and the table read the deposit's F0 fronts
and the committed motif sets, and are closed wherever the deposit is.

### E3: UTR-aware design (running)

Four fronts per protein and UTR pair, all re-scored on the whole molecule by one folder: the
CDS designed alone with the UTRs pasted on afterwards, the UTRs inside the DP with CDS–UTR pairs
forbidden (`--forbid 5-cds,cds-3`; 5′–3′ pairs allowed), mRNAfold with the UTRs in its energy
model, and the modular front. Instances: the short band's proteins
(`inputs/e3/`, with their source checksums), each with up to five UTR pairs used in practice
(`inputs/utrs/utrs.fasta` and its provenance); the middle and long bands wait on E2's
selections. The stop codon is a recorded dimension of every cell. The driver is
`drivers/e3/drive_e3.py`, and its design was pre-registered in a plan document that is not
published here. The grid ran against a development build per arm rather than one commit, so
it has no pin and no rule here; what is published is the figure bundle below.

The figure bundle copied from the host (`paper-e4/`, 664 files, checked against
`pins/results-e3.sha256`) is read by two figures, the main one and a supplement one, which are rules:

| Rule | Output | Content |
|---|---|---|
| `e3_hull_input` | `e3/hullcheck/<cell>.in.jsonl` | one cell's concatenation hull vertices as a derna result file |
| `e3_refold_deposit` | `e3/hullcheck/<cell>.out.jsonl` | the deposited refold of one cell, which is what the figures read; needs no binary |
| `e3_refold` | the same file | `derna check` of those vertices, which carries the positions of each design's CDS–UTR pairs; runs the `derna` built at `derna2.pin`. The deposit is preferred, so recomputing one is `--force-run e3_refold` |
| `e3_cand_input` | `e3/cand/<cell>.{conc,mrf}.in.jsonl` | the main figure's instance (`e3.story`): the best concatenation design and the best mRNAfold design at the CAI floor, as derna result files |
| `e3_cand_refold` | `e3/cand/<cell>.{conc,mrf}.out.jsonl` | `derna check` of those, the whole-mRNA structure and CDS–UTR pairs panel (e) draws |
| `e3_fig_story` | `figures/fig_e3_story.pdf` | the main figure: UTR elements violated, modular designs per method, DERNA2's windows by state, the instance's fronts, and its three designs |
| `e3_fig_story_verify` | `stamps/e3_story.verified` | the figure against `pins/e3_story.outputs` (its bytes depend on matplotlib and TeX) |
| `e3_story_stats` | `e3/story_stats.txt` | the numbers the paper's UTRs subsection quotes, one `key: value` line each: modular shares and counts per method and UTR pair, DERNA2's windows by state, and the main figure's instance, with the comparison over all instances at the CAI floor |
| `e3_story_stats_verify` | `stamps/e3_story_stats.verified` | that file against `pins/e3_story_stats.outputs` |
| `e3_contacts_data` | `e3/supp_contacts.json` | per UTR pair, the share of the non-modular concatenation designs with a CDS–UTR contact endpoint in each bin of the 5′ UTR, the CDS and the 3′ UTR |
| `e3_fig_supp_contacts` | `figures/fig_supp_contacts.pdf` | the supplement figure: one row per UTR pair, drawn through usetex (a TeX install) |
| `e3_supp_contacts_verify` | `stamps/e3_supp_contacts.verified` | the histograms against `pins/e3_supp_contacts.outputs` |

`e3_figures` is both figures and the statistics file. Link the bundle's `paper-e4/` at `data/results/e3` and give
E2's results at `data/results/e2`, since the concatenation arm points at E2's fronts. Both are
checked: `pins/results-e3.sha256` and, since the cut of 2026-11-05, `pins/results-e2.sha256`.

### Figure 1: GLPE_HUMAN with the hbb UTRs

One real instance (P15421, 78 aa; hbb UTRs, 50 nt 5′ and 132 nt 3′) carrying both
contributions: the exact CDS-only front, what attaching the UTRs does to its designs, a
constrained front with no CDS–UTR pairs, and the exact modular front under that constraint
(free MFE = constrained MFE). `drivers/fig1/FIG1.md` describes the method.

| Rule | Output | Content |
|---|---|---|
| `fig1_so_front` | `fig1/so_front.jsonl` | the modular front reassembled from its window reports; fails unless the windows settle the whole log-CAI range and the front equals the deposited one |
| `fig1_plot` | `figures/fig1_glpe_hbb.pdf` | panels a–d |
| `fig1_figures` | | the two above |
| `fig1_overview_designs` | `fig1/overview_designs.json` | the overview figure's three designs at CAI ≈ 0.87: mRNAfold's design (CDS–UTR pairs discouraged, from the `fig1_fronts` deposit) with the forbidden XhoI site, a hull vertex of the concatenation front (what a weighted sum returns) and the nearest modular design, refolded with the pinned `fold_seq`; fails unless the refold matches the deposited energy with no CDS–UTR pair and no modular design contains XhoI |
| `fig1_overview` | `figures/fig1_overview.pdf` | the paper's overview figure; needs a LaTeX install (Times) for its text |
| `fig1_front` | `runs/fig1/{arm}/front.jsonl` | the CDS-only (`cdsonly`) or symmetric-constrained (`symmetric`) front, re-walked |
| `fig1_check_front` | `runs/fig1/{arm}/check.txt` | that front against the deposited one, energy at every log-CAI level |
| `fig1_runs` | | both fronts, re-walked and checked |
| `fig1_overview2_prep_hbb` | `fig1/overview2_hbb/` | the revised overview's inputs for the hbb pair from the fig1 and fig1_fronts deposits, hull vertices computed, and its three designs picked as `fig1_overview_designs` does; the modular design refolded with the pinned `fold_seq` and checked against its recorded energy, no CDS–UTR pair and no XhoI site |
| `fig1_overview2_prep_hba` | `fig1/overview2_hba/` | the same for the hba pair, from `fig1_hba/` |
| `fig1_overview2_plot` | `figures/fig1_overview2_{pair}.pdf` | the revised overview figure: (a) the problem, with the AlphaFold model of the protein (`inputs/structures/`); (b) the fronts with a table of what each front guarantees (UTRs, motifs, complete, modular) and a one-row marker legend (square modular, circle not; filled on the hull, open off it); (c)–(e) the three designs' RNApuzzler structures; layout settings fixed per pair in `config.yaml` (`fig1.overview2`) |
| `fig1_overview2` | | the revised overview figure for both pairs; needs ViennaRNA and a LaTeX install with Times and the `pifont` ZapfDingbats font for the table's marks (on EL7, `texlive-times` and `texlive-zapfding`); the PDF carries no creation date, so a rebuild on the same TeX fonts is byte for byte the same |

### Figure 1's three fronts

The same instance with three fronts drawn together: the exact plain front (114 points, 23 hull
vertices); the exact front under each of three CDS–UTR constraints (`5-cds,cds-3`: 117 points,
none modular; `5-cds,5-3`: 119 points, 42 modular; `5-5,5-cds,5-3`:
119 points, none, since the hbb 5′ UTR pairs with itself); and, per constrained hull vertex,
the first modular design a lazy k-best walk finds (21 of 22 vertices under
`5-cds,cds-3`, 22 of 22 under `5-cds,5-3`). mRNAfold is run at each vertex's weight. The deposit
is derna2's `experiment/glpe-hbb-fronts-data` at `f67a04a6` (merged as derna2 PR #398), built
from `1b3739b`; `drivers/fig1-fronts/README.md` describes the runs. The scripts there carry the
paths of the host they ran on and are not pipeline steps.

| Rule | Output | Content |
|---|---|---|
| `fig1_fronts_tables` | `fig1_fronts/tables.md` | per front: points, hull vertices, modular points and wall; per front-3 walk: λ, the vertex's energy, rank found, the design's energy and the walk's cost; mRNAfold against each vertex |

### The proteins

One supplement table spans the experiments: every protein the paper designs, with the
experiments run on each. It reads the corpus of the Pareto-front comparison, E2's per-protein
AUC table for the proteins finished at every motif level, E3's cells, and the FASTA files under
`inputs/`. A protein run under a label that names another UniProt entry is renamed here, at the
one place a reader sees an accession; the label still joins a result to its protein, because the
frozen deposits carry it in their file names.

| Rule | Output | Content |
|---|---|---|
| `supp_proteins` | `figures/tab_supp_proteins.tex` | one row per protein: accession, name, residues, and a mark per experiment (front, motif levels, UTRs), ordered by length |
| `supp_proteins_verify` | `stamps/supp_proteins.verified` | that table against `pins/supp_proteins.outputs` |

## Tests

```bash
.venv/bin/pytest
```

The pipeline tests dry-run every target the tables above name, and run `run.sh` with no
arguments.

## Open items

Each of these can be closed, and says what closes it.

- **The whole-front experiment (E1) is in progress**, and its deposit has no checksum
  manifest until the grid finishes. `e1.pin` is now `8603d37d`, the pin every row of the
  grid's sweep carries; it selects rows and drives no build, which is why it may name a
  commit the public history does not carry. What remains open is the grid itself: on the
  sweep of 2026-10-31, the last pass to have run, the cells not finished per arm are DERNA2
  2, mRNAfold 19, LinearDesign 21, LinearCDSfold 50, DERNA1 57, of 85 proteins. The deposit
  freezes when those close.

  Until 2026-11-05 `data/results/e1` pointed at the first host's pilot instead: 437 rows and
  9 instances, all of them `P0DTC2_pNN` ladder rungs stamped `acfaae60`, last written
  2026-10-16, while the grid's 258,083 rows over 85 proteins sat unread beside it. Those
  rungs are prefixes of one protein and not instances, so they are not part of the grid, and
  `e1_figures` is in the `all` target — a default run would have drawn E1 from them. The link
  now names the grid's output. The three figures had never been built, because `e1.pin` read
  `pending` and no row carries that, and the overview cannot be drawn from the pilot at all
  (no finished walk of P47211 under `acfaae60`).
- **The overview figure reads a wall time the result format does not carry.**
  `e1_plot_overview` takes the pruned hull search's wall (`sweep_s`) from
  `hull/<protein>.walk.json`, which the sweep rows do not carry, and the walk's rescored front
  points where `sweep.jsonl` has them (30 of the 83 walks are not rescored and use the walk's
  own energies). `derna front --report` carries the phase times, so this is a matter of reading
  them from there rather than from a file the dev binary wrote — and the rerun that produces
  those reports is the grid above, so this item closes with that one.
- **The remaining deposits have not been scrubbed of this host's paths.** Figure 1's three
  fronts were, on 2026-11-05: four substitutions over 97 of its 191 files, all of them under
  `raw/`, with nothing a rule reads touched and the manifest regenerated
  (`pins/fig1_fronts-PROVENANCE.md`). E1's is partly affected and not frozen yet — one of the
  seven files under `hull/` names the host, and its 216 KB sweep names it nowhere. E2's cut
  needed no scrub: the 205 root-level files that named the host and the development tooling
  are run scripts, backups and session notes that no reader can reach, and the cut leaves
  them behind rather than rewriting them (`pins/e2-PROVENANCE.md`); no file in the deposit
  names either. `scrub_host_paths.py` reports before it writes and
  exits non-zero on any path its substitutions miss, so this is one command and one reading per
  deposit, at the point each is frozen.
- **E2's deposit is cut and packed, and waits only on publication.** It was cut on 2026-11-05:
  `data/results/e2` is now 2,514 files and 2.099 GiB, against the 132 GiB run tree it was
  cut from, and `pins/results-e2.sha256` names every one of them. All eleven rules that
  reach the tree now declare `source_ok("e2")`, so the digests are a precondition of every
  E2 figure rather than a record kept beside it; before the cut there was no
  `pins/results-e2.sha256` for `check_results` to read, and nothing stopped a figure being
  drawn from a tree that had moved underneath it. `pins/e2-PROVENANCE.md` records how the
  cut was chosen and what it was checked against: the thirteen outputs the four E2 and E3
  pins name were deleted and rebuilt from the cut, and all thirteen came back
  byte-identical. It packs and reconstitutes: `package_deposit` writes
  `deposit/e2.tar` from the manifest alone, verifying every file against its digest on the
  way in, and the unpacked tar verifies against the manifest again — 2,514 files both times.
  What remains is publishing it, which is where every deposit here stands: no experiment's
  tar is committed, since `package_deposit` builds one on demand from the manifest. So this
  is not a debt particular to E2, and it closes for all of them together, when the deposits
  are published.
- **Two of Figure 3's inputs come from outside any deposit**, `e2.e1_sweep` and `e2.e1_extra`,
  and are set in `host.yaml`. The sweep is E1's own result, read for one thing: the set of
  proteins the screen is scored over, which `composite_data.py` then asserts equal to the
  committed `inputs/e2/corpus.json`. That assertion is why the sweep cannot be replaced by a
  list derived from it — two derivations that agree until they do not is the point of it, and
  a list compared against itself agrees always. It closes when E1's deposit lands, since the
  path currently names a snapshot of that sweep rather than the sweep. `e2.e1_extra` is
  separate and small: two screen files read whole, and eight front files read for their first
  line alone.

  `e2.e1_sweep` names the grid's final sweep (2026-10-31) since 2026-11-05; it had named the
  10-27 snapshot. The sweep is read for one thing only, the set of its rows' `protein`
  values, so the change is a no-op and was checked as one: all eleven outputs of
  `pins/e2_composite.outputs` and `pins/e2_supp_motif.outputs` rebuild byte-identical from
  it. One consequence is left deliberately. `draft_panel_a.py` unions `P22301` and `P41159`
  into E1's instance list by name, because they post-dated the 10-27 snapshot; both are in
  the 10-31 sweep, so the union is now inert. It stays because `composite_data.py` records
  that script's own sha256 in its provenance, so deleting a line that changes nothing would
  still change `composite-data.json` and force a pin re-cut. It goes when a re-cut is
  wanted for another reason.
- **License**: none chosen yet.

## Known limits

These are properties of how the results were made, not work that is pending. Each names a
step that a published checkout cannot repeat, and what stands in for repeating it.

- **Part of Figure 1 is checked rather than reproduced.** Its two fronts rerun through
  `derna front` and check against the deposit (`fig1_runs`), and the overview figures refold
  their modular design with `mfe_fold`. Panel (b)'s scoring and the window search do not
  rerun: those came from the scripts in `drivers/fig1/`, and the window rounds from a binary
  that is not part of the command line. The pipeline checks their output (`fig1_so_front`)
  without reproducing it.
- **E2's grid has no rule**, because `derna front --motif` is not in a pinned build yet — not
  because its fronts are unavailable: they are deposited and checksummed
  (`pins/results-e2.sha256`). The hull probe ran on a binary that predates the command line,
  so what is under `drivers/e2/` for it is a record of the run rather than a route to repeat
  it. The union rows of the forced-motif census came
  from code the plan does not include, and are not checked.
- **E3's plate is not a run target.** Its figures, its statistics and the refolds of
  concatenation's hull vertices all rerun — `derna check` at `derna2.pin` rebuilds
  `e3_story_stats` and `e3_supp_contacts` byte-identical to their pins — but the census those
  read is the deposit's, produced by `drivers/e3/drive_e3.py`. The build the deposited refolds
  were made with is not a commit of the public repository and is identified by its bytes in
  `pins/e3-PROVENANCE.md`.
- **The comparator builds** (LinearDesign, LinearCDSfold, CDSfold, mRNAfold) are frozen once
  and identified by their checksums rather than rebuilt; `pins/comparators-PROVENANCE.md`
  records how each was made and why one of the two mRNAfold builds is the parallel one.
