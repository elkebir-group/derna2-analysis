# derna2-analysis

This repository holds the code that reproduces the analyses and figures of the
[derna2](https://github.com/elkebir-group/derna2) paper. The paper covers exact Pareto-optimal
mRNA codon design, its cost against other design tools, and constrained design with UTRs.
There are five parts.

- **`derna2_analysis/`** is the shared library. It reads result rows (`io.py`), computes Pareto
  fronts, hull vertices and hypervolume (`fronts.py`), and checks files against `pins/`
  (`provenance.py`). The per-experiment modules are `e1.py` and `fig1.py`, and the figure style
  lives in `figures.py`.
- **`pipelines/`** is a Snakemake pipeline that draws every figure. See "Running the pipeline"
  below.
- **`drivers/`** holds the benchmark drivers that produced the results, one directory per
  experiment. `drivers/README.md` names their derna2 commits.
- **`pins/`** holds checksums of the binaries and result files every figure depends on, with
  their provenance.
- **`inputs/`** holds the protein sequences, UTRs, the codon table and E2's motif sets.

## Install

Python 3.11–3.13.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.lock
.venv/bin/pip install -e '.[dev,pipeline]' --no-deps
```

`requirements.lock` pins the package versions the published figures were made with.

## Data

The figures are drawn from the results of the paper's benchmark runs, which include the fronts,
timings and refolded structures the drivers produced. Regenerating them takes hours to days, so
they are archived as a dataset instead. That dataset is not published yet. Once it is, unpack it
under `data/results/`, or point `results_root` in `pipelines/config/host.yaml` at it. Every file
is checked against `pins/results-<experiment>.sha256` before a figure is drawn.

```
data/results/
├── e1/            # E1 sweep (sweep.jsonl) and each walk's front and hull vertices (hull/)
├── e1_extra/      # unconstrained fronts of the E1 instances E2 did not screen
├── e1_fig4/       # Figure 4's frozen snapshot of e1
├── e2/            # E2 hull probes, the short/middle/long band fronts, the binding screen
├── e3/            # E3 figure bundle (paper-e4/)
├── fig1/          # Figure 1's CDS-only, constrained and modular fronts
├── fig1_hba/      # the overview figure for the hba UTR pair (not archived yet)
└── fig1_fronts/   # Figure 1's three fronts, per CDS–UTR constraint, and mRNAfold
```

Everything the pipeline writes goes to `data/derived/`. Nothing under `data/` is committed.

## The figures

| # | Experiment | Target |
|---|---|---|
| 1 | GLPE_HUMAN with the hbb UTRs | `fig1_figures`, `fig1_overview2` |
| 2 | E1, the whole front | `e1_figures` |
| 3 | E2, the price of motif-freedom | `e2_figures`, `e2_fig_composite` |
| 4 | E3, UTR-aware design | `e3_figures` |
| S1–S4 | supplement | `e1_memory_figures`, `e2_fig_supp_motif`, `e3_fig_supp_contacts`, `supp_proteins` |

Every front panel has x = CAI and y = MFE (kcal/mol), with hull vertices filled and off-hull
points open.

## Running the pipeline

`pipelines/run.sh` runs `pipelines/Snakefile` with `pipelines/config/config.yaml`. It then loads
`pipelines/config/host.yaml` if that file exists (copy `host.example.yaml`). A target is a rule
name or an output path relative to the repository root.

```bash
pipelines/run.sh --dry-run                                     # what would run
pipelines/run.sh                                               # every figure, from the results
pipelines/run.sh -- fig1_figures                               # one experiment's figures
pipelines/run.sh -- data/derived/figures/e1_fronts.pdf         # one output
pipelines/run.sh -- e1_runs                                    # regenerate E1 from binaries
```

`run.sh` passes `--rerun-triggers mtime`. An output is rebuilt when one of its input files is
newer, but not when code or config changes. To rebuild one, delete it or add `--forcerun <rule>`.
Output paths in the tables below are relative to `data/derived/`.

There are two kinds of target. **Figure targets** read the archived results and need no
binaries, and `all` is these. **Run targets** regenerate the results with `derna` and `mfe_fold`
built at `derna2.pin`, or taken from `derna2.bin_dir`. They take hours to days and run only when
named. Set `figures_from: runs` to draw the figures from your own runs instead. The comparator
tools are described in `pins/comparators-PROVENANCE.md`.

The `*_verify` targets check a built figure against the checksums in `pins/`. They need the
archived results, so CI does not run them.

### Provenance

These rules check the archived results against their manifests, build DERNA2 at the pinned
commit, and package one experiment's results for release.

| Rule | Output | Content |
|---|---|---|
| `check_results` | `stamps/results_{exp}.ok` | one experiment's files against `pins/results-{exp}.sha256` |
| `build_derna2` | `bin/derna2-{pin}/` | `derna` and `mfe_fold` built at `derna2.pin` |
| `package_deposit` | `deposit/{exp}.tar` | the files of one experiment's manifest, each one checked |

### E1, the whole front

E1 reports the time and memory each method needs to recover the Pareto front between CAI and
MFE. It runs on 85 proteins of up to 1,629 aa, which are UniProt entries, twelve longer native
proteins, the SARS-CoV-2 spike (P0DTC2) and five Cas orthologs. DERNA2 enumerates the whole
front exactly, on one core and on eight. The other four methods are DERNA1, LinearDesign,
LinearCDSfold and mRNAfold. Each of them optimizes a weighted sum of the two objectives, so each
is run at every hull weight λ that DERNA2 found, and can reach hull vertices at best. A cell is
limited to 2 h below 1,000 aa and 6–12 h above. Outputs are under `data/derived/e1/`,
`runs/e1/` and `figures/`.

| Rule | Output | Content |
|---|---|---|
| `e1_plot_fronts` | `figures/e1_fronts.pdf`, `e1/cells.tsv` | wall time and Pareto points recovered per cell against length |
| `e1_fig2_data` | `e1/fig2.json` | Figure 2's data, which is wall time, peak RSS and designs per arm, and the share of each front that no weight reaches |
| `e1_figures` | | the two above |
| `e1_config` | `runs/e1/config.json` | the driver config with this host's paths and the pin |
| `e1_verify` | `runs/e1/verify.done` | each comparator's λ map against derna2 |
| `e1_grid` | `runs/e1/grid.done` | the four timed passes into `runs/e1/sweep.jsonl` (needs a reserved socket, resumable) |
| `e1_summarize` | `runs/e1/summary.txt` | fronts, recovery and anytime curves |

#### Figure 4

Figure 4 summarizes E1. It gives each method's exit and wall time per protein, and how much of
DERNA2's front lies between hull vertices, where no weight reaches. It is drawn from a frozen
snapshot of the sweep in `data/results/e1_fig4`. The protein of panel (a) is `e1_fig4.focus` in
`config.yaml`. Outputs are under `data/derived/e1/fig4/` and `figures/`.

| Rule | Output | Content |
|---|---|---|
| `e1_fig4_cells` | `e1/fig4/cells.csv` | one row per (arm, protein) |
| `e1_fig4_per_protein` | `e1/fig4/per_protein.csv` | each arm's exit and wall time per protein |
| `e1_fig4_steps` | `e1/fig4/steps.npz` | every front, its hull vertices and the steps between neighbouring designs |
| `e1_fig4` | `figures/fig4_e1_overview.pdf` | the figure |
| `e1_fig4_verify` | `stamps/e1_fig4.verified` | the tables and figure against `pins/e1_fig4.outputs` |
| `e1_fig4_figures` | | the two above |
| `e1_runs` | | the run chain above |

### E2, the price of motif-freedom

E2 reports what forbidding sequence motifs costs in CAI and MFE. For each protein, DERNA2
computes the exact front over all designs, called F0, and then over the designs that avoid each
of four cumulative motif sets, called F1 to F4. The four sets are manufacturing sites,
processing signals, repeats and internal initiation, 24 motifs in all, listed under
`inputs/motifs/`. The proteins come from E1 and are grouped into short, middle and long bands.
There are no UTRs here. No other method can express this constraint, so there are no
comparators. Outputs are under `data/derived/e2/`.

| Rule | Output | Content |
|---|---|---|
| `e2_forced_motifs` | `e2/forced_motifs.jsonl` | motifs forced by a peptide, checked against the recorded census |
| `e2_probe_table` | `e2/probe.tsv`, `e2/probe_fit.txt` | wall time, peak RSS and vertices per rung of the hull probe, and the peak-RSS fit |
| `e2_tables` | | the two above |
| `e2_fig3_cells` | `e2/fig3_cells.json` | state, wall time, peak RSS and points for every cell of the three bands |
| `e2_fig3_panels` | `e2/fig3_panels.json` | coverage per protein and motif level, and peak-RSS fits |
| `e2_fig3_planned` | `e2/fig3_planned.json` | ΔMFE at matched CAI against the unconstrained front, per short-band cell |
| `e2_fig3_strip` | `e2/fig3_strip.json` | per motif, proteins admitting it against proteins whose front spells it |
| `e2_figures` | | Figure 3's data, the three above |

#### Figure 3 (`e2_fig_composite`)

Figure 3 has seven panels. Panel (a) shows which motifs each protein's unconstrained designs
contain, panel (b) what each motif level costs in MFE and CAI, and panels (c) to (g) one
protein's nested fronts and designs. Panel (a) also reads E1's sweep and the E1 extra screen,
which are named by `e2.e1_sweep` and `e2.e1_extra` in `host.yaml`. Outputs are under
`data/derived/e2/composite/` and `figures/`.

| Rule | Output | Content |
|---|---|---|
| `e2_composite_auc` | `e2/composite/auc01_fin.json` | the finished cells and each one's AUC per motif level |
| `e2_composite_data` | `e2/composite/composite-data.json` | the screen catalog of panel (a), and panel (b)'s per-level cost at four CAI floors |
| `e2_composite_dcai` | `e2/composite/dcai-at-mfe.json` | the CAI given up at four MFE floors, for panel (b) |
| `e2_composite_example` | `e2/composite/example-*.json` | one protein's nested fronts and designs for panels (c)–(g), needs ViennaRNA |
| `e2_fig_composite` | `figures/fig3_e2_composite.pdf`, `.png` | the figure |
| `e2_fig_composite_verify` | `stamps/e2_composite.verified` | the figure against `pins/e2_composite.outputs` |

#### Motif supplement (`e2_fig_supp_motif`)

The supplement shows three things. For each motif and protein, whether the protein can avoid
that motif, must contain it, or cannot contain it at all. Then how often each motif set still
appears on the front as levels are added. Then what each level costs in wall time, memory and
front size, relative to F0. Outputs are under `data/derived/e2/supp/` and `figures/`.

| Rule | Output | Content |
|---|---|---|
| `e2_supp_admits` | `e2/supp/admits.json` | per protein and motif, whether any synonymous sequence spells it |
| `e2_supp_presence` | `e2/supp/presence.json` | the share of front designs spelling each motif set, per protein and level |
| `e2_supp_cost` | `e2/supp/cost.json` | wall, peak RSS and front size per level, relative to the unconstrained run |
| `e2_fig_supp_motif` | `figures/fig_supp_motif.pdf`, `.png` | the figure |
| `e2_supp_table` | `figures/tab_supp_motifs.tex` | the 24 motifs by set |
| `e2_fig_supp_motif_verify` | `stamps/e2_supp_motif.verified` | the figure against `pins/e2_supp_motif.outputs` |

Both E2 figures refuse to save if any label overlaps another label or panel.

### E3, UTR-aware design

E3 reports what designing the CDS together with its UTRs delivers. It runs on 150 instances,
which are the 30 proteins under `inputs/e3/` paired with five UTR pairs used in practice. Those
pairs are hbb, hba, BNT162b2, mRNA-1273 and CureVac, listed in `inputs/utrs/`. Each instance has
four fronts, all re-scored on the whole mRNA by one folder. DERNA2 designs the CDS alone and the
UTRs are attached afterwards, DERNA2 designs with the UTRs in the model and CDS–UTR pairs
forbidden, DERNA2 computes the modular front, and mRNAfold runs with the UTRs in its energy
model. The figures read the E3 bundle at `data/results/e3` and E2's fronts at `data/results/e2`.
Outputs are under `data/derived/e3/` and `figures/`.

| Rule | Output | Content |
|---|---|---|
| `e3_hull_input` | `e3/hullcheck/<cell>.in.jsonl` | one cell's hull vertices as a derna result file |
| `e3_refold_deposit` | `e3/hullcheck/<cell>.out.jsonl` | the archived refold of one cell, which needs no binary |
| `e3_refold` | the same file | `derna check` of those vertices, run with `--force-run e3_refold` to recompute one |
| `e3_cand_input` | `e3/cand/<cell>.{conc,mrf}.in.jsonl` | the main figure's designs at the CAI floor |
| `e3_cand_refold` | `e3/cand/<cell>.{conc,mrf}.out.jsonl` | `derna check` of those |
| `e3_fig_story` | `figures/fig_e3_story.pdf` | the main figure |
| `e3_fig_story_verify` | `stamps/e3_story.verified` | the figure against `pins/e3_story.outputs` |
| `e3_story_stats` | `e3/story_stats.txt` | the numbers quoted in the paper's UTR section |
| `e3_story_stats_verify` | `stamps/e3_story_stats.verified` | that file against `pins/e3_story_stats.outputs` |
| `e3_contacts_data` | `e3/supp_contacts.json` | where CDS–UTR contacts fall, per UTR pair |
| `e3_fig_supp_contacts` | `figures/fig_supp_contacts.pdf` | the supplement figure |
| `e3_supp_contacts_verify` | `stamps/e3_supp_contacts.verified` | the histograms against `pins/e3_supp_contacts.outputs` |

`e3_figures` builds both figures and the statistics file.

### Figure 1, GLPE_HUMAN with the hbb UTRs

Figure 1 shows both of the paper's contributions on one instance, GLPE_HUMAN (P15421, 78 aa)
with the hbb UTRs. It draws DERNA2's exact CDS-only front, what attaching the UTRs does to those
designs, the front with CDS–UTR pairs forbidden, and the exact modular front under that
constraint. The overview figure compares three designs at CAI ≈ 0.87, one from mRNAfold, one a
hull vertex, and the nearest modular design. Its revised version is drawn for the hbb and hba
pairs. `drivers/fig1/FIG1.md` describes the method. Outputs are under `data/derived/fig1/`,
`runs/fig1/` and `figures/`.

| Rule | Output | Content |
|---|---|---|
| `fig1_so_front` | `fig1/so_front.jsonl` | the modular front from its window reports, checked against the archived one |
| `fig1_plot` | `figures/fig1_glpe_hbb.pdf` | panels a–d |
| `fig1_figures` | | the two above |
| `fig1_overview_designs` | `fig1/overview_designs.json` | the overview figure's three designs at CAI ≈ 0.87 |
| `fig1_overview` | `figures/fig1_overview.pdf` | the overview figure, needs LaTeX |
| `fig1_front` | `runs/fig1/{arm}/front.jsonl` | the CDS-only or constrained front, re-walked |
| `fig1_check_front` | `runs/fig1/{arm}/check.txt` | that front against the archived one |
| `fig1_runs` | | both fronts, re-walked and checked |
| `fig1_overview2_prep_hbb` | `fig1/overview2_hbb/` | the revised overview's inputs for the hbb pair |
| `fig1_overview2_prep_hba` | `fig1/overview2_hba/` | the same for the hba pair |
| `fig1_overview2_plot` | `figures/fig1_overview2_{pair}.pdf` | the revised overview figure, per pair |
| `fig1_overview2` | | both pairs, needs ViennaRNA and LaTeX with Times and `pifont` |

### Figure 1's three fronts

On the same instance, this compares DERNA2's plain front with its exact fronts under three
CDS–UTR constraints. For each constrained hull vertex it also finds the first modular design,
and mRNAfold is run at that vertex's weight for comparison. `drivers/fig1-fronts/README.md`
describes the runs. Outputs are under `data/derived/fig1_fronts/`.

| Rule | Output | Content |
|---|---|---|
| `fig1_fronts_tables` | `fig1_fronts/tables.md` | points, hull vertices, modular points and wall time per front, and mRNAfold per vertex |

### The proteins

The supplement table lists every protein in the paper with the experiments run on it. Outputs
are under `figures/`.

| Rule | Output | Content |
|---|---|---|
| `supp_proteins` | `figures/tab_supp_proteins.tex` | every protein in the paper and the experiments run on it |
| `supp_proteins_verify` | `stamps/supp_proteins.verified` | that table against `pins/supp_proteins.outputs` |

Each script under `pipelines/scripts/` also runs by hand. Its docstring gives the arguments.

## Tests

```bash
.venv/bin/pytest
```

The pipeline tests dry-run every target named in the tables above.

## Copyright and data credits

Copyright (c) 2026, El-Kebir Group. All rights reserved.

The input sequences under `inputs/` are redistributed unmodified from UniProtKB/Swiss-Prot,
NCBI RefSeq, a CureVac patent sequence listing, GSRS and the AlphaFold DB. `NOTICE` gives each
source and its terms.

DERNA2 and the refolds here use the Turner 2004 nearest-neighbor parameters from the NNDB.

- Mittal, A., Turner, D. H. & Mathews, D. H. NNDB: an expanded database of nearest neighbor
  parameters for predicting stability of nucleic acid secondary structures. *Journal of
  Molecular Biology* 436, 168549 (2024). doi:10.1016/j.jmb.2024.168549
- Turner, D. H. & Mathews, D. H. NNDB: the nearest neighbor parameter database for predicting
  stability of nucleic acid secondary structure. *Nucleic Acids Research* 38, D280–D282
  (2010). doi:10.1093/nar/gkp892
- Mathews, D. H. et al. Incorporating chemical modification constraints into a dynamic
  programming algorithm for prediction of RNA secondary structure. *Proceedings of the
  National Academy of Sciences* 101, 7287–7292 (2004). doi:10.1073/pnas.0401799101

## Related repositories

- [`elkebir-group/derna2`](https://github.com/elkebir-group/derna2), the derna2 method.
