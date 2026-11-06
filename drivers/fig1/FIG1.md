# GLPE_HUMAN with hbb UTRs: the data behind the overview figure

One real instance carries both contributions: the exact CDS-only Pareto front, what attaching
UTRs does to its designs, a guarded front that forbids CDS–UTR pairs, and the exact
**structurally optimal front** under that guard.

## Instance

GLPE_HUMAN (P15421), all 78 residues
([fixture](../scope-25-beam-runtime-quality/fixtures/GLPE_HUMAN_78.fasta)), with the hbb flanks
from [`scope-19-utr-support/m5`](../scope-19-utr-support/m5/): 50 nt 5′, 234 nt CDS, 132 nt 3′,
416 nt in all. Turner 2004, dangles 0, internal-loop cap 30. Codon weights:
`packages/derna-codon-dp/data/codon_usage_human.csv`. `cai_log` is Σ ln w over the 78 codons.

Guard: **symmetric** (`--region 5-cds,cds-3`, `--table symmetric`): no pair between the CDS and
either UTR; 5′–3′ UTR pairs allowed. A design is **structurally optimal** when its free MFE equals
its guarded MFE (`F = E`).

## Results

| file | contents |
|---|---|
| `results/a_cdsonly_front.jsonl` | exact CDS-only front, 145 points; `hull_vertex` marks the 21 vertices of the hull in (`cai_log`, E) |
| `results/b_cdsonly_with_utrs.scored.jsonl` | every CDS-only design: fold alone (`cds_mfe_kcal_mol`, `structure`), free fold with the UTRs (`free_mfe_kcal_mol`, `free_structure`), and guarded energies |
| `results/cd_symmetric_guarded_front.scored.jsonl` | exact symmetric-guarded front, 117 points: guarded `structure` and `free_structure` per point |
| `results/d_so_front.jsonl` | **exact structurally optimal front, 114 points** |

- Attaching the UTRs lowers each CDS-only design's free MFE by 33.0–42.5 kcal/mol, and none of the
  145 is structurally optimal under the guard.
- The structurally optimal front spans E −180.3 … −123.9 kcal/mol and `cai_log` −24.86 … −0.58. It
  lies 0–10.3 kcal/mol above the guarded front (median 3.0).
- No design with `cai_log` above −0.5816 is structurally optimal: the most adapted structurally
  optimal design gives up 0.58 log-CAI against the maximum.

## Method

**Fronts.** `scope35_walk_front`, CDS only, unguarded with UTRs, and `--region 5-cds,cds-3`. The
binary includes the scan's interior-loop and multiloop leading-run fixes; the fronts reproduce the
energies of earlier runs on the same instance at every point (16 of 145 CDS-only points are spelled
differently at equal energy).

**Structurally optimal front.** Seeds: 21 structurally optimal designs from one lazy k-best walk
per guarded hull vertex ([`glpe-hbb-fronts`](../glpe-hbb-fronts/)), 16 of them non-dominated
(`so_windows/seed_stair.jsonl`). Between each pair of adjacent seeds, a `scope37_window` walk
capped by the upper seed (`--lo-seq`, `--incumbent`, `--slope band --front <guarded front>`) finds
every structurally optimal design in the stretch that the cap does not dominate. Below the
least-adapted seed, the same walk runs from the least possible `cai_log`.

- Every window over `cai_log` [−77.64, −0.5816] **completed** (`complete: true`). Windows at the
  adapted end span a wide energy range and completed only with the energy bound
  `--energy-slopes 0,0.25x,0.5x,2x,4x` (one: 5.7 min, 1.4 GB with it; unfinished at 1 h, 12 GB
  without). One window stopped at its memory cap; its finds bounded a narrower second window, which
  completed (`so_windows/r2_w8/`).
- Above the most adapted seed, no design caps the stretch. A ceiling walk proved that no structurally
  optimal design there has E ≤ −115 at `cai_log` > −0.4945, and then stopped at its memory cap. The
  stretch was then settled **exhaustively**: all 156,698 designs with `cai_log` ≥ −0.5816 were
  folded free (`fold_seq`) and under the guard (`guard_fold`), and only the seed at −0.5816 is
  structurally optimal (`top_exhaustive/`).
- Window reports are under `so_windows/` by round (`r1` plain, `r2_eb`/`r3`/`r3b` with the energy
  bound, `r2_w8` the follow-up window, `top` the ceiling walk).

**Checks.** All 114 points refold with free MFE = guarded MFE = the stated energy, and `cai_log`
recomputed from the sequence matches. `fold_seq` and `guard_fold` were first checked against known
values: 20 structurally optimal designs read `E = F`, and 12 sampled guarded-front points reproduce
their energies.

## Reproduce

Frozen binaries and run logs: `run/fig1-glpe-hbb/` (checksums in
`results/BINARIES.sha256`). The scripts in `scripts/` were run from a working directory given by
`FIG1_WORK` (binaries under its `bin/`), and carry absolute paths to a derna2 checkout; set both
before re-running.

```text
scripts/run.sh                                     # the four fronts
scripts/windows.py STAIR.jsonl OUT JOBS [k ...]    # capped windows (EXTRA="--energy-slopes ..." for the bound)
scripts/top_ceiling.sh C                           # ceiling walk above the most adapted seed
scripts/enumerate_top.py; scripts/top_free_folds.sh; scripts/post.py   # exhaustive top stretch
scripts/assemble.py                                # merge seeds + window finds into the front
```
