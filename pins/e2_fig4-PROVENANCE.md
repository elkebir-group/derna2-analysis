# Figure 4's deposit: a frozen snapshot of the E2 sweep, verified 2026-10-30

## The deposit

`data/results/e2_fig4` holds two things from the E2 grid's output directory on the second
host (`/home/melkebir/projects/derna2/e2/out` on tamazirt):

- `sweep.jsonl`: the file saved there as `sweep.jsonl.before-add-2026-10-30`, the sweep as it
  stood before rows were added to it on 2026-10-30. 256,214 rows, sha256 `1a40f9f9…`. The
  live `sweep.jsonl` keeps growing, so it is not the deposit.
- `hull/`: the walk's per-protein outputs, of which the figure reads `<P>.front.jsonl` (the
  Pareto front as the walk wrote it) and `<P>.walk.json` (its hull sweep: vertices and
  `sweep_s`) for the 79 proteins of the corpus.

`pins/results-e2_fig4.sha256` names those 159 files. The directory may hold more; nothing
else under it is read.

## What the figure reads from it, and what that reproduces

- **Run table (panel d).** Cell rows stamped with the replay pin `8603d37d`, one per
  (arm, protein), the delivered and then faster attempt where a cell has two
  (`fig4_cells.py`). Byte-identical to the table the figure was first drawn from.
- **Per-protein table (panels a and e).** Each arm's exit and wall time from the later cell row
  per (arm, protein) of any pin, and the dichotomic hull search with Carrillo-Lipman pruning
  (the walk's `sweep_s`) and without it (projected as 2V−1 exact DP runs at derna2-hull's wall
  per λ) (`fig4_per_protein.py`). Equal, value for value, to those columns of E2's analysis
  table `out/analysis/e2_per_protein.csv` (sha256 `a52bba30…`, 2026-10-27). Two proteins of
  equal length appear there in set-iteration order; here they are sorted by name.
- **Front steps (panels a, b and c).** Each front Pareto-filtered on raw (energy, log CAI),
  which gives 153,199 points over the 79 proteins, and the sweep's 7,571 hull vertices
  matched into it (`fig4_steps.py`).

## Corpus

79 proteins: those with a derna2 walk cell, less P05160, Q06187, Q8TAA9 and Q9UNN5, whose
length band was stopped partway and never attempted for mRNAfold or LinearDesign.

## Reproducing the bytes

`pins/e2_fig4.outputs` lists the digests of the three tables and the figure;
`e2_fig4_verify` checks a build against them. The tables need only `requirements.lock`
(Python 3.12, numpy 2.2.6, pandas 2.3.2, matplotlib 3.10.9). The figure typesets its text
with LaTeX (`text.usetex`, Times via `mathptmx`), so its bytes also depend on the TeX
installation: the pinned digests were made with pdfTeX 3.1415926-2.5-1.40.14 (TeX Live
2013). Neither file carries a creation date, and every jitter is seeded, so the same inputs
and toolchain give the same bytes.
