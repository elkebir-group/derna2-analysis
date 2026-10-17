# Figure 4's deposit: a frozen snapshot of the E1 sweep, verified 2026-10-31

## The deposit

`data/results/e1_fig4` holds two things from the E1 grid's output directory on the second
host:

- `sweep.jsonl`: the sweep as it stood when the grid's last pass ended (2026-10-31 00:04),
  after the rows added on 2026-10-30/31: P22301 and P41159 on every arm, and the middle band's
  cells that ran. 258,083 rows, every one stamped `8603d37d`, sha256 `9d2bd7e6…`. Its first
  256,214 rows are byte for byte the earlier deposit (`1a40f9f9…`, saved on that host as
  `sweep.jsonl.before-add-2026-10-30`).
- `hull/`: the walk's per-protein outputs, of which the figure reads `<P>.front.jsonl` (the
  Pareto front as the walk wrote it) and `<P>.walk.json` (its hull sweep: vertices and
  `sweep_s`) for the 85 proteins of the corpus. The 158 files of the earlier deposit's 79
  proteins are unchanged.

`pins/results-e1_fig4.sha256` names those 171 files. The directory may hold more; nothing
else under it is read.

## What the figure reads from it, and what that reproduces

- **Run table (panel d).** Cell rows stamped with the replay pin `8603d37d`, one per
  (arm, protein), the delivered and then faster attempt where a cell has two
  (`fig4_cells.py`). P41159's `derna1-own` cell has a crash row from a duplicate driver
  beside the delivered re-run; the delivered one is kept.
- **Per-protein table (panels a and e).** Each arm's exit and wall time from the later cell row
  per (arm, protein) of any pin, and the dichotomic hull search with Carrillo-Lipman pruning
  (the walk's `sweep_s`) and without it (projected as 2V−1 exact DP runs at derna2-hull's wall
  per λ) (`fig4_per_protein.py`). Equal, value for value, to those columns of E1's analysis
  table `out/analysis/e1_per_protein.csv` as regenerated over all 85 proteins on 2026-10-31.
  Two proteins of equal length appear there in set-iteration order; here they are sorted by
  name.
- **Front steps (panels a, b and c).** Each front Pareto-filtered on raw (energy, log CAI),
  which gives 165,887 points over the 85 proteins, and the sweep's 8,246 hull vertices
  matched into it (`fig4_steps.py`).

## Corpus

85 proteins: every protein with a derna2 walk cell. Four of them, P05160, Q06187, Q8TAA9 and
Q9UNN5 (524-661 aa), belong to a length band that was stopped partway, and seven of their
cells were never run: derna2-hull on Q06187, and mRNAfold and LinearDesign on Q9UNN5, Q06187
and P05160. Timing one λ of each and multiplying by its planned λ count puts every one of
the seven over the 2 h budget (2.3-4.2 h), and the four cells run that way all timed out,
taking 12-22% longer per λ than the one-λ timing. Panel (d) counts a cell that never ran
among the runs a method did not complete, and draws no point for it.

## Reproducing the bytes

`pins/e1_fig4.outputs` lists the digests of the three tables and the figure;
`e1_fig4_verify` checks a build against them. The tables need only `requirements.lock`
(Python 3.12, numpy 2.2.6, pandas 2.3.2, matplotlib 3.10.9). The figure typesets its text
with LaTeX (`text.usetex`, Times via `mathptmx`), so its bytes also depend on the TeX
installation: the pinned digests were made with pdfTeX 3.1415926-2.5-1.40.14 (TeX Live
2013). Neither file carries a creation date, and every jitter is seeded, so the same inputs
and toolchain give the same bytes.
