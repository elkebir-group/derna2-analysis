# E1's deposit: the whole-front comparison, frozen 2026-11-05

## The deposit

`data/results/e1` is **368 files, 1.364 GiB**: the grid's `sweep.jsonl` and its `hull/`
directory, taken from the grid's output tree on the second host
(`e2-out-1031`). `pins/results-e1.sha256` names every one of them and all 368 verify; the
manifest's own sha256 is
`9109ed05c0de96f8fa714ff34913743dfb61395af22831dc812fb8a40f0835d6`.

    sweep.jsonl   258,083 rows, every one stamped 8603d37d
                  sha256 9d2bd7e64991d88d2c4fff1aeb4655a409ef1727479ee598efabc65dc883a9c7
    hull/         367 files: 85 each of walk.json, vertices.json and front.jsonl, and 56
                  each of derna2-walk.walk.json and derna2-walk.front.jsonl

`sweep.jsonl` is byte-identical to `data/results/e1_fig4/sweep.jsonl`: Figure 4's deposit is
a frozen snapshot of this same sweep, pinning 171 of these files
(`pins/e1_fig4-PROVENANCE.md`). The two are not alternatives — Figure 4 pins the subset it
draws, and this manifest covers everything E1's own figures read.

The sweep carries **89** distinct `protein` values, not 85. Four are `P0DTC2_pNN` ladder
rungs, prefixes of one protein rather than instances, and appear under one arm only; the
corpus is the other 85 (`inputs/e2/corpus.json`). A reader that counts proteins by taking
the sweep's distinct values gets the wrong population, which is why `draft_panel_a.py`
filters the `P0DTC2_p` prefix before intersecting.

`e1.pin` is `8603d37d`, and the figures read only rows carrying it: a sweep file accumulates
across re-pins, and two pins' rows must not be mixed. It selects rows and drives no build,
so it may name a commit the public history does not carry. Every row here carries it, so the
filter is currently a no-op — but it is not redundant, because it is what keeps a later
re-pin's rows out of these figures.

## What the grid completed, and what it did not

The grid's last pass ended **2026-10-31 00:04**. It is deposited as finished, with its
timeouts as results rather than as work outstanding. The cells not finished, per arm, of 85
proteins:

    DERNA2          2
    mRNAfold       19
    LinearDesign   21
    LinearCDSfold  50
    DERNA1         57

Those are the figures' own counts, printed by `e1_plot_overview` as `cells not finished per
arm`. The overview's pruning speed-up is accordingly measured on **68** of the 85, and the
figure says so rather than reporting a mean over an unstated subset. A reader comparing arms
must read these counts with the arm: an arm that finished 28 cells and one that finished 83
are not being compared over the same instances, and the longest instances are exactly the
ones the slower arms dropped.

## Host paths

None. The deposit was checked before freezing: no file under `hull/` and no row of the
258,083 names `/home` or `/scratch`. Unlike Figure 1's three fronts, this deposit needed no
scrub — the run wrote its reports without recording the filesystem it ran on.

The earlier note that "one of the seven files under `hull/` names the host" described a
different tree: until 2026-11-05 `data/results/e1` pointed at the first host's **pilot**,
437 rows over nine `P0DTC2_pNN` ladder rungs stamped `acfaae60`. That pilot is not part of
the grid and is not deposited.

## What the freeze is checked against

All eleven rules that read a deposit here now declare its stamp, E1's three included, so the
digests are a precondition of every E1 figure rather than a record kept beside it. Before
this, `e1_plot_fronts`, `e1_plot_overview` and `e1_fig2_data` passed `[]` as their `ok` input
whenever `figures_from` was `deposit` — they waited on nothing at all, and a figure could be
drawn from a tree that had moved underneath it. They call `source_ok("e1")` now, like E3's
and Figure 1's.

The freeze was gated on a rebuild: the six outputs of `e1_figures` were deleted and remade
from the frozen deposit. `fig2.json`, `cells.tsv` and both PNGs came back **byte-identical**;
the two PDFs are identical once their `/CreationDate` is stripped, which is the only respect
in which a matplotlib PDF is not reproducible. There is no byte pin for the PDFs for that
reason.
