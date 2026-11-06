# E2's deposit: what it holds, how it was cut, and what the cut is checked against

## The deposit

`data/results/e2` is a cut of the grid's run tree: **2,514 files, 2.099 GiB**, against the
tree's own 132.458 GiB. `pins/results-e2.sha256` names every one of them and all 2,514
verify; the manifest's own sha256 is
`f1f65539f439ee19991a278b61f129588b3340b290e70a0717c594c355ba6486`.

The run tree is not a checkout of anything — it holds no `.git` — so no commit identifies
it. What identifies the deposit is the manifest above, and what identifies the *build* the
rows were produced by is the same thing that identifies E3's: a checksum, in
`pins/e3-PROVENANCE.md`. The grid ran on two hosts (`e2.hosts`), and its early rounds
recorded which only through a home directory in each line's environment.

Before this cut the figures read the run tree directly and `e2.pin` was `null`, so E2 was
the one experiment whose figures had no deposit and no checksums to wait on: `check_results`
had no `pins/results-e2.sha256` to read, and none of the eleven rules that reach the tree
declared a stamp. They all declare `source_ok("e2")` now, so the deposit's digests are a
precondition of every E2 figure rather than a record kept beside it.

## How the cut was chosen

A `sys.addaudithook` hook recorded every `open` of every process in a full run of the E2 and
E3 chains, as the path the caller passed and the path it resolved to. That audit names
**930 files, 2.086 GiB** under the run tree.

**An opens audit cannot see a directory listing, and these readers depend on one.** They do
not merely locate files with a glob, they take a glob's matches as a *set*:
`composite_data.py` builds its cell names from `{E}/{band}/*_F[1-4].jsonl`, `fig2.py` and
`auc01_fin.py` iterate one, `fig3_planned.py` intersects two, and `supp_common.py` calls
`os.listdir` on each `screen/<arm>`. A file that was globbed but never opened would be
absent from an opens-only cut, those sets would silently shrink, and the figures would
change — which reads as a result rather than as an error. So the unit copied is the
**directory**: every file in every directory that contributed a read. Below the root that
costs 0.013 GiB over the opens alone, so the robustness is very nearly free.

**The root is the exception, and there it is not free.** Every glob that reaches the run
tree's own root spells a *directory* wildcard — `{E}/*/{cell}_{nest}.r*.jsonl` in
`e3common.py`, `example_coords.py` and `pick_example3.py` — so what the root supplies a
reader is its subdirectory entries, never its files. `best_round` is always called on a band
directory, and `supp_common.py` lists only `screen/<arm>`. Of the root's 206 files exactly
one, `FINISHED.tsv`, matches any pattern any reader can apply there, and that one is opened
by name. The other 205 are the run's shell scripts, pid files, editor backups and session
notes: 1.3 MiB that no reader can reach, and that name the host and the development tooling.
Taking the root whole would have published them. It is cut to what was opened.

The fifteen directories that contributed a read, with what was opened in each against what
the deposit holds:

    .                  1 / 1      (the root: opens only; see above)
    e3-probe           7 / 10
    e3-probe2          8 / 12
    grid             102 / 254
    long               1 / 6
    long/cells        36 / 248
    mid              435 / 1183
    motifs             4 / 5
    screen             7 / 13
    screen/f0         61 / 62
    screen/f0-five     5 / 5
    screen/f0-l       19 / 33
    screen/f0-xl      18 / 25
    short2           130 / 356
    short3            96 / 301

## What the cut is checked against, and what the check can and cannot show

A superset argument is not a proof. Taking a directory whole answers the glob that selects
members of a directory the audit entered; it cannot say whether some reader walks a
directory the audit never entered at all. Only rebuilding can, so the gate is a rebuild:
with `data/results/e2` repointed at the cut, the thirteen outputs the four pins name were
deleted and remade, and every one came back **byte-identical** to the cut it was published
with — the five of `pins/e2_composite.outputs`, the six of `pins/e2_supp_motif.outputs`, and
E3's `pins/e3_story_stats.outputs` and `pins/e3_supp_contacts.outputs`, which read E2's tree
through `e3_hull_input`. A missing directory listing is the one defect that test fails.

E3's chain had to be forced to be audited at all. `e3_refold_deposit` is preferred over
`e3_refold`, so a plain run executes it 150 times and **never runs `e3_hull_input`** — and
that is the rule that relocates E3's recorded paths into E2's tree. A run that reported "no
reads of E2" was therefore reporting which rule it had skipped, not what the chain reads.
`e3_hull_input` was forced over all 150 cells under the hook; it read `screen/f0/*.jsonl`,
and those files were already in the 930.

Two reaches outside this deposit remain, both in the composite plate's panel (a): `e1_sweep`
(the proteins the screen is scored over) and `e1_extra` (an auxiliary tree of E1 instances
screened for this figure). Both are `null` in the committed config and set in `host.yaml`,
so that chain reproduces on a host that has both and nowhere else. Depositing E1 is what
closes it, and it must be closed by depositing E1's sweep — not by depositing a list derived
from it, which would make the figure's population unauditable.

## What is not recorded

No row of this deposit carries a build stamp. Attribution runs through this file and
`pins/e3-PROVENANCE.md`, not through the rows, so rows from this build cannot be told apart
from another's by inspection alone.

## The scrub of 2026-11-06

This deposit named the hosts it ran on: 610 of its 2,514 files carried a `/home/<user>` or
`/scratch/<user>` path. Those paths were run directories, worktrees of the checkout, and **another user's home directory** (`/home/melkebir/projects/derna2`, 96 occurrences), so the deposit disclosed a third party's account and layout as well as this host's. Substitutions, applied with
`pipelines/scripts/scrub_host_paths.py --apply`, which reports before it writes and exits
non-zero on any path its substitutions miss:

    /scratch/xinyugu4/paper-e3        -> results/e2   (this deposit's own tree)
    /home/melkebir/projects/derna2    -> derna2       (the second host's checkout)
    /scratch/xinyugu4/wt-f6-carry     -> derna2       (worktrees of that checkout)
    /scratch/xinyugu4/wt-midslice-secs-> derna2
    /scratch/xinyugu4/wt-carry15      -> derna2
    /scratch/xinyugu4/wt-levers2      -> derna2

Afterwards no file here contains either string, checked independently of the tool, and
`pins/results-e2.sha256` was regenerated: exactly 610 of its 2,514 digests moved, which is
the same count as the files that carried a path.

**The paths were not only cosmetic, and that is what made this more than a rewrite.**
`FINISHED.tsv`'s `round_file` column is a path the pipeline **opens**, not a note about where a file once was. Spelled absolutely, it resolves only on the machine that wrote it, so a reader who unpacked the published tar anywhere else could not have run `e2_composite_auc`, `e2_supp_cost` or `e2_supp_presence` at all -- the scrub turned that latent defect into three failing rules, which is how it was found. The column is now relative to the deposit root, and `auc01_fin.py` and `supp_common.py` resolve it against this deposit; both accept the absolute spelling too, so a tree cut before today still loads. `supp_common.py` had already noticed half of this -- its digest over the column uses the basename, precisely because an absolute path 'moves when the deposit is mounted anywhere else' -- and kept opening the path.

What the rewrite is checked against: every figure and table this deposit feeds rebuilt
**byte-identical** -- fig3_e2_composite.pdf and .png, fig_supp_motif.pdf and .png, tab_supp_motifs.tex, tab_supp_proteins.tex, and the caches dcai-at-mfe.json and example-Q6PF06-c80-none.json. Two pins were re-cut, and their headers say what moved: `e2_composite.outputs` for composite-data.json, whose 135 recorded front paths and 49 front digests are provenance of which file each measurement came from, and `e2_supp_motif.outputs` for one leaf in each of three caches, the sha256 of the script that gained the resolver. No measured value in any of them changed.
