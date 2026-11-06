# E3: deposit, binary and refolds, verified 2026-11-04

## The deposit

`data/results/e3` points at the figure bundle's `paper-e4/`, wherever `results_root` puts it. `pins/results-e3.sha256` names its 664 files and `check_results` verifies them.

The bundle was assembled by running the plate with `open` and the stat calls wrapped and
recording every path the build actually touched, so it holds what the figures read rather than
the whole tree on the host; the bundle's own `README.md`, one level above `paper-e4/`, describes
that method and is not covered by the manifest. It was first cut for the 22-protein pool and
later extended with the eight middle-band proteins and with the refolds.

The census rows carry, per design, the gap between the constrained and the free fold, but not
the positions of the CDS–UTR base pairs. Those come from the refolds, which are deposited under
`hullcheck/` for that reason: a host with the deposit draws both figures without the binary.

## The binary

The deposited refolds were made by a `derna` built before the command line was published,
from a commit of the development repository that the public repository's history does not
carry. A commit that cannot be fetched is not a pin, so that build is identified here by
its bytes and nowhere else:

    41c91ffbd477571a3425de418b40b70de0fb2d595471e6869f697dfd20ef5646  derna

A rerun uses the `derna` the pipeline builds at `derna2.pin`, and its rows carry that
commit instead. The rows it writes are the deposited rows: `e3_story_stats` and
`e3_supp_contacts` rebuild byte-identical to `pins/e3_story_stats.outputs` and
`pins/e3_supp_contacts.outputs` from refolds made with a `derna` built from the public
repository (checked 2026-11-04). So nothing here needs the older build; it is recorded
because the deposited rows name it.

## The refolds

`data/derived/e3/hullcheck/<cell>.out.jsonl`, one per cell, 150 in all. They are in the deposit
and the manifest; `e3_refold_deposit` copies them out of it and `e3_refold` recomputes one from
the binary. A rebuild of `e3_figures` from the deposit alone reproduces the pinned figure
(`41ff2d8e`) and statistics file (`b7c895c5`) byte for byte.

Of the 150:

- **110 short-band cells** are the run the earlier `pins/e3_story.outputs` was cut on, copied on
  2026-11-04 from that run's own `data/derived/e4/hullcheck/` on the second host. Their
  `.in.jsonl` are byte-identical to the ones `e3_hull_input` writes here.
- **40 middle-band cells** were computed on the second host on 2026-11-03/04 with the binary
  above, two threads per cell.

31 of the 110 were independently recomputed on the second host during that run. Over their 3,183 points the
two copies agree on every design's `mfe`, `gap` and `violations`; the only fields that differ
anywhere are the header's `wall_s` and `peak_rss_kb`. That is the evidence that a refold does not
depend on the host it ran on.

They were deposited on 2026-11-04. Before that they lived only where they had been computed,
which is how the 110 short-band ones came to exist on one machine and nowhere else, and why the
figures could be rebuilt in exactly one place.

## The scrub of 2026-11-06

This deposit named the hosts it ran on: 311 of its 664 files carried a `/home/<user>` or
`/scratch/<user>` path. Those paths were run directories and **another user's home directory** (`/home/xinyu/e4`, 24 occurrences). Substitutions, applied with
`pipelines/scripts/scrub_host_paths.py --apply`, which reports before it writes and exits
non-zero on any path its substitutions miss:

    /scratch/xinyugu4/paper-e4 -> results/e3   (this deposit's own tree)
    /scratch/xinyugu4/paper-e3 -> results/e2   (E2's fronts, which this reads)
    /home/xinyu/e4             -> run/e3       (the same experiment on another
                                                host; kept distinct from the
                                                deposit rather than conflated)

Afterwards no file here contains either string, checked independently of the tool, and
`pins/results-e3.sha256` was regenerated: exactly 311 of its 664 digests moved, which is
the same count as the files that carried a path.

**The paths were not only cosmetic, and that is what made this more than a rewrite.**
The `source`, `front_from`, `census` and `from` fields of the `a0.jsonl` rows are paths the readers **open**, and they point into E2's tree as well as this one. The repository already had the mechanism for that -- `e3common.relocate` rewrites a recorded prefix to wherever the tree is mounted -- but the prefixes it was given, `e3.recorded_prefixes` in the config, were themselves this host's directories. They are now the deposits' own names, so the committed configuration no longer spells anyone's filesystem either. `run/e3` is deliberately not relocatable: it appears only in `derived_from` fields, which record what a row came from and are never opened.

What the rewrite is checked against: every figure and table this deposit feeds rebuilt
**byte-identical** -- fig_e3_story.pdf, story_stats.txt and supp_contacts.json. No pin of E3's moved.
