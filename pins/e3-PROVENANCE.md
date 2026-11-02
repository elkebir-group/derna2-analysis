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

`derna` at the pin `db5982d8` of config `e3.derna_pin`. The deposited refolds were made with
that build; a rerun uses the `derna` the pipeline builds, so its rows carry its own commit.
The build the deposit was made with hashes to:

    41c91ffbd477571a3425de418b40b70de0fb2d595471e6869f697dfd20ef5646  derna

It is not pinned by a manifest here, so this file is its only record. A commit does not identify
a binary, so a rebuild at the same pin in another working tree may differ.

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
