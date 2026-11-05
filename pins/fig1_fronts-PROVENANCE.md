# Figure 1's three fronts: deposit and binaries, verified 2026-10-25, rescrubbed 2026-11-05

## The deposit

`data/results/fig1_fronts` is `experiments/glpe-hbb-fronts/results` of derna2
`f67a04a6acd62d52116a9889670d7040f4ebf882` (2026-10-12, *experiment(glpe-hbb): raw run
outputs for audit*), scrubbed of the host's paths as described below. It was byte-identical
to that tree when it was deposited, checked file by file, and
`pins/results-fig1_fronts.sha256` names the same 191 files throughout. The commit is one of
the development repository, which the published history does not carry, so it identifies the
originals rather than offering them.

## The scrub of 2026-11-05

The run's reports and logs named the filesystem they ran on. Four substitutions, applied
with `pipelines/scripts/scrub_host_paths.py --apply`, in this order:

    <checkout>/                        -> derna2/
    <cargo target directory>/          -> target/
    <run directory>/                   -> run/
    derna2/<worktree of the checkout>/ -> derna2/

The fourth exists because the first three leave the path that a worktree of the checkout
reached it by, so a fasta that the repository holds at
`experiments/scope-25-beam-runtime-quality/fixtures/GLPE_HUMAN_78.fasta` was still named
through that indirection. Collapsing it is what makes the recorded path the repository's own.

It rewrote **97 of the 191 files, every one of them under `raw/`**: the per-walk JSON
reports, the per-walk logs, the three arms' `stderr.log` and the three build and test logs.
No file that a rule reads changed — `front3.jsonl`, `front3_designs.scored.jsonl`, the
`v<N>.load` files and `mrnafold/compare.jsonl` are byte-identical to the deposited originals,
and everything outside `raw/` is identical as a tree. `fig1_fronts_tables` accordingly
rebuilds to the same bytes it did before the scrub
(`fd52c465731ff6ac6c42d982f4eef2ecb3864e974b1ea0ac6a4679bbb0e0b5a6`).

`pins/results-fig1_fronts.sha256` was regenerated afterwards and all 191 files match it. The
digests of the originals are the revision of that file immediately before this one in this
repository's history, whose own sha256 is
`36c922c7e886abba2581f05acfb0acbc76cb013b4343881756f2eda3bdeb9485`; comparing the two
revisions names the 97. The deposit is writable by nobody once deposited, and was unfrozen
for the rewrite and frozen again.

## The binaries

Built from derna2 `1b3739b7aab3fae889dd7a40bcd27abc981b5c83` (2026-10-06,
*experiment(glpe-hbb): 5'-unpaired guard, rescorer and drivers for three fronts*), which
the build directory records as `experiment/glpe-hbb-fronts = origin/xg 7bb14ec + step 0`;
`7bb14ecc6a7081f2d7f4ca44c570762a0887dd5f` is the merge of PR #394. They live in that run's
own build directory and are **not** pinned by a manifest here, so this file is their only
record.

    2e89e0e51daa20f00c2f1c014837c5f459120dc96ee0d2f60ef9fb67310048f4  glpe_hbb_rescore
    c9a808fba2022aae0adf16826cb738651061c6c9546cddc7b1907ed4952021fa  lazy_kbest_probe
    de0d7b3f00de29bb0ef449967bbe325d6823468d88824b756bbb2e1e8f60e7be  scope35_walk_front
    6e1fe849e463b9cf388e61799d680e11acf832a5ae2d142d8973ab3f4f26abd7  BUILT_FROM.txt

**`scope35_walk_front` is not the `scope35_walk_front` of `pins/fig1.sha256`.** That one is
`8a36ca7a…`, this one is `de0d7b3f…`: the two halves of Figure 1 ran on different builds of
a binary with the same name, and a panel that puts their rows on one axis is comparing
across builds. The name does not identify a build; the digest does.

## What is not recorded

No row of this deposit carries a build stamp — the arm reports hold the configuration
(region, threads, slacks, residues) and the timings, and `build_s` is seconds of build
time, not a build. Attribution runs through this file, not through the rows, so rows from
this deposit cannot be told apart from another build's by inspection alone.
