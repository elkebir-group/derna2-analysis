# Figure 1's data panels: binaries and what identifies them, verified 2026-10-25

## The binaries

`pins/fig1.sha256` pins four binaries, which sit in the directory `fig1.bin_dir` names and
match that manifest exactly on the host that holds them. **Their source commit was never recorded**: the
build directory holds no `BUILT_FROM.txt`, so the figure is reproducible from the freeze
but not rebuildable from source until a candidate commit is identified.

## What ties the deposited rows to those binaries

The window reports under `data/results/fig1/so_windows` carry `binary_fnv64`: 19 report
rows across the round directories, each repeated in the matching `.log`, and one distinct
value between them, `1ff256ead61c92ba`. That value is FNV-1a over the whole
binary file — `derna2_analysis.provenance.fnv1a64` — and it is `scope37_window` of the
manifest above and of no other pinned binary:

    scope37_window      1ff256ead61c92ba   <- the stamp on every window report
    scope35_walk_front  15e4032c88c1a27c
    fold_seq            403ef0a281ee489b
    guard_fold          362f759135fad5c8

So the rows are tied to a pinned file even though the commit is unknown, and the same
function tests a candidate commit: build it, hash the binary, compare. A candidate whose
`scope37_window` hashes to `1ff256ead61c92ba` is the build; any other is not.

What carries no stamp: the four front files (`a_cdsonly_front.jsonl` and its siblings),
and the 18 design rows in `so_windows/seed_stair.jsonl` and `w8b_stair.jsonl`, which hold
a sequence and its two coordinates and nothing about the run. Only the reports are
attributable by inspection; everything else is attributable through this file and
`pins/results-fig1.sha256`.

## A name is not a build

`scope35_walk_front` here is `8a36ca7a…`; the binary of the same name behind Figure 1's
three fronts is `de0d7b3f…` (see `fig1_fronts-PROVENANCE.md`). Rows from the two deposits
are not comparable on a shared axis without saying so.
