# E1's auxiliary screen: deposit and provenance, frozen 2026-11-05

## What this is, and why it is separate from E1

`data/results/e1_extra` is **22 files, 169.0 MiB**, cut from the screening tree at
`e3-claims-out/e2-extra`. `pins/results-e1_extra.sha256` names every one and all 22 verify;
the manifest's own sha256 is
`56736dbcb96eb0155ab12851239be0fd10fdf67bdfc1302dbe6f64a2f3bc6a6d`.

It holds the eight E1 instances that **E2's own motif screen never covered** — the five Cas
orthologs `PRT-03` to `PRT-07`, the SARS-CoV-2 spike `P0DTC2`, and the two proteins added to
E1 after its first snapshot, `P22301` and `P41159` — screened here from E1's own
unconstrained fronts rather than from E2's grid. That is why it is a deposit of its own
rather than part of either: the rows are E1 instances, but the measurement is E2's screen,
and neither experiment's manifest is the honest home for it.

    f0/<cell>.jsonl       8 files: the unconstrained front of each instance (the bulk of the
                          169 MiB; PRT-07 alone is 43 MiB)
    corpus/<cell>.fasta   8 files: the sequence each front was designed for
    binding-e2extra.jsonl the screen's per-protein, per-set motif occurrences
    lattice-e2extra.jsonl the codon-lattice census over the same instances
    raw/                  4 files: two gzipped walk logs and two fetch stderr logs

## What reads it, and why the whole of `f0/` is deposited

`draft_panel_a.py`, which `composite_data.py` execs, takes a **set of cell names** from
`glob(f"{XTRA}/f0/*.jsonl")` and reads `binding-e2extra.jsonl` and `lattice-e2extra.jsonl`
by name; `supp_common.py` appends `f0/` to the directories it calls `os.listdir` on. Both
therefore depend on which files that directory *contains*, not merely on the contents of the
ones they open. A cut that kept only the files some run happened to open would silently
shrink those sets and change panel (a) without failing, so `f0/` is deposited whole.

`raw/` is read by nothing. It is kept because four small files are not worth a second
manifest, and because the two walk logs are the only record of how these fronts were
produced.

## Host paths

None: no file here names `/home` or `/scratch`. Unlike E2's and E3's deposits, this one
needed no scrub.

## What the freeze is checked against

Before this, `e2.e1_extra` was a path set in `host.yaml`, so panel (a) reproduced on a host
that happened to have the tree and nowhere else, and nothing checked that the tree had not
moved. The rules now resolve it through `source("e1_extra")` like any other deposit, and the
six rules that read it — `e2_composite_data` and the five supplement rules whose scripts
import `supp_common` — declare `source_ok("e1_extra")`, so its digests are a precondition of
the figures.

Gated on a rebuild: with both this deposit and E1's resolving through `source`, the eleven
outputs of `pins/e2_composite.outputs` and `pins/e2_supp_motif.outputs` were deleted,
`auc01_fin.json` with them, and all eleven came back byte-identical.
