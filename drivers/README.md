# Drivers

The programs that produced the results, so that the code beside a result is the code that
made it. They are standard-library Python and bash, and the pipeline renders their
configurations with this host's paths and calls them.

They came from the derna2 repository, and the table names the commit each was taken from.
They are **not** byte-identical to those copies any more: the calls into derna2 were
rewritten to the `derna` command line, because the binaries the original runs invoked are
not part of the published tool. What each one does, and the configuration it does it under,
is unchanged. A reader comparing a driver with its source commit should expect to find the
difference in the argument lists and in how a result file is read — `derna design` for a
design at one weight, `derna front` for a walk, `derna check` for a verdict, `mfe_fold` for
the free fold of one sequence.

| Directory | Source in derna2 | Content |
|---|---|---|
| `e1/` | `experiments/paper-benchmark/` on `experiment/e2-tamazirt` at `42be727e`, over the shared driver from `chore/paper-e2-plan` at `2483e9b9` | `drive_e1.py` (the whole-front grid, importing the `drive.py` beside it), `drive.py` and `measure.py` (the shared grid and one measured run), `isolate.sh` (moves this account's processes off the reserved socket), `run_e1_grid.sh` (the four passes behind an idle-core check), `pilot_concurrency.py` and `pilot_report.py` (how many single-core cells run at once without slowing each other), the DERNA v1 helpers (`build_derna1.sh`, `derna1_table.py`, `verify_derna1_map.py`) and `e1_cost_table.py` |
| `e2/` | `experiments/paper-benchmark/` on `chore/paper-e3-plan` at `5e098add` | `E2-PIN.md` (the pin and why); the grid's runners (`grid_short.sh`, `grid_threaded.sh`, `grid_resume.sh`, `pass2_rounds.sh` with `jobs-pass2.txt`, `stop_at_gap.sh`, `flank_arm.sh`); the binding screen and band selection (`screen_f0*.sh`, `screen_lattice.py`, `screen_occurs.py`, `a_prime_counts.py`, `select_long_band.py`, `select_middle_band.py`, with their records in `screen/` and `calib/`); the analysis scripts (`front_identity.py`, `convert_e1_front.py`, `check_hull_recompute.py`, `price_per_lambda.py`, `band_closeout.py`, `fix_summary.py`, `icount.py`); the hull and memory probes (`e2_hull_probe*.sh`, `fit_probe.py`, `mem_probe.sh`, `mem_ladder.sh`, `disc1082/`); and `forced_motifs.py`. The motif sets are under `inputs/motifs/` |
| `e3/` | `experiments/paper-benchmark/` on `chore/paper-e4-plan` at `1dc30fdb` | `drive_e3.py` (plan, run, score and summarize) and `measure.py`; its proteins are under `inputs/e3/` (with the source's `PIN` and `SHA256SUMS`) and its UTR pairs under `inputs/utrs/` |
| `fig1/` | `experiments/fig1-glpe-hbb/` on `experiment/fig1-glpe-hbb` at `0d56c280` | the scripts behind the overview figure, and `FIG1.md` (that directory's README) |
| `fig1-fronts/` | `experiments/glpe-hbb-fronts/` on `experiment/glpe-hbb-fronts-data` at `f67a04a6` | the scripts behind Figure 1's three fronts (`run_fronts.sh`, `run_front3.sh`, `run_mrnafold.py`, the checks, `make_tables.py`, `build_plot.py`), and that directory's README |

## The configurations

`e1/e1.config.json` is `config.e2.tamazirt.json` from `42be727e`, the configuration of the
run on the second host, with its paths turned into `{bin}`, `{derna2_bin}`,
`{lineardesign_dir}`, `{lcds_dir}`, `{inputs}`, `{out}` and `{root}` and its pin into
`{pin}`. `{derna2_bin}` is new: derna2's two binaries are built from source at a commit and
the comparators are not, so they no longer come from one directory. Its protein directory
and codon table, which that file names inside a derna2 checkout, are the identical copies
under `inputs/`, and the twelve 511–962 aa proteins it adds are under `inputs/e1-natives/`.

`e3/e3.config.json` is `config.e4.example.json` from `1dc30fdb` with the same placeholders.
It is the example configuration: the run's own, with A4's builds and memory cap, stays on
its host.

## What is a record rather than a route

The `e2/` runners, screens and probes are the record of what ran and carry the absolute
paths of the hosts they ran on; the pipeline does not call them. The hull probe in
particular ran a binary that predates the command line, so it cannot be repeated from a
published checkout at all. `forced_motifs.py` is the exception: it is run unchanged by
`pipelines/scripts/e2/forced_motifs.py`, which points it at `inputs/` instead of a derna2
checkout.

The `fig1/` scripts are kept as the record of what was run, not as pipeline steps. The
working directory they ran in, which the derna2 copies spell as an absolute path, is read
here from `FIG1_WORK` (default `work`). They still carry absolute paths to a derna2
checkout, and the window search among them ran a binary the command line does not have.

The `fig1-fronts/` README differs from its source in one sentence: its link to a hosted copy
of the interactive plot names the deposit's own copy
(`results/raw/plot/glpe-hbb-fronts.html`) instead. `make_tables.py` is the one script the
pipeline runs (`fig1_fronts_tables`).

`e1/drive.py` and `e1/measure.py` descend from the copies the **time** pass ran (and from
`chore/paper-e1-plan` at `24eade73`), plus a sync to disk after every row and a parallel
re-scorer. The **work** pass ran an earlier `drive.py`, from `c7bdbd32`, before the time
pass's concurrency, contamination and retry changes; its `measure.py` was the same.
`e3/measure.py` and `e1/isolate.sh` are the same files as `e1/`'s.

## Host paths, and why these files are no longer byte-identical

These drivers are vendored as the record of how each run was invoked, and the rule above is
that they stay byte-identical to the commit that produced the deposited results. On
2026-11-06 that rule was broken deliberately, in 33 of them, for one reason: they spelled
absolute paths on the machines that ran them, including **another contributor's home
directory**, 21 times. A published repository should not carry someone's account name and
directory layout, and no reading of byte-identity outweighs that.

Three substitutions, the same shape as the deposits' scrub (`pins/e2-PROVENANCE.md`):

    <contributor's checkout>/derna2 -> derna2
    <this host's checkout>/derna2   -> derna2
    <run root>/                     -> run/

Only the leading path changed. Every tree name the drivers referred to — `paper-e3`,
`glpe-hbb-fronts`, `paper-e2-bin` and the rest — is kept, because which run directory a
driver read is provenance and the user it belonged to is not. Nothing else in these files
was touched: no flag, no value, no comment, no commit stamp.

Nothing in the pipeline reads a changed file. The two drivers a rule does read,
`e1/drive_e1.py` and `e1/e1.config.json`, carried no host path and are untouched; the
rendered-config check (`tests/test_pipeline.py`) passes as before, and all 14 manifests and
output pins verify unchanged.
