# Drivers

The programs that produced the results, copied from the derna2 repository so that the code
beside a result is the code that made it. They are standard-library Python and bash; the
pipeline renders their configurations (`e1/e1.config.json`, `e2/e2.config.json`) with this
host's paths and calls them. Each experiment has its own copy of the shared driver
(`drive.py`, `measure.py`, `isolate.sh`), because the copies differ and each is the one its
experiment ran.

| Directory | Source in derna2 | Content |
|---|---|---|
| `e1/` | `experiments/paper-benchmark/` on `chore/paper-e2-plan` at `2483e9b9` | `drive.py` (E1 grid), `measure.py` (one measured run), `isolate.sh` (moves this account's processes off the reserved socket), `BENCHMARK.md` (that directory's README) |
| `e2/` | `experiments/paper-benchmark/` on `experiment/e2-tamazirt` at `42be727e` | `drive_e2.py` (E2 grid, imports the `drive.py` beside it), `run_e2_grid.sh` (the four passes behind an idle-core check), `pilot_concurrency.py` (how many single-core cells run at once without slowing each other), `pilot_report.py`, the DERNA v1 helpers, `e2_cost_table.py` and `E2-PLAN.md` |
| `e3/` | `experiments/paper-benchmark/` on `chore/paper-e3-plan` at `5e098add` | `E3-PLAN.md`, `E3-PIN.md` (the pin and why), `E3-GRID-PASS1.md`, `LONG-BAND-SUPPLY.md`, `PIN-ada.txt`; the grid's runners (`grid_short.sh`, `grid_threaded.sh`, `grid_resume.sh`, `pass2_rounds.sh` with `jobs-pass2.txt`, `stop_at_gap.sh`, `flank_arm.sh`); the binding screen and band selection (`screen_f0*.sh`, `screen_lattice.py`, `screen_occurs.py`, `a_prime_counts.py`, `select_long_band.py`, `select_middle_band.py`, with their records in `screen/` and `calib/`); the analysis scripts (`front_identity.py`, `convert_e2_front.py`, `check_hull_recompute.py`, `price_per_lambda.py`, `band_closeout.py`, `fix_summary.py`, `icount.py`); the hull and memory probes (`e3_hull_probe*.sh`, `fit_probe.py`, `mem_probe.sh`, `mem_ladder.sh`, `disc1082/`); and `forced_motifs.py`. The motif sets are under `inputs/motifs/` |
| `e4/` | `experiments/paper-benchmark/` on `chore/paper-e4-plan` at `1dc30fdb` | `drive_e4.py` (plan, run, score and summarize), `measure.py`, `E4-PLAN.md`; its proteins are under `inputs/e4/` (with the source's `PIN` and `SHA256SUMS`) and its UTR pairs under `inputs/utrs/` |
| `e5/` | `experiments/paper-benchmark/` on `chore/paper-e5-plan` at `10dda4c2` | `E5-PLAN.md` only: E5 has no driver yet |
| `fig1/` | `experiments/fig1-glpe-hbb/` on `experiment/fig1-glpe-hbb` at `0d56c280` | the scripts behind the overview figure, and `FIG1.md` (that directory's README) |

`e1/e1.config.json` is the configuration the deposited E1 timing run used
(`data/results/e1/config.json`), with its host paths turned into `{bin}`, `{inputs}` and
`{out}`. `{w}` is the driver's own beam-width field. `e2/e2.config.json` is
`config.e2.tamazirt.json` from `42be727e`, the configuration of the second host's run, with its
paths turned into `{bin}`, `{lineardesign_dir}`, `{lcds_dir}`, `{inputs}`, `{out}` and `{root}`
and its pin into `{pin}`; its protein directory and codon table, which that file names inside a
derna2 checkout, are the identical copies under `inputs/`, and the twelve 511–962 aa proteins it
adds are under `inputs/e2-natives/`. `e4/e4.config.json` is `config.e4.example.json` from
`1dc30fdb` with the same placeholders. It is the example configuration: the run's own, with A4's
builds and memory cap, stays on its host.

`BENCHMARK.md` differs from its source in one sentence, which pointed at a planning note outside
both repositories. `e1/drive.py` and `e1/measure.py` are byte-identical to the copies the E1
**time** pass ran (and to `chore/paper-e1-plan` at `24eade73`). The E1 **work** pass ran an
earlier `drive.py`, from `c7bdbd32`, before the time pass's concurrency, contamination and retry
changes. Its `measure.py` was the same. `e2/drive.py` is the E1 time-pass copy plus a sync to
disk after every row and a parallel re-scorer; `e2/measure.py`, `e4/measure.py` and
`e2/isolate.sh` equal E1's.

The `e3/` runners, screens and probes are the record of what ran and carry the absolute paths of
the hosts they ran on; the pipeline does not call them. `forced_motifs.py` is run unchanged
by `pipelines/scripts/e3/forced_motifs.py`, which points it at `inputs/` instead of the derna2
checkout. The plans refer to files by their derna2 paths (`experiments/paper-benchmark/…`).

The `fig1/` scripts are kept as the record of what was run, not as pipeline steps. They differ
from the derna2 copies in one respect: the working directory they ran in, which those copies
spell as an absolute path, is read here from `FIG1_WORK` (default `work`). They still carry
absolute paths to a derna2 checkout.
