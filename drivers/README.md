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
| `e2/` | `experiments/paper-benchmark/` on `experiment/e2-tamazirt` at `dd57114b` | `drive_e2.py` (E2 grid, imports the `drive.py` beside it), `run_e2_grid.sh` (the four passes behind an idle-core check), `pilot_concurrency.py` (how many single-core cells run at once without slowing each other), the DERNA v1 helpers, `e2_cost_table.py` and `E2-PLAN.md` |
| `fig1/` | `experiments/fig1-glpe-hbb/` on `experiment/fig1-glpe-hbb` at `0d56c280` | the scripts behind the overview figure, and `FIG1.md` (that directory's README) |

`e1/e1.config.json` is the configuration the deposited E1 timing run used
(`data/results/e1/config.json`), with its host paths turned into `{bin}`, `{inputs}` and
`{out}`. `{w}` is the driver's own beam-width field. `e2/e2.config.json` is
`config.e2.tamazirt.json` from `dd57114b`, the configuration of the second host's run, with its
paths turned into `{bin}`, `{lineardesign_dir}`, `{lcds_dir}`, `{inputs}`, `{out}` and `{root}`
and its pin into `{pin}`; its protein directory and codon table, which that file names inside a
derna2 checkout, are the identical copies under `inputs/`.

`BENCHMARK.md` differs from its source in one sentence, which pointed at a planning note outside
both repositories. `e1/drive.py` and `e1/measure.py` are byte-identical to the copies the E1
**time** pass ran (and to `chore/paper-e1-plan` at `24eade73`). The E1 **work** pass ran an
earlier `drive.py`, from `c7bdbd32`, before the time pass's concurrency, contamination and retry
changes. Its `measure.py` was the same. `e2/drive.py` is the E1 time-pass copy plus a sync to
disk after every row; `e2/measure.py` and `e2/isolate.sh` equal E1's.

The `fig1/` scripts are kept as the record of what was run, not as pipeline steps. They differ
from the derna2 copies in one respect: the working directory they ran in, which those copies
spell as an absolute path, is read here from `FIG1_WORK` (default `work`). They still carry
absolute paths to a derna2 checkout.
