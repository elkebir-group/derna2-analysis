# Drivers

The programs that produced the deposited results, copied from the derna2 repository so
that the code beside a result is the code that made it. They are standard-library Python and
bash; the pipeline renders their configurations (`benchmark/e1.config.json`,
`benchmark/e2.config.json`) with this host's paths and calls them.

| Directory | Source in derna2 | Content |
|---|---|---|
| `benchmark/` | `experiments/paper-benchmark/` on `chore/paper-e2-plan` at `2483e9b9` | `drive.py` (E1 grid), `drive_e2.py` (E2 grid, imports `drive.py`), `measure.py` (one measured run), `isolate.sh` (moves this account's processes off the reserved socket), the DERNA v1 helpers, `BENCHMARK.md` (that directory's README) and `E2-PLAN.md` |
| `fig1/` | `experiments/fig1-glpe-hbb/` on `experiment/fig1-glpe-hbb` at `0d56c280` | the scripts behind the overview figure, and `FIG1.md` (that directory's README) |

`benchmark/e1.config.json` is the configuration the deposited E1 timing run used
(`data/results/e1/config.json`), with its host paths turned into `{bin}`, `{inputs}` and
`{out}`. `{w}` is the driver's own beam-width field. `benchmark/e2.config.json` is
`config.e2.example.json` from the same commit, with the same substitution.

`BENCHMARK.md` differs from its source in one sentence, which pointed at a planning note outside
both repositories. `drive.py` and `measure.py` here are byte-identical to the copies the E1 **time** pass ran
(and to `chore/paper-e1-plan` at `24eade73`). The E1 **work** pass ran an earlier `drive.py`,
from `c7bdbd32`, before the time pass's concurrency, contamination and retry changes. Its
`measure.py` was the same.

The `fig1/` scripts are kept as the record of what was run, not as pipeline steps. They differ
from the derna2 copies in one respect: the working directory they ran in, which those copies
spell as an absolute path, is read here from `FIG1_WORK` (default `work`). They still carry
absolute paths to a derna2 checkout.
