# P26436: disturbed by an external signal, not by the budget

**2026-10-19.** While stopping the E2 grid, another session killed processes by
binary path under `e2/bin/`. The F0 screen runs `e2/bin/derna`, so the cell in
flight took a SIGTERM at 19:05:46 — PIDs 995265 (`timeout`) and 995266
(`derna`), started 18:42:41 on `P26436`.

The screen was already SIGSTOPped at 18:50:00 for a reserved timing window, so
**the signal did not kill anything; it queued.** Read off `/proc/<pid>/status`
at 19:07:45, with both processes in state `T`:

    995265 timeout  ShdPnd 0000000000004000   SigIgn ...0001  SigCgt ...010000
    995266 derna    ShdPnd 0000000000004000   SigIgn ...1001  SigCgt ...00440

`ShdPnd` bit 15 is SIGTERM pending on the shared queue of both. `derna` ignores
only SIGHUP and SIGPIPE and catches only 7, 11 and 33, so SIGTERM has its
default action there: the cell would have died the moment the screen resumed.
The runner itself (880881) had only SIGCHLD pending, so the loop survives.

**So "it was still alive afterwards, so it seems not to exit on TERM" is
explained by the stop, not by any signal handling in `derna`.** A stopped
process does not act on a pending SIGTERM, and SIGTERM does not resume it.

Both PIDs were then killed deliberately, rather than left to die on resume:
the outcome was already fixed, and killing them returned 1.2 GB of resident
memory during the window instead of holding it to no purpose.

**Consequence for the table.** The runner records a non-zero exit as a censored
cell, which in `f0/summary.csv` otherwise means *this protein exceeded the
3600 s budget*. P26436's row means nothing of the kind, and at 188-200 aa it
was nowhere near the budget. The cell is recomputed on the screen's second
pass, and `screen_occurs.py` treats a protein that later succeeded as
successful rather than censored, so the stale row cannot reach a result. This
file is the provenance for that row.

## Suspended time is inside the recorded walls (2026-10-20)

Two proteins were suspended mid-run so that another session's timed arms had a
constant L2 neighbour; cluster {24-27} was free of this screen from 11:45:14 to
12:25:57 and from 13:15:22 to 13:48:01. The runner records elapsed wall, so those
intervals are *inside* the `wall_s` column:

| protein | recorded `wall_s` | suspended | compute |
|---|---|---|---|
| Q9UBM8 (478 aa) | 7371 s | 4418 s | ~2953 s |

Nothing downstream reads these numbers — A′ publishes occurrence counts, which are
exact integers that neither load nor suspension can perturb, and no wall from it
enters a table. The row is here so that a later reader comparing `wall_s` across
proteins does not treat 7371 s as a cost.

**The hazard worth carrying elsewhere: a SIGSTOP counts against the suspended
process's own `timeout`.** `timeout` measures wall clock and does not care that its
child is stopped. Q9UBM8 ran under `timeout -k 30 10800` only because the cap had
been raised from 3600 s at 05:54 for an unrelated reason — the tail's cost exponent.
Under the original cap, a protein whose compute was comfortably inside it would have
been killed by a courtesy to a peer. Any budgeted job suspended for someone else
needs its budget raised by at least the hold, or the hold spends the job's budget.

## A6NKF1: a complete cell with no summary row (found 2026-10-20)

`screen/f0/A6NKF1.jsonl` was written on 2026-10-20 at 04:51 and holds a complete
front — 1719 points, 87 on the hull, 86 gaps, `wall_s` 1908.94 in its own header —
yet `f0/summary.csv` carried no row for it, and neither does the surviving
`f0.run.log`, whose 61 cells run from P15421 to P47974 without it. The runner
truncates that log on each relaunch, so the pass that produced the cell left no
trace of itself.

**The row has been reconstructed from the front record, not from the log**:
`A6NKF1,404,0,1909,1719,87,86`. Every field but the exit code is read off the header
derna wrote; `points` and `gaps` are the header's own, `hull` is a count of the
point records' `on_hull`, and `wall_s` is the header's `wall_s` rounded, which is
derna's own measure of the walk rather than the runner's outer clock. The exit code
is recorded as 0 on the strength of a header and 1719 point records being present at
all: a cell cut by `timeout` leaves no output.

**The harness hazard this exposes.** The runner skips a protein when
`[ -s "$OUT/$acc.jsonl" ]`, and writes the summary row only after the front returns.
A pass interrupted between those two — the window is one `python3` invocation wide —
leaves output that every later pass skips, so the row is never written and the cell
is **silently absent from the summary for good** while appearing complete in the
counts, which are built from the `.jsonl` files directly. The skip condition should
require a summary row as well as output.

That fix is **not applied to `screen_f0_long.sh` while the band is running**, bash
reading a script by offset. `screen_f0.sh` is idle and can take it; the band's runner
takes it at the next relaunch.

**What it did not affect.** The counts are unharmed: `screen_occurs.py` enumerates
`f0/*.jsonl` and never consulted the summary for them, so A6NKF1 was in
`binding.jsonl` throughout and the pool and per-set tallies are as published. What
read the summary was the tail fit behind the band's per-protein cap, and the missing
404 aa cell moved it in the **safe** direction: with A6NKF1 the exponent above 400 aa
is 3.53 rather than 3.79, which puts the band's longest protein at 4.0 h at +1 sd
instead of 4.5 and leaves `PER=28800` at 2.0× rather than 1.8×.
