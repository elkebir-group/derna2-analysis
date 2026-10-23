# E3 pass 1: the short band's motif-exclusion grid, as measured

Six proteins at 265-478 aa by four cumulative motif-exclusion depths. F0 is
free; F1 excludes set I, F2 adds set II, and so on to 25 motifs at F4. No UTRs,
no `--forbid`, no `--structurally-optimal`; every cell carries `stop: null` and a
nucleotide length of exactly 3 x aa, so no cell appends a stop codon.

`complete` is whether the front is whole (gaps == vertices - 1), which is a
different question from whether an external budget cut the run off. Most cells
here are budget-stopped at 4 h and are therefore ANYTIME fronts, not whole ones;
`open` counts the boxes still unresolved. Peak RSS is read from each front's own
header rather than from the `.rss` sidecar beside it, for the reason in
MOTIF-FRONT-SPEEDUP-FINDINGS.md section 8b.

| cell | aa | depth | motifs | points | hull gaps | open | complete | stopped_by | wall s | peak GiB | threads | pin |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P26436_F1 | 265 | I manufacturing (17) | 17 | 627 | 38 | 3 | False | null | 9360 | 1.59 | 1 | `-` |
| P26436_F2 | 265 | II +processing (3) | 20 | 642 | 37 | 4 | False | null | 9908 | 1.59 | 1 | `-` |
| P26436_F3 | 265 | III +repeats (4) | 24 | 189 | 17 | 3 | False | null | 8001 | 1.18 | 1 | `-` |
| P26436_F4 | 265 | IV +initiation (1) | 25 | 189 | 17 | 3 | False | null | 8070 | 1.18 | 1 | `-` |
| Q9HC38_F1 | 313 | I manufacturing (17) | 17 | 58 | 10 | 4 | False | null | 4674 | 1.49 | 1 | `-` |
| Q9HC38_F2 | 313 | II +processing (3) | 20 | 52 | 8 | 4 | False | null | 4397 | 1.49 | 1 | `-` |
| Q9HC38_F3 | 313 | III +repeats (4) | 24 | 9 | 1 | 6 | False | null | 1379 | 1.40 | 1 | `-` |
| Q9HC38_F4 | 313 | IV +initiation (1) | 25 | 63 | 9 | 12 | False | null | 4255 | 5.14 | 4 | `549074e9` |
| P47211_F1 | 349 | I manufacturing (17) | 17 | 1283 | 88 | 0 | True | null | 19345 | 3.24 | 1 | `-` |
| P47211_F2 | 349 | II +processing (3) | 20 | 1283 | 88 | 0 | True | null | 19926 | 3.25 | 1 | `-` |
| P47211_F3 | 349 | III +repeats (4) | 24 | 1252 | 85 | 4 | False | time | 28801 | 3.10 | 1 | `-` |
| P47211_F4 | 349 | IV +initiation (1) | 25 | 1048 | 59 | 2 | False | time | 28800 | 2.47 | 1 | `-` |
| P15086_F1 | 417 | I manufacturing (17) | 17 | 110 | 17 | 9 | False | time | 14654 | 7.25 | 4 | `549074e9` |
| P15086_F2 | 417 | II +processing (3) | 20 | 130 | 19 | 9 | False | time | 14812 | 7.27 | 4 | `549074e9` |
| P15086_F3 | 417 | III +repeats (4) | 24 | 138 | 19 | 8 | False | time | 14606 | 7.47 | 4 | `06605630` |
| P15086_F4 | 417 | IV +initiation (1) | 25 | 123 | 18 | 9 | False | time | 14815 | 7.47 | 4 | `06605630` |
| Q9NWZ3_F1 | 460 | I manufacturing (17) | 17 | 868 | 49 | 8 | False | time | 14401 | 8.78 | 4 | `549074e9` |
| Q9NWZ3_F2 | 460 | II +processing (3) | 20 | 639 | 46 | 8 | False | time | 14551 | 8.88 | 4 | `549074e9` |
| Q9NWZ3_F3 | 460 | III +repeats (4) | 24 | 161 | 19 | 11 | False | null | 14567 | 9.35 | 4 | `06605630` |
| Q9NWZ3_F4 | 460 | IV +initiation (1) | 25 | 103 | 15 | 9 | False | null | 14746 | 9.85 | 4 | `06605630` |
| Q9UBM8_F1 | 478 | I manufacturing (17) | 17 | 537 | 37 | 6 | False | null | 22317 | 10.75 | 4 | `549074e9` |
| Q9UBM8_F2 | 478 | II +processing (3) | 20 | 491 | 34 | 6 | False | null | 22209 | 10.74 | 4 | `549074e9` |
| Q9UBM8_F3 | 478 | III +repeats (4) | 24 | 276 | 24 | 8 | False | null | 29424 | 9.76 | 4 | `549074e9` |
| Q9UBM8_F4 | 478 | IV +initiation (1) | 25 | 202 | 20 | 10 | False | null | 27298 | 10.11 | 4 | `549074e9` |

24 cells, 2 complete fronts, 22 anytime. Pass 1 is closed: no E3 cell is running.

## `stopped_by` and `complete` discriminate only as a pair

Eleven cells read `stopped_by: null` with `complete: false`. That was read here
as proof of an **external kill**, on the reasoning that `stopped_by` records the
budget that cut a run off, so a run with no budget hit and an unfinished front
must have been stopped from outside. **That inference is withdrawn: the field is
not always written when the walk stops itself.** A deadline that fires inside a
gap fold does not set the stop reason, so a walk that times out there and has no
other task waiting ends `stopped_by: null`, `complete: false` — the same reading
as a kill. The defect is live on this grid's pin and on `xg`; the memory budget
records `"memory"` on the same path, which is what made the omission look like a
kill rather than like a field that was never filled in.

So the eleven split by **where the wall sits relative to the cell's own budget**,
and only the first group is a kill:

- Well under the budget: a genuine kill at a round boundary to free cores —
  P26436 and Q9HC38 at 1,379-9,908 s against 28,800 s.
- At or past the budget: the defect, not a kill — Q9NWZ3 F3 and F4, a few
  hundred seconds past their own 14,400 s, previously attributed to the wrapper's
  hard stop. Q9UBM8 at 22,209-29,424 s against 64,800 s is the ambiguous middle
  and is not attributed either way here.

Each file is intact — header plus one line per point, all
parseable — so the fronts are valid anytime fronts; it is the *reason they stopped*
that the header does not carry.

Neither field alone separates the three outcomes. `complete: true` with
`stopped_by: null` is a whole front (P47211 F1 and F2). `complete: false` with
`stopped_by: time` is a budget-stopped front. `complete: false` with
`stopped_by: null` is an external kill, a walk that timed out inside a gap fold
without recording it, or what a **live** front file looks like when its header is
read mid-run, because a front on disk is a checkpoint. Three causes for one
reading, so the pair narrows the outcome and does not fix it. Read the pair, and for a live file check the process before reading
the header at all: the two F4 rows in this table were first recorded from running
checkpoints and both moved afterwards — P15086 F4 from 107 points to 123 and
Q9NWZ3 F4 from 79 to 103.

## What the grid is for, and what it does not yet support

The price of motif-freedom is a comparison ACROSS depths at one protein, so a
depth pair is only comparable when both its cells resolved the same region of the
front. Budget-stopped cells at different open-box counts are not that, which is
why `open` is a column rather than a footnote: a cell that stopped with 13 boxes
open has not priced the same front as its sibling that stopped with 5.

A count over a set also hides substitution. P47211 F1 and F2 agree on points,
hull size and gaps while 122 of 1283 designs differ, so equal columns are not an
unchanged front and the per-depth comparison has to be made point by point.

Generated from the fronts on disk at 2026-10-22 10:45, after the last cell exited.
