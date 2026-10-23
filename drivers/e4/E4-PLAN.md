# E4 — UTR-aware design: four fronts per instance

> Status: **plan, not run.** Written 2026-10-19. Nothing below has been measured
> except where a line names an earlier scope's result and says so.

E4 is the experimental side of the paper's second contribution (UTR-PORD,
`prob:utr-pord`). It answers one question per instance:

**what does a design practice actually deliver once the UTRs are attached?**

Four fronts on the same `(MFE(u5·v·u3), CAI(v))` axes, for one protein and one
UTR pair. The MFE axis is always the **re-scored free MFE of the whole molecule**,
computed by one independent scorer (`mfe_fold`) for every arm, so the four curves
are comparable by construction and no arm is plotted on its own reported number.

| arm | what it is | command |
|---|---|---|
| `A1` unconstrained | the UTRs inside the DP and nothing forbidden: the same molecule and energy model as `A2` without the constraint, so the two differ by the constraint alone | `derna front --protein P --utr5 U5 --utr3 U3 --forbid none`, then each design re-scored under `5-cds,cds-3` |
| `A2` constrained | the UTRs inside the DP, every candidate forming a forbidden pair rejected | `derna front --protein P --utr5 U5 --utr3 U3 --forbid 5-cds,cds-3` |
| `A3` mRNAfold | the comparator, UTRs in its energy model, nothing rejecting CDS–UTR pairs, run at `A2`'s hull λs | `mrnafold_fold_codon_graph` with `utr_5p`/`utr_3p`, λ' = (1−λ)/λ |
| `A4` structurally optimal | the answer the problem asks for; partial under a budget, each point `certified` or `provisional` | `derna front … --structurally-optimal --front-from A2.jsonl --max-secs 43200 --max-rss-gb 48` |

`A4 ⊆ A2` as design sets, by construction: structural optimality is an extra test
applied to the constrained front. `A1` and `A3` are **not** fronts of structurally
optimal designs and are not claimed to be; what they contribute is the headline
number, the share of each arm's points whose own free fold respects the constraint.

## The constraint is the two crossing classes, not three

The paper's problem statement explicitly allows 5′–3′ pairs (the two ends of an mRNA
are intrinsically close and do pair). **E4's headline constraint is therefore
`--forbid 5-cds,cds-3`, passed explicitly on every `A2` and `A4` cell.**

`derna`'s own default has since moved to exactly that pair of classes
([ADR 0044](../../docs/decisions/0044-derna-forbids-only-the-classes-that-reach-the-coding-sequence.md),
2026-10-19, superseding 0042), on the same argument this section makes below and on
0035's measurement that `5-3` left 132 of 132 designs byte-identical. **The flag stays
written out on every cell regardless**, and is now load-bearing for a different reason
than when it was written: the default differs *between binaries*, so a cell whose
command line omits it is not reproducible without also knowing which `derna` commit ran
it. The pin of step 1 records the commit and the header echoes `forbid`, but a command
line that states the constraint cannot be misread by either.

This is not a formality, and it is the reason `A1` and `A2` are different
experiments:

- **Under the three-class constraint, `A2` is `A1` shifted by a constant.** At
  `dangles = 0` (derna2 folds nowhere else) no pair leaves its region, the outermost
  loop scores its branches independently, and the constrained minimum is
  `MFE(u5) + MFE(v) + MFE(u3)`. Scope 37 proved this and measured the per-pair shift
  (hbb −31.40, gapdh −61.30, hbb25 −2.50, gu-poor 0.00 cost units); a smoke test on
  2026-10-19 reproduced it (`MWWKFYCLMW`, −6.1 alone against −13.6 under the three-class constraint). The
  argmax set is then **identical** to the no-UTR front's, so a three-class `A2`
  carries no design `A1` does not.
- **Allowing `5-3` breaks that separability.** A 5′–3′ pair encloses the CDS, so the
  CDS's structures become branches of that loop and pay its branch terms; the CDS's
  optimal substructure is no longer its free optimum. `A2`'s designs can then differ
  from `A1`'s, which is what makes the arm worth running.

So the three-class setting stays in the grid as a **control**, one instance only, and
after ADR 0044 it is reached by writing `--forbid 5-cds,cds-3,5-3` rather than by
omitting the flag: it
must reproduce the constant shift exactly, which verifies the whole UTR plumbing
against a closed-form answer. If it does not, no other cell is readable.

Both settings are accepted by the shipped CLI (`derna` at `e2/bin`; the header echoes
`"forbid":["5-cds","cds-3"]`), and the two-class constraint is the library's
`symmetric_isolation` — exactly `{Utr5Cds, CdsUtr3}`, since `RegionConstraints`
derives its default and forbids nothing. It is **already covered by the repo's own
oracle tests**: `packages/derna-lib/tests/guarded_front.rs` runs its enumeration
oracle, its independent pair classifier and its binding case over three constraint
tables, `symmetric` among them, on 4 plans × 1/4 workers. The standing caveat is
length, not the class set: those cases are 10-residue position sequences, the same
standing every exact claim in this project has.

## Instances — three length bands, matching E3

The axis is **length**, swept at a single UTR pair, in the three bands E3 uses: under
500 aa, 500–1000, and 1000 and above. Ten instances per band was the target; the short
band comes in at 22 because matching E3 means taking its pool whole rather than
sampling it.

**Short band: E3's selection, plus the companion proteins.** E3's short band is the 18
proteins that bind its motif sets I–III, which it takes entire rather than as a sample.
Two of those 18 are already companion proteins here (`P09919`, `RBD`), so adding the
remaining four gives **22 instances at 163–478 aa**:

| | |
|---|---|
| E3's band | `P54852` 163, `P09919` 207, `Q9HAE3` 211, `RBD` 223, `P98173` 230, `P26436` 265, `P30281` 292, `Q9HC38` 313, `Q6PF06` 316, `P78382` 337, `P47211` 349, `Q8TDG2` 376, `P01860` 377, `Q92734` 400, `P15086` 417, `Q9NWZ3` 460, `Q6B0B8` 471, `Q9UBM8` 478 |
| added companions | `P01563` 188, `P01588` 193, `P01241` 217, `eGFP` 239 |

**What matching buys, and what it does not.** It buys reuse: the same instances carry
E3's screens and its measured baseline walls, so a wall here can be read beside one
there. It does **not** import a criterion — E3 selects on whether a protein's designs
spell its motifs, which is E3's constraint and says nothing about whether a 5′ UTR pairs
into a coding sequence. The band's length distribution is therefore a property of E3's
pool, not of anything E4 measures, and it clumps at 200–249 aa where the four
therapeutic natives sit. No instance is selected on any quantity E4 reports.

**Middle and long bands: E3 has not selected them yet, and E4 waits.** The middle band
is ten of the 32 mirror proteins at 500–1000 aa, chosen once E3's 501–700 binding screen
reports; only the lattice census for those 19 exists today, and the candidate list runs
to 962 aa through five fronts inherited from E2. The long band is ten of 15 mirror
proteins plus the four Cas orthologs at 1053–1629 aa, chosen by lattice census — and no
long-band census has been run. Cas9 at 984 aa falls in the **middle** band, as in E3.

So the two upper bands are a dependency on another experiment's schedule rather than on
anything here. What that costs is stated plainly: **the short band is the whole of what
E4 can run until E3 selects**, which is 22 of its 46 cells, and it is also the band where
a complete front is reachable. The supply is on disk either way — 33 candidates at
500–1000 aa and 19 above 1000 — so if waiting becomes the binding constraint, ten per
band spaced evenly by length is a one-line change and needs no screen.

Native length throughout; no truncation rungs. At 1000 aa and above most cells will be
partial, which E3 settled by publishing anytime curves rather than fronts, and E4 does
the same: a partial cell is a result, reported as partial, never as a front.

## UTR pairs — what is actually used in practice

Five pairs. Three are the UTRs of real constructs, one is the classic IVT pair the
paper's intro already uses, and one is a deliberately unstructured minimal pair that
functions as the control.

| key | 5′ | 3′ | why it is in the grid |
|---|---|---|---|
| `hbb` | human β-globin, 50 nt | human β-globin, 134 nt | the paper's own anchor (intro, Fig. 1); the standard IVT pair in kits and preclinical constructs; characterized by Scope 37 |
| `hba` | human α-globin, 37 nt | human α-globin, 111 nt | the other globin pair in clinical use; `HBA1` MANE Select |
| `bnt` | α-globin-derived with an optimized Kozak, 53 nt | the AES + mtRNR1 ("FI") element, 295 nt | the UTRs of an **approved** mRNA vaccine (BNT162b2); the most-deployed pair in existence |
| `mrna1273` | the construct's own 5′ UTR, 57 nt | its α-globin-derived element, 111 nt | the UTRs of a second approved construct, chosen independently of `bnt` |
| `curevac` | HSD17B4-derived, 77 nt | PSMB3-derived with a histone stem-loop, 95 nt | the UTRs of **CV2CoV**, a clinically tested construct, from two genes unrelated to the above — the longest block it shares with any other pinned record is 7 nt, so this is the one pair carrying independent material on both sides |

Sequences and boundaries: [`utrs/PROVENANCE.md`](utrs/PROVENANCE.md), pinned from the
source documents on 2026-10-19.

**The axis is narrower than five pairs suggests, and the measurement says so.** Longest
common substring over the pinned records: `bnt_5p` shares a **32-nt** block with
`hba_5p` (32 of α-globin's 37 nt — BNT162b2's 5′ UTR *is* α-globin plus a Kozak), and
`mrna1273_3p` shares **91 of 111 nt** with `hba_3p`, being the α-globin 3′ UTR with the
two substitutions that remove its in-UTR `ATG`. So three of those four pairs are
α-globin on at least one side, at GC 0.57–0.66. The independent material is HBB (both
sides, GC 0.38–0.44), the 295-nt FI element, mRNA-1273's synthetic 5′ UTR, and
`curevac`. Results are read with that in mind: a pattern shared by `hba`, `bnt` and
`mrna1273` is weaker evidence than three pairs would imply, and `min` was rejected for
this slot because a low-structure control adds no new material.

**No sequence is quoted in this plan**; `utrs/utrs.fasta` holds them, every record
extracted programmatically from its source document and carrying its own length and
SHA-256. A pair whose sequence cannot be sourced to a public document is **dropped, not
reconstructed** — which is what happened to the originally planned fifth pair, a
minimal-UTR construct: the published design (`minRNA`, Mamaghani et al. 2024) has a
4-nt 5′ UTR and **no 3′ UTR at all** by design, so there was no pair to source.

Two boundaries are worth knowing before reading any result. `bnt`'s 5′ UTR span is
**stated by its FDA registration** (`Site Range 1_1-1_53`); `mrna1273` has no such
annotation anywhere in its record, so its 3′ UTR boundary rests on secondary sources
and is pinned at the framing under which the element aligns 1:1 to the native α-globin
3′ UTR. And the `hbb` 3′ UTR pinned in this repo since Scope 19 is its MANE source
**truncated by 2 nt**; E4 runs the 134-nt source sequence, so its `hbb` numbers are not
interchangeable with Scope 37's.

Two decisions about the molecule's edges, fixed now:

- **No poly(A) tail in the headline cells.** The tail is added enzymatically, its
  length is not a design choice, and including 100 nt of A would put a large
  U-complementary partner in the 3′ region — a real effect, but one that would be read
  as a property of the UTR pair. One **cross-check cell** appends `A100` to `hbb`'s 3′
  UTR on one companion-panel instance, reported separately.
- **Stop codon included** (`--stop UAA`), coding to the constraint and outside CAI, as
  the CLI already defines it: the real molecule has one, and it sits exactly at the
  CDS–3′ boundary where the constraint bites.

**The length grid runs at `hbb` alone, and the pair axis is a second panel.** Sweeping
both axes at once would multiply a 30-instance grid by five pairs for a question neither
axis answers. So: every band instance runs at `hbb`, and the six therapeutic natives run
at the other four pairs as well. The panel is a property of the **cell**, not of the
instance, which is what keeps the two from overlapping — a companion protein's `hbb`
cell is run once and serves as the pair panel's `hbb` column. The two panels are
reported separately and never pooled: the length panel carries any claim about length,
the pair panel any claim about which pair.

Grid size today: 22 short-band instances × `hbb` = **22 length cells**, plus 6
companions × 4 further pairs = **24 pair cells**, for **46 cells**. The middle and long
bands add ten instances each at `hbb` when E3 selects them, taking the length panel to
42 and the grid to 66. `A1` is **per protein, not per instance** (one no-UTR front,
re-scored against each pair), so it is paid 22 times rather than 46.

## Run order

Each step's output gates the next; the gate is written here so it cannot drift.

1. **Pin.** One commit, one frozen bin directory, `SHA256SUMS` + `PIN`, as E2 does.
   `derna` must be at or past `0c6f46f4` (the merge of #428), which is where
   `--front-from` and the as-emitted `--out` arrive and where the default constraint
   becomes the two classes; the seeded walk (`8e588452`) is below that. The pin records
   the commit, and since the default moved on 2026-10-19 the pin is the only thing that
   says what a command line without `--forbid` would have meant — one more reason every
   cell writes it. `mrnafold_fold_codon_graph` must be the build whose `utr_5p` /
   `utr_3p` path is exercised — **mRNAfold's UTR mode has not been run in this
   project before**, so step 3 is where it is first trusted.
2. **Two-class exactness** — `just rust-test-name guarded_front`, whose `symmetric`
   table is this constraint: the walk against whole-molecule enumeration, the
   classifier that no emitted point carries a `5-cds` or `cds-3` pair while `5-3`
   pairs may appear, and the constructed flank on which the constraint is known to
   bind. **Gate: all green.** No cell runs before it passes.
3. **mRNAfold selftest with UTRs.** On a 30-residue protein with the `hbb` pair: its
   design translates back to the protein, its UTRs come back unmodified, and its
   reported MFE is recovered as `Score + λ'·Σ ln w` and agrees with `mfe_fold` on the
   emitted molecule. **`EFE` is not the MFE** and is never plotted. Its partition
   function runs after the search and is inside its wall; that is noted, not corrected.
4. **The three-class control.** One companion instance, both settings: the three-class
   `A2` front must equal `A1`'s front shifted by `MFE(u5) + MFE(u3)`, point for point,
   with byte-identical designs. **Gate: exact equality.**
5. **Companion panel, all four arms**, longest-first.
6. **Long panel, all four arms**, longest-first, one cell at a time.
7. **Summarize** — row counts against the plan first, then exit categories,
   completeness, and the verification table below. No number is read before the row
   counts match.

Per cell, in order, so a stopped cell still yields the arms that depend on it:

`A3` takes **the hull λs that actually completed**, and `A4` runs on whatever
constrained front `A2` produced, with its own 12 h budget on top. A cell therefore
reports a quadruple of completeness verdicts, not a single pass/fail.

`A3` runs at **every** one of `A2`'s hull λs in the short band (`a3_max_lambdas: null`),
not the first 24 in bisection order. The comparator is a scalarized solver, so a weight
can reach only a hull vertex, and one weight per hull interval is everything it can be
asked; 24 of them covered 24 of 46 intervals at 163 aa and 24 of 93–117 at 417–478 aa.
Measured cost at one weight: 17 s at 163 aa, 97–163 s at 417–478 aa, so a whole hull is
13 min to ~4.7 h a cell. Cells first run under the 24-weight cap continue from weight 24,
since the order is the same. Whether the long band keeps every weight is set from its
own per-weight cost.

**Superseded at `549074e9`: the walk now takes a budget.** `derna front` there takes
`--max-secs` and `--max-rss-gb` on a plain front ("stop the walk after S seconds and
return the partial front, naming the boxes left open") and `--checkpoint` (rewrite
`--out` after each settled pair of hull points). The driver passes `--max-secs
walk_budget_s --checkpoint` on `A1` and `A2`, with the outer kill `emit_s` beyond it, so
a walk that does not finish yields a front file marked `complete: false` instead of no
file, and `--resume FRONT.jsonl` can later walk only the boxes it left open. That is the
whole of the gain. Tested on `P54852__hbb` (1 thread, core 16, 2026-10-21):

- `A4`'s `--front-from` **refuses** an incomplete front: *"its walk was stopped with 14
  box(es) open, and --front-from takes a whole front"*, exit 1. `A3` reads hull weights
  from the same front. So a cell's `A3` and `A4` still need `A2` to finish.
- A walk stopped before its grid ends holds **no points**. The deadline skips only
  probes and gap folds. The credit, the ends and the weight grid always run, and points
  are settled only after the grid, by the gaps between its vertices. At `--max-secs 20`
  and `--max-secs 200` both stopped walks ran the credit, an end and all 16 grid passes
  (342 s in the run with `--report`), then skipped everything else: 15 hull vertices,
  0 points. The unbudgeted walk ran the same opening (grid 248 s), then 77 probes and
  45 gap folds, for 495 points in 307 s.
- So a binding budget overshoots by **the rest of the grid**, and the outer kill's
  `emit_s` headroom has to cover a whole grid at the cell's length. Up to the deadline
  the walk is the unbudgeted one task for task, so a budget does not change what is
  found, only where it stops. The summary line's "1 sweep fold(s)" on a stopped walk is
  a counter reading (`sweep_s` 0.0), not an extra probe. The paragraphs below describe the build this plan
was written against (`e31c5dbe`) and are kept as the record of why the driver changed.

**`A2` is the one arm with no anytime curve, which this plan originally assumed it had.**
It was written as running with `--stream`, on the model of E2's `replay_partial_hull`.
That flag is `scope35_walk_front`'s, not `derna front`'s, and the distinction is not
cosmetic: a plain `front` takes no `--max-secs` at all — ADR-era `--max-secs` and
`--max-rss-gb` budget the *structurally optimal* search and are refused without it — so
the only limit on `A2` is the driver's outer kill, and a walk killed there leaves **no
file**. Measured while shaking the driver down: a 300-residue protein behind the `hbb`
pair, killed at 30 s, wrote no `--out` at all.

The consequence is recorded rather than engineered around, because at 984–1629 aa it
decides what the long panel can report: **a cell whose `A2` does not finish has no `A3`
and no `A4` either**, and the driver records all three as missing rather than as zero.
A partial `A2` is not available at any price through the shipped command; the only
mechanism that yields one is `scope35_walk_front --stream`, whose output cannot feed
`--front-from` (that flag reads a `derna front` result, and hand-building its header
would be fabricating one). If the long panel comes back mostly empty for this reason,
the fallback is to take **`A3`'s λs alone** from a streamed partial hull, which needs
no format conversion — and to say so, rather than reporting a grid of absences.

`A4` consumes `A2`'s front literally, not conceptually: `--front-from A2.jsonl` reads
it in place of walking again (#428). This is the single biggest change to the grid's
cost — the constrained walk is most of what a budgeted structurally optimal run spends,
and it is paid **once per cell instead of twice**. The flag is refused unless the
file's header names the same molecule down to protein, flanks, stop codon, constraint
and codon table, and says which field differs, so `A2`'s own output is the only thing
that can feed `A4` and a mismatched pairing cannot silently answer a different
question. `--front-out` is the same front written from inside a structurally optimal
run; E4 does not need it, because `A2` is a cell in its own right and already writes
one.

### The build follows the fastest source, re-checked between cells

The pin is the **newest** `origin/xg` source, not a frozen commit: E4 reports what the
shipped method costs, so a merged speedup belongs in its numbers. That makes the pin a
standing question rather than a one-time choice, and it is re-asked at every cell
boundary, before each band, and before any result is read as final:

1. `git fetch origin`, then `git log --merges --oneline <pin>..origin/xg` and
   `git log --oneline <pin>..origin/xg -- packages/`. Open PRs into `xg` are read too
   (`just gh-pr-list`), so that a speedup about to land can be waited for instead of
   switched to mid-band.
2. For each merge, read what it touches and sort it: **default-on** changes on a path an
   arm runs (`derna front`, the constrained walk, the structurally optimal search,
   `check`, the census, `mfe_fold`) can move a row; opt-in flags E4 does not pass, and
   paths no arm reaches, cannot. Only the first kind is a candidate.
3. A candidate is built into a **new** bin directory named for its commit
   (`bin-<sha>/`, with its own `PIN` and `SHA256SUMS`), never over the running one, and
   A/B'd on one finished cell against the current build: `A2`'s and `A4`'s emitted
   fronts compared as point sets and designs, and the wall. **Identical output and a
   lower wall** means switch. The A/B is a **short fixed window**, not a whole arm: both
   builds run the arm on the same input under the same small `--max-secs` (600 s), and
   the faster build is the one that certifies more points in the window, or settles a
   region sooner, with every point both certify agreeing. A whole-cell `A4` takes 9 h
   and more at 163 aa, and four candidate builds were superseded on 2026-10-21 while a
   whole-cell gate waited. **Different output** is a different experiment, not a
   speedup: it is reported, and the build does not move until that is decided.
   Set equality is the test only where both runs report `complete: true`. A run that
   stopped on its budget spends the same budget further along under a faster build, so
   there the test is that every point both runs return carries the same energy and CAI.
   A change that also tightens a check can legitimately drop a point. Each dropped point
   must then show the new check's free MFE below the old one's, which is a false pass
   caught, not a regression.
4. A switch happens at a cell boundary — never under a running cell, whose binary must
   not change — and points the config at the new bin directory **and a new `out_dir`**.
   The driver stamps each row with the pin and the binary digests and counts a row as
   done only under its own stamp, so a switch into the old `out_dir` would both re-run
   every finished cell and overwrite its outputs.
5. Cells finished under the old build are re-run under the new one where that is cheap
   enough that the band reads as one build; otherwise they stay, and every table that
   spans the switch says which build each row ran. The decision, the A/B and the commit
   range examined are recorded in this section either way, including a check that found
   nothing — a pin that was re-examined and kept looks the same afterwards as one nobody
   looked at.

Record of checks (pin `549074e9`, built 2026-10-21 on ada):

| when | xg tip | merges since last check | on an arm's path? | decision |
|---|---|---|---|---|
| 2026-10-21 04:09 | `549074e9` | — | — | pin |
| 2026-10-21 05:13 | `ce45436d` | #443 (`556295d2`, `999b0011`): gap-fold threads in `motif_front.rs`, help text in `front.rs` | no — E4 passes no `--motif` | keep `549074e9` |
| 2026-10-21 06:00 | `add77e48` | #444: certification's free fold is a Zuker fold instead of LinearFold (no 30-nt multiloop leading-run cap), constrained energy from the point's own structure; 6.7× per check at 650 nt | yes — A4 and `score` | candidate; built into `bin-add77e48/`, workers stop claiming; one-cell A4 A/B on `P54852__hbb` from 07:07 (core 12) |
| 2026-10-21 07:20 | `ef5c7337` | #445 (`846889a7`..`bc450b51`): budget read inside gap folds, bounded first probe, energy cut on every table, all in `motif_front.rs` | no — E4 passes no `--motif` | candidate stays `add77e48` |
| 2026-10-21 07:40 | `6f891cb9` | #446: the check's free fold is the sparse Zuker fold (same MFE, 2.9–3.8×) | yes — A4 and `score` | candidate replaces `add77e48`; built into `bin-6f891cb9/`; the `add77e48` A/B run is stopped unfinished at 1 h 35 min and the A/B restarted on `6f891cb9` (08:45, core 12). Stack samples then put all of A4's time in the check's constrained scan |
| 2026-10-21 12:30 | `20c60faf` | #447 (motif walk only; `fb593793` shrinks a front's capacity, output unchanged); #448: the top band refutes a design by the sparse free and constrained Zuker folds before the scan (refuted designs 4.2–6.0×, passing 0.93×; the scan's energy still decides), and reads `--max-secs` before each design instead of once per 1,024-design batch | yes — A4 | candidate replaces `6f891cb9`; built into `bin-20c60faf/`; the `6f891cb9` A/B run is stopped unfinished at 4 h 36 min and the A/B restarted on `20c60faf` (13:19, core 12) |
| 2026-10-21 13:48 | `6d42ff3f` | #449: emission takes a seed's structure from its own check instead of re-folding it with the scan, which ran outside `--max-secs` (a 900 s search had spent another 1,500 s emitting); walk-found points still fold. Only `structure` can differ, and only between co-optimal traces | yes — A4 emission | candidate replaces `20c60faf`; built into `bin-6d42ff3f/`; the `20c60faf` A/B run is stopped after ~30 min and the A/B restarted on `6d42ff3f` (core 12). The A/B now requires every field but `structure` to match, and counts structure differences apart for a re-score |
| 2026-10-21 18:00 | `762d509e` | #450: the structurally optimal search checkpoints (`--checkpoint`) and resumes (`--resume-search`), records its phases in the header, takes `--anchors-from`, searches an uncapped region by `--top-search band\|window`, walks its rounds side by side; emission and `derna check` fast paths; a tie in both axes is broken by the sequence. A shared fold cut (`20602804`, unmerged) is undecided. A first whole-run A/B (one cell, 3 h: 533 points and 2 open regions against the base's 535 and 1, at 1.4 GB more peak) was read as a rejection and is retracted: its two arms ran on different sockets, one saturated, and the cut arm's walked stretches summed 3,330 s against the base's 6,599 s. A same-core sequential A/B decides it | yes — A4 and `score` | built into `bin-762d509e/`. `derna check` on `P54852__hbb`'s A2 front: all 495 rows identical in every field to `6d42ff3f`'s, 115 s against 666 s at one thread; used to count seeds on every A2 front. A4's 600 s window and a checkpoint kill-at-300 s / resume A/B against an uninterrupted run follow on a cell with seeds |

## Budget and host

- **12 h wall per cell**, 48 GB per cell under `systemd-run --user --scope -p
  MemoryMax=48G -p MemorySwapMax=0` (the 52 GB of E2's note, less headroom for the
  re-scoring that follows each cell). A cell over either limit records `partial` with
  its open regions, or `oom`; neither is censored.
- **`--max-secs` bounds the search, not the whole run.** A structurally optimal run has
  three stages — the constrained walk, the search, then emitting the points it found —
  and that budget governs only the middle one. So an outer `timeout` set equal to
  `--max-secs` can kill a cell during emission, after the answer was found and before
  it was written. E4 avoids this by construction rather than by a margin: `A2` is its
  own cell, and `A4` starts from `--front-from`, so the only stage left under the
  budget is the search. The driver's outer limit is `--max-secs` plus an emission
  allowance, and a cell that trips the outer limit rather than its own budget is
  recorded as a distinct exit category, because the two mean different things.
- **A killed `A4` leaves a usable prefix.** Points are written to `--out` as they are
  emitted, so an outer kill loses only what had not been reached. Reading one needs the
  check the guide names: the header's `points` is how many the front holds, so a file
  with fewer point lines than that is a prefix. The driver applies this to every `A4`
  file, not only to cells it believes were killed — a complete file passes the same
  check, and a count that disagrees unnoticed is the one way this arm could report a
  truncated front as a whole one.
- **`--threads 1` for everything timed.** The long cells are memory-bound before they
  are core-bound, and E2 measured the multi-core arms' interference; a threaded arm
  here would buy wall at the cost of the one resource that decides these cells.
- **The tamazirt reservation rules below are retired by the host move**, not
  overridden: E4 runs on another host, so E2's walks, the walk predicate, the E-core
  rescoring note and the P-core wait are no longer constraints on it. They are kept
  because they record *why* no figure here was ever measured under contention, and
  because anything still run on tamazirt — a correctness check, a quick census — is
  still subject to them. **Every number in this section was derived on tamazirt and
  none of it transfers.** The core lists, the memory cap and the 12 h budget are all
  properties of a 61 GB, 16+12-core host under three other sessions; on the new host
  they are re-derived from the first cell's own wall and peak RSS before any band runs.
- **Memory, not cores, is the binding reservation.** E2's three longest walks are
  projected to 53.4 GB on a 61 GB host, so anything of E4's above ~1 GB RSS would make
  the host OOM killer, not the cap, decide the walk's fate. Before launching such a
  step, check `pgrep -af 'scope35_walk_front.*--n (1273|1368|1629)($| )'` and hold while
  it prints. **No clock time is recorded here on purpose.** E2 revised its projected
  windows three times on 2026-10-19 alone, later and then earlier by up to half a day
  each way, because a per-cell cost it has not measured yet sets them. A window written
  into this plan would be wrong by the time anyone read it, and wrong in the direction
  that matters — a launch believing it was between walks. The walks in question are
  `PRT-06` (1368 aa), `PRT-07` (1629 aa, the 53 GB one) and `P0DTC2_p1273`, which is
  what the predicate matches on; for current timing, ask E2's session.
  Steps 2–4 are all well under 1 GB and ran before any walk started.
- **The E-cores are not idle between the walks.** E2's grid runs 10 rescoring folds on
  16–27, ~2 GB in total, so an E4 step sharing those cores runs slower then. That
  affects nothing recorded here — E4's cells are timed on their own reservation, and
  steps 2–4 are correctness checks — but a wall figure read off an E-core step during
  that period is not a measurement.
- **E4 runs on another host** (decided 2026-10-20, after E3 moved for the same
  reason). The wait for E2's `par` phase was the alternative and is no longer the plan.
  What the move costs is one measurement: the first cell on the new host reports its own
  wall and peak RSS, and the cap, the budget and the concurrency are set from it rather
  than carried over.

At 12 h × (46 `A2` + 46 `A4`) the grid is weeks of wall in the worst case. That ceiling
does not move — a budget is a budget — but what the `A4` budget now buys is entirely
search, since `--front-from` means no cell walks the constrained front twice. The
cells expected to *finish* should therefore finish sooner than this plan first assumed,
and the long cells that exhaust the budget will have spent all of it on the stage that
decides whether they finish. The same release also seeds `front`'s sweep from the pass
grid and cuts its folds under the grid's chords, measured at 2.7–5.3× less wall for the
same front; `A2`'s projections here predate that and are now loose upper bounds. That
is accepted on purpose: the cells that matter most are the ones that do not finish,
and their partiality is the measurement. The run order is longest-first inside each
panel so the information arrives before the budget does.

## How it runs

`drive_e4.py`, beside this plan, with a config of the shape
`config.e4.example.json` (validated: 22 short-band instances at `hbb` plus the six
companions at four further pairs = 46 cells, every protein and every pinned UTR record
resolving).

```bash
python3 experiments/paper-benchmark/drive_e4.py plan      --config CONFIG
python3 experiments/paper-benchmark/drive_e4.py run       --config CONFIG [--cell ID__pair]
python3 experiments/paper-benchmark/drive_e4.py score     --config CONFIG
python3 experiments/paper-benchmark/drive_e4.py summarize --config CONFIG
```

**`run` and `score` are separate commands on purpose.** Every fold that is not the
measured search — the census of `A1` and `A3`, the certification of `A2` and `A4` —
happens under `score`, on `score_cores`, after the measurement has ended. No arm's wall
carries the cost of reading it.

Three things the driver does that the plan asks for and would otherwise be left to
whoever runs it:

- **Each arm's rows carry both its exit category and the result's own verdict**, which
  are different facts. A structurally optimal run that reaches `--max-secs` returns what
  it found and exits **0**, so `exit: ok` beside `complete: false` with a `stopped_by`
  of `time` is the normal partial case, not a contradiction. Confirmed on the shakedown
  cell, which came back exactly that way with one open region carrying an `mfe_floor`
  and a `slope`.
- **Every `A4` file is checked for being a prefix** — the header's `points` against the
  point lines actually present — whether or not the cell is believed to have been
  killed. A complete file passes the same check, and a count that disagrees unnoticed is
  the one way this arm could report a truncated front as a whole one.
- **`summarize` reports whether the constraint bound at all**: how many of `A2`'s
  designs are byte-identical to `A1`'s, and the range of the per-point energy shift
  between them. A price against a constraint that forbade nothing is not a price, and
  that column is what says which it was.

## Gate results — run 2026-10-19, before any cell

All three prerequisites passed, on the E-cores (16–27) beside the other sessions'
jobs; these are correctness checks, so contention does not affect them.

**Step 2 — two-class exactness.** `guarded_front` 6/6 green in 149 s, the enumeration
case taking the bulk of it; re-run after merging the 2026-10-19 release and green
again in 151 s, which is what says the default change did not move the constrained
walk itself. The two-class constraint is covered by the repo's oracle
tests as `symmetric`, so this step is a re-run, not new work. The same file's
`the_fully_isolated_front_is_the_coding_front_shifted_by_the_flanks` is the
separability result that made the two-class constraint the headline.

**Step 3 — mRNAfold with UTRs**, its first run in this project. 30-residue protein,
the `hbb` pair, `lambda 1.0`, custom codon table:

| check | result |
|---|---|
| 5′ and 3′ UTR byte-identical in the emitted molecule | yes |
| emitted molecule is the whole transcript, UTRs included | yes, 5′ UTR at offset 0 |
| CDS back-translates to the target protein | yes |
| MFE recovered as `Score + Σ ln w`, against `mfe_fold` on the molecule | −97.000020 against −97.0000, **2 × 10⁻⁵ kcal/mol** |
| its emitted structure against `mfe_fold`'s | byte-identical |
| `EFE` against the MFE | −101.2041 against −97.00, **4.20 kcal/mol apart** |

Two things follow for the driver. `Σ ln w` is computed **from the codon table and the
emitted CDS**, not from the printed `CAI`: the printed value is 4 dp, which alone puts
a 1.4 × 10⁻³ kcal/mol floor on the recovery, while the table-exact sum closes it to
2 × 10⁻⁵. And `EFE` is 4.20 kcal/mol below the MFE on this cell, so plotting it would
have shown mRNAfold beating every arm at a quantity no other arm reports.

**Step 4 — the three-class control**, `hbb`, same protein, against the closed form:

| check | result |
|---|---|
| points, three-class against no-UTR | 45 and 45 |
| energies equal the no-UTR front shifted by `MFE(u5) + MFE(u3)` = −31.40 | **45/45** |
| CAI identical | 45/45 |
| CDS **byte-identical** | 45/45 |

The shift is −31.40 kcal/mol, which is exactly the figure Scope 37 measured for `hbb`.
So the three-class setting carries no design the concatenation arm does not, on the
real UTR pair and through the shipped CLI — the reason E4's headline constraint is the
two classes. It was the CLI default when this control was run, and is not one now;
the control's command line therefore names all three classes explicitly, and the
figures above stand unchanged because the constraint, not the default, produced them.

Re-run on the **134-nt** `hbb` 3′ UTR that E4 actually uses: identical in every column —
45/45 at the shifted front, 45/45 byte-identical, shift −31.40. The terminal `AA` moves
the 3′ UTR's own MFE by **0.00 kcal/mol**, two unpaired adenosines being free at
`dangles = 0`. So the 132/134 difference is inert for *this* quantity, and Scope 37's
shift transfers. It does **not** follow that the structural-optimality verdicts
transfer: those two adenosines are available to a **free** fold, where an `A` can pair
with a `U` in the coding sequence, which is exactly what arms `A1` and `A4` measure. A
verdict comparison back to Scope 37 names its sequence.

**Step 5 — the driver, end to end on a toy cell.** A 20-residue protein behind the
`hbb` pair, all four arms, `run` then `score` then `summarize`, with the pinned
binaries. It is a shakedown of the instrument, not a result about the problem: at this
length the numbers are what the CLI guide's own worked example already reports. What it
establishes is that each arm's plumbing holds and that the columns are legible —
`A1` 0/3, `A2` 0/3, `A3` 0/2 structurally optimal, `A4` partial with 0 points at
`--max-secs 300`; the census's reported shift for this pair is −31.40, Scope 37's
figure, reached through a different code path than step 4's.

Two things it measured that are worth carrying, because both are properties of the
constraint rather than of the toy length:

- **The two-class constraint did not bind on the design here.** All 3 of `A2`'s coding
  sequences are byte-identical to `A1`'s, in order. So on this cell `A2` carries no
  design `A1` does not, and the arm contributes only its share — the situation the
  three-class control was built to detect, now reachable under the two-class setting
  too, and the reason `summarize` carries that column on every cell rather than on the
  control alone.
- **Separability still broke, by a little, without changing the argmax.** The per-point
  shifts are **−31.90, −31.50, −31.90**, not one constant and not the separable
  −31.40: the two-class constrained fold gains up to 0.50 kcal/mol from a 5′–3′ pair,
  which this constraint leaves allowed. So the closed form of step 4 is genuinely a
  three-class statement, and a two-class cell is a different measurement — exactly what
  arms `A1` and `A2` are separate for — but at 20 residues the gain was too small to
  move any design. Whether it stays too small at 984–1629 aa is an open question this
  grid answers, and the `shift_min..shift_max` column is where the answer shows up.

The failure path is checked as well as the success path: a cell forced to time out
records `exit: timeout` for `A1` and `A2`, `missing-a2` for `A3` and `A4`, no fabricated
rows, and `score` and `summarize` both degrade to `-` rather than to zero.

## What each cell records

Per point: `mfe` (re-scored, `mfe_fold`), `cai`, `on_hull`, `sequence`, `structure`,
the structural-optimality verdict and its **gap** (`free_check`: constrained energy
minus free MFE, positive when the molecule prefers a fold the constraint forbids),
and for `A4` the `certified` / `provisional` label.

Per cell: arm, instance ID, UTR key, constraint classes, points, hull vertices, gaps solved and
gaps left open, wall, peak RSS, `VmHWM`, exit category, completeness verdict, and the
open regions of a partial answer.

## The numbers the paper takes from it

1. **The failure rate of concatenation** — the share of `A1` points that are
   structurally optimal, per UTR pair. This is the generalization of the intro's
   claim (all 145 designs form UTR structures) and of Scope 37's 527 of 1,527, now on
   pre-registered practical UTR pairs rather than development contexts.
2. **Constraining the DP is not enough** — the same share for `A2`. Scope 37's census
   already shows it is far below 1 (hbb certified 186 of 258 points at n = 100 and
   none at 200–500); E4 reports it per practical pair at native length.
3. **The comparator** — the same share for `A3`, plus its MFE/CAI distance from `A2`
   at matched λ. mRNAfold puts the UTRs in its energy model but has nothing that
   rejects a CDS–UTR pair, which is precisely the gap the problem statement names.
4. **The price of structural optimality** — the MFE distance from `A2`'s front to
   `A4`'s at matched CAI, and the CAI distance at matched MFE. Under a partial `A4`
   this is a **bound**, signed and labelled as one: the open regions could only hold
   points that make the price smaller. Every open region has a floor to report it
   with, because `A2`'s front is kept: a missing point at CAI `A` cannot cost less than
   the constrained front's least energy at CAI `A` or above, so a region the run left
   without its own `mfe_floor` is still bounded. A partial `A4` therefore yields an
   interval for the price rather than one-sided silence.
5. **The isolation gap distribution** — the per-design gap of arm `A1`, which is what
   "the molecule does not keep the structure you designed" costs in kcal/mol.

A per-pair table and one four-curve figure per selected instance; the figure's
instances are chosen **before** the fronts are read, by the same rule as the panels.

## Verification

| check | what it establishes |
|---|---|
| every arm's MFE re-scored by `mfe_fold` on the emitted molecule | the four curves share one scale; no arm is plotted on its own report |
| back-translation of every design to the target protein | no arm silently changed the protein |
| UTRs byte-identical in every emitted molecule | no arm modified a fixed flank (mRNAfold in particular) |
| classifier over emitted structures: no `5-cds` or `cds-3` pair in any `A2`/`A4` point | the constraint bound; built from region bounds alone, sharing no code with the DP |
| `A4`'s design set ⊆ `A2`'s | the nesting the plan claims |
| every `A4` point's free fold respects the constraint, re-checked independently | structural optimality is a property of the molecule, not of the walk |
| three-class control equals the shifted no-UTR front | the UTR plumbing, against a closed-form answer |
| `A3`'s MFE recovered as `Score + λ'·Σ ln w` | mRNAfold read correctly, `EFE` not mistaken for the MFE |
| ε = 2⁻²³ at both corners, and each exact arm unchanged at ε/2 | the corners are the lexicographic optima (E2's `epscheck`) |

## Risks, and what each one would cost

- **No long cell finishes any arm.** Then the long panel reports completeness and
  cost only, the price of structural optimality comes from the companion panel, and
  the paper's practical-runtime claim is a frontier length, as E3's ladder reports one.
  The experiment still answers questions 1–3 and 5, which need no complete front.
- **The two-class constraint turns out not to be exact at the lengths E4 runs.** Step 2
  establishes it at 10 residues, which is where the oracle can enumerate; nothing
  establishes it at the 478 residues the grid now reaches, let alone at the long band's. This is the standing caveat on every exact claim here, and
  the mitigation is the same: the classifier runs on **every** emitted point at every
  length, so a point that crossed a boundary would be caught in the cell that emitted
  it rather than inferred from a test. If the walk were found inexact, E4 falls back to
  the three-class constraint — and then `A2` is `A1` shifted, the arm is a verification,
  and the paper's problem statement and its experiment have diverged, which is itself
  the finding to report.
- **mRNAfold's UTR mode misbehaves** (modifies a flank, or its λ mapping does not hold
  with flanks present). Step 3 catches it; `A3` is then dropped with its evidence,
  and the paper states that no comparator takes fixed UTRs correctly.
- **A UTR pair forbids nothing** — a pair with no complementarity to any design on the
  front. Then its `A1` failure rate is near zero and the price is zero. Reported as
  non-binding, which is a real property of that pair, not a null result; `min` is in
  the grid precisely to make the non-binding end of the axis visible.
- **`A4`'s provisional points are read as certified.** Every table carries the two
  counts separately, and a price taken from a partial `A4` is labelled a bound.

## Dependencies

None in the code, and less than when this plan was written. `derna front
--structurally-optimal --max-secs --max-rss-gb` and `--forbid` all ship, `--stream`
for partial cells exists on E2's branch (`scope35_walk_front --stream`,
`run_glpe78 --stream`), and the merge of 2026-10-19 added the two things this plan had
listed as its own work: `--front-from`, so `A4` reuses `A2`'s front instead of the
driver re-walking it, and an `--out` written as points are emitted, so a killed `A4`
needs no separate partial-recording mechanism. The `utrs.fasta` pin and the two-class
exactness check of step 2 are both done.

**The middle and long bands wait on E3**, which has selected neither: its 501–700
binding screen is still running and no long-band lattice census has been run. That is a
dependency on another experiment's schedule, and it is the reason only 22 of 46 cells
are runnable today. The alternative, ten per band spaced evenly by length from the 33
and 19 candidates already on disk, needs no screen and is a one-line change if waiting
becomes the binding constraint.

**The driver is written**: `drive_e4.py` beside this plan, with `measure.py` copied
unchanged from E2's branch at the same path, so the two copies merge trivially. It has
run all four arms end to end on a toy cell and on a cell forced to fail (step 5).
Nothing in the code is now outstanding. What E4 waits on is a host of its own and E3's
two remaining band selections.

One thing it needed that the shipped command does not provide, and which is recorded
above rather than built around: a partial `A2`. `derna front` has no budget for a plain
front and writes nothing when killed.

## Run decisions on ada, 2026-10-21

**A4 starts only where it has seeds.** On `6d42ff3f`, `P54852__hbb`'s A4 certified its
495 constrained points in ~70 s (1 thread) and found **0 seeds**: no point of `A2`'s front
is structurally optimal. With no seed the whole range is one uncapped region, and the
search can only enumerate designs in adaptation order from the top band and check each;
530 s of band found 0 points, and at 163 aa a band slice can hold astronomically many
designs, so a longer budget is not expected to change that. That is an algorithmic
question (what bounds the band when nothing caps it), not a speed-up, and it is open. So
each cell's seeds are counted first, by certifying its `A2` front, and `A4` runs only on
cells with at least one seed; zero-seed cells are held, with their certification recorded,
until that question is settled. `549074e9`'s whole-cell A4 on `P54852__hbb` (10 h, still in
certification under LinearFold) was stopped for this.

**`A4` has a 6 h budget** (`budget_s: 21600`, set 2026-10-21), half the 12 h first planned; the walks and the comparator keep 12 h (`walk_budget_s`, `a3_budget_s`, `mrnafold_budget_s`). With the 4 h emission allowance the outer limit on an `A4` is 10 h.

**`A4` runs with 4 threads, `A1`–`A3` with one** (`a4_threads: 4`). `A4` checks designs on
every thread it has and runs few cells for hours; certification of `P54852__hbb`'s front
took 70.3 s on 1 thread and 26.1 s on 3 (2.7×, 283 % CPU). The walks stay single-threaded
while more cells are queued than cores, since 2.05× at 8 threads is less throughput per
core than one cell per core; the queue's tail takes the spare cores as threads.

**`A1` designs with the UTRs** (set 2026-10-21 ~15:20): `--utr5 --utr3 --forbid none`, where it had designed the coding sequence alone and pasted the UTRs on. Earlier `A1` rows lack `flanks_designed` and count as missing, so every cell's `A1` is re-run; designing with the flanks also makes `A1` a per-pair arm, so nothing is shared across a protein's pairs.

**`A4` runs on its own build, `d1aa9eb8`, on the seeded cells** (set 2026-10-21 ~22:30).
Seeds by `derna check` (`762d509e`, rows identical to `6d42ff3f`'s check) over every `A2`
front: **17 of 46 cells have seeds, 7–672 each; 29 have none** — every `hbb` cell and
`eGFP`/`hba`, `P01563`/`hba`, `eGFP`/`mrna1273`, `RBD`/`mrna1273`, `P01588`/`mrna1273`,
`P01563`/`curevac`, `P01563`/`mrna1273`. Whether a cell
has seeds follows its UTR pair more than its protein: `eGFP` has 0 under `hbb` and 211 under
`bnt`. On `eGFP__bnt` (2 threads, `--max-secs 600`) builds through `73084c00` walked all 52
windows between seeds serially at a ~39 s setup floor each, ran ~2,100 s and found nothing
beyond the seeds; `24f5dbfb` (#455: the uncapped region first, a wave whose share is below
its setup deferred) stopped at 601 s with **228 certified points, 17 beyond the seeds**, the
211 common ones identical in every field, and `d1aa9eb8` (#456) returned the same 228 at
600 s. `d1aa9eb8` reproduces `A2` byte-identically (`eGFP__bnt`, 586 rows, 154 s against
171 s), so `A1`–`A3` keep their `6d42ff3f` rows and only `A4` moves: the config's
`arm_builds` gives `A4` its build and its rows that build's stamp. Measurements in
`/scratch/xinyugu4/paper-e4/ab-short/egfp_bnt/FINDINGS.md`. Seedless cells stay held.
From 2026-10-22 ~01:15 `A4` moves on to `f2112884` (#457: the check's fold cut built once
and shared) for the 15 cells not yet started: on `eGFP__bnt` at 600 s it returned 225
certified points against `d1aa9eb8`'s 228, 223 identical, and none of the 7 others
dominates a point of the other front — the budget completed different windows — while each
window both walked ran 1.3-2.5x faster. `P09919__bnt` and `P09919__mrna1273`, started on
`d1aa9eb8`, finish there; under the `f2112884` stamp their rows count as missing until that
is decided.

**`A2` (and `A1`) first, for every cell, before any more `A3` or `A4`.** The four
`549074e9` band cells' `A4` runs were stopped at 5.9–6.8 h with nothing written (their rows
read `crash`, an operator stop), keeping their `A1`, `A2` and first 24 `A3` weights. `A1`
and `A2` for the other 41 cells run on `6d42ff3f` (`--arms A1,A2`), `hbb` cells first.

## Where this stands, and the host move

Everything above is written, checked and committed, and nothing of E4's has ever run as a
measurement. The move therefore carries inputs, code and a plan, and there is no result to
port.

**What the branch carries.** `drive_e4.py`, `measure.py`, `config.e4.example.json`, the 22
instance fastas under `inputs/` with their `PIN` and `SHA256SUMS`, the pinned UTR records
under `utrs/`, and the codon table under `packages/derna-codon-dp/data/`. The fastas are in
git because of where they came from rather than for convenience: 16 of the 22 were read out
of the C++ mirror under `docs/reference/_cpp/`, which is gitignored and restored by `just
sync-cpp` — a re-clone of an external repository, which need not return the revision this
panel was read at. Re-deriving them on another host would give the same instance ids over
sequences nobody had checked. `PIN` names the mirror commit and the E3 commit the six
companions were copied from, and `sha256sum -c SHA256SUMS` in `inputs/` is the check that
the copies are the sequences these rows are about.

**What the branch does not carry: the binaries.** The release binaries are copied out of
the worktree's `target/release` into `/home/melkebir/projects/derna2/e4/bin`, beside a `PIN`
naming the commit they were built at and a `SHA256SUMS` over all four. The cargo target
directory is shared between every worktree on a host, so a rebuild in any of them replaces
the file the config points at; a copy is what makes the binary an input rather than a side
effect of whoever built last. That reasoning is why another host repeats the freeze instead
of pointing a config at `target/release`: rebuild at the commit `PIN` names, copy the four
out, write a fresh `SHA256SUMS`. Three of the four are this repository's; `mfe_fold`'s
digest differs from E2's copy, as it must, being a different build of the same source.

**`A3`'s comparator is not built from this repository at all, and nothing in git records
how it was made.** `mrnafold_fold_codon_graph` is mRNAfold's `fold_codon_graph` target,
built out of a clone of `https://github.com/maxhwardg/mRNAfold` at
`d5544e4f86efb58b5307224e1fad6744ce4d387c` that sits beside E2's data rather than inside
any branch, so it travels no better than the binaries do. It is a CMake out-of-source
**Release** build (`-O3 -DNDEBUG`) of `cpp_src/`, compiled through a one-line wrapper that
is the whole of why the arm is serial: `g++ -D_GLIBCXX_USE_TBB_PAR_BACKEND=0`. Without
that flag libstdc++ routes parallel algorithms to TBB, and a binary built the obvious way
is not the single-threaded comparator every `A3` row claims; E2's `bin` keeps the threaded
build separately, as `mrnafold_fold_codon_graph_tbb`. The tamazirt copy digests to
`2003f475…` and matches E2's record, **and that match is a same-host check a rebuild
elsewhere will not reproduce**: a different g++ emits a different binary from the same
source. So on another host the identity that carries is the source commit and that
compiler flag, not the digest; record the new digest beside them, and say in the row that
`A3`'s comparator is a different build of the same revision. The repository ships a
`Dockerfile` if a build independent of the host's toolchain is wanted instead.

All four binaries were exercised after the copy here, and should be again there.

**The grid.** `config.e4.example.json` resolves against the frozen directory and against
every protein and every pinned UTR record: `plan` reports **46 cells, 4 arms, 138 timed
arms owed**, which is the three timed arms over the full grid with nothing yet done. That
line is the cheapest check that an input list or a path has not moved, because a dropped
fasta or a stale mirror path shows up there and nowhere else.

**The first cell is the shortest one, not a band.** `run` takes `--cell ID__pair`, and the
grid's own order is longest-first, which is right for a sweep and wrong for a first cell.
`P54852__hbb`, 163 aa, is the cheapest cell in the grid, and its wall and `maxrss_kb` are
what the cap, the 12 h budget and the concurrency are re-derived from — every one of those
numbers in "Budget and host" above is a property of tamazirt and none of it transfers.
`score` and `summarize` that cell too, so all four arms are exercised end to end before a
band is launched. Then the short band at `hbb`, longest-first.

**Two host properties to establish before relying on them.** `systemd-run --user --scope -p
MemoryMax=…` needs a systemd user session with delegation, and the split between cores
reserved for a measurement and cores for re-scoring — `score_cores` against the timed arm's
affinity — is a property of the host's topology rather than of this plan.

### The launch gate was never a reservation

The predicate `pgrep -af 'scope35_walk_front.*--n (1273|1368|1629)($| )'` was written as a
hold-and-resume: wait while E2's walk prints, launch when it does not. That is the right
check before a step of a few hundred megabytes, and it does not reach a cell. A cell is 12 h
under a 48 GB cap, and 48 GB cannot be fitted beside the 53.4 GB E2 reserves on a 61 GB
host — not between two walks, not at any clock time, because the two reservations exceed the
host whatever the gap. The predicate going silent means a walk is not running this second;
it does not mean a 12-hour window is open. That arithmetic is what moved E3's long band off
tamazirt on 2026-10-20 and E4 with it, so the gate binds nothing E4 now does. It is kept
because it records why no figure in this plan was ever measured under contention, and
because a correctness check or a census still run on tamazirt is still subject to it.

Its second half transfers to any host. A predicate going silent is not the host becoming
free: with E2's walks stopped, three other sessions' jobs held about 28 GB between them and
one of them sat on a P-core, and a 48 GB cell cannot be launched against 25 GB available no
matter what a predicate says. So the reading to take before a cell anywhere is `free -g`,
and `taskset -pc <pid>` on whatever is large.

### The cap fit is owed by the long band, not by the grid as it stands

`mem_cap_gb` is 48, derived as E2's 52 GB note less re-scoring headroom, and that derivation
never used E4's own molecule lengths. Both halves of that have now moved. The cap is one of
the numbers the first cell on the new host replaces, and the lengths that prompted the
question are no longer in the grid: the 5038–5238 nt cells were the two longest Cas
orthologs of the earlier eleven-instance panel, which the regrid by length dropped until
E3 selects a long band, while E4's longest molecule is now about 1.6 kb — 478 aa behind `hbb`'s 184 nt of flank — against the 4890
nt of coding sequence that carries E2's 53.4 GB projection. Nothing in the short band or the
pair panel is near any plausible cap.

The obligation moves to the long band rather than disappearing. A fold's tables grow with
the square of the length and `A2` is all-or-nothing, so a long cell above the cap spends its
whole budget and returns nothing — no partial front, no anytime curve. When E3's selection
arrives, measure peak RSS on three cells spanning the band, fit, and read the cap off the
fit before launching the longest; if the fit lands above what the host can reserve, the
choice is a raised cap or those cells reported as out of reach. Do not discover it by
launching the longest cell.

**Not verified:** `just health-check` could not run in this worktree — the CodeScene CLI is
absent — so `drive_e4.py` and `measure.py` have not been scored against the boy-scout rule.
`ruff` and `pyright` are likewise not on PATH here; the file compiles and carries no editor
diagnostics, and that is the whole of what was checked. No PR is open for
`chore/paper-e4-plan`.

## A4's first cells, 2026-10-22

**The first two A4 rows.** Both on `d1aa9eb8`, 4 threads, the 6 h budget binding to within
29 s:

| cell | wall s | points (certified) | seeds | beyond the seeds | open regions | peak RSS |
|---|---|---|---|---|---|---|
| `P09919__bnt` | 21,629 | 693 (693) | 672 | 21 | 2 | 35.6 GB |
| `P09919__mrna1273` | 21,601 | 674 (674) | 666 | 8 | 1 | 18.6 GB |

Neither front is complete, and what each has left is adjacent to its most adapted end —
`P09919__bnt`'s two open regions span log CAI (−0.83, −0.82] and (−0.73, −0.37], against a
front reaching −0.37. That is the uncapped band, the region a budget does not finish, and
the same place `eGFP__bnt`'s band spent its whole allowance without a find. So a 6 h A4 row
is to be read as a certified front over the adapted range it reached, with the top of the
range open, rather than as the whole structurally optimal front.

**`A4` checkpoints its search** (`--checkpoint`, set 2026-10-22 ~04:00). A4 runs six hours a
cell and wrote `a4.jsonl` only at the end, so a run in progress showed nothing and a killed
one lost its hours. Each check, find and finished band slice now appends to
`a4.jsonl.search.jsonl` as it settles, which `--resume-search` continues from; the A/B on
`eGFP__bnt` returned a resumed front identical in every field to an uninterrupted one, and
the resume also saves the certification phase (61 s → 0.4 s). The two rows above predate the
flag. `P09919__curevac` was restarted 8.5 min in to take it, before any row was written.

**The A4 cap is 120 GB, and it is a host-sharing bound rather than a fit.** `mem_cap_gb`
was null for A4, and `P01241__mrna1273` reached a peak of **248 GB** resident inside its
first two hours, settling back to ~80 GB — against the 35.6 GB of the largest finished cell
and the 1–6 GB of the short A/B runs. One such peak fits on a 503 GB host; three A4 workers
peaking together does not, and the host carries another experiment's grid. 120 GB per run
bounds three concurrent workers at 360 GB with ~35 GB held elsewhere. A capped run returns
its partial front rather than racing the host, which is the same anytime reading the budget
already gives. The peak is not yet explained: it is far above what a 207-residue molecule's
fold tables account for, so the cells queued behind it are to be watched rather than trusted
to stay under the cap.

### A zero-seed result is a statement about one stop codon

An external report, on `P54852__hbb` with everything else identical and each arm computing
its own constrained front: `--stop UAA` gives 0 seeds and 0 structurally optimal points,
`UAG` the same, **`UGA` gives 303 seeds and 311 points, all certified**, and no stop at all
gives 1 and 1. An inside-outside pass attributes it: about 60 % of the near-optimal pairs
crossing a flank boundary anchor at one of the two coding positions with no synonymous
choice, the initiator and the stop, and only the stop is a choice at all.

E4 runs `--stop UAA` on every arm, and every front header records it. So the 29 cells with
no seeds are **"none at UAA"**, not "no structurally optimal design exists", and every A4
row states its stop. The zero is a property of a constraint the grid imposes, and a reader
who takes it for a property of the molecule would conclude the opposite of what one codon
substitution shows. Whether the grid should instead run the front over all three stops —
an opt-in `--free-stop` exists, each point naming its own stop, the budget split three ways
and certification only where every arm ran whole — is a question about what the paper claims,
and is open.

**A second cell refutes the rescue, and the wording is what survives.** `P01563__hbb`,
3,600 s on the search, each arm computing its own constrained front: `UAA` 506 constrained
points, **`UGA` 510 — and 0 seeds in both**, with `unchecked: 0` on each arm. Every point of
the constrained front was checked, so the seed set is certified empty at `UGA` rather than
merely unfound, and that verdict comes from a stage that finished in 62 s. So `UGA` does not
rescue a zero-seed cell; one cell was not a mechanism. What stands is only the wording: a
0-point row says **none at UAA**, and claims nothing about the molecule.

The pair is sharper than that for whoever writes up the open problem. `P54852__hbb` and
`P01563__hbb` have **md5-identical `utr5` and `utr3`** and differ only in protein, 163 aa
against 188 aa. At a suffix beginning `UGA` the same flank pair gives 311 tight points on
one protein and a certified-empty seed set on the other — so emptiness is not decidable
from the flank pair alone, which is what the conjecture on that problem assumes. At `UAA`
the ten `hbb` instances were consistent with it; at `UGA` it is false.

### What the peak is made of, read off the fronts already in hand

The `stretches[]` the SO header already records carry `cai_log_lo`, `cai_log_hi`, `walked`,
`complete`, `points`, `designs` and `wall_s`, so the two finished cells can be decomposed
without a new instrument. Both ran at `threads: 4` by their own headers — read from there
rather than from the config, since a peak normalised by a factor measured on one arm is how
E3's § 8a decomposition went wrong — so lanes is 2 in both and the peaks compare directly.

| | `P09919__bnt` | `P09919__mrna1273` |
|---|---|---|
| walked wall | 22,494 s | 22,202 s |
| the uncapped band | 50.6 %, 9 visits, 178,890 designs, 3 points | 93.3 %, 10 visits, 645,469 designs, **0 points** |
| widest capped window walked | 0.770, one visit, 9,044 s, 8 points, incomplete | 0.136 |
| capped windows walked | 15 of 674 | 12 of 665 |
| peak RSS | 35.6 GB | 18.6 GB |

Two readings, and they point the same way. The band is where the **time** goes — on
`mrna1273` it took 93 % of the search and returned nothing for 645,469 designs — but it is
not where the **memory** goes: that cell has the lower peak of the two. What orders the
peaks is the widest capped window the run actually walked, 0.770 against 0.136 for 35.6 GB
against 18.6 GB, which is the B&B node expansion the walk carries and whose depth the
window's width buys. `P01241__mrna1273`, with 145 unseeded points against 49 and 50, is
consistent at 248 GB, and it is the same length as both.

The sizing rule that follows is a **dynamic** quantity, not a static one: the risk is the
widest window a run *walks*, which the scheduler chooses within the budget, not the widest
window the cell *has*. Both cells above hold a 3.598-wide window neither ever reached. So a
cell cannot be declared safe from its A2 front alone, and a whole-run RSS cap is the wrong
control — it fires once, at 1 s resolution, after the allocation that caused it. The right
control is per-window, which is re-armed by construction on the next window, and the field
that would confirm all of this is peak RSS per stretch, which the header does not yet carry.

### The memory cap is usable, and A4 runs at 120 GiB

Four arms on `eGFP__bnt`, f2112884, 600 s, 2 threads, sequential on two cores
(`ab-short/egfp_bnt/FINDINGS.md` § 9). `so_front.rs` computes `limit = gb * 1024³`,
so the flag is GiB; peaks below are `ru_maxrss`, in GiB throughout — the same field
quoted earlier in this plan as 18.6 and 17.7 GB for one cell, which was a decimal
figure and a binary one for one measurement.

| arm | cap | wall | peak GiB | points | stopped by |
|---|---|---|---|---|---|
| nocap | — | 10:17.30 | 5.35 | 225 | time |
| cap120 | 120 GiB | 10:16.18 | 5.36 | 244 | time |
| cap4 | 4 GiB | 2:51.89 | 4.61 | 211 | **memory** |
| nocap2 | — | 10:07.21 | 5.36 | 244 | time |

- The flag's **presence** costs −0.2 % wall and +0.3 % peak, and `nocap2` returns a
  front identical to `cap120`'s in every field, so it is neutral in what comes back
  as well. E3 measured +72 % from the same flag on the motif walk; that does not
  transfer to this code path. The 19-point spread between the first arm and the
  other three is a first-arm effect, which is also the noise floor for any
  point-count reading on this cell: **8.4 %**.
- The cap **bounds**: `cap4` stopped at 4.61 GiB against 4 GiB. The overshoot is
  absolute rather than fractional — about (allocation rate × the 1 s poll) — so
  0.613 GiB measures a ~0.6 GiB/s ramp and must not be carried to another cap as a
  percentage. Peak minus cap, in GiB, is the figure that generalises.
- A capped front is a **valid prefix**: 211 common with `nocap`, 0 points only in
  `cap4`, no certified point dominated.
- The watcher is one-shot and polls at 1 s, and `so_front.rs:189` documents that it
  never stops the constrained walk — only the structurally optimal front's walks.
  Every E4 A4 reads A2's front with `--front-from`, so no E4 row has an unguarded
  walk phase; a row that walked its own front would have 500–1,000 s of one.

So A4 carries **`rss_cap_gb: 120`**, a key of its own, and no cell is held back for
memory. It is not `mem_cap_gb`: that key wraps the whole command in
`systemd-run --user --scope -p MemoryMax=…`, a cgroup cap the kernel enforces, and
`systemd-run --user` has no session bus on a host reached over ssh without one.
Setting it emptied the queue — twelve cells claimed and failed in six seconds, each
writing a row with `exit: "crash"` and `wall_s: 0` on `Failed to create bus
connection`, and `done_keys` does not read `exit`, so all twelve then counted as
done. The rows were dropped (`results.jsonl.bak-10-48` holds them), the claims
released, and the cells re-queued. `capped()` now refuses to start when `mem_cap_gb`
is set on a host where `systemd-run --user` fails, rather than turning a queue into
crash rows.

The lesson is narrower than "test the flag": the flag **was** tested, in four arms,
by calling the binary directly under `taskset`. What went untested was the *driver's*
path to it, and `mem_cap_gb` turned out to mean a different mechanism there than the
one measured. A config key named for a quantity can select a mechanism, and a cap
validated on the binary says nothing about the harness that sets it. Two
capped runs plus a few GiB of overshoot each fit beside the ~35 GiB another
experiment holds on a 503 GB host; three do not, which is one reason A4 runs two
workers. A cell that fires the cap returns its seeds and its finds with
`stopped_by: "memory"` in the header, which is a reportable row rather than a hole.

**Read `complete` and `stopped_by` as a pair.** `complete: false` with
`stopped_by: null` is not a run that stopped itself: it is an external kill, or a
`--checkpoint` file read while its run was still going, which A4 now produces for
every cell. Gate "ran whole" on `complete: true`, and "stopped by its budget" on
`stopped_by: "time"` with `wall_s` near `max_secs`. Every A4 row and A/B arm recorded
so far carries `stopped_by` of `"time"` or `"memory"`, so none of them is a kill.

### What drives A4's peak, and what does not

| cell | build | aa | nt | widest window walked | band | peak GiB |
|---|---|---|---|---|---|---|
| `P09919__bnt` | d1aa9eb8 | 207 | 972 | 0.770 | 50.6 % of walked wall | 33.95 |
| `P09919__mrna1273` | d1aa9eb8 | 207 | 792 | 0.136 | 93.3 %, 645,469 designs, 0 points | 17.71 |
| `P09919__curevac` | f2112884 | 207 | 796 | 0.112 | 99 %, 625,165 designs | 3.47 |
| `P01241__mrna1273` | f2112884 | 217 | 822 | — | — | 247.9 (VmHWM) |

- **The band is not the allocator.** One cell spent 99 % of six hours in it, for
  625,165 designs, and peaked at 3.47 GiB. The band is where the *time* goes and it
  returns almost nothing: 0 points for 645,469 designs on the cell above it.
- **Width orders the peaks within a build and across cells**, which is what the
  reservation of the widest cells rested on before the cap made reservations
  unnecessary. It does **not** order them across builds: `P09919__curevac` and
  `P09919__mrna1273` are the same protein at 796 against 792 nt, the same 3.598
  widest candidate window, near-identical band work and 0.112 against 0.136 widest
  walked — and the newer build uses **5.1× less memory**. So f2112884 is not a memory
  regression in this regime, and the 248 GiB peak is a property of its cell rather
  than of its build.
- The **below-front** region — from the least adaptation the codon table allows up to
  the front's lowest point, about 0.63 units per residue — does not order the peaks
  either: it is widest (150.8) on the cells with the lowest peaks. On d1aa9eb8 it is
  never walked at all (`cai_log_lo: null`, `walked: false` on both finished cells),
  and f2112884 did not walk it on `P09919__curevac` either.

### The 248 GiB cell: one 11.6-wide window, walked ten times, for nothing

`P01241__mrna1273` finished on f2112884 at 10:56: **572 points, all certified**, 2
regions open, `stopped_by: "time"`, `wall_s` 21,718 against a 21,600 s budget, peak
247.9 GiB (VmHWM). 59 walked visits over 35 regions, 42,716 s of walk wall across
lanes. Decomposed:

| region | visits | wall | share | designs | points | complete |
|---|---|---|---|---|---|---|
| one window, width **11.635** | 10 | 20,420 s | 47.8 % | 0 | **0** | no |
| the band | 10 | 19,892 s | 46.6 % | 289,452 | **0** | no |
| 33 windows, 0.087–1.025 wide | 39 | 2,404 s | 5.6 % | 0 | **65** | 31 of 33 |

- **The widest window walked is 11.635, which is exactly the widest candidate window
  computed from the seed census before the run** — the largest gap between
  consecutive seeds in log CAI. So it is one of the cell's own windows, not a region
  the build invented: the below-front region was **not walked at all**, and no walked
  window spans anything the seed set does not define. Both mechanisms proposed for
  the earlier cross-build memory observation on another cell are therefore absent here, and
  width predicts the peak from static data, on this build, to three decimals.
- **94.4 % of the search bought nothing.** The wide window took 20,420 s over ten
  visits and never completed or found a point; the band took 19,892 s and 289,452
  designs for none. Every one of the 65 points beyond the seeds came from the 5.6 %
  spent on narrow windows, which complete in about 60 s a visit.

That reframes a per-window width limit as a **throughput** change and not only a
memory control: on this cell, refusing any window wider than ~1 unit would have freed
roughly half the budget, and the only region that converts time into points is the
narrow one. It also predicts that the queued cells with wide candidate windows —
`P01241__bnt` at 35.999, `RBD__hba` and `RBD__curevac` at 26.576, `P01588__hba` at
18.417 — will spend most of a 6 h budget on a window that returns nothing, which is a
reason to want the limit well before any of them peaks.

### What a seedless cell can report

Measured outside the grid (`ab-short/seedless/`), so no A4 row is committed while
what to report is open: A4's own flags on f2112884, 1 thread, 1,800 s, on the two
`hbb` cells an external session had already measured from the other side.

| | `P54852__hbb` | `P01563__hbb` |
|---|---|---|
| protein | 163 aa | 188 aa |
| constrained points, all checked | 495 | 506 |
| seeds | **0** | **0** |
| certification | 36.9 s | 46.3 s |
| structurally optimal points | **0** | **0** |
| open region, log CAI | one, the whole range [−153.48, 0] | one, [−171.03, 0] |
| `mfe_floor_constrained` | **−357.9** | **−402.7** |
| `mfe_floor` | −640.37 | −897.09 |
| region slope | 2.277 | 3.399 |
| band | 21,755 designs in 1,749 s | 18,322 designs in 1,735 s |

So a seedless cell is not an empty table entry. It reports: no structurally optimal
design **at UAA**; its whole constrained front checked, so the seed set is certified
empty rather than merely unfound; one open region covering the entire adaptation
range, with the energy floor no design in it can beat; and a band that enumerated
~20,000 designs without a find. The two floors are 282 and 494 kcal/mol apart, which
is the quantity that says why nothing is structurally optimal there — a design's free
fold can run that much below what the constraint allows it, so the gap is far from
zero across the region.

The band's rate here is **12.4 designs/s at one thread** on a 163 aa cell, against the
~4/s measured on a 239 aa one, so the enumeration is not the obstacle at these
lengths — the absence of a find is.

This is the measurement behind option 2 of the seedless-cell question. It does not
settle what the paper reports, which is open: a 0-point row with these four numbers
beside it, or no A4 entry for those 29 cells at all.

### Why those cells are seedless: the initiator pairing into the 3′ UTR

`7937f01c` (`scope/37-obstruction-report`, three commits on f2112884, no PR) records
`gap_kcal_mol` and the **blocking pairs** on every structural optimality verdict, so
one 300 s run over a cell's constrained front is an obstruction census. Run on the
two seedless `hbb` cells, 1 thread, `--stop UAA`, `bin-7937f01c` (diagnostic only,
never a grid build; `--out` under `ab-short/seedless/`).

| | `P54852__hbb` | `P01563__hbb` |
|---|---|---|
| verdicts with a gap | 495 | 506 |
| gap kcal/mol, min / median / max | 0.2 / 2.9 / 6.0 | 0.6 / 5.0 / 12.7 |
| blocking pairs in all | 5,864 | 3,239 |
| distinct codons carrying them | 20 of 164 | 47 of 189 |
| designs within 1.0 kcal/mol | 151 | 86 |
| designs blocked by exactly one pair | 1 | **86** |

**The closest designs are blocked by one unrepairable pair.** Every single-pair design
on both cells is blocked by the same shape: the **first base of the initiator** paired
into the **3′ UTR** — molecule position 637 on `P54852__hbb`, 712 on `P01563__hbb`.
Neither end can move. `AUG` has no synonymous alternative, and the 3′ UTR is given
input. So `P01563__hbb` has **86 designs sitting 0.6 kcal/mol from structural
optimality that no codon choice can repair**, which is a structural statement about
the cell rather than a statement about the budget. It also shows why a different stop
codon need not help: the obstruction is anchored at the *initiator*, which `--free-stop`
does not touch.

**A figure of the same shape from another instrument, not to be pooled with this
one.** An external session reported ~60 % of near-optimal crossing pairs anchoring at
the two coding positions with no synonymous choice. That figure counts *pairs* within
a max-excess of a constrained minimum taken under **full mutual isolation**, over 12
front points; the census below counts blocking pairs over all 495 and 506 verdicts
under this grid's own `5-cds,cds-3`. Different populations, denominators and
baselines, so the agreement below is **not** a reproduction and the two must not be
pooled — they are recorded side by side with their instruments named. What the census
establishes on its own terms is that the share depends entirely on the gap band:

| subset | `P54852__hbb` | `P01563__hbb` |
|---|---|---|
| gap < 1.0 | **54.9 %** | **100 %** |
| gap < 2.0 | 55.4 % | 84.9 % |
| every checked design | 17.4 % | 20.6 % |

So within this one instrument the share is 55–100 % among near-optimal designs and
17–21 % over all of them: the word "near-optimal" is worth a factor of three, and
neither number may be quoted without the band it was taken in.

**The fixed positions explain the closest designs, not the cell.** 80 % of all
blocking pairs sit at interior codons that do have synonymous choice, so "the
forbidden pairs are concentrated where there is no choice" holds only inside the
near-optimal band. Outside it, most obstructions are at slots the search is free to
respell, and why it still finds no seed is unexplained.

**What the census does not establish.** Classifying pairs by whether an endpoint is
coding cannot come out any other way here: `--forbid 5-cds,cds-3` forbids exactly the
pairs with one coding endpoint, so a pair with none could not have been reported, and
"0 of 495 designs blocked only by unrepairable-by-position pairs" is vacuous rather
than informative. The informative columns are the gap, the pair count per design, and
*which* coding position the pair anchors at. Note also that 80 % of all blocking pairs
sit at interior codons that **do** have synonymous choice, so unfixable pairs do not
explain why the search finds no seed in general — only why its closest designs stay
closed.

### Which artifact decides that a cell has no seeds, and what certifying costs

The seed determination is `derna check` (762d509e) over each cell's whole A2 front:
certify only, no search, **no budget**. It folds every point and exits, so a row's
`points` equals the cell's A2 row count and `seeds: 0` means every point was folded
and none was structurally optimal. All 45 rows are `rc 0` with `points` matching the
front — including the four longest, at 1,627 to 2,329 points and walls of 1,210 to
2,385 s. So all 29 seedless cells are certified empty at UAA with nothing unchecked.

Beside those files sit `seeds/*__hbb*.w900.t2/` and `.t3/`: exploratory
`derna front --structurally-optimal` runs at `--max-secs 900` from before A4's
policy was set. At 900 s the certify stage does not finish past ~2,000 points, so
`Q9UBM8__hbb` carries 75 unchecked, `Q9NWZ3__hbb` 365 and `Q6B0B8__hbb` 520, and
their `seeds: 0` is a silence rather than a certified empty. Nothing cites them. An
external session read them as the determination and retracted the conclusion;
`seeds/README.md` now says which artifact decides, since the hazard is the directory
and not the reader.

**Certify cost, since `derna check`'s wall is certify cost.** Over the 21 `hbb`
cells, per-point cost runs 0.151 s at 188 aa (529 points, 80 s) to 1.086 s at 471 aa
(2,004 points, 2,177 s) — roughly **linear in length** over 188–478 aa rather than
N³, so much of it is per-point overhead rather than the fold. Totals: 16,440 s for
all 21. The consequence for any budgeted run is that a front of ~2,000 points owes
2,000–4,000 s of certification **before its search starts**, which is what makes a
900 s budget on a long cell return a silence. A4's own 6 h budget absorbs this
comfortably — 0.02 to 0.66 h — but a short arm does not.

### The stop codon is a measured dimension, not a confound

An external census, one `--free-stop` run or pinned arms per cell, every arm gated
on `unchecked == 0`:

| cell | aa | UAA | UAG | UGA |
|---|---|---|---|---|
| `P54852__hbb` | 163 | 0 | 0 | **303** |
| `P01563__hbb` | 188 | 0 | 0 | 0 |
| `P01588__hbb` | 193 | 0 | **6** | 0 |
| `Q9HAE3__hbb` | 211 | 0 | 0 | 0 |
| `RBD__hbb` | 223 | 0 | 0 | 0 |
| `P01563__hba` | 188 | 0 | 0 | 0 |

Four of six are certified empty at **all three** stops, and the two rescues come from
different suffixes — `UGA` on one cell, `UAG` on another. So neither "UGA rescues" nor
"no stop rescues" is right, and E4's wording for a 0-point row becomes **none at any
stop on the cells measured, with the exceptions named**, which is stronger than
hedging on UAA alone. 22 of E4's 29 seedless cells are `hbb`, the pair that census
covers, so its rows are, cell for cell, three quarters of what E4 must report. The two
instruments stay labelled rather than pooled: that census reads the SO front's own
certify stage inside `derna front`, at a budget; E4's reads `derna check`, unbudgeted.

### The obstruction class of a design, and a rule that the second cell breaks

`7937f01c`'s blocking pairs can be classed by which flank the non-coding end sits in,
`5-cds` or `cds-3`, which is the finer reading of the census above. On the two cells
measured at UAA, 1 thread, over the whole constrained front:

| | `P54852__hbb` | `P01588__hbb` |
|---|---|---|
| protein | 163 aa | 193 aa |
| designs | 495 | 529 |
| blocking pairs | 5,864 | 3,046 |
| `cds-3` pairs | 68.3 % | 68.5 % |
| `5-cds` pairs | 31.7 % | 31.5 % |
| designs carrying both | 280 | 194 |
| designs `5-cds` only | **130** | **0** |
| designs `cds-3` only | 85 | 335 |

**The pair share does not determine the design composition.** The two cells agree on
the per-pair split to within 0.2 points and differ completely one level up: 130 designs
blocked only at the 5′ end on one, none at all on the other. So any argument that runs
on the design-level class cannot be made from the pair-level share, and an earlier
prediction in this plan's favour — `P01563__hbb` at 90.3 % `cds-3` implying "rescued by
nothing" — was right for the wrong reason, since it used the share.

**What this tests.** An external census finds `P54852__hbb` rescued by `UGA` (303
seeds) and `P01588__hbb` by `UAG` (6). Matching that session's 311 `UGA`-tight designs
to this cell's front by coding sequence gave 114 counterparts, and **all 114 are
`5-cds` only** — 114 of these 130, against 0 of the 365 designs carrying any `cds-3`
pair. At a population share of 26.3 % that is not chance. The candidate rule it
suggests is that respelling the stop rescues a design when its obstruction is confined
to the 5′ end, which is a selection rule and not pair removal: the same test found **0
of 114** blocked by a pair touching the stop, so "the edit removed the blocking pair"
is refuted.

But `P01588__hbb` has **no** `5-cds`-only design, and 6 of its designs are rescued by
`UAG`. So the rule predicts zero rescues there and cannot be softened — 0⁶ admits no
sampling story. Either one of those 6 has a counterpart on this cell's front and the
rule is refuted outright, or none matches and the test is vacuous, as it was for the
197 unmatched designs on the other cell. That intersection is running externally; the
class table is written for it at
`ab-short/seedless/P01588__hbb.obstruct/blocking_class_by_cds.tsv`.

What stands either way is narrower and worth keeping: on one cell and one suffix, every
rescued design was blocked only at the 5′ end and no design carrying a `cds-3`
obstruction was rescued — with both directions budget-limited, since that arm stopped
on time with three regions open. And the near-miss designs of both cells are `cds-3`
blocked — the initiator's first base paired into the 3′ UTR — so the local certificate
of 86 designs on `P01563__hbb` covers designs that no stop can rescue, which is
consistent with that cell being empty at all three.

### An unfixed stop codon for A1 and A2 needs a code change, and waits for one

Wanted: A1 and A2 over designs whose stop codon is **not** fixed, the UTRs being
included. The current `derna` cannot express that for a walk.

- `--free-stop` exists but **refuses** any run without `--structurally-optimal`, by
  design: it chooses the stop by whether the design it ends is structurally optimal,
  which a walk has no notion of. The refusal is a hard error, checked.
- For a walk, `--stop` takes exactly one of `UAA|UAG|UGA|none` — a *fixed* codon,
  "coding sequence to the constraint and outside CAI". So the stop can anchor a
  `cds-3` pair while contributing no adaptation, and there is no optimise-over-stops
  mode. All three codons are individually accepted; none of them is a choice.

Three pinned walks plus a Pareto union **is** exactly the unfixed-stop front — the
stop being outside CAI puts all three arms on one adaptation scale, so the union is
the lowest energy at each adaptation level across the stops, ties kept where designs
differ only in their stop. But that union is post-processing, not an artifact the code
emits, and nothing downstream consumes it: `--front-from` validates the molecule
**including the stop codon** and refuses an incomplete walk, so A4 cannot read a
mixed-stop front; and A3 takes its weights from a front file's `on_hull` flags, which
a union would not carry without recomputing the hull outside the binary.

Cost if it were done by three arms: A1+A2 over the 46 cells measured **48,334 s at 3
threads**, so three stops is **40.3 h** of 3-thread work — about 20 h wall on E4's
eight cores, two cells at a time, 201 s for the cheapest cell and 3,230 s for the
dearest. An external session's own reading is that a stop-as-slot DP would be ~3×
cheaper than three arms on certification, which this plan's own curve shows is the
dominant term at length, and that one 3-block slot among 160–480 codons barely moves
front size.

**Decision taken 2026-10-22 ~13:30:** the rerun waits for a code change, A4 continues
on the eight cores in the meantime, and another session will say when the rerun can
run. One cell is being measured first — `P54852__hbb`, A2's own flags at all three
stops, 1 thread — to establish whether an unfixed stop changes the front at all before
anyone is asked to build it. Its UAA arm returned 495 points against the grid's 495,
so f2112884 reproduces 6d42ff3f on that cell.

### f2112884 is not reproducible in memory or walk structure, and every peak here is one draw

An external session replicated its own cross-build memory A/B and retracted the ratio
it had reported. On `eGFP__curevac`, same binary, same flags, interleaved in one
script: `6d42ff3f` gave 0.88 GB and 32 finite-walked stretches **twice, to four
decimals**, while `f2112884` gave **3.93 and 4.97 GB** with **14 and 86** walked
stretches, and neither replicate reproduced that arm's original 9.48 GB. So the
treatment spans 3.93–9.48 GB over three runs of one binary on one cell at one budget,
the cross-build figure is a **range of 4.5–10.8×** with an unstable top rather than a
point estimate, and a sizing argument must take the maximum observed.

The finding underneath is larger than the ratio: **the older build is deterministic in
memory and walk structure and the newer one is not.** That matches what E4 measured on
`eGFP__bnt` from the other direction — three identical f2112884 arms returning
225/244/244 points and 14/14/72 walked stretches — which had been written off as that
cell's peculiarity precisely because the external session's *control* replicated
bit-identically. A control that replicates licenses nothing about the treatment arm.
In every replicate the widest walked region was identical, so what varies is how many
*other* stretches are visited, not which one is the big one; the band hand-back making
region visits load-dependent is the suspect, and `--threads 1` would test it.

**Consequence for this plan's own numbers.** Every A4 peak recorded here — 3.47, 6.45,
17.71, 33.95 and 247.9 GiB — is a single draw from a build with no demonstrated
reproducibility in that column, and so is every walked-stretch count. The point counts
and the fronts are not in question: those were checked by domination across arms and
were consistent. But the width-orders-peaks reading rests on one measurement per cell,
and until an A4 cell is run twice, the honest statement is that the ordering has never
been shown to survive a repeat. Re-running one cell at the same settings is the cheap
test and is worth more than another cell at a new width.

### Lower energy is not bought at the junction: a proposed mechanism, tested and refuted

An external session proposed a mechanism that would unify E4's stop-codon energy result
with its own seed census: if a stop lowers MFE by pairing into the 3′ UTR, it creates
exactly the pairs `5-cds,cds-3` forbids, so a lower-energy stop should yield **fewer**
structurally optimal designs. It fits the inputs — `UAG` beats `UAA` on 30 of 30 fixed
designs by a median 2.49 kcal/mol and yields **0** seeds, while `UGA` is an energetic
wash (median +0.00) and yields **303**.

Tested directly, by folding the same 30 designs under each stop and counting the pairs
of the **free** fold that cross a flank boundary (means over the designs):

| stop | median ΔE vs UAA | 5′UTR–CDS | CDS–3′UTR | total crossing |
|---|---|---|---|---|
| UAA | +0.00 | 3.53 | 7.17 | 10.70 |
| **UAG** | **−2.50** | **1.07** | **3.77** | **4.83** |
| UGA | +0.00 | 2.60 | 6.87 | 9.47 |

**The mechanism is refuted, and the sign is the opposite of the proposal.** `UAG` buys
its 2.5 kcal/mol while **halving** the crossing pairs, 4.83 against `UAA`'s 10.70, and
cutting the `CDS–3′UTR` class specifically from 7.17 to 3.77. So its energy does not
come from reaching across the junction; it comes from elsewhere, and it leaves the
junction *freer*. Fewer forbidden pairs in the free fold should make structural
optimality more likely, not less, which leaves `UAG`'s zero seeds more puzzling than
before rather than explained.

What the test does not cover: these are designs from the `UAA` front with the stop
swapped, so it measures what a stop swap does to a **fixed** design, not what each
arm's own constrained front contains. The seed asymmetry is a statement about the
latter, and remains unexplained. Anything built on the junction-pairing story should
stop here.

### The stop that wins the constrained front is the one that destroys structural optimality

Measured on the population that matters — each arm's **own** constrained front, not one
arm's designs rebuilt under another's stop. 60 points per arm spanning the front of
`P54852__hbb`, each folded free by `mfe_fold`, pairs classed as in the obstruction
census, and `gap` taken as the point's constrained energy minus its own free energy:

| arm | 5′UTR–CDS | CDS–3′UTR | total crossing | mean gap kcal/mol | seeds |
|---|---|---|---|---|---|
| UAA | 3.63 | 8.07 | 11.70 | 3.25 | 0 |
| UAG | 0.15 | 4.03 | **4.18** | **5.79** | 0 |
| UGA | 1.30 | 3.37 | **4.67** | **1.29** | **303** |

**The number of forbidden pairs does not explain the rescue.** `UAG` and `UGA` carry
almost the same crossing-pair count, 4.18 against 4.67, and one yields no structurally
optimal design while the other yields 303. So a pair-counting story fails on the right
population as well as on the wrong one.

**The gap does explain it, and it is the energy those pairs are worth.** Structural
optimality is `gap = 0`, and the arms' mean gaps are 5.79, 3.25 and 1.29 — `UGA`'s front
sits closest to tightness by a factor of four and is the arm that rescues. The
quantitative match is tight: `UAG`'s gap exceeds `UAA`'s by **2.54** kcal/mol, against
the **2.50** median by which `UAG` lowered the free fold on 30 fixed designs. So the
free energy a stop buys is very nearly the gap it adds.

The reading that follows is that **a stop which lowers the free fold's energy raises the
gap**, and `UAG` is the case: −2.50 kcal/mol on 30 fixed designs, +2.54 kcal/mol of gap,
no seeds.

**But there is no tension between A2 and A4 here, and an earlier version of this section
claimed one.** Recomputing the union's composition: of its 465 points, **`UGA` supplies
293, `UAG` 172 and `UAA` none**. So `UGA` — the arm with the smallest gap and all 303
seeds — also supplies 63 % of the union. It is the better stop on *both* objectives on
this cell, and the two arms agree rather than conflict. The claim that the union "would
improve A2 and empty A4" was wrong, and so was the reading that `UAG` dominates the
union: `UAG` wins on *designs taken from `UAA`'s front and re-stopped*, which is a third
population again, and that does not carry to each arm's own front. A two-arm `UAA`/`UAG`
shortcut would therefore be doubly wrong — it drops the arm supplying most of the union
*and* every seed.

A free-stop structurally optimal front also cannot be emptied by whichever stop the plain
front prefers, as the scope-37 successor pointed out: structural optimality is a property
of a whole sequence, so every `UGA` design that is tight stays tight in the union, and the
union's SO set contains all three arms' SO designs. The cost of a free stop is three
searches, not a worse answer. The one real hazard is a **pipeline** one: if A2's union
front were used to seed a single pinned-stop A4, the seeds would be points from an arm
whose designs are not tight, and A4 would find fewer. That is a wiring question about
`--front-from`, not a property of the free-stop front.

Caveats: one cell, 60 points an arm, and the gap comparison is across three design sets
rather than within one. An external session's pair-count version of this mechanism was
refuted separately and should not be revived; what survives is the energy version, and
it predicts that any stop beating `UAA` on the plain front will do worse on seeds.

### Two cells on the stop codon: what generalises and what does not

Both cells, A2's own flags, f2112884, 1 thread, three pinned walks each, with the UAA
arm reproducing the grid's `6d42ff3f` front exactly in both cases (495 of 495 and 586
of 586, nothing either side — the build control for everything below). Gaps and
crossing-pair counts are per arm over its **own** front points, 40–60 a cell.

| | `P54852__hbb` (163 aa, hbb) | `eGFP__bnt` (239 aa, bnt) |
|---|---|---|
| UAA / UAG / UGA front size | 495 / 482 / 477 | 586 / 586 / 591 |
| UAG against UAA | 153 of 495 shared | **identical, 586 of 586** |
| distinct union points | 465 | 591 |
| reachable at UGA | 293 of 465 | **591 of 591** |
| reachable at UAA | **0** of 465 | 329 of 591 |
| UAA points dominated by another arm | 495 of 495 | 257 of 586 |
| mean gap, UAA / UAG / UGA | 3.25 / 5.79 / **1.29** | 3.29 / 7.15 / 4.86 |

**What holds on both.** `UGA`'s front contains the whole union on one cell and 63 % of
it on the other, so it is the best single stop for the plain front twice over, and
pinning `UAA` loses points on both — all 495 on one cell, 257 of 586 on the other.
`UAG` has the **largest gap on both cells**, 5.79 and 7.15 against `UAA`'s 3.25 and
3.29, so it is consistently the worst stop for structural optimality. And `UAA`'s own
gap is nearly identical across two unrelated cells and flank pairs, 3.25 and 3.29.

**What does not hold.** "`UAA` contributes nothing to the union" is false on the second
cell: it reaches 329 of 591 points there, where it reached none on the first. `UAG`'s
relation to `UAA` is not stable either — the two fronts are **identical** on `eGFP__bnt`
and share only 153 of 495 on `P54852__hbb`. And `UGA` is best on gap on one cell (1.29)
and middle on the other (4.86 against `UAA`'s 3.29), so the gap ordering is not fixed
beyond `UAG` being worst.

So the defensible cross-cell statement is narrow: **a `UAA`-pinned plain front is
dominated in part or in whole on both cells measured, and `UGA` recovers it; which stop
sits closest to tightness is cell-specific except that `UAG` is furthest.** Anything
stronger needs more than two cells, and the two here disagree on most of the detail.

### Two cells, two rescue mechanisms, and the obstruction census predicts which

An external session's third obstruction census refutes a generalisation this plan made
and separates two things that had been pooled. Verified independently from E4's own
census files, which reproduce their numbers exactly:

| cell | designs blocked by one pair | at the **initiator** | at the **stop** |
|---|---|---|---|
| `P54852__hbb` | 1 | 1 | 0 |
| `P01563__hbb` | 86 | 86 | 0 |
| `P01588__hbb` | **281** | 107 | **174** |

So "every single-pair design is blocked by the initiator paired into the 3′ UTR" holds
on the two cells it was measured on and **fails on the third**, where the majority
anchor at the stop. The narrow claim about `P01563__hbb`'s 86 stands; the
generalisation does not. E4 had the census for the third cell and did not recompute the
split on it, which is how the claim survived a cell that refutes it.

**The initiator and the stop must not be pooled.** `AUG` has no synonymous alternative,
so an initiator-anchored block cannot be repaired by any codon choice — a property of
the molecule. The stop has three alternatives and reads as fixed *only because every run
pins `--stop UAA`* — a property of the run. Pooling them mixes a molecular fact with a
flag, and the two diverge by cell:

| cell | share at gap < 1.0, initiator only | initiator **or** stop |
|---|---|---|
| `P54852__hbb` | 51.7 % | 54.9 % |
| `P01563__hbb` | 100 % | 100 % |
| `P01588__hbb` | **20.3 %** | **55.7 %** |

On the third cell the pooled figure is 2.7× the initiator-only one. Every share this
plan quotes for a fixed-codon position should be read as the split pair.

**The rescue mechanisms differ by cell, and the census predicts which.** `P01588__hbb`
is rescued by `UAG` (6 points) and is the cell whose single-pair blocks are mostly
stop-anchored; `P54852__hbb` is rescued by `UGA` (303 points) and has **zero**
stop-anchored single-pair designs. Testing the first directly: each of the 6
`UAG`-rescued designs was rebuilt with `UAA` and folded free, and **all six carry
exactly one crossing pair, anchored at the stop**. So on that cell the rescue *is* pair
removal — changing the stop deletes the one pair that blocks the design — which is the
mechanism proposed and then refuted on `P54852__hbb`, where the rescued designs carried
only `5-cds` blocks and none touched the stop. Both results stand; they are different
cells with different obstructions.

A set intersection could not have shown this. Matching the 6 rescued designs against
the 174 stop-anchored ones returns **0 of 174** — but vacuously, because none of the 6
is on `UAA`'s front at all, the same gap as the 197 unmatched designs on the other cell.
The answer came from rebuilding the designs, not from intersecting the sets.

### Neither composition nor the constraint's energy cost predicts whether a flank pair seeds

Two candidate explanations for `hbb`'s 0 of 22 are now refuted against matched controls,
both computed by an external session from E4's existing A1/A2 fronts with no new runs.

**Composition.** 3′ purine fraction is **0.43 for both `bnt` and `hbb`**, at 100 % and
0 % seeded. So the AU-rich story is refuted outright rather than left unsupported, and
with a matched control rather than a correlation.

**The constraint's front-level energy cost.** With A1 and A2 being the same molecule,
protein and stop and differing only in `forbid`, the cost at adaptation level L is
`min{A2.mfe : cai_log ≥ L} − min{A1.mfe : cai_log ≥ L}`. The prediction was that a pair
whose constraint costs more energy should seed less. The result reverses it:

| pair | seeded | median of cell medians |
|---|---|---|
| `bnt` | 6 of 6 | 13.60 |
| `curevac` | 5 of 6 | 11.70 |
| `hba` | 4 of 6 | 19.40 |
| `mrna1273` | 2 of 6 | 15.00 |
| `hbb` | **0 of 22** | **8.00** |

`hbb` has the **cheapest** constraint of all five pairs and yields nothing, and on the
matched companion panel it is the cheapest column at 5 of the 6 proteins. Note the
quantity is a **front-level** cost and not a per-design gap: the two come apart, and
`P01241__bnt` has 11 seeds at a cost floor of 7.20, so a nonzero front cost does not
imply that no design is tight.

**Why it cannot be a predictor even where it correlates.** Within a pair, cost and seeds
run broadly opposite as predicted (`bnt`: 1.00 → 672 seeds, 25.60 → 37); across pairs
the ordering inverts. A quantity that orders cells inside a group and anti-orders the
groups is unusable for an unseen pair, which is the only use the generalisation claim
would have for it.

**Five candidate properties of a flank pair are now refuted against matched controls**,
by an external session over E4's artifacts: composition, the front-level constraint cost,
the initiator obstruction as an explanation of emptiness, the obstruction's
concentration, and its **span**. Three of the five fail with the two matched protein rows
pointing in **opposite** directions, which rules them out as predictors rather than
merely failing to establish them. On span: median `cds-3` span is 561 on the empty cell
against 107 on the seeded one for one protein, and 74 against 690 — the reverse — for the
other. And `P01588__bnt` carries **15,600** `cds-3` pairs at median span 690, an
obstruction across the whole coding sequence, and still yields 37 tight designs, so
long-range `cds-3` blocking is fully compatible with seeding.

**So the method statement for the generalisation claim is that an unseen pair must be
tested, not modelled** — and that is affordable rather than a dead end: one exact
constrained front plus one unbudgeted `derna check`, 60–285 s over the companion cells by
this plan's own certify curve. "The screen is cheap enough that no model of the pair is
needed" is a usable claim, and it is also what gives the held-out pair its shape.

**Consequence for the generalisation claim.** If no compositional or energetic summary
of a flank pair predicts whether it seeds, then screening an unseen pair may require
**folding designs on it** rather than inspecting the pair. That is affordable rather
than a dead end: by this plan's own certify curve, `derna check` over a whole A2 front
takes **60–285 s** across the 30 companion cells (the external estimate of 66–212 s is
narrower than the measured range). So the screen for a new pair is one A2 walk plus one
check, not a model of the UTRs.

### The matched control: the fixed-codon obstruction is absent on the seeded pair

The obstruction censuses had only ever run on `hbb` cells, which are 0 of 22 — a
mechanism measured only where the outcome is uniform cannot explain the outcome. The
control is the **same proteins under `bnt`**, which is 6 of 6 seeded: same protein, same
length, same `--forbid`, same `--stop UAA`, same binary, same 300 s budget, the `run()`
body copied from the `hbb` script so both sides are one instrument. Every census is
**complete** — `checked` equals the cell's front size on all four — so nothing here is a
sample. Failing designs are those not structurally optimal and carrying blocking pairs;
a tight design has `gap_kcal_mol` of 0 and must be excluded or it inflates the failing
count.

| cell | tight | failing | median gap | blocked by one pair | at initiator | at stop | interior | failing within 1.0 |
|---|---|---|---|---|---|---|---|---|
| `P01563__hbb` | 0 | 506 | 5.00 | **86** | 86 | 0 | 0 | **86** |
| `P01563__bnt` | **78** | 429 | 2.90 | **0** | 0 | 0 | 0 | **0** |
| `P01588__hbb` | 0 | 529 | 0.60 | **281** | 107 | 174 | 0 | **294** |
| `P01588__bnt` | **37** | 493 | 4.10 | **7** | 0 | 0 | **7** | **0** |

**The pre-registered refuting outcome DID occur, read on the quantity the rule names.**
The rule was stated over the initiator shape's *share*, and the share that matters is
design-level — the fraction of failing designs carrying at least one initiator-anchored
pair. That is **83.9 % and 86.0 % on the seeded `bnt` cells**, against 90.9 % and 28.4 %
on the empty `hbb` ones: comparable on one protein and **three times larger** on the
other. The shape is at least as common where seeds exist. So **the initiator obstruction
is withdrawn as an explanation of emptiness**, surviving only as a description of
particular cells' near-miss designs, with `P01563__hbb`'s certificate standing on its own
terms. That is the third candidate refuted against a matched control, after composition
and the front-level constraint cost.

An earlier version of this section read the refutation off the **single-pair count** (0
and 7 on `bnt` against 86 and 281) and the **pair-level** share (19.2 % and 4.3 %), and
concluded the outcome had not occurred. Those are different quantities from the one the
rule fixed, and reading a pre-registered rule on a quantity it did not name is how a
registration gets defeated by its own author.

**And no failing design on either `bnt` cell is within 1.0 kcal/mol of tight** — 0 and 0,
against 86 and 294 on the `hbb` cells. So the seeded pair's designs either succeed (78
and 37 tight) or fail by a margin, while the empty pair holds hundreds of designs a
fraction of a kcal/mol short, each blocked by a single pair at a position with no
synonymous choice.

**The real finding is the gap distribution, and it inverts the natural reading.** The
empty pair **crowds the threshold and never crosses it**: minimum gaps of 0.60 and 0.50,
with 86 and 294 failing designs inside 1.0 kcal/mol, each blocked by a single pair. The
seeded pair is **bimodal** — 78 and 37 designs reach gap 0 exactly, and its failing
designs start at **1.60 and 2.60**, so *none* lies in the gap < 1.0 band at all. Hence:

> A population of near-miss designs each blocked by one pair is the signature of a cell
> that is **empty**, not of one that nearly succeeded. Where seeds exist, the
> almost-tight designs *are* the seeds; nothing is left sitting just below the line.

So `hbb`'s 86 and 281 single-pair near-misses are not evidence it almost worked — they
are what empty looks like.

Two confounds, stated rather than buried. The `bnt` 3′ UTR is **295 nt against `hbb`'s
134**, so pairs per failing design (8.75 and 32.48 against 6.40 and 5.76) is unsafe as a
pair-to-pair comparison; it does not touch the gap distribution, which is an energy, or
the single-pair counts. It also cuts *against* an obstruction story, since the pair with
twice the 3′ UTR and four times the blocking pairs is the pair that seeds. And the class
split does not order seeding: `P01563__bnt` is 66.2 % `5-cds` while `P01588__bnt` is
97.4 % `cds-3`, two seeded cells of one pair with opposite composition.

`P01588__bnt`'s 7 single-pair designs anchor at neither no-choice position but at
**interior codon 8**, repairable in principle unlike the initiator or a pinned stop.

Median gap among failing designs does not separate the pairs: 5.00 against 2.90 on one
protein and 0.60 against 4.10 on the other, in opposite directions. So the discriminating
quantities are the single-pair count, where it anchors, and how many failing designs sit
within a kcal/mol — not the gap's central tendency.

### Pending: rerun A1 and A2 with the stop as a variable, when the DP lands

**The decision, 2026-10-22:** A1 and A2 are to be rerun over designs whose stop codon is
**not fixed**, the UTRs being included. This waits on a code change and is the first
thing to run when it arrives.

What it waits for is the stop as a **lattice variable** — one walk over a design space
that includes the stop slot, one certification — and not the merged `--free-stop`, which
merges three structurally optimal searches and so pays three walks and three
certifications, refuses `--report`, `--front-from` and `--checkpoint`, and cannot be
expressed for a plain walk at all. The two differ by 3× in cost and by everything in
drivability:

| | the merged three-arm flag | the stop as a lattice slot |
|---|---|---|
| walks per cell | 3 | 1 |
| certifications | 3 | 1 |
| A1+A2 over the 46 cells | **40.3 h** of 3-thread work | **~13.4 h**, i.e. what A1+A2 already cost |
| works on a plain walk | no, refused | yes, that is the point |
| A3 can read its weights | no — `--report` refused | yes |
| A4 can seed from it | no — `--front-from` refused, header stop differs | yes, if A4 also runs free |

Why the cost barely rises: every stop codon weighs zero in CAI, so the adaptation axis is
unchanged, and front size tracks synonymous freedom — one 3-nucleotide slot among 160–480
codons hardly moves it.

**The validation targets already exist.** Three pinned walks and a Pareto union is
mathematically the same front the DP should return, and that is measured on two cells:
`P54852__hbb` (union 465 distinct points, `UGA` 293, `UAG` 172, `UAA` 0) and `eGFP__bnt`
(union 591 distinct points, all reachable at `UGA`, `UAA` and `UAG` identical at 586
each). On both, the `UAA` arm reproduces the grid's own `6d42ff3f` front exactly — 495 of
495 and 586 of 586, nothing either side — which is the bit-for-bit control the DP needs
before anything else is read off it.

**What has to be decided when it lands, not before:** a free-stop A2 changes A3's hull
weights and A4's seeds, and `--front-from` matches the header's stop, so a free-stop A2
requires a free-stop A4 rather than silently seeding a pinned one. That is a grid-wide
choice — rerun A3 (2,888 rows) and A4 (17 seeded cells at 6 h) against the new front, or
carry the free-stop arms beside the pinned ones.

Status at the time of writing: `scope/37-stop-dp` is off `scope/38-so-pools`, with the
21st codon family in and type-checking and the CLI half — pinned/free/absent and the CDS
window over the slot — still to come. Watcher armed on `origin` for that branch moving
and for `xg` gaining stop-slot commits.

## The obstruction census already existed over the whole grid

The eight budgeted obstruction runs on `bin-7937f01c` re-derive data the grid already holds.
`seeds/<cell>.check.762d509e.jsonl` carries, per point, `gap`, `structurally_optimal` and
`violations` — the forbidden pairs the free fold reinstates — for all **46** cells, unbudgeted,
every point folded. `7937f01c` added `gap_kcal_mol` to the *SO front's* verdicts; `check` had been
emitting the equivalent all along. The whole-grid table therefore supersedes the eight-cell one as
the figure of record, and the eight runs are kept only as an equivalence control, read by
intersecting designs rather than comparing counts (each is a 300 s prefix of its cell's A2 front).

Two format facts, both traps at the flank boundary where these obstructions live:

- **`violations` is 1-based; `blocking` is 0-based.** `packages/derna/src/result.rs`
  `verdict_fields()` maps `gap.reinstated_pairs` through `json!([i+1, j+1])`, while `--report`
  emits that same field unshifted. Classifying positions on `violations` with a layout derived for
  `blocking` reclassifies a pair at the first coding base as `5-cds` when it is `cds-3`. The shift
  is in the serializer, so it holds identically on every cell and does not need a per-cell
  agreement test to establish.
- **A verdict row's `sequence` and `free_structure` span the whole mRNA**, not the CDS
  (`free_fold.rs:191`). `P01563__hbb` reads 751 nt = 50 + 567 + 134. Reading that 751 as the CDS
  and concatenating the UTRs a second time produces a length mismatch and 37 % pair
  complementarity — a null manufactured by the wrong layout.

## Designs blocked solely by one pair: the separator that survives

Over all 46 cells, designs whose obstruction is **exactly one** pair number **7,008 across the 22
hbb cells against 22 across the other 24** — about 300× per cell. Within hbb that single pair
anchors at the initiator 5,666 times and at the stop 1,342.

This restates the initiator finding in the only form the measurements support. *Carrying* an
initiator-anchored pair is common on both flank pairs (83–91 % of designs) and was already
withdrawn as an explanation of emptiness. Being blocked by **nothing else** — and so being
unrepairable by any single re-spelling — is hbb-specific. The two quantities are different and had
been collapsed.

The strong cross-pair form is not quoted: "every initiator-anchored single pair is on hbb" rests on
0 of **22** designs, since the non-hbb cells have only 22 single-pair designs in total. The
load-bearing number is the 7,008-against-22 rate, which has a real denominator.

Two hbb cells break the separator — `eGFP__hbb` and `Q9HAE3__hbb` have zero single-pair designs and
no dense band. They are the same two cells that break the band display, which makes one cell
property more likely than two independent exceptions.

A companion 2×2 display, "seeded cells have nothing in the gap < 1.0 band", is **withdrawn**: 6 of
17 seeded cells do have a dense band (`RBD__curevac` 411 tight with 167 band designs,
`P01241__hba` 7/73, `RBD__bnt` 483/38, `P01241__mrna1273` 507/33, `eGFP__curevac` 163/14,
`P09919__curevac` 631/13), and 8 of 29 empty cells have none. What survives is weaker and
quantitative: 72 % of empty cells have dense bands against 35 % of seeded ones, and the magnitudes
separate — empty cells reach 1,030, 1,011, 790 and 786 band designs while no seeded cell exceeds
167.

## The obstruction locus: one A–U pair at one position

The single-pair obstruction resolves to **one base pair at one position**, which is why five bulk
statistics read as five independent nulls when they were one fact seen five times.

Reproduced independently off the 22 `seeds/*__hbb.check.762d509e.jsonl` files, classifying the 5′
endpoint by `i == n5` and the 3′ endpoint as `j − (n5 + ncds)` after the serializer's −1:

| quantity | value |
|---|---|
| single-pair designs, 22 hbb cells | 7,008 |
| initiator-anchored | 5,666, over 17 cells |
| distinct 3′ offsets among those | **{95}** |
| stop-anchored | 1,342, distinct offsets **{14}** |
| other 5′ endpoints (CDS interior, 3′ UTR) | 0 |
| pair letters | `('A','U')` |
| 3′ context | `GCAUC[U]GGAUU`, one string across all 17 cells |
| 5′ context | `ACACC[A]UG` + second codon, 18 spellings |

The invariant part of the 5′ context is the Kozak-proximal `ACACC` and the `AUG`; all variation is
downstream of the paired base.

**Base availability is refuted as the mechanism, inside hbb.** hbb's 3′ UTR is 134 nt and carries
**48 uridines**, every one as complementary to the initiator `A` as offset 95. **One of the 48 ever
appears as a sole obstruction.** Length, protein and context are held fixed, so this is a tighter
control than any cross-pair comparison: the instrument had 48 chances and took one. What decides it
is the loop weights of the surrounding context, not whether a complementary base exists.

**Offset 95 is not a cross-pair coordinate.** The base there is `U` in hbb, `A` in bnt, `G` in hba,
`G` in mrna1273 — and curevac's 3′ UTR is **exactly 95 nt**, so the offset does not exist in it
(noted so the absence is not later read as a bounds error). The supportable cross-pair claim is the
simpler one: **zero initiator-anchored single-pair designs at any offset** across all 18
bnt/curevac/hba cells.

**The stop endpoint supplies the cross-pair control.** Offset 14 is a `U` in curevac as well as in
hbb, and curevac's six cells hold 7 single-pair designs with **0 stop-anchored**. So on the stop
side a complementary base is present in a pair that yields no obstruction — insufficiency shown
rather than inferred. bnt, hba and mrna1273 all carry `G` at offset 14, making curevac the only
pair that can serve as this control. mrna1273's 8 stop-anchored designs sit at offset 22, so the
stop-side locus is pair-specific.

**Scope.** This is a sufficient obstruction for 5,666 identified designs — **5,666 of 25,906 hbb
constrained points, 21.9 %, on 17 of 22 cells** — and not an explanation of hbb's emptiness: all 22
cells are certified empty, the 5 cells it never touches included, and multi-pair designs are
untouched by it. It is also not predictive, for the reason above. The share is stated so a later
reader cannot promote it into an explanation.

## The memory cap's cost appears as coverage, not only as early termination

`RBD__curevac` is the grid's first memory-stopped A4 row and behaves as the hold-out decision rule
assumes: 457 points, **all 457 `status == "certified"`**, zero with a nonzero gap, 46 beyond its
411 seeds, `stopped_by: "memory"`, peak 127.85 GiB. `complete: false` here is the benign kind — the
pair `complete: false` + `stopped_by: "memory"` is a self-stop, not the `stopped_by: null` that
means an external kill or a checkpoint read mid-flight. The cell was verified exited, not live,
before the row was read.

The cap ended it at **4.06 h with 1.95 h of its 6 h budget unspent**.

**Coverage is where the cost shows, and it must be normalized.** Raw covered log-CAI span confounds
the cell, since span depends on how wide that cell's front is and `RBD` is the only `RBD` row in
the set. Against each cell's own available A2 span:

| cell | A4 covered | A2 available | share | peak |
|---|---|---|---|---|
| P09919__bnt | 65.36 | 65.69 | 99.5 % | 33.95 GiB |
| P09919__mrna1273 | 64.26 | 65.69 | 97.8 % | 17.71 GiB |
| P01241__mrna1273 | 68.36 | 71.04 | 96.2 % | 247.93 GiB |
| P09919__curevac | 62.96 | 65.69 | 95.8 % | 3.47 GiB |
| P09919__hba | 62.96 | 65.69 | 95.8 % | 6.45 GiB |
| **RBD__curevac** | **27.53** | **56.37** | **48.8 %** | 127.85 GiB |

The confound does not explain it away — the front is 56.37 wide and the row reached under half of
it, against 96–100 % for every time-stopped row. **Report the share, not the raw span.**

The missing coverage is at the **low-adaptation end**: A4 spans `[−29.79, −2.26]` against A2's
`[−56.37, 0.00]`, so 26.58 units are missing below and 2.26 above, where `P09919__bnt` is
`[−65.69, −0.34]` against `[−65.69, 0.00]` and misses nothing below.

**Held at n = 1.** One capped row against five uncapped ones, and the protein is still confounded
with the cap. Recorded as a registered expectation — coverage share is the first column to read on
the next capped row — not as an established property of the cap. Peak and coverage do not order
together: 247.93 GiB at 96.2 % against 127.85 GiB at 48.8 %.

**A partiality column must not sum a below-front region with in-front windows.** Dividing total
open-region width by the covered span printed an "open share" of 558.5 % on this row, because its
open regions include the below-front region `(−181.47, −29.79]` of width 151.68 against a covered
span of 27.53. A below-front region's width runs to the least adaptation the codon table allows, so
summing it is meaningless. Split the column three ways: in-front (clipped to the covered range),
above the most adapted point reached — the one that matters, since that is where a budget runs out
— and below-front counted rather than summed.

## Amendment: "initiator" is wrong for two cells, and the locus is a property of the position

Two corrections to "The obstruction locus", both tightening it.

**The word.** 20 of 22 hbb cells begin with `AUG`; `P01860__hbb` begins with `GCC` (Ala) and
`RBD__hbb` with `CGC`/`CGG` (Arg). On those two the first coding position has synonymous spellings,
so "the initiator has no alternative" does not apply and the label "initiator" is simply wrong.
Read `cds_start` as *the head of the coding sequence* and reserve "initiator" for the 20 cells
where it is one. The 5,666 figure is unaffected: every sole obstruction in the census reads `A–U`
at offset 95, so all of them sit at a genuine `AUG`.

**The locus is about the 3′ position, not the initiator's `A`.** Measured at `cds_start` on the two
non-`AUG` cells:

| cell | first codon | pair at `cds_start` | 3′ offset | count |
|---|---|---|---|---|
| `P01860__hbb` | `GCC` | **G–U** | **95** | 42 |
| `RBD__hbb` | `CGC`/`CGG` | C–G | 96 | 295 |

`P01860` reaches **the same offset 95** through a `G–U` wobble where the 17 `AUG` cells make a
Watson–Crick `A–U`. So offset 95 attracts whatever occupies the head of the coding sequence, across
two base identities and two pairing geometries — a stronger and more positional claim than "the
initiator's `A` finds a `U`".

**Reachability is pooled, and thin per molecule.** The pooled set of offsets the head reaches is
{87, 95, 96, 102}, but **18 of 22 cells only ever use 95**; `P01563__hbb` and `P78382__hbb` use
{87, 95}, `Q6PF06__hbb` and `Q9HAE3__hbb` use {95, 102}, and `RBD__hbb` uses only {96}. No cell has
more than two reachable offsets, so "95 out-competes three rivals" is true pooled and false per
molecule, where it out-competes at most one. Same shape as quoting a pooled count as a per-instance
one.

**The designable-first-position contingency is one cell, not two.** `P01860`'s first codon is `GCC`
on **every design in its front** — never respelled — so "designable, therefore repairable" is
demonstrated only on `RBD`, which does vary `CGC`/`CGG`. State it as one cell demonstrably
respelled and carrying no sole obstruction, and one respellable in principle but never respelled in
practice.

**So the five cells the locus never touches split three ways:** `P01860__hbb` and `RBD__hbb` by
first-codon identity, and `Q92734__hbb`, `Q9HAE3__hbb`, `eGFP__hbb` unexplained — the last two
being the same pair that breaks the band table. Two cells failing two independent tables together
is the unit worth investigating next.

## The omitted `--stop` flag was a silent-change hazard, now closed

With the stop codon becoming a slot of the codon lattice, **omitting `--stop` with a flank present
changes meaning**: it used to select no stop codon and now means *design the stop*. A run that
meant "no stop" would compute a free-stop molecule on a newer build and emit an ordinary-looking
row.

The grid was never exposed — all 25 E4 configs carry `"stop": "UAA"`, `drive_e4.py` builds the flag
from `cfg["stop"]`, and the running cell's command line shows `--stop UAA` — but the driver held a
branch that omitted the flag whenever a config said `none`, and that branch is the hazard rather
than the config. `--stop` accepts `none` on every build the grid uses (and `none` is its default),
so the driver now spells the stop out unconditionally and the molecule is named on every command
line, where a header can be checked against it. A config that drops the `stop` key still raises.

## The two cells that broke two tables are one cause with three symptoms

`eGFP__hbb` and `Q9HAE3__hbb` kept appearing together — no single-pair designs, no dense band, no
sole obstruction at the locus. They are not a separate phenomenon. Over all 22 hbb cells, minimum
gap against the fewest pairs blocking any design:

| fewest pairs blocking any design | cells | minimum gap |
|---|---|---|
| **1** | **20** | 0.10 – 0.60 |
| 3 | 1 (`eGFP__hbb`) | 2.70 |
| 4 | 1 (`Q9HAE3__hbb`) | 1.40 |

Every cell with a single-pair design has a small minimum gap and every cell without one has a large
one, with no exceptions: a design blocked by a single pair carries that pair's energy and nothing
else. So the large minimum gap, the empty band and the missing single-pair designs are one cause
with three symptoms, and the named exception dissolves.

All 22 hbb cells are empty, and they split by **how**: **20 crowd the threshold**, their repairs
blocked by one immovable pair — the group the offset-95 locus speaks to — and **2 are comfortably
empty**, every design multiply obstructed, needing no explanation beyond having no near-miss
designs at all. That is a tidier account than "17 explained, 5 exceptions".

**Report the hole, not a threshold.** 20 cells land in 0.10–0.60 and 2 at 1.40 and 2.70, so **no
cell has a minimum gap between 0.60 and 1.40** and every cutoff in that interval separates
identically. Quoting a threshold of 1.0 implies a fitted boundary where there is an empty interval,
and it would make an incoming cell at 0.9 read as a near-miss on a rule rather than as the
interesting case it would be.

**The separation is perfect and its margin is one design on two cells.** `P54852__hbb` has exactly
**1** single-pair design and `P01860__hbb` has **2**; both join the 20 by that margin. The
mechanism predicts cells at the margin, so this does not weaken it, but 20/2 with no exceptions
reads as more robust than it is.

**Two limits on the relation itself.** It separates 1 from >1 and does **not** order within >1 —
three blocks gives 2.70 while four gives 1.40 — so "k pairs carries k pair-energies" is a loose
heuristic, not additive. And no cell has a fewest-block count of 2, which is a hole at n = 22
rather than a finding.

## A4 moves to 7bbdf525, and the floor arm is what licensed it

The eGFP__bnt A/B at 600 s, 2 threads, cores 18,22, with a same-binary floor arm
on the **treatment** build:

| run | points | certified | peak GiB | walked stretches | stretches |
|---|---|---|---|---|---|
| f2112884 | 230 | 230 | 5.423 | 65 | 224 |
| 7bbdf525 r1 | 234 | 234 | 5.423 | 90 | 249 |
| 7bbdf525 r2 | 238 | 238 | 5.432 | 14 | 173 |

The registered gate passes: no certified point of any arm is dominated by any point of another, and
the three fronts are **nested** — 230 ⊂ 234 ⊂ 238, with nothing exclusive to an earlier run.

**The floor arm is the reason the +4 is not reported as an effect.** The cross-build difference is
+4 points; the treatment build against *itself* is also +4. So the A/B detects no difference in
front content beyond the build's own run-to-run spread, and `7bbdf525` is recorded as **not worse**,
not as faster or better. Walked-stretch counts moved **90 → 14** between two identical runs, a
6.4× swing that confirms those counts can carry no claim; `stretches` moved 249 → 173 with them.

**A correction to the earlier note on this pin.** The delta `f2112884 → 7bbdf525` is *not* "one
error-message string" — that describes `7937f01c → 7bbdf525` (11 lines in `front.rs`, a doc
change). `f2112884 → 7bbdf525` is four commits and the whole of PR #458: 781 insertions over 7
files, including a new `could_dominate` predicate **in the certification path** and
`gap_kcal_mol`/`blocking` on the verdict event. So the switch crosses a real change to how points
are certified, which is why a floor arm was owed rather than optional.

The move is taken because `7bbdf525` is the xg tip and because `could_dominate` is documented sound
and strictly more conservative than the search's own test ("fires at least as often, and a point it
leaves uncontested really is uncontested"), so staying on `f2112884` would run the hold-out on a
superseded build for no measured gain.

## The hold-out block's hold was described but not on disk

The pre-registration states the three held-out cells are kept from the workers by a `RESERVED` note
in each of `claims-f2112884-A4/{RBD__hba,P01588__hba,P01241__hba}`. **Those directories did not
exist.** Nothing mechanical was keeping a worker off them; what happened to be protecting them was
the xg pin guard having stopped the only idle worker. `RBD__hba` sits fourth in the remaining queue
order, so an ordinary worker start would have claimed it.

Two further holes in the staged switch, both of which would have been silent:

- The new worker reads a **fresh** claims directory, so every already-run cell would have been
  re-claimed — including `RBD__bnt` while it was still running, which is two writers on one
  `--out`. The five claimed cells are now carried forward with a `CARRIED` note.
- A fresh claims directory also discards the hold-out notes, so the hold would have evaporated
  exactly at the moment of the switch it was most needed for.

The notes are now written on **both** pins, and they say what releasing them requires rather than
merely that the cell is reserved: the block's build must be settled first.

**The hold-out block runs without an xg guard, deliberately.** Every other A4 worker stops claiming
when `origin/xg` moves, which is right when the question is whether a cell should start on a
superseded build. For this block the commitment is the opposite — all three cells on one pin, a new
tip waits — so a guard would stall the block halfway and leave two cells on one build and one on
another, the exact outcome the registration exists to prevent. `run_hba_holdout.sh` carries no
fetch and no tip comparison, and says so at the top. The settled pin is recorded in
`HBA_HOLDOUT_PIN.SETTLED` with the A/B and floor-arm numbers that chose it.

## Certifying the union costs one certification, not three

`P54852__hbb`, 163 aa, hbb flanks, one binary (`scope/37-stop-dp` @ `dedac53f`), `--threads 1`,
`taskset -c 2`, five runs back to back, `--structurally-optimal --max-secs 300`:

| arm | stop | certify_s | constrained | ms/pt | seeds | points | unchecked |
|---|---|---|---|---|---|---|---|
| UAA | UAA | 175.65 | 495 | 354.9 | 0 | 0 | 0 |
| UAG | UAG | 130.57 | 482 | 270.9 | 0 | 0 | 0 |
| UGA | UGA | 209.01 | 477 | 438.2 | 303 | 314 | 0 |
| **free** | **free** | **182.02** | **465** | 391.5 | 293 | 293 | 0 |
| UAA2 | UAA | 175.12 | 495 | 354.0 | 0 | 0 | 0 |

**Count** (load-immune): free 465 against a pinned sum of 1,454 and a largest arm of 495 —
**0.320 of the sum, 0.939 of the largest arm**. **Wall:** `certify_s` 182.02 against a pinned sum of
515.24 — **0.353**, and 0.87× the largest single arm. Both inside the predictions registered before
the run (count 0.9–1.3× of one arm; wall 0.9–1.5× of one arm, i.e. 0.3–0.5 of the sum).

`unchecked == 0` on every arm, so each `certify_s` is the cost of certifying that arm's **whole**
front and the sum is meaningful. An arm stopped partway through certification would report a
truncated reading that looks ordinary; the reader refuses to print a ratio if any arm shows
`unchecked > 0`.

**The free arm's 465 points reproduce, by an independent route, the Pareto union of the three
pinned arms computed earlier on a different build** (465 from 1,454 candidates, UGA 293, UAG 172,
UAA 0). Two unrelated computations of the same front agree exactly.

**The floor pair reclassifies the spread from noise to signal.** UAA ran first and again last:
175.65 and 175.12 s, a **0.3 %** spread. So the instrument is reproducible to 0.3 % on this cell at
one thread, and the **1.60×** spread among the three pinned arms (130.57–209.01 s over fronts
differing by 3.8 %) is therefore a **real stop-dependent effect, not run-to-run variation**. An
earlier reading of that 1.60× as the instrument's noise floor was wrong in the conservative
direction: the error bar on the ratio is ~0.3 %, not 1.60×, so the ~2.8× saving is far better
supported than the cautious version claimed.

**Certification cost is a property of the cell, not of the method**, so the ratio is the quotable
quantity and the absolute seconds are not. Per constrained point, across the grid's six existing A4
rows on `f2112884`, all `unchecked == 0`: 28.69, 30.16, 30.57, 31.13, 31.99 and 57.75 ms/pt — a
**2.01×** spread at comparable front size. `P54852__hbb` runs 271–438 ms/pt, about **10×** the
seeded cells.

**Why a seed-rich front certifies slower, and it is not what it first looks like.** `certify` is not
folding every point: it carries a 64-entry pool of free structures taken from points already found
**not** to be seeds, and rules later points out *unfolded* by scoring them against that pool. The
pool is populated only by failures, so a seedless front fills it immediately and skips most folds
while a seed-rich front keeps it thin and pays a full fold nearly every time. The front is visited
in **bit-reversed** order (`so_front.rs:spread_order`), so the pool holds structures spread across
the whole front rather than a local neighbourhood — which is why a handful of folds can rule out
most of the rest.

So the claim is **not** "certification cost tracks the seed count" but "it tracks how many points
must actually be folded", and the falsifiable form is that `certify_s / (points − ruled_out)`
should be flat across arms where `certify_s / points` spans 1.60×. `DERNA_SO_CENSUS=1` prints the
ruled-out count to stderr; census re-runs are queued. The pool is a `Mutex<Vec<String>>` shared
across workers and each worker takes every *n*-th element of the visit order, so **the ruled-out
count is deterministic only at `--threads 1`** and is not comparable across thread counts.

**Dimensions of the claim**, stated as part of it: one cell, one protein, one flank pair, one
length, one thread, **seedless on two of three pinned stops**, on a branch whose structurally
optimal search differs from `xg`'s. A second cell (`eGFP__bnt`, 239 aa, bnt flanks, seeded on every
arm) is queued with the census on.
