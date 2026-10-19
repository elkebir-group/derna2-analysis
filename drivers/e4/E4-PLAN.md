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
| `A1` concatenation | design the CDS alone, paste the UTRs on afterwards — the practice the paper's intro describes | `derna front --protein P` (no UTRs), then re-score each `v` as `u5·v·u3` |
| `A2` constrained | the UTRs inside the DP, every candidate forming a forbidden pair rejected | `derna front --protein P --utr5 U5 --utr3 U3 --forbid 5-cds,cds-3` |
| `A3` mRNAfold | the comparator, UTRs in its energy model, nothing rejecting CDS–UTR pairs, run at `A2`'s hull λs | `mrnafold_fold_codon_graph` with `utr_5p`/`utr_3p`, λ' = (1−λ)/λ |
| `A4` structurally optimal | the answer the problem asks for; partial under a budget, each point `certified` or `provisional` | `derna front … --structurally-optimal --max-secs 43200 --max-rss-gb 48` |

`A4 ⊆ A2` as design sets, by construction: structural optimality is an extra test
applied to the constrained front. `A1` and `A3` are **not** fronts of structurally
optimal designs and are not claimed to be; what they contribute is the headline
number, the share of each arm's points whose own free fold respects the constraint.

## The constraint is the two crossing classes, not three

`derna front`'s default `--forbid` is *every class whose flanks are given*, which
includes `5-3`. The paper's problem statement explicitly allows 5′–3′ pairs (the two
ends of an mRNA are intrinsically close and do pair). **E4's headline constraint is
therefore `--forbid 5-cds,cds-3`, passed explicitly on every `A2` and `A4` cell.**

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

So the three-class setting stays in the grid as a **control**, one instance only: it
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

## Instances

**The long instances of E3, at native length** — addressed by ID as in E2 and E3, so
no ortholog name reaches a row:

| ID | length (aa) |
|---|---:|
| `PRT-03` … `PRT-07` | 984–1629 |

Native length means each protein's own full coding sequence; no truncation rungs. The
consequence is stated up front rather than discovered: **at 984–1629 aa most cells
will be partial.** E3's hull probe projects ≈15 GB at 984 aa and ≈40–45 GB at 1629 aa
serial for the motif hull, and the *front* is the hull plus its off-hull gap solves,
so a complete exact front at the top of this range is not expected inside 12 h × 48 GB.
Every arm therefore has to degrade gracefully, which is what the run order below is
built around; a partial cell is a result, reported as partial, never as a front.

A short **companion panel** is run alongside at the same four arms so the paper has
at least one figure of four *complete* curves: the six therapeutic natives already
screened in E3 (EPO `P01588`, G-CSF `P09919`, hGH `P01241`, IFN-α2 `P01563`, eGFP,
SARS-CoV-2 RBD; 188–239 aa, inside the lengths Scope 37 ran). The long panel carries
the claim about real transcripts; the companion panel carries the claim about the
complete front. They are reported as two panels and never pooled.

## UTR pairs — what is actually used in practice

Five pairs. Three are the UTRs of real constructs, one is the classic IVT pair the
paper's intro already uses, and one is a deliberately unstructured minimal pair that
functions as the control.

| key | 5′ | 3′ | why it is in the grid |
|---|---|---|---|
| `hbb` | human β-globin 5′ UTR | human β-globin 3′ UTR | the paper's own anchor (intro, Fig. 1); the standard IVT pair in kits and preclinical constructs; already pinned in the repo (`experiments/scope-19-utr-support/m5/hbb5.txt`, `hbb3.txt`) and characterized by Scope 37 |
| `hba` | human α-globin 5′ UTR | human α-globin 3′ UTR | the other globin pair in clinical use; `HBA1` MANE Select |
| `bnt` | α-globin-derived 5′ UTR with an optimized Kozak | the AES + mtRNR1 ("F") 3′ element | the UTRs of an **approved** mRNA vaccine (BNT162b2); the most-deployed pair in existence |
| `mrna1273` | the construct's 5′ UTR | its β-globin-derived 3′ element | the UTRs of a second approved construct, chosen independently of `bnt` |
| `min` | a short low-structure 5′ UTR ending in a Kozak | a short low-structure 3′ UTR | used in practice as a minimal UTR, and here the control: a pair with little structure of its own should forbid least and leave `A4` closest to `A2` |

**Sequences are not quoted in this plan.** Each is fetched from a cited public source
(MANE Select RefSeq records for `hbb`/`hba`; the published regulatory/INN sequence
filings for `bnt` and `mrna1273`), written into one `utrs.fasta` with a per-record
source line, and pinned by SHA-256 in `PROVENANCE.md` before any run. A pair whose
sequence cannot be sourced to a public document is **dropped, not reconstructed**.

Two decisions about the molecule's edges, fixed now:

- **No poly(A) tail in the headline cells.** The tail is added enzymatically, its
  length is not a design choice, and including 100 nt of A would put a large
  U-complementary partner in the 3′ region — a real effect, but one that would be read
  as a property of the UTR pair. One **cross-check cell** appends `A100` to `hbb`'s 3′
  UTR on one companion-panel instance, reported separately.
- **Stop codon included** (`--stop UAA`), coding to the constraint and outside CAI, as
  the CLI already defines it: the real molecule has one, and it sits exactly at the
  CDS–3′ boundary where the constraint bites.

Grid size: 5 long + 6 companion = 11 proteins × 5 UTR pairs = 55 instances. `A1` is
**per protein, not per instance** (one no-UTR front, re-scored against each pair), so
the walks cost 11 + 55 (`A2`) + 55 (`A4`) cells, with `A3` cheap per λ.

## Run order

Each step's output gates the next; the gate is written here so it cannot drift.

1. **Pin.** One commit, one frozen bin directory, `SHA256SUMS` + `PIN`, as E2 does.
   `derna` must be at or past the commit whose `front` is the seeded walk
   (`8e588452`), and `mrnafold_fold_codon_graph` must be the build whose `utr_5p` /
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

`A2` runs with `--stream`, so a walk stopped by the clock or the cap still leaves its
hull and its finished gaps on disk (E2's `replay_partial_hull` is the same mechanism).
`A3` then takes **the hull λs that actually completed**, and `A4` runs on whatever
constrained front `A2` produced, with its own 12 h budget on top. A cell therefore
reports a quadruple of completeness verdicts, not a single pass/fail.

## Budget and host

- **12 h wall per cell**, 48 GB per cell under `systemd-run --user --scope -p
  MemoryMax=48G -p MemorySwapMax=0` (the 52 GB of E2's note, less headroom for the
  re-scoring that follows each cell). A cell over either limit records `partial` with
  its open regions, or `oom`; neither is censored.
- **`--threads 1` for everything timed.** The long cells are memory-bound before they
  are core-bound, and E2 measured the multi-core arms' interference; a threaded arm
  here would buy wall at the cost of the one resource that decides these cells.
- **E4 cannot start while E2 is timing on this host.** E2 holds tamazirt (61 GB, P-cores
  0–15, E-cores 16–27) and its `par` phase uses 8 cores. E4 runs after E2's `par`
  phase finishes, or on the collaborator's host. Until then only steps 2–4 run, which
  are small and can go on CPUs 16–27.

At 12 h × (55 `A2` + 55 `A4`) the long panel is weeks of wall in the worst case. That
is accepted on purpose: the cells that matter most are the ones that do not finish,
and their partiality is the measurement. The run order is longest-first inside each
panel so the information arrives before the budget does.

## Gate results — run 2026-10-19, before any cell

All three prerequisites passed, on the E-cores (16–27) beside the other sessions'
jobs; these are correctness checks, so contention does not affect them.

**Step 2 — two-class exactness.** `guarded_front` 6/6 green in 149 s, the enumeration
case taking the bulk of it. The two-class constraint is covered by the repo's oracle
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
So the three-class default carries no design the concatenation arm does not, on the
real UTR pair and through the shipped CLI — the reason E4's headline constraint is
`--forbid 5-cds,cds-3` and not the default.


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
   points that make the price smaller.
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
  establishes it at 1629. This is the standing caveat on every exact claim here, and
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

None in the code: `derna front --structurally-optimal --max-secs --max-rss-gb` and
`--forbid` all ship today, and `--stream` for partial cells exists on E2's branch
(`scope35_walk_front --stream`, `run_glpe78 --stream`). What E4 needs that does not
exist yet is the driver (modelled on `drive_e2.py`), the `utrs.fasta` pin, and the
two-class exactness check of step 2.
