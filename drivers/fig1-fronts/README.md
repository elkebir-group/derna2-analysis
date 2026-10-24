# GLPE_HUMAN with hbb UTRs: three fronts and mRNAfold

## The question

UTR inclusion and the exact Pareto front have so far been measured apart. This
puts both on one real instance and draws them in one picture:

- **Front 1**: the exact unguarded front of the whole molecule, UTRs attached.
- **Front 2**: the exact guarded front, where a region guard forbids classes of
  pairs in the *designed* structure. One per guard.
- **Front 3**: structurally optimal designs. A design is structurally optimal
  when its free MFE equals its guarded energy (`F = E`), i.e. the molecule's
  own fold respects the guard. No finished method computes the whole
  structurally optimal front, so front 3 is one lazy k-best walk per **on-hull
  vertex of front 2**, at the midpoint of that vertex's weight interval, to the
  first structurally optimal design (2 h cap per walk).
- **mRNAfold**, the closest competitor (exact DP, fixed UTRs, a CAI weight, soft
  penalties on the designed structure), at the same weights in four arms.

## Instance

GLPE_HUMAN (P15421), first 78 residues
([fixture](../scope-25-beam-runtime-quality/fixtures/GLPE_HUMAN_78.fasta)), with
the hbb flanks from [`scope-19-utr-support/m5`](../scope-19-utr-support/m5/):
50 nt 5′, 132 nt 3′, 416 nt in all. Positions (0-based): 5′ UTR 0–49, CDS
50–283, 3′ UTR 284–415. Turner 2004, dangles 0, internal-loop cap 30.

| arm | `scope35_walk_front --region` | `lazy_kbest_probe --table` |
|---|---|---|
| symmetric | `5-cds,cds-3` | `symmetric` |
| 5′ no_clamp | `5-cds,5-3` | `five-prime+no-clamp` |
| 5′ unpaired | `5-5,5-cds,5-3` | `five-prime+unpaired` |

The `five-prime+unpaired` table (`RegionConstraints::five_prime_unpaired`) was
added for this, with a walk test that its designs never pair a 5′ position
while the unguarded walk on the same case does.

The objective is `λ·E − (1 − λ)·L` with `E` in kcal/mol and `L = Σ ln w`, the
accumulated log-CAI (the walk reports it ×100). Hulls are taken in `L`, not CAI.

## Results

### Tables

#### Fronts

| front | guard | points | hull vertices | SO points | SO hull vertices (rank 1) | wall s |
|---|---|---|---|---|---|---|
| 1 | unguarded (`none`) | 114 | 23 | – | – | 252 |
| 2 | symmetric (`5-cds,cds-3`) | 117 | 22 | 0 | 0 | 178 |
| 2 | no_clamp (`5-cds,5-3`) | 119 | 22 | 42 | 7 | 223 |
| 2 | unpaired (`5-5,5-cds,5-3`) | 119 | 22 | 0 | 0 | 219 |

#### Front 3

**symmetric**

| vertex | λ | vertex E | found | rank | E = F | ΔE vs vertex | walk s | peak GB | start load |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.9378 | -180.60 | yes | 336 | -180.30 | +0.30 | 6.8 | 0.5 | 6.13 |
| 1 | 0.7607 | -180.20 | yes | 94 | -179.90 | +0.30 | 3.8 | 0.5 | 6.13 |
| 2 | 0.6345 | -179.40 | yes | 209 | -179.10 | +0.30 | 5.4 | 0.6 | 6.13 |
| 3 | 0.5745 | -178.00 | yes | 185 | -177.00 | +1.00 | 5.1 | 0.6 | 6.80 |
| 4 | 0.5205 | -176.60 | yes | 233 | -174.20 | +2.40 | 6.3 | 0.5 | 6.80 |
| 5 | 0.5099 | -175.30 | yes | 258 | -174.20 | +1.10 | 6.8 | 0.5 | 6.80 |
| 6 | 0.4444 | -171.60 | yes | 1779 | -172.80 | -1.20 | 33.2 | 0.7 | 7.06 |
| 7 | 0.3766 | -166.60 | yes | 54 | -166.60 | +0.00 | 3.2 | 0.6 | 7.29 |
| 8 | 0.3591 | -165.70 | yes | 35 | -166.60 | -0.90 | 2.6 | 0.6 | 7.29 |
| 9 | 0.344 | -164.50 | yes | 70 | -165.40 | -0.90 | 3.7 | 0.6 | 7.59 |
| 10 | 0.2983 | -161.20 | yes | 1060 | -165.40 | -4.20 | 22.8 | 0.7 | 7.78 |
| 11 | 0.2533 | -159.40 | yes | 23492 | -163.60 | -4.20 | 464.1 | 0.7 | 7.88 |
| 12 | 0.199 | -157.10 | yes | 95035 | -152.70 | +4.40 | 1992.6 | 1.7 | 7.64 |
| 13 | 0.1322 | -153.30 | yes | 32689 | -150.70 | +2.60 | 710.7 | 1.9 | 7.64 |
| 14 | 0.1131 | -151.10 | yes | 31024 | -142.90 | +8.20 | 710.0 | 2.8 | 7.69 |
| 15 | 0.112 | -139.20 | yes | 31195 | -142.90 | -3.70 | 745.8 | 2.8 | 6.96 |
| 16 | 0.09556 | -136.80 | yes | 30144 | -134.60 | +2.20 | 761.6 | 2.8 | 10.57 |
| 17 | 0.06392 | -132.60 | yes | 22099 | -130.40 | +2.20 | 568.3 | 3.8 | 14.34 |
| 18 | 0.03844 | -130.90 | yes | 35044 | -128.70 | +2.20 | 995.5 | 11.0 | 10.12 |
| 19 | 0.02729 | -124.90 | yes | 61921 | -123.90 | +1.00 | 1838.1 | 24.1 | 8.02 |
| 20 | 0.0197 | -123.10 | yes | 66128 | -123.90 | -0.80 | 2051.2 | 40.3 | 7.92 |
| 21 | 0.006557 | -119.40 | no (seconds) | 19057 | – | – | 7200.0 | 139.4 | 16.95 |

**no_clamp**

| vertex | λ | vertex E | found | rank | E = F | ΔE vs vertex | walk s | peak GB | start load |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.9378 | -191.80 | yes | 1 | -191.80 | +0.00 | 1.9 | 0.6 | 20.15 |
| 1 | 0.7607 | -191.40 | yes | 1 | -191.40 | +0.00 | 2.1 | 0.6 | 20.15 |
| 2 | 0.6345 | -190.60 | yes | 1 | -190.60 | +0.00 | 2.0 | 0.6 | 20.57 |
| 3 | 0.5745 | -189.20 | yes | 1 | -189.20 | +0.00 | 2.1 | 0.6 | 20.57 |
| 4 | 0.5023 | -187.80 | yes | 1 | -187.80 | +0.00 | 1.9 | 0.6 | 20.76 |
| 5 | 0.4577 | -186.40 | yes | 1 | -186.40 | +0.00 | 1.8 | 0.5 | 20.76 |
| 6 | 0.3956 | -180.20 | yes | 1 | -180.20 | +0.00 | 2.4 | 0.7 | 22.07 |
| 7 | 0.3529 | -176.70 | yes | 3 | -180.20 | -3.50 | 2.5 | 0.7 | 22.07 |
| 8 | 0.3073 | -175.50 | yes | 63 | -179.00 | -3.50 | 3.5 | 0.6 | 22.30 |
| 9 | 0.2637 | -174.70 | yes | 2224 | -179.00 | -4.30 | 47.2 | 0.5 | 22.30 |
| 10 | 0.2332 | -166.50 | yes | 5890 | -167.00 | -0.50 | 124.1 | 0.8 | 22.10 |
| 11 | 0.2022 | -162.20 | yes | 9089 | -158.40 | +3.80 | 209.1 | 0.8 | 22.23 |
| 12 | 0.1835 | -160.70 | yes | 5544 | -158.40 | +2.30 | 124.5 | 0.8 | 20.03 |
| 13 | 0.1548 | -157.90 | yes | 8320 | -158.40 | -0.50 | 187.1 | 1.0 | 20.37 |
| 14 | 0.1372 | -155.10 | yes | 26970 | -158.40 | -3.30 | 726.4 | 1.4 | 20.37 |
| 15 | 0.1213 | -150.30 | yes | 28976 | -146.10 | +4.20 | 860.5 | 1.6 | 25.71 |
| 16 | 0.1016 | -147.90 | yes | 16241 | -143.70 | +4.20 | 490.7 | 1.5 | 44.50 |
| 17 | 0.08517 | -145.90 | yes | 16932 | -143.70 | +2.20 | 433.0 | 1.7 | 63.93 |
| 18 | 0.07527 | -141.70 | yes | 19578 | -139.50 | +2.20 | 498.3 | 2.8 | 27.14 |
| 19 | 0.06207 | -137.70 | yes | 21664 | -134.90 | +2.80 | 574.3 | 3.0 | 24.63 |
| 20 | 0.03843 | -133.70 | yes | 8690 | -126.20 | +7.50 | 256.6 | 7.1 | 33.91 |
| 21 | 0.0121 | -128.20 | yes | 1209 | -126.20 | +2.00 | 235.5 | 40.8 | 29.30 |

**unpaired**

| vertex | λ | vertex E | found | rank | E = F | ΔE vs vertex | walk s | peak GB | start load |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.9378 | -182.40 | no (seconds) | 436544 | – | – | 7200.2 | 1.7 | 29.50 |

#### mRNAfold

| arm | counterpart | runs | objective = ours | same sequence | max Δobj | SO symmetric | SO no_clamp | SO 5′ unpaired | median s |
|---|---|---|---|---|---|---|---|---|---|
| plain | unguarded | 22 | 22 | 21 | 0.000 | 0 | 6 | 0 | 3.9 |
| guard | symmetric | 22 | 22 | 18 | 0.000 | 0 | 5 | 0 | 4.0 |
| unpaired5 | unpaired | 22 | 22 | 21 | 0.000 | 0 | 9 | 0 | 4.4 |
| unpaired_both | unpaired | 22 | 3 | 3 | 6.244 | 0 | 5 | 0 | 4.1 |

### Fronts

Every check passes on every front: no structure on a guarded front has a pair
its guard forbids; every guarded front differs from front 1 (the symmetric and
5′-unpaired fronts share no point with it, no_clamp shares 34); the hull
recomputed here in log-CAI matches the walk's own; and every hull vertex refolds
to exactly its front energy, under its guard or free.

The guard is cheap where adaptation is free and dear where it is not. Under the
symmetric guard the most stable design loses 11.2 kcal/mol (−191.8 → −180.6) and
the most adapted 10.4 (−129.8 → −119.4). Under no_clamp the stable half of the
front is front 1's own: the first six hull vertices are front 1's vertices, the
same designs at the same energy, and with the seventh they are already
structurally optimal, which is why those seven walks answer at rank 1. Further
along, a design can be a vertex of both fronts at different energies: the guard
moves its energy up by the 1.6 kcal/mol its free fold gains from a forbidden
pair.

### Front 3

The walk's rank-1 design is the hull vertex at every weight (checked: each
walk's optimum equals the vertex objective to 1e-4), and every design it
returned refolds as structurally optimal with the objective it reported
([`verify_front3.py`](verify_front3.py)).

- **Symmetric**: no hull vertex is structurally optimal (gaps 1.0–10.4 kcal/mol),
  and 21 of 22 walks found one, at ranks 35–95 035 and 2.6–1 993 s. The
  structurally optimal design is often *more* stable than the vertex and less
  adapted (ΔE −4.2 to +8.2 kcal/mol). The 22nd, vertex 21 at λ = 0.0066 (the
  most adapted end), stopped at the 2 h cap after 19 057 designs and 139 GB:
  every structurally optimal design there has objective ≥ −34.99 against the
  vertex's −78.30. Of its 7 200 s, 7 187 went to the walk's own next-derivation
  step and 0.2 s to checking designs, so at that end the walk, not the check,
  is the limit.
- **5′ no_clamp**: 22 of 22 found; 7 at rank 1, the rest within 29 000 designs and
  15 min.
- **5′ unpaired: empty by construction.** hbb's 5′ UTR folds a −9.4 kcal/mol
  hairpin inside positions 0–49. A guarded structure has no pair touching the
  5′ UTR, so that hairpin can always be added to it without crossing a pair, and
  at dangles 0 the exterior loop is additive: every sequence has
  `F ≤ E − 9.4`, and none is structurally optimal under this guard. The guarded
  hull vertices' gaps are 9.4 or 11.0, never less. One control walk was run at
  vertex 0: it visited 436 544 designs in 2 h, the windowed check rejected
  every one without a full fold, and none was structurally optimal. It covered
  301 objective units above the optimum, of the ≥ 0.9378 × 940 ≈ 882 the bound
  says lie below any structurally optimal design at that weight, so the walk
  could not have answered at any cap a run would afford.

Walk cost rises steeply toward low λ, where adaptation dominates and many
near-optimal designs differ only in codons the free fold reorganizes: ranks go
from tens to ~10⁵ and resident memory from 0.5 GB to 40+ GB.

### mRNAfold

Run at the symmetric front's 22 hull weights with `λ_mf = (1 − λ)/λ`,
one trace, our codon weights as `codon_freq_path`. Arms:

- **plain**: no constraint;
- **guard emulation**: `discourage_pair i,j` on every 5′↔CDS and CDS↔3′ pair, both
  orientations (a check on a 12-residue case shows indices are 0-based and only
  `i < j` acts: `discourage_pair 16,21` removes the pair, `21,16` alone does
  not);
- **5′ unpaired**: `encourage_unpaired 0,49`;
- **both UTRs unpaired**: `encourage_unpaired 0,49` and `encourage_unpaired 284,415`.

Each penalty is +100 kcal/mol per pair on the designed structure.

mRNAfold reproduces our optimum at all 22 weights in every arm with a
counterpart: plain = our unguarded front, guard emulation = our symmetric
front, 5′ unpaired = our 5′-unpaired front (objective equal to 1e-3; sequences
identical except at ties). No designed structure in a guarded arm carries a
forbidden pair, and the energy of each designed structure equals our refold
under the matching guard. The both-UTRs arm has no counterpart; its designs sit
22–31 kcal/mol above our 5′-unpaired refold, the price of the 3′ UTR's own
pairs.

**No mRNAfold design is structurally optimal under the symmetric or the
5′-unpaired guard.** Its soft penalties act on the designed structure only; the
molecule's own fold is not consulted. That is the gap front 3 measures.

Cost: ~4 s per mRNAfold weight; 3–4 min per whole exact front (a busy host).

## Plot

Interactive: `results/raw/plot/glpe-hbb-fronts.html` (data embedded). Static:
[`results/plot/fronts.png`](results/plot/fronts.png). Front 1 circles,
front 2 diamonds (5′ unpaired: squares), front 3 triangles (5′ unpaired: ▽),
hollow off-hull or unfinished; mRNAfold `+` plain, `×` guard, hollow square 5′
unpaired, hollow circle both UTRs. Each panel uses three colour slots (front 1,
2, 3) so the set passes an all-pairs colour-vision check; arms are told apart by
shape. Each design is drawn at the energy its method optimized.

## Caveats

- Walls are from a shared host. Each walk's 1-minute load at start is in
  `front3/v<k>.load` and in the table; loads ran 6–64. Ranks, objectives and
  designs do not depend on load.
- Front 3 is one walk per hull vertex, not the structurally optimal front:
  it says how far below each guarded vertex the first structurally optimal
  design lies, not which structurally optimal designs are Pareto-optimal among
  themselves.

## Reproduce

Binaries from this branch, frozen before any run:
`scope35_walk_front`, `lazy_kbest_probe`, `glpe_hbb_rescore`.

```text
run_fronts.sh ROOT [ARM...]                 # fronts 1-2
analyze_fronts.py --root ROOT                # checks, hull, vertices.jsonl
glpe_hbb_rescore --utr5-nt 50 --cds-nt 234 --in X.jsonl --out X.scored.jsonl
summarize_vertices.py ROOT                   # rank-1 structural optimality
run_front3.sh ROOT ARM [JOBS]                # front 3 (ONLY="k ..." to pick vertices)
collect_front3.py ROOT; verify_front3.py ROOT ARM...
run_mrnafold.py --precheck ...; run_mrnafold.py --vertices ROOT/symmetric/vertices.jsonl ...
compare_mrnafold.py --runs ... --fronts ROOT --out compare.jsonl
build_plot.py ROOT compare.jsonl data.json fronts.png page.html
make_tables.py ROOT compare.jsonl
```

Small results are in [`results/`](results/).
