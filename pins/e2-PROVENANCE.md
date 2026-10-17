# Frozen binaries for the front-enumeration benchmark (E2), 2026-10-16

Every derna2 binary in this directory is built from one commit. The comparator
binaries are not rebuilt: they are the E1 freeze at `/scratch/xinyugu4/paper-e1-bin`,
whose hashes are re-verified below and match that directory's own PROVENANCE.md.

## derna2, pinned at `acfaae60` (origin/xg tip)

Source worktree `/scratch/xinyugu4/paper-e2-src` (detached at `acfaae60`),
cargo 1.95.0 release, `CARGO_TARGET_DIR=/scratch/xinyugu4/paper-e2-target`.

| binary | crate | target |
|---|---|---|
| `scope35_walk_front` | `packages/derna-beam` | `--bin scope35_walk_front` |
| `run_design_per_lambda` | `packages/derna-beam` | `--bin run_design_per_lambda` |
| `mfe_fold` | `packages/rna-mfe` | `--example mfe_fold` |

    c9c287049a46a0e1bd83f4f990568ff34c5cb3d3058fdbf7432e059e4189451e  scope35_walk_front
    4b1d23da0f8e7ce0ff65a8a49450c07c3fdff67c1cd4549eab27ae031a928584  run_design_per_lambda
    88425013dfc959293226c27b157925eb83369e382af8facd8ebc1b8d3695a60c  mfe_fold

The pin moved from `d3aa667a` to `acfaae60` when `experiment/walk-at-length`
merged into xg (#414): that merge carries the two walk levers the earlier pin
could not include, so the plan's own re-pin rule applied. `mfe_fold` rebuilt
byte-identical across the two, so everything re-scored under either pin is
comparable; the other two binaries changed and every row records which pin
produced it.

E1 pinned derna2 at `17ee9d8b`, so E2's derna2 is a later commit than E1's and
the two experiments' derna2 numbers are not interchangeable.

### What the pin does and does not contain

`--fold-by-width` and `--slice-width 250` (`f54cb1e4`, measured at 0.69x
instructions at n = 1000) are in this pin, and are the walk's **defaults** --
`--no-fold-by-width` and `--no-slice-width` opt out. The arm therefore runs the
walk as it actually ships, not a slower stand-in.

Still absent: the hull-walk anchor schedule (`ad963f76`,
`experiment/motif-hull-anchors`), which is unmerged. It only makes the walk
cheaper, so the arm remains an upper bound by that one lever. If it merges
before the grid runs, re-pin and rebuild rather than mixing commits.

## DERNA v1

Upstream C++ at `docs/reference/_cpp/derna-main`. Its CMake wants cmake >= 3.1
and C++17; this host has cmake 2.8.12 and g++ 4.8.5, so it is compiled directly
over CMakeLists' source list with
`/scratch/xinyugu4/cxx-env/bin/x86_64-conda-linux-gnu-g++ -std=c++17 -O3`
(120 warnings, no errors) by `build_derna1.sh`.

    940ff6f0c494a95484373eb859e22bbdd5529d14c0666005fc002f7a4a1a103c  derna1
    8b5925445412f5e6bbae503eb1a662bde30ef5fadcf25d6098ad51982a117127  derna1_codon_usage_human.csv

The CSV is the human codon table converted into v1's two-line-per-amino-acid
`-c` format. v1's shipped `sample_codon_usage.csv` is not a human table, and
passing it would have had v1 optimizing a different objective than every other
arm. v1 normalizes per amino acid (`codon_cai = log(usage / max_c)`), so only
the within-amino-acid ratios matter.

## Comparators, from the E1 freeze (hashes re-verified 2026-10-16)

    16a5a46cf3a761dbf1158cb1c357a28af7daa082bd4b13b192ce874913dee5dc  LinearDesign/bin/LinearDesign_2D
    0c25256f09809495933e82d02c392973cb0c511f0b3e4f4b1c17ba17f2d1e64d  LinearDesign/src/Utils/libraries/LinearDesign_linux64_old.so
    e96a8e2ebeba05a91c0887f0b9bb6132906fdd828d48d8f3d8bd70c10126ad2e  LinearCDSfold/LinearCDSfold
    ba06279c21961c018cc4036a36aef546fb83db9c5a08fced0b3f6bab7bc6efd9  CDSfold/CDSfold
    77cbb78732ca153e99ac6cca52c5b2cf0e19c2e81e3195a25e2c3f43356dbbfc  mrnafold_fold_codon_graph
    05b4a223f49c94ce62707db645e51adbfaa1ca35f025c2e9c17197382ea78dbf  mrnafold_fold_codon_graph_tbb

## The weight these binaries are handed is not the weight they are given

DERNA v1 and LinearCDSfold score `lambda*E_dcal - (1 - lambda)*sum(ln w)` where
derna2 carries a factor 100 on the adaptation term, so at one lambda they weight
CAI 100x less. Each is handed `1/(1 + 100*(1 - lambda)/lambda)`; LinearDesign and
mRNAfold take plain `(1 - lambda)/lambda`. Verified against derna2 at interior
lambdas, not at corners: at a corner the second objective only breaks ties, so a
corner check cannot see this at all.
