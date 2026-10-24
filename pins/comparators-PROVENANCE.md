# The comparator binaries: where each came from, 2026-10-11

Every method derna2 is measured against, as built for the benchmark. These are frozen
once and never rebuilt, so a row identifies a binary rather than a recipe to repeat; the
checksums are of the files the runs used. derna2's own two binaries are not here — they
are built from source at the commit `derna2.pin` names.

| tool | source | build |
|---|---|---|
| LinearDesign_2D | LinearDesignSoftware/LinearDesign tag v1.0.0 (`5a635d1`) | system g++ 4.8.5, links the bundled `LinearDesign_linux64_old.so` (the primary `.so` needs glibc 2.29; the host has 2.17). Run with cwd = `LinearDesign/` |
| LinearCDSfold | LinearCDSfold `de1db7e` | prebuilt binary copied as found |
| CDSfold | gterai/CDSfold `06f3ee8` | rebuilt from source (the repository's shipped binary discarded), g++ 4.8.5 -O3, against ViennaRNA 2.1.1 (git tag v2.1.1) static `libRNA.a`, MAXLOOP 30 |
| mRNAfold | maxhwardg/mRNAfold `d5544e4`, `cpp_src/exe/src/fold_codon_graph.cpp` | conda `x86_64-conda-linux-gnu-g++ -std=c++20 -O3` (system gcc 4.8 is too old) |

    16a5a46cf3a761dbf1158cb1c357a28af7daa082bd4b13b192ce874913dee5dc  LinearDesign/bin/LinearDesign_2D
    0c25256f09809495933e82d02c392973cb0c511f0b3e4f4b1c17ba17f2d1e64d  LinearDesign/src/Utils/libraries/LinearDesign_linux64_old.so
    e96a8e2ebeba05a91c0887f0b9bb6132906fdd828d48d8f3d8bd70c10126ad2e  LinearCDSfold/LinearCDSfold
    ba06279c21961c018cc4036a36aef546fb83db9c5a08fced0b3f6bab7bc6efd9  CDSfold/CDSfold
    77cbb78732ca153e99ac6cca52c5b2cf0e19c2e81e3195a25e2c3f43356dbbfc  mrnafold_fold_codon_graph
    05b4a223f49c94ce62707db645e51adbfaa1ca35f025c2e9c17197382ea78dbf  mrnafold_fold_codon_graph_tbb

**Two mRNAfold builds, and only one of them is parallel.** `std::execution::par_unseq` is
serial in libstdc++ without TBB, so the plain build ignores its own `parallel true`. The
`_tbb` build links libtbb and is the one whose threaded arm actually runs threaded; an arm
that means to measure mRNAfold on several cores has to use it.

## The two binaries that produced the deposit

The deposited rows of the whole-front experiment were made by binaries that predate the
`derna` command line, and are recorded here because no commit of the public repository
rebuilds them:

    577bd177acb5fa0e0233fe58ff2dc1ab70c8f93c1937cad7145235118b3640c5  run_design_per_lambda
    dea62e117c588a3019309b5ce63d77327e20901a4913e65ff288777d89c388fb  mfe_fold

`run_design_per_lambda` designed one weight at a time, which `derna design` does; the
`mfe_fold` beside it is the same free fold of one sequence the pipeline still runs. A
re-run of the experiment uses the command line, so its rows carry a different build stamp
(`sha_derna`) from the deposit's.
