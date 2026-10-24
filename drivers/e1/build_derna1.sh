#!/bin/bash
# Build DERNA v1 into the E1 frozen binary directory.
#
# The repo ships a CMakeLists requiring cmake >= 3.1 and C++17; this host has
# cmake 2.8.12 and g++ 4.8.5, so neither works. The project is a single
# add_executable over a fixed source list, so it is compiled directly with the
# same conda g++ E1 used for mRNAfold.
set -euo pipefail

SRC=/home/xinyugu4/derna2/docs/reference/_cpp/derna-main
OUT=/scratch/xinyugu4/paper-e2-bin
CXX=/scratch/xinyugu4/cxx-env/bin/x86_64-conda-linux-gnu-g++

mkdir -p "$OUT"

"$CXX" --version | head -1

# The source list is CMakeLists' add_executable, headers dropped.
"$CXX" -std=c++17 -O3 -Wall -Wextra -Wno-uninitialized \
    -o "$OUT/derna1" \
    "$SRC/src/main.cpp" \
    "$SRC/src/utils.cpp" \
    "$SRC/src/Nussinov.cpp" \
    "$SRC/src/NussinovAlgorithm.cpp" \
    "$SRC/src/Zuker.cpp" \
    "$SRC/src/default.cpp" \
    "$SRC/src/ZukerAlgorithm.cpp"

echo "BUILD OK"
ls -l "$OUT/derna1"
sha256sum "$OUT/derna1"
"$OUT/derna1" -h 2>&1 | head -20 || true
