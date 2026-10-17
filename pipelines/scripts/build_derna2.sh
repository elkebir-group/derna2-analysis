#!/usr/bin/env bash
# Build derna2's benchmark binaries at one commit into a directory of their own.
#
#   build_derna2.sh REPO_URL PIN OUT_DIR
#
# The source is a fresh clone checked out detached at PIN, and the cargo target directory
# sits beside it, so nothing is shared with a working checkout whose HEAD may move.
# Writes OUT_DIR/SHA256SUMS and OUT_DIR/COMMIT; the caller compares them with pins/.
set -euo pipefail

repo=$1
pin=$2
out=$3
build=$(dirname "$out")/derna2-src-$pin

if [ ! -d "$build/.git" ]; then
    git clone --quiet "$repo" "$build"
fi
git -C "$build" fetch --quiet origin
git -C "$build" checkout --quiet --detach "$pin"

export CARGO_TARGET_DIR=$build/target
export PATH=$HOME/.cargo/bin:$PATH
cargo build --release --quiet --manifest-path "$build/packages/derna-beam/Cargo.toml" \
    --bin run_design_per_lambda --bin scope35_walk_front
cargo build --release --quiet --manifest-path "$build/packages/rna-mfe/Cargo.toml" \
    --example mfe_fold

mkdir -p "$out"
cp "$CARGO_TARGET_DIR/release/run_design_per_lambda" "$out/"
cp "$CARGO_TARGET_DIR/release/scope35_walk_front" "$out/"
cp "$CARGO_TARGET_DIR/release/examples/mfe_fold" "$out/"
git -C "$build" rev-parse HEAD > "$out/COMMIT"
(cd "$out" && sha256sum run_design_per_lambda scope35_walk_front mfe_fold > SHA256SUMS)
cat "$out/SHA256SUMS"
