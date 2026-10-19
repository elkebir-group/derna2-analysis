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

# The algorithm crate is packages/derna-beam up to derna2's 27ffb226 and packages/derna-lib
# from it on; build whichever the pinned commit has.
crate=$build/packages/derna-lib
[ -d "$crate" ] || crate=$build/packages/derna-beam

export CARGO_TARGET_DIR=$build/target
export PATH=$HOME/.cargo/bin:$PATH
cargo build --release --quiet --manifest-path "$crate/Cargo.toml" \
    --bin run_design_per_lambda --bin scope35_walk_front
cargo build --release --quiet --manifest-path "$build/packages/rna-mfe/Cargo.toml" \
    --example mfe_fold
# E2's DERNA v1 arms run the Rust port.
cargo build --release --quiet --manifest-path "$build/packages/derna-codon-dp/Cargo.toml" \
    --bin run_brd_per_lambda --bin run_glpe78

bins=(run_design_per_lambda scope35_walk_front run_brd_per_lambda run_glpe78)
mkdir -p "$out"
for b in "${bins[@]}"; do
    cp "$CARGO_TARGET_DIR/release/$b" "$out/"
done
cp "$CARGO_TARGET_DIR/release/examples/mfe_fold" "$out/"
git -C "$build" rev-parse HEAD > "$out/COMMIT"
(cd "$out" && sha256sum "${bins[@]}" mfe_fold > SHA256SUMS)
cat "$out/SHA256SUMS"
