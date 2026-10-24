#!/usr/bin/env bash
# Archive a run's raw outputs into results/raw/: every front-3 walk's full
# output line, start load and time log; the full walk reports and logs; the
# raw mRNAfold outputs and precheck, and its configs as one compressed tar; the
# gate and launcher logs; the interactive page. Frozen binaries are left out:
# they rebuild from the commit recorded in bin/BUILT_FROM.txt.
# Usage: copy_raw.sh RUN_ROOT   (e.g. /scratch/xinyugu4/glpe-hbb-fronts)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN="$1"
OUT="$HERE/results/raw"
mkdir -p "$OUT"
for arm in unguarded symmetric no_clamp unpaired; do
  mkdir -p "$OUT/fronts/$arm"
  cp "$RUN/fronts/$arm/report.json" "$RUN/fronts/$arm/front.jsonl" "$RUN/fronts/$arm/stderr.log" \
     "$RUN/fronts/$arm/vertices.jsonl" "$OUT/fronts/$arm/"
  cp "$RUN/fronts/$arm.driver.log" "$OUT/fronts/"
done
for arm in symmetric no_clamp unpaired; do
  mkdir -p "$OUT/fronts/$arm/front3"
  cp "$RUN/fronts/$arm/front3/"v*.json "$RUN/fronts/$arm/front3/"v*.load "$RUN/fronts/$arm/front3/"v*.log \
     "$OUT/fronts/$arm/front3/"
  cp "$RUN/fronts/$arm/front3_designs.jsonl" "$OUT/fronts/$arm/" 2>/dev/null || true
done
mkdir -p "$OUT/mrnafold"
cp "$RUN/mrnafold/runs.jsonl" "$RUN/mrnafold/runs.scored.jsonl" "$RUN/mrnafold/driver.log" "$OUT/mrnafold/"
tar -C "$RUN/mrnafold" -czf "$OUT/mrnafold/configs.tar.gz" work precheck
mkdir -p "$OUT/sanity" "$OUT/logs" "$OUT/plot" "$OUT/bin"
cp "$RUN/sanity/"*.json "$OUT/sanity/"
cp "$RUN/front3.log" "$RUN/front3.log.gate" "$RUN/rust-test.log" "$RUN/build.log" "$OUT/logs/"
cp "$RUN/plot/glpe-hbb-fronts.html" "$OUT/plot/"
cp "$RUN/bin/BUILT_FROM.txt" "$OUT/bin/"
sha256sum "$RUN/bin/"glpe_hbb_rescore "$RUN/bin/"lazy_kbest_probe "$RUN/bin/"scope35_walk_front \
  | sed "s#$RUN/##" > "$OUT/bin/SHA256SUMS"
du -sh "$OUT"
