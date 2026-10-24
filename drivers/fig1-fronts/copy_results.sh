#!/usr/bin/env bash
# Copy the small results from a run directory into results/ for commit.
# Usage: copy_results.sh RUN_ROOT   (e.g. /scratch/xinyugu4/glpe-hbb-fronts)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN="$1"
OUT="$HERE/results"
for arm in unguarded symmetric no_clamp unpaired; do
  mkdir -p "$OUT/$arm"
  cp "$RUN/fronts/$arm/front.scored.jsonl" "$RUN/fronts/$arm/vertices.scored.jsonl" "$OUT/$arm/"
  # The report without its per-task timeline.
  python3 -c 'import json,sys; r=json.load(open(sys.argv[1])); r.pop("timeline",None); json.dump(r,open(sys.argv[2],"w"),indent=1)' \
    "$RUN/fronts/$arm/report.json" "$OUT/$arm/report.json"
done
for arm in symmetric no_clamp unpaired; do
  cp "$RUN/fronts/$arm/front3.jsonl" "$OUT/$arm/"
  [[ -s "$RUN/fronts/$arm/front3_designs.scored.jsonl" ]] && cp "$RUN/fronts/$arm/front3_designs.scored.jsonl" "$OUT/$arm/"
done
mkdir -p "$OUT/mrnafold" "$OUT/plot"
cp "$RUN/mrnafold/compare.jsonl" "$OUT/mrnafold/"
cp "$RUN/plot/plot_data.json" "$RUN/plot/fronts.png" "$OUT/plot/"
du -sh "$OUT"
