#!/usr/bin/env bash
# Front 3: for each on-hull vertex of a guarded front, a lazy k-best walk at
# the midpoint of the vertex's weight interval to the first structurally
# optimal design under that arm's guard. One JSON line per vertex.
#
# Usage: [ONLY="k ..."] run_front3.sh OUT_ROOT ARM [JOBS]   (arms: symmetric no_clamp unpaired)
# Reads OUT_ROOT/ARM/vertices.jsonl (from analyze_fronts.py); writes
# OUT_ROOT/ARM/front3/v<k>.json and v<k>.log. A vertex with a finished
# v<k>.json is skipped, so a stopped run resumes.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
ROOT="$1"
ARM="$2"
JOBS="${3:-1}"
MAX_SECS="${MAX_SECS:-7200}"
BIN="${BIN:-/scratch/xinyugu4/glpe-hbb-fronts/bin}"
FASTA="$REPO/experiments/scope-25-beam-runtime-quality/fixtures/GLPE_HUMAN_78.fasta"
UTR5="$(cat "$REPO/experiments/scope-19-utr-support/m5/hbb5.txt")"
UTR3="$(cat "$REPO/experiments/scope-19-utr-support/m5/hbb3.txt")"

case "$ARM" in
  symmetric) TABLE=symmetric ;;
  no_clamp) TABLE=five-prime+no-clamp ;;
  unpaired) TABLE=five-prime+unpaired ;;
  *) echo "unknown arm $ARM" >&2; exit 2 ;;
esac

OUT="$ROOT/$ARM/front3"
mkdir -p "$OUT"

walk() {
  local k="$1" lam="$2"
  [[ -s "$OUT/v$k.json" ]] && return 0
  rm -f "$OUT/v$k.json.part"
  echo "[$(date +%T)] $ARM v$k λ=$lam start load=$(cut -d' ' -f1-3 /proc/loadavg)" >&2
  cut -d' ' -f1-3 /proc/loadavg > "$OUT/v$k.load"
  /usr/bin/time -v "$BIN/lazy_kbest_probe" --fasta "$FASTA" --utr5 "$UTR5" --utr3 "$UTR3" \
    --production --table "$TABLE" --lambda "$lam" --max-secs "$MAX_SECS" \
    --out "$OUT/v$k.json.part" 2> "$OUT/v$k.log"
  mv "$OUT/v$k.json.part" "$OUT/v$k.json"
  echo "[$(date +%T)] $ARM v$k done" >&2
}
export -f walk
export OUT ARM TABLE BIN FASTA UTR5 UTR3 MAX_SECS

# ONLY="0 3" restricts the run to those vertices.
python3 -c '
import json, sys
only = set(sys.argv[2].split())
for line in open(sys.argv[1]):
    v = json.loads(line)
    if not only or str(v["vertex"]) in only:
        print(v["vertex"], repr(v["lambda_mid"]))
' "$ROOT/$ARM/vertices.jsonl" "${ONLY:-}" | xargs -P "$JOBS" -n 2 bash -c 'walk "$0" "$1"'
