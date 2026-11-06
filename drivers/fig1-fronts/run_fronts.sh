#!/usr/bin/env bash
# Fronts 1 and 2: the exact Pareto front of GLPE_HUMAN n = 78 with the hbb
# UTRs, unguarded and under each guard, one walk per arm.
#
# Usage: run_fronts.sh OUT_ROOT [ARM...]   (arms: unguarded symmetric no_clamp unpaired)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
ROOT="$1"
shift
ARMS=("${@:-unguarded symmetric no_clamp unpaired}")
BIN="${BIN:-run/glpe-hbb-fronts/bin}"
FASTA="$REPO/experiments/scope-25-beam-runtime-quality/fixtures/GLPE_HUMAN_78.fasta"
UTR5="$(cat "$REPO/experiments/scope-19-utr-support/m5/hbb5.txt")"
UTR3="$(cat "$REPO/experiments/scope-19-utr-support/m5/hbb3.txt")"

region() {
  case "$1" in
    unguarded) echo "" ;;
    symmetric) echo "5-cds,cds-3" ;;
    no_clamp) echo "5-cds,5-3" ;;
    unpaired) echo "5-5,5-cds,5-3" ;;
    *) echo "unknown arm $1" >&2; exit 2 ;;
  esac
}

for arm in ${ARMS[*]}; do
  out="$ROOT/$arm"
  mkdir -p "$out"
  classes="$(region "$arm")"
  args=(--fasta "$FASTA" --n 78 --utr5 "$UTR5" --utr3 "$UTR3"
        --out "$out/report.json" --front "$out/front.jsonl")
  [[ -n "$classes" ]] && args+=(--region "$classes")
  echo "[$(date +%T)] $arm start" >&2
  /usr/bin/time -v "$BIN/scope35_walk_front" "${args[@]}" 2> "$out/stderr.log"
  echo "[$(date +%T)] $arm done: $(grep '^n=' "$out/stderr.log")" >&2
done
