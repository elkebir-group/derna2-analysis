#!/usr/bin/env bash
# E2 hull probe, round 2. Three fixes over round 1:
#   - icount wraps timeout, not the reverse, so a censored rung still records
#     its wall and peak RSS;
#   - n=200 is re-run first as a load control, since round 1's n=300 launched
#     at load 11.9 against 3.9 for n=200 and a bound cannot be read without it;
#   - intermediate rungs 150 and 250, so the fit gets its third and fourth
#     points without another hour-plus run.
set -u
BIN=/home/melkebir/projects/derna2/motif-work/bin/tip          # commit 5d12df5f
BASE=/home/melkebir/projects/derna2/e3-plan/experiments/paper-benchmark
MOTIFS=$BASE/motifs/probe-all-sets.bare.txt
ICOUNT=/home/melkebir/projects/derna2/motif-work/tools/icount.py
OUT=$BASE/e3-probe2
CORES=16-27
MEMMAX=24G

mkdir -p "$OUT"
printf 'label,n,exit,wall_s,peak_rss_kb,vertices,load_at_launch\n' > "$OUT/summary.csv"

run() {  # label n fasta budget
  local label=$1 n=$2 fa=$3 budget=$4
  local load; load=$(cut -d' ' -f1 /proc/loadavg)
  echo "=== $label n=$n at $(date -Is), load $load, budget ${budget}s" >&2
  systemd-run --user --scope -q -p MemoryMax=$MEMMAX -p MemorySwapMax=0 -- \
    taskset -c "$CORES" \
    python3 "$ICOUNT" "$OUT/$label.ic.json" -- \
      timeout -k 30 "$budget" \
      "$BIN" --fasta "$fa" --motifs "$MOTIFS" --arms walk --workers 1 --epsilon 0 \
      > "$OUT/$label.jsonl" 2> "$OUT/$label.err"
  local code=$?
  local wall rss verts
  wall=$(python3 -c "import json;print(json.load(open('$OUT/$label.ic.json')).get('wall_s',''))" 2>/dev/null)
  rss=$(python3 -c "import json;print(json.load(open('$OUT/$label.ic.json')).get('maxrss_kb',''))" 2>/dev/null)
  verts=$(grep -c '"walk-vertex"' "$OUT/$label.jsonl" 2>/dev/null || echo 0)
  printf '%s,%s,%s,%s,%s,%s,%s\n' "$label" "$n" "$code" "$wall" "$rss" "$verts" "$load" >> "$OUT/summary.csv"
  echo "=== $label exit $code wall ${wall}s rss ${rss}kB vertices $verts" >&2
}

FA=/home/melkebir/projects/derna2/motif-work/inputs/fa
run control200 200 "$FA/P31645_200.fasta"        1800
run n150       150 "$OUT/../e3-probe/fa/P31645_150.fasta" 1800
run n250       250 "$OUT/../e3-probe/fa/P31645_250.fasta" 3600
run n300       300 "$FA/P31645_300.fasta"       10800
echo "=== done $(date -Is)" >&2
