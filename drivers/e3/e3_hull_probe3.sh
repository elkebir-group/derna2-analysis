#!/usr/bin/env bash
# E3 hull probe, round 3: the discriminator between length and AC-axis density.
# The same rungs as round 2, under the 8-pattern p8 set instead of the 25
# cumulative patterns. If the 200 -> 250 per-fold jump appears under p8 too, it
# is length; if it does not, it is the size of the blacklist.
#
# Waits for round 2 to finish first: two walks on the same cores would spoil
# both walls.
set -u
BASE=/home/melkebir/projects/derna2/e3-plan/experiments/paper-benchmark
until grep -q '^=== done' "$BASE/e3-probe2.run.log"; do sleep 60; done

BIN=/home/melkebir/projects/derna2/motif-work/bin/tip          # commit 5d12df5f
MOTIFS=/home/melkebir/projects/derna2/motif-work/inputs/p8.txt
ICOUNT=/home/melkebir/projects/derna2/motif-work/tools/icount.py
FA=/home/melkebir/projects/derna2/motif-work/inputs/fa
OUT=$BASE/e3-probe3
CORES=16-27
MEMMAX=24G

mkdir -p "$OUT"
printf 'label,n,motifs,exit,wall_s,peak_rss_kb,vertices,load_at_launch\n' > "$OUT/summary.csv"

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
  printf '%s,%s,p8,%s,%s,%s,%s,%s\n' "$label" "$n" "$code" "$wall" "$rss" "$verts" "$load" >> "$OUT/summary.csv"
  echo "=== $label exit $code wall ${wall}s rss ${rss}kB vertices $verts" >&2
}

run p8_200 200 "$FA/P31645_200.fasta"                 1800
run p8_250 250 "$BASE/e3-probe/fa/P31645_250.fasta"   3600
echo "=== done $(date -Is)" >&2
