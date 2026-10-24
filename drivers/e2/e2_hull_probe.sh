#!/usr/bin/env bash
# E2 hull probe: the exact motif hull at a length ladder, as a lower bound on
# the motif front's cost. Pre-registered in E2-PLAN.md § 3a before this ran.
#
# Primary readout is peak RSS. Walls are load-qualified lower bounds: the
# E-cores carry another session's builds.
set -u

BIN=/home/melkebir/projects/derna2/motif-work/bin/tip          # commit 5d12df5f
MOTIFS=/home/melkebir/projects/derna2/e3-plan/experiments/paper-benchmark/motifs/probe-all-sets.bare.txt
FA=/home/melkebir/projects/derna2/motif-work/inputs/fa
ICOUNT=/home/melkebir/projects/derna2/motif-work/tools/icount.py
OUT=/home/melkebir/projects/derna2/e3-plan/experiments/paper-benchmark/e3-probe
CORES=16-27          # E-cores; the P-cores carry E1 and the motif-speed timing
BUDGET=3600          # per rung, seconds
MEMMAX=24G

mkdir -p "$OUT"
printf 'rung,exit,wall_s,peak_rss_kb,vertices,load_at_launch\n' > "$OUT/summary.csv"

for n in 100 200 300 400 600; do
  fa="$FA/P31645_${n}.fasta"
  [ -f "$fa" ] || { echo "no rung $n" >&2; continue; }
  load=$(cut -d' ' -f1 /proc/loadavg)
  echo "=== rung n=$n at $(date -Is), load $load" >&2

  systemd-run --user --scope -q -p MemoryMax=$MEMMAX -p MemorySwapMax=0 -- \
    timeout -k 30 "$BUDGET" \
    taskset -c "$CORES" \
    python3 "$ICOUNT" "$OUT/n${n}.ic.json" -- \
      "$BIN" --fasta "$fa" --motifs "$MOTIFS" --arms walk --workers 1 --epsilon 0 \
      > "$OUT/n${n}.jsonl" 2> "$OUT/n${n}.err"
  code=$?

  wall=$(python3 -c "import json,sys;d=json.load(open('$OUT/n${n}.ic.json'));print(d.get('wall_s',''))" 2>/dev/null)
  rss=$(python3 -c "import json,sys;d=json.load(open('$OUT/n${n}.ic.json'));print(d.get('maxrss_kb',''))" 2>/dev/null)
  verts=$(grep -c '"walk-vertex"' "$OUT/n${n}.jsonl" 2>/dev/null || echo 0)
  printf '%s,%s,%s,%s,%s,%s\n' "$n" "$code" "$wall" "$rss" "$verts" "$load" >> "$OUT/summary.csv"
  echo "=== rung n=$n exit $code wall ${wall}s rss ${rss}kB" >&2

  # A rung that blew the budget or the cap settles the ladder: longer rungs cost more.
  [ "$code" -eq 0 ] || { echo "stopping: rung $n hit a limit (exit $code)" >&2; break; }
done
echo "=== done $(date -Is)" >&2
