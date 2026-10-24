#!/usr/bin/env bash
# Phase B's remaining 20 cells, at --threads 4.
#
# Why width 4 and no --fold-threads, measured on P54852 set I, complete runs,
# fronts identical at every arm:
#
#   cores   pure --threads        mixed
#     2     1.76x  (88%)
#     4     3.07x  (77%)          t2 x ft2  2.35x (59%)   ft4  1.56x (39%)
#     8     6.00x  (75%)          t4 x ft2  4.29x (54%)   t2 x ft4  2.84x (36%)
#
# Across-gap parallelism wins at every budget because these cells have 44-88
# independent gaps; moving threads into folds costs 28% at 8 cores and 53% at
# 4. Throughput is flat between widths 4 and 8 (3 cells x 3.07x = 9.2 against
# 1 x 6.00x + 1 x 3.07x = 9.07 on twelve cores), so the width is chosen for
# latency: at 4 the dearest cell lands near 13.7 h instead of 42.
#
# --deterministic is NOT passed: it clamps --threads to 1. The user approved
# dropping byte-identical output, and in every arm measured the front came out
# identical anyway (sequence, energy and cai per point, in emission order) --
# but the guarantee is gone, so a cell's row records its width.
#
# Longest-first, so the dearest cells start earliest and the tail is shortest.
set -u
BIN=${BIN:?set BIN to the pinned binary}
PINNAME=${PINNAME:?set PINNAME to the pinned commit}
CORPUS=/scratch/xinyugu4/paper-e3/screen/corpus
SETS=/scratch/xinyugu4/paper-e3/motifs   # durable copy, NOT the branch-switched repo path
OUT=/scratch/xinyugu4/paper-e3/grid
THREADS=4
BUDGET=64800              # 18 h: the dearest projection is 13.7 h, which is a
PER=$(( BUDGET + 600 ))   # forecast rather than a runaway, so the budget covers it
SUM=$OUT/summary.csv
JOBS=$OUT/jobs-threaded-ordered.txt
CLAIM=$OUT/.claims
SLOTLIST="1 3 5 7 9 11 13 15 17 19 21 23"
mkdir -p "$CLAIM"

SETFILE_1=set1-manufacturing
SETFILE_2=set2-processing
SETFILE_3=set3-repeats
SETFILE_4=set4-initiation

# A core is free when its claim is absent OR its owner has exited. Checking only
# existence leaves a dead holder's claim blocking a core for good, which cost a
# harvester several minutes earlier tonight.
core_free() {
  local f=$CLAIM/core.$1 pid
  [ -e "$f" ] || return 0
  pid=$(cat "$f" 2>/dev/null)
  [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null && return 1
  return 0
}

# Claim THREADS cores at once, or none: a cell pinned to fewer cores than it has
# threads oversubscribes and its wall stops meaning anything.
claim_set() {
  local want=$1 got="" c
  for c in $SLOTLIST; do
    [ "$(echo $got | wc -w)" -ge "$want" ] && break
    core_free "$c" && got="$got $c"
  done
  [ "$(echo $got | wc -w)" -lt "$want" ] && { echo ""; return 1; }
  for c in $got; do echo $$ > "$CLAIM/core.$c"; done
  echo $got
}

cell() {  # acc nest mask cores...
  local acc=$1 nest=$2 mask=$3; shift 3
  local cores="$*"
  local fa=$CORPUS/$acc.fasta label=${acc}_F${nest}
  local aa load start code wall rss pts hull comp openb k args=() names=()
  for ((k=1;k<=nest;k++)); do
    eval "local f=\$SETFILE_$k"
    args+=(--motif "$SETS/$f.txt"); names+=("${f%%-*}")
  done
  local setlist; setlist=$(IFS=+; echo "${names[*]}")
  aa=$(grep -v '^>' "$fa" | tr -d '\n\r ' | wc -c)
  load=$(cut -d' ' -f1 /proc/loadavg)
  start=$(date +%s)
  /usr/bin/time -f '%M' -o "$OUT/$label.rss" \
    taskset -c "$mask" timeout -k 30 "$PER" \
      "$BIN" front --protein "$fa" --threads "$THREADS" \
        --max-secs "$BUDGET" --checkpoint \
        "${args[@]}" --out "$OUT/$label.jsonl" --force \
      > "$OUT/$label.log" 2>&1
  code=$?
  wall=$(( $(date +%s) - start ))
  rss=$(tail -1 "$OUT/$label.rss" 2>/dev/null)
  pts=; hull=; comp=; openb=
  if [ -s "$OUT/$label.jsonl" ]; then
    read -r pts hull comp openb <<<"$(python3 -I -c "
import json
recs=[json.loads(l) for l in open('$OUT/$label.jsonl')]
pt=[r for r in recs if 'sequence' in r]
hdr=next((r for r in recs if 'points' in r), {})
print(hdr.get('points',len(pt)), sum(1 for r in pt if r.get('on_hull')),
      hdr.get('complete','?'), len(hdr.get('open',[])))
")"
  fi
  printf '%s,%s,F%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$acc" "$aa" "$nest" "$setlist" "$PINNAME" "$THREADS" "$mask" \
    "$code" "$wall" "${rss:-}" "$pts" "$hull" "$comp" "$openb" "$load" >> "$SUM"
  echo "=== $label cores $mask exit $code ${wall}s ${rss}kB $pts point(s) complete=$comp open=$openb" >&2
  for c in $cores; do rm -f "$CLAIM/core.$c"; done
}

done_already() {  # label -- only a COMPLETE front retires a cell
  [ -s "$OUT/$1.jsonl" ] || return 1
  python3 -I -c "
import json,sys
recs=[json.loads(l) for l in open('$OUT/$1.jsonl')]
hdr=next((r for r in recs if 'points' in r), {})
sys.exit(0 if hdr.get('complete', '?') in (True, '?') else 1)
" 2>/dev/null
}

while read -r acc nest; do
  [ -n "$acc" ] || continue
  label=${acc}_F${nest}
  if done_already "$label"; then
    echo "=== $label already complete, skipping" >&2
    continue
  fi
  got=""
  while [ -z "$got" ]; do
    got=$(claim_set "$THREADS") || got=""
    [ -z "$got" ] && sleep 20
  done
  mask=$(echo $got | tr ' ' ',')
  echo "=== $label -> cores $mask at $(date +%H:%M:%S)" >&2
  ( cell "$acc" "$nest" "$mask" $got ) &
done < "$JOBS"
wait
echo "=== done $(date -Is)" >&2
