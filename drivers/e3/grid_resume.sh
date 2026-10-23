#!/usr/bin/env bash
# Phase B in 4 h rounds, resuming rather than restarting.
#
# The 18 h driver is retired. Every cell now runs for at most 4 h and is stopped
# at a gap boundary, and the next round continues from the open boxes the last one
# left. Resume is not an assumption here: measured over five arms today, a resumed
# walk reproduces an uninterrupted one bit-identically in the gap phase at one, two
# and four threads.
#
# Why rounds rather than one long run. A 4 h cell is a cell whose cost is known
# before it starts, where an 18 h cell is a forecast -- and the forecasts were
# wrong in the same direction every time, because per-gap cost CLIMBS as a walk
# proceeds. Q9UBM8_F1 was projected at 13.9 h from its first 24 gaps and was
# running at 26.7 h by its 35th. Rounds also mean the host can be reclaimed at any
# boundary without losing more than the gap in flight.
#
# The round convention. `--resume A --out B` reads A and writes B, so each round
# needs its own file: L.jsonl for round 0, then L.r1.jsonl, L.r2.jsonl, ... The
# current state of a cell is its HIGHEST round file, and that is what the next
# round resumes from. Rounds are never merged -- a resumed front already returns
# its own points together with the ones it inherited, so the newest file is the
# whole front and the earlier ones are history.
#
# A cell is retired only when its front reports `complete`, never when a file
# merely exists: a budget-stopped cell leaves a perfectly valid short file, and an
# existence test would retire it forever.
set -u
BIN=${BIN:?set BIN to the pinned binary}
PINNAME=${PINNAME:?set PINNAME to the pinned commit}
CORPUS=/scratch/xinyugu4/paper-e3/screen/corpus
SETS=/scratch/xinyugu4/paper-e3/motifs
OUT=/scratch/xinyugu4/paper-e3/grid
THREADS=${THREADS:-4}
# PASS1=1 runs only cells with NO front yet and never resumes one. Pass 1 buys
# breadth -- a 4 h front for every cell -- and pass 2 returns to the unfinished
# ones. Resuming inside pass 1 would spend the breadth budget on depth for
# whichever cell happened to start first.
PASS1=${PASS1:-0}
BUDGET=14400              # 4 h, the limit itself rather than a forecast
PER=$(( BUDGET + 600 ))   # the outer kill sits above it so the graceful stop wins
SUM=$OUT/summary.csv
JOBS=${JOBS:-$OUT/jobs-threaded-ordered.txt}   # nest-major list passed in by the caller
CLAIM=$OUT/.claims
# Node 0 is the EVEN cores and node 1 the odd ones; an arm per socket has
# decided an A/B here before, so a run names one node and stays in it.
SLOTLIST=${SLOTLIST:-"1 3 5 7 9 11 13 15 17 19 21 23"}
mkdir -p "$CLAIM"

SETFILE_1=set1-manufacturing
SETFILE_2=set2-processing
SETFILE_3=set3-repeats
SETFILE_4=set4-initiation

core_free() {  # free when the claim is absent OR its owner has exited
  local f=$CLAIM/core.$1 pid
  [ -e "$f" ] || return 0
  pid=$(cat "$f" 2>/dev/null)
  [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null && return 1
  return 0
}

busy_cores() {  # cores... -> "core(pct)" for each at or above 20%, one window
  local cs="$*" c out=""
  local -A b0 i0
  for c in $cs; do
    read -r x y < <(awk -v k="cpu$c" '$1==k{print $2+$3+$4+$6+$7+$8, $5}' /proc/stat)
    b0[$c]=$x; i0[$c]=$y
  done
  sleep 10
  for c in $cs; do
    read -r x y < <(awk -v k="cpu$c" '$1==k{print $2+$3+$4+$6+$7+$8, $5}' /proc/stat)
    local db=$(( x - b0[$c] )) di=$(( y - i0[$c] ))
    local pct=$(awk -v b="$db" -v i="$di" 'BEGIN{printf "%d", (b+i)?100*b/(b+i):100}')
    [ "$pct" -ge 20 ] && out="$out $c($pct%)"
  done
  echo "$out"
}

# Claim THREADS cores: claim first, then measure, then give back what is busy.
# Checking before claiming loses every race against another harvester, and an
# unclaimed-but-busy core is a real state that outlives a run by about a minute.
# Each write is read back, because a free-list scan followed by a write is a
# read-then-write and the last writer takes the file while both sides proceed as
# owner -- which put two four-thread cells on one set of four cores earlier today.
# The candidate order ROTATES by attempt. A fixed order starves: with one core
# per cell the first core is claimed, measured busy, handed back and tried again
# every round, so a persistently loaded first core means no later core is ever
# reached and the driver waits forever while the rest of the node is idle.
ROTATE=0
claim_set() {
  local want=$1 got="" c busy list
  list=$(python3 -I -c "
cs='$SLOTLIST'.split()
n=$ROTATE % len(cs)
print(' '.join(cs[n:]+cs[:n]))")
  ROTATE=$(( ROTATE + 1 ))
  for c in $list; do
    [ "$(echo $got | wc -w)" -ge "$want" ] && break
    core_free "$c" || continue
    echo $$ > "$CLAIM/core.$c"
    if [ "$(cat "$CLAIM/core.$c" 2>/dev/null)" = "$$" ]; then
      got="$got $c"
    else
      echo "=== lost the race for core $c to $(cat "$CLAIM/core.$c" 2>/dev/null)" >&2
    fi
  done
  if [ "$(echo $got | wc -w)" -lt "$want" ]; then
    for c in $got; do
      [ "$(cat "$CLAIM/core.$c" 2>/dev/null)" = "$$" ] && command rm -f -- "$CLAIM/core.$c"
    done
    echo ""; return 1
  fi
  busy=$(busy_cores $got)
  if [ -n "$busy" ]; then
    echo "=== claimed$got but BUSY:$busy -- giving them back" >&2
    for c in $got; do
      [ "$(cat "$CLAIM/core.$c" 2>/dev/null)" = "$$" ] && command rm -f -- "$CLAIM/core.$c"
    done
    echo ""; return 1
  fi
  echo $got
}

latest_round() {  # label -> the newest round file, or empty if none has content
  local lbl=$1 best="" n=0 f
  [ -s "$OUT/$lbl.jsonl" ] && best=$OUT/$lbl.jsonl
  while :; do
    n=$(( n + 1 ))
    f=$OUT/$lbl.r$n.jsonl
    [ -s "$f" ] || break
    best=$f
  done
  echo "$best"
}

next_round_path() {  # label -> the file this round should write
  local lbl=$1 n=0
  [ -s "$OUT/$lbl.jsonl" ] || { echo "$OUT/$lbl.jsonl"; return; }
  while :; do
    n=$(( n + 1 ))
    [ -s "$OUT/$lbl.r$n.jsonl" ] || { echo "$OUT/$lbl.r$n.jsonl"; return; }
  done
}

settled() {  # jsonl -> gaps, open, complete
  python3 -I -c "
import json
try:
    h=next(r for r in (json.loads(l) for l in open('$1')) if 'points' in r)
    print(h.get('gaps') or 0, len(h.get('open',[])), h.get('complete'))
except Exception:
    print(-1, -1, 'unreadable')" 2>/dev/null
}

is_complete() {  # label -> 0 if its newest round reports complete
  local cur; cur=$(latest_round "$1")
  [ -n "$cur" ] || return 1
  python3 -I -c "
import json,sys
h=next(r for r in (json.loads(l) for l in open('$cur')) if 'points' in r)
sys.exit(0 if h.get('complete') is True else 1)" 2>/dev/null
}

cell() {  # acc nest mask cores...
  local acc=$1 nest=$2 mask=$3; shift 3
  local cores="$*"
  local fa=$CORPUS/$acc.fasta label=${acc}_F${nest}
  local aa load start code wall rss pts hull comp openb k args=() names=()
  local from to
  for ((k=1;k<=nest;k++)); do
    eval "local f=\$SETFILE_$k"
    args+=(--motif "$SETS/$f.txt"); names+=("${f%%-*}")
  done
  local setlist; setlist=$(IFS=+; echo "${names[*]}")
  from=$(latest_round "$label")
  to=$(next_round_path "$label")
  if [ "$PASS1" = "1" ]; then
    from=""                      # pass 1 never resumes
  fi
  if [ -n "$from" ]; then
    args+=(--resume "$from")
    echo "=== $label round -> $(basename "$to"), resuming $(basename "$from") at $(settled "$from") gaps/open/complete" >&2
  fi
  # A front composed across rounds can be composed across BUILDS, which is sound
  # -- both builds compute the exact front, so inherited points and newly walked
  # boxes are each exact -- but the row's single build column would then name one
  # build for a front produced by two. So each round appends its own build to a
  # per-cell provenance file, and that file is the authority on what made a front.
  printf '%s %s %s %s\n' "$(date -Is)" "$(basename "$to")" "$PINNAME" \
    "resumed=$([ -n "$from" ] && basename "$from" || echo none)" >> "$OUT/$label.builds"
  aa=$(grep -v '^>' "$fa" | tr -d '\n\r ' | wc -c)
  load=$(cut -d' ' -f1 /proc/loadavg)
  start=$(date +%s)
  /usr/bin/time -f '%M' -o "$OUT/$label.rss" \
    taskset -c "$mask" timeout -k 30 "$PER" \
      "$BIN" front --protein "$fa" --threads "$THREADS" \
        --max-secs "$BUDGET" --checkpoint \
        "${args[@]}" --out "$to" --force \
      > "$OUT/$label.log" 2>&1
  code=$?
  wall=$(( $(date +%s) - start ))
  rss=$(tail -1 "$OUT/$label.rss" 2>/dev/null)
  pts=; hull=; comp=; openb=
  if [ -s "$to" ]; then
    read -r pts hull comp openb <<<"$(python3 -I -c "
import json
recs=[json.loads(l) for l in open('$to')]
pt=[r for r in recs if 'sequence' in r]
hdr=next((r for r in recs if 'points' in r), {})
print(hdr.get('points',len(pt)), sum(1 for r in pt if r.get('on_hull')),
      hdr.get('complete','?'), len(hdr.get('open',[])))
")"
  fi
  # the core mask is QUOTED: unquoted it occupies one field per core and shifts
  # every later column, which read as exit 11 and a 13 second wall on one row
  printf '%s,%s,F%s,%s,%s,%s,"%s",%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$acc" "$aa" "$nest" "$setlist" "$PINNAME" "$THREADS" "$mask" \
    "$code" "$wall" "${rss:-}" "$pts" "$hull" "$comp" "$openb" "$load" >> "$SUM"
  echo "=== $label $(basename "$to") exit $code ${wall}s $pts point(s) complete=$comp open=$openb" >&2
  for c in $cores; do
    [ "$(cat "$CLAIM/core.$c" 2>/dev/null)" = "$$" ] && command rm -f -- "$CLAIM/core.$c"
  done
}

# Refuse to start while any cell of this grid is already running. Killing a
# driver does NOT kill the cells it forked: they are reparented to init and
# keep writing, so a second driver re-derives the same job list and starts the
# same cell again, two processes truncating and writing one --out path. That
# happened once and cost nothing only because both were killed before either
# emitted; a merged short file is indistinguishable from a small front.
live=$(ps -eo args | grep "[d]erna.* front .*--out $OUT/" | wc -l)
if [ "$live" -gt 0 ]; then
  echo "=== REFUSING: $live cell(s) of this grid are already running:" >&2
  ps -eo pid,etime,args | grep "[d]erna.* front .*--out $OUT/" \
    | sed "s/.*--out //;s/ .*//" | sort | uniq -c >&2
  echo "=== kill them by PID (the cells, not just the driver) before relaunching" >&2
  exit 1
fi
echo "=== grid_resume: $BUDGET s rounds, $THREADS threads, pin $PINNAME ($(date -Is))"
while read -r acc nest; do
  [ -n "$acc" ] || continue
  label=${acc}_F${nest}
  if is_complete "$label"; then
    echo "=== $label complete, retired" >&2
    continue
  fi
  if [ "$PASS1" = "1" ] && [ -n "$(latest_round "$label")" ]; then
    echo "=== $label already has a front, deferred to pass 2" >&2
    continue
  fi
  got=""
  while [ -z "$got" ]; do
    got=$(claim_set "$THREADS") || got=""
    [ -z "$got" ] && sleep 20
  done
  mask=$(echo $got | tr ' ' ',')
  ( cell "$acc" "$nest" "$mask" $got ) &
done < "$JOBS"
wait
echo "=== grid_resume done $(date -Is)"
