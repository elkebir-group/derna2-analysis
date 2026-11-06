#!/usr/bin/env bash
# The short band's four-nest block: the six all-four binders at F1-F4, 24 cells.
#
# Ordered six-binders-first (decided 2026-10-20) because these are the only
# instances where set IV has a price at all, so they carry the complete nested
# F0-F4 table the paper leads with. The twelve three-nest members follow.
#
# The sets are NESTED and cumulative: F1 forbids set I, F2 sets I-II, F3 I-III,
# F4 all four. Each cell therefore passes one --motif per set up to its nest.
#
# Twelve cells at a time, one core each, from socket 1 — threads are a sublinear
# lever (2.05x at 8 on record, and that a free fold's) while concurrency is
# near-linear, so the throughput choice needs no thread curve. `--deterministic`
# makes the width result-neutral, so a one-thread cell here is comparable with
# any other width. Cores 25,27 are left for the flank arm.
#
# No UTRs and no --forbid: the headline grid is coding sequence alone.
set -u
BIN=${BIN:?set BIN to the pinned binary}
PINNAME=${PINNAME:?set PINNAME to the pinned commit}   # no apostrophe: inside
# ${VAR:?word} an unescaped ' opens a quote that swallows the rest of the file,
# which cost three failed launches reported as unbound variables far below.
CORPUS=run/paper-e3/screen/corpus
SETS=derna2/experiments/paper-benchmark/motifs
OUT=run/paper-e3/grid
BUDGET=28800       # the walk's own budget: at the limit it RETURNS the partial front
PER=$(( BUDGET + 600 ))   # the outer kill, above BUDGET so the graceful stop always fires first
# No --max-rss-gb. The flag is not a dormant ceiling: `budget.max_rss_gb.is_some()`
# alone puts the walk in residency mode, capping held passes at 2*threads+1 -- 3
# at one thread, against a 4-weight pass grid at every length here, so one pass is
# dropped and RE-FOLDED on each miss -- and `Grid::prefetch` returns false
# outright, so no worker ever takes a pass ahead of need. A cell peaks at 1.6 GB
# and twelve of them at about 20 GB against 419 GB available, so the ceiling could
# never bind and all of that was paid for nothing. It also sat on the constrained
# side of the constrained-to-free comparison and not on the free side, where the
# screen ran with max_rss_gb None, which made part of that ratio the cap's cost
# rather than the constraint's.

# Why the walk carries its own budget rather than only an outer `timeout`:
# a SIGTERM mid-walk writes NO output at all -- measured, exit 124 and no --out
# file -- so an outer kill alone makes an overrunning cell a total loss. The two
# saves are not redundant. `--max-secs` stops deliberately and returns every
# point plus the boxes it was mid-way through; `--checkpoint` rewrites --out
# after each settled hull pair, which is coarser (a kill at 200 s kept 243
# points against the budget's 278) but is the only thing that survives a kill
# --max-secs cannot intercept: OOM, a reboot, a stray pkill.
#
# A budgeted cell is NOT reproducible: --deterministic fixes the output at any
# thread width, but where a wall-clock budget falls depends on host speed and
# load. A row with complete=true is a publishable artifact; complete=false is a
# resumable intermediate and must be carried to completion with --resume before
# it enters a price table. The open_boxes column is what names the work left.
mkdir -p "$OUT"
SUM=$OUT/summary.csv
[ -f "$SUM" ] || printf 'protein,aa,nest,sets,build,threads,core,exit,wall_s,peak_rss_kb,points,hull,complete,open_boxes,load_at_launch\n' > "$SUM"

# Plain word lists rather than arrays. Not required -- the arrays were fine, and
# the failures blamed on them were the apostrophe above -- but a flat list needs
# no quoting care, which is the fault this file already paid for once.
BINDERS="P26436 Q9HC38 P47211 P15086 Q9NWZ3 Q9UBM8"   # ascending length, 265-478 aa
NESTS="1 2 3 4"
SETFILE_1=set1-manufacturing
SETFILE_2=set2-processing
SETFILE_3=set3-repeats
SETFILE_4=set4-initiation
SLOTLIST="1 3 5 7 9 11 13 15 17 19 21 23"   # concurrency IS this list's length

cell() {  # acc nest core
  local acc=$1 nest=$2 core=$3
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
    taskset -c "$core" timeout -k 30 "$PER" \
      "$BIN" front --protein "$fa" --threads 1 --deterministic \
        --max-secs "$BUDGET" --checkpoint \
        "${args[@]}" --out "$OUT/$label.jsonl" --force \
      > "$OUT/$label.log" 2>&1
  code=$?
  wall=$(( $(date +%s) - start ))
  rss=$(tail -1 "$OUT/$label.rss" 2>/dev/null)
  pts=; hull=; comp=; openb=
  # Read whatever is on disk whichever way the cell ended: a budgeted stop exits
  # 0 with complete=false, and a checkpoint file left by a kill is readable too.
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
  printf '%s,%s,F%s,%s,%s,1,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$acc" "$aa" "$nest" "$setlist" "$PINNAME" "$core" \
    "$code" "$wall" "${rss:-}" "$pts" "$hull" "$comp" "$openb" "$load" >> "$SUM"
  echo "=== $label core $core exit $code ${wall}s ${rss}kB $pts point(s) complete=$comp open=$openb" >&2
}

# One job per slot, each slot a core. The job list is a file so the set of cells
# is inspectable after the fact rather than only in the shell's memory.
#
# A cell is skipped only when BOTH its output and its summary row exist -- output
# alone would skip a cell no row was ever written for, which is how the screen
# lost A6NKF1.
JOBS=$OUT/jobs.txt
: > "$JOBS"
for acc in P26436 Q9HC38 P47211 P15086 Q9NWZ3 Q9UBM8; do
  for nest in 1 2 3 4; do
    printf '%s %s\n' "$acc" "$nest" >> "$JOBS"
  done
done

# A cell is skipped only when its front is COMPLETE and a row exists. A budgeted
# cell now exits 0 with a partial front, so the old "output exists" test would
# retire an unfinished cell as done -- the quiet way a price table ends up built
# on fronts that were still being walked. A partial is re-walked from scratch
# here; carrying it on with --resume is the cheaper follow-up and is deliberately
# not done by this driver, which has never exercised that path.
done_already() {  # label
  [ -s "$OUT/$1.jsonl" ] || return 1
  python3 -I -c "
import json,sys
recs=[json.loads(l) for l in open('$OUT/$1.jsonl')]
hdr=next((r for r in recs if 'points' in r), {})
sys.exit(0 if hdr.get('complete', '?') in (True, '?') else 1)
" 2>/dev/null
}

# One cell per listed core, claimed by a marker file. A modular index does NOT
# give this: with slot launching the freed slot is not the one `i % NSLOT` names,
# so cell 12 would take the first core whether or not cell 0 still held it --
# two cells on one core while another idles.
#
# A claim records its owner's pid, so a core is free when the claim is absent OR
# its owner is gone. That makes a killed run self-healing without deleting
# anything at startup: a wipe of core.* would be a glob delete run before the
# driver knows whether another copy of itself is live, and two drivers sharing
# one core list is the fault the claims exist to prevent.
CLAIM=$OUT/.claims
mkdir -p "$CLAIM"

core_free() {  # core -- absent claim, or one whose owner has exited
  local f=$CLAIM/core.$1 pid
  [ -e "$f" ] || return 0
  pid=$(cat "$f" 2>/dev/null)
  [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null && return 1
  return 0
}

while read -r acc nest; do
  [ -n "$acc" ] || continue
  label=${acc}_F${nest}
  if done_already "$label" && grep -q "^$acc,[0-9]*,F$nest," "$SUM"; then
    echo "=== $label already complete, skipping" >&2
    continue
  fi
  pick=
  while [ -z "$pick" ]; do
    for c in $SLOTLIST; do
      if core_free "$c"; then pick=$c; break; fi
    done
    [ -z "$pick" ] && sleep 10
  done
  # The claim is written BEFORE the fork and is never empty, in that order for
  # two reasons: if the subshell wrote its own claim, the parent could scan and
  # pick the same core again before it got there; and core_free reads an empty
  # file as free, so a zero-byte placeholder would not reserve anything. The
  # parent's pid holds the core for the instant between here and the handover.
  echo $$ > "$CLAIM/core.$pick"
  ( cell "$acc" "$nest" "$pick"; rm -f "$CLAIM/core.$pick" ) &
  echo $! > "$CLAIM/core.$pick"
done < "$JOBS"
wait
echo "=== done $(date -Is)" >&2
