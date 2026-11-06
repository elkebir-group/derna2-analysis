#!/usr/bin/env bash
# Stage 1 of the flank cross-check: the hbb F0 baseline, which is the price of
# the arm and the gate on the remaining eight cells.
#
# The arm asks whether E2's headline prices, all measured with no flanks, survive
# a realistic 5' context. P47211 is the instance, fixed by a median rule before
# any price was visible. Two flank pairs, differing only in the 5' context:
#
#   hbb   F1-F3 flanked vs unflanked. F4 needs no run: it is already REFUSED,
#         and that refusal is the arm's strongest row rather than a gap in it.
#         RCCATGG is the Kozak consensus, so it straddles the 5' junction
#         whenever the UTR ends in RCC -- which is what a strong initiation
#         context is -- and the coding side closes it, residue 1 being M with a
#         unique codon and P47211's residue 2 being E, every codon of which
#         begins with G. Set IV and a real initiation context are not expensive
#         together, they are empty.
#   alb   all four nests, so the row licensing unqualified headline prices still
#         exists. alb carries no Kozak context, which is exactly the price of
#         having one.
#
# alb is paired with hbb's 3' UTR rather than a shorter one, so the two arms
# differ only in the variable under test and their folds stay comparable at 1229
# against 1220 nt rather than 1229 against 1103.
#
# Why only F0 here. The arm is unpriced and a flanked fold is 1229 nt against
# 1047 for coding sequence alone, so its cost is not known from any cell E2 has
# run. The baseline is needed whatever the rest of the arm does, so it is run
# first and alone, and the eight remaining cells are committed only once it
# reports a wall.
#
# --threads 2 rather than 4 because 25,27 are the arm's reserved pair. Width does
# not affect the front -- measured identical at 1, 2, 4 and 8 threads on complete
# runs -- and this arm publishes a front comparison, not a wall, so the narrower
# width costs it nothing it claims.
set -u
P=run/paper-e3
CLAIM=$P/grid/.claims
OUT=$P/flank/out
BIN=$P/bin/derna_549074e9
FA=$P/screen/corpus/P47211.fasta
U5=$P/flank/hbb5.txt
U3=$P/flank/hbb3.txt
CORES="25 27"
THREADS=2
BUDGET=${BUDGET:-21600}          # 6 h: generous, since the flanked cost is unmeasured
PER=$(( BUDGET + 600 ))
mkdir -p "$OUT"

for f in "$FA" "$U5" "$U3" "$BIN"; do
  [ -s "$f" ] || { echo "missing $f"; exit 1; }
done

# A core is taken only when it is BOTH unclaimed and measurably idle.
#
# Claim absence alone is not enough and this arm is where that bites: 25 and 27
# are lent to another experiment that is using them without writing a claim, so
# the claim file says free while core 25 runs at 100%. Launching on the claim
# alone would put this arm on top of theirs and corrupt both -- the same
# collision that cost a peer half a core per cell earlier, from the same mistake
# of trusting a bookkeeping record over a measurement.
#
# Idleness is read from /proc/stat tick deltas, not from ps: %CPU there is a
# lifetime average and psr is the last CPU a thread ran on, neither of which says
# what a core is doing now.
idle_pct() {  # core -- busy percentage over a 10 s window
  local c=$1 b0 i0 b1 i1
  read -r b0 i0 < <(awk -v k="cpu$c" '$1==k{print $2+$3+$4+$6+$7+$8, $5}' /proc/stat)
  sleep 10
  read -r b1 i1 < <(awk -v k="cpu$c" '$1==k{print $2+$3+$4+$6+$7+$8, $5}' /proc/stat)
  awk -v b="$((b1-b0))" -v i="$((i1-i0))" 'BEGIN{printf "%d", (b+i)?100*b/(b+i):100}'
}

# The arm's own pair, not harvested from the grid's list. Still claimed through
# the same protocol so a third party sees them held, and still checked for a live
# owner rather than for mere existence.
sleep 36000 & HOLD=$!
MINE=""
# Release ONLY claims this script still owns -- never by name. Deleting by name
# released cores the grid driver had re-claimed and put two four-thread cells on
# one set of four, which is the defect this arm must not repeat.
release() {
  local c f owner
  for c in $MINE; do
    f=$CLAIM/core.$c
    owner=$(cat "$f" 2>/dev/null)
    if [ "$owner" = "$HOLD" ]; then
      command rm -f -- "$f"
    else
      echo "=== core $c now owned by ${owner:-nobody}, left alone" >&2
    fi
  done
  kill $HOLD 2>/dev/null
}
trap release EXIT

echo "=== waiting for cores $CORES ($(date +%H:%M:%S))"
while :; do
  free=""
  for c in $CORES; do
    f=$CLAIM/core.$c
    if [ -e "$f" ]; then
      pid=$(cat "$f" 2>/dev/null)
      [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null && continue
    fi
    free="$free $c"
  done
  if [ "$(echo $free | wc -w)" -eq "$(echo $CORES | wc -w)" ]; then
    busy=""
    for c in $free; do
      pct=$(idle_pct "$c")
      [ "$pct" -ge 20 ] && busy="$busy $c($pct%)"
    done
    if [ -z "$busy" ]; then
      MINE=$free
      for c in $MINE; do echo "$HOLD" > "$CLAIM/core.$c"; done
      echo "=== claimed$MINE at $(date +%H:%M:%S), all measurably idle"
      break
    fi
    echo "=== unclaimed but BUSY, so not taken:$busy ($(date +%H:%M:%S))"
  fi
  sleep 30
done
MASK=$(echo $MINE | tr ' ' ',')

u5=$(tr -d '\n\r ' < "$U5"); u3=$(tr -d '\n\r ' < "$U3")
cds=$(( $(grep -v '^>' "$FA" | tr -d '\n\r ' | wc -c) * 3 ))
echo "=== hbb F0: 5' ${#u5} nt + CDS $cds nt + 3' ${#u3} nt = $(( ${#u5} + cds + ${#u3} )) nt"

start=$(date +%s)
/usr/bin/time -f '%M' -o "$OUT/hbb_F0.rss" \
  taskset -c "$MASK" timeout -k 30 "$PER" \
    "$BIN" front --protein "$FA" --utr5 "$u5" --utr3 "$u3" \
      --threads "$THREADS" --max-secs "$BUDGET" --checkpoint \
      --out "$OUT/hbb_F0.jsonl" --force > "$OUT/hbb_F0.log" 2>&1
code=$?
wall=$(( $(date +%s) - start ))
echo "=== hbb F0 exit $code, ${wall}s, $(tail -1 "$OUT/hbb_F0.rss" 2>/dev/null)kB peak"

python3 -I -c "
import json, os
def load(p):
    rs=[json.loads(l) for l in open(p)]
    h=next(r for r in rs if 'points' in r)
    return h,[r for r in rs if 'sequence' in r]
h,pt = load('$OUT/hbb_F0.jsonl')
print()
print('hbb F0 : %d points, %d hull, %s gaps, complete=%s, %d open, wall %.0f s'
      % (len(pt), sum(1 for r in pt if r.get('on_hull')), h.get('gaps'),
         h.get('complete'), len(h.get('open',[])), h.get('wall_s') or 0))
cds='$P/screen/f0/P47211.jsonl'
if os.path.exists(cds):
    hc,pc = load(cds)
    print('CDS F0 : %d points, %d hull, %s gaps, complete=%s, wall %.0f s'
          % (len(pc), sum(1 for r in pc if r.get('on_hull')), hc.get('gaps'),
             hc.get('complete'), hc.get('wall_s') or 0))
    print()
    print('The two are NOT comparable point for point -- a flanked fold scores a')
    print('different molecule, so the energies live on different scales. What the')
    print('arm compares is the SHAPE of the price as the nests deepen, which is')
    print('why F1-F3 are needed and why this cell only reports a cost.')
    if h.get('complete') and hc.get('complete'):
        print()
        print('cost of the flank on the free front: %.2fx wall, %+d points'
              % ((h.get('wall_s') or 0)/(hc.get('wall_s') or 1), len(pt)-len(pc)))
print()
if h.get('complete'):
    print('GATE: PASSES -- commit the remaining eight cells (hbb F1-F3, alb F1-F4,')
    print('      plus the recorded hbb F4 refusal row, which needs no run).')
    print('      Budget each at this wall times a margin, not at this wall.')
else:
    print('GATE: the baseline did not complete in %d s (%d open box(es)).' % ($BUDGET, len(h.get('open',[]))))
    print('      Eight more flanked cells at this cost is the decision to take')
    print('      before committing any of them.')
" 2>&1 | tee "$OUT/hbb_F0.verdict.txt"
echo "=== releasing $MASK"
