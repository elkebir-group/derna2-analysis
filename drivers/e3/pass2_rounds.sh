#!/usr/bin/env bash
# Pass 2 as a sequence of 4 h rounds, until every cell reports `complete`.
#
# grid_resume.sh walks its job list ONCE and exits, so one invocation is one
# round. This drives round after round, which is what pass 2 needs: a cell is
# retired only when its front reports complete, and the open regions that make
# E3's depth pairs incomparable close a few boxes per round.
#
# It also picks up the cells a previous round could not start. Q9UBM8_F4 was
# skipped in round 1 because an orphaned cell from a killed driver shared its
# --out path and both writers had to be killed; the next round simply finds it
# unfinished and runs it.
#
# Between rounds it waits for the driver AND for every cell the driver forked:
# killing or exiting a driver does not reap its cells, and starting a new round
# while an old cell still writes puts two writers on one path. grid_resume.sh's
# own preflight refuses that, so a round that starts early would abort the whole
# round rather than corrupt a front -- this wait is what keeps it from happening.
set -u
P=/scratch/xinyugu4/paper-e3
OUT=$P/grid
BIN=$P/bin/derna_06605630
PINNAME=06605630
JOBS=${JOBS:-$OUT/jobs-pass2.txt}
MAXROUNDS=${MAXROUNDS:-12}
SLOTS=${SLOTS:-"1 3 5 7 9 11 13 15 17 21 23 25 27"}

cells_live() { ps -eo args | grep -c "[d]erna.* front .*--out $OUT/"; }

unfinished() {  # how many cells in the job list are not yet complete
  python3 -I -c "
import json, os
n=0
for line in open('$JOBS'):
    line=line.split()
    if len(line)<2: continue
    label='%s_F%s' % (line[0], line[1])
    best=None
    p='$OUT/'+label+'.jsonl'
    if os.path.exists(p) and os.path.getsize(p): best=p
    k=0
    while True:
        k+=1
        p='$OUT/%s.r%d.jsonl' % (label, k)
        if os.path.exists(p) and os.path.getsize(p): best=p
        else: break
    if best is None: n+=1; continue
    try:
        h=next(r for r in (json.loads(l) for l in open(best)) if 'points' in r)
        if h.get('complete') is not True: n+=1
    except Exception:
        n+=1
print(n)"
}

for r in $(seq 1 "$MAXROUNDS"); do
  while [ "$(cells_live)" -gt 0 ]; do sleep 60; done
  left=$(unfinished)
  echo "=== round $r: $left cell(s) unfinished, $(date -Is)"
  [ "$left" -eq 0 ] && { echo "=== pass 2 complete: every cell reports complete"; break; }
  BIN=$BIN PINNAME=$PINNAME THREADS=1 SLOTLIST="$SLOTS" JOBS=$JOBS \
    bash $P/grid_resume.sh >> $OUT/pass2-round$r.out 2>&1
  echo "=== round $r driver exited $(date -Is); $(unfinished) still unfinished"
done
echo "=== pass2_rounds done $(date -Is), $(unfinished) unfinished"
