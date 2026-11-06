#!/usr/bin/env bash
# Which exponent is peak RSS on? One short run at band length settles it.
#
# Two fits of the same quantity disagree: n^2.00 over nineteen complete,
# unconstrained, one-thread cells at 511-694 aa, and n^3.63 over six F1 cells at
# 265-478 aa. At 1082 aa they predict 20.5 GB and 81.7 GB at one thread, a 4x
# separation, and at 2828 aa they differ by 18x -- the difference between the
# band's top instance fitting a 503 GB host and needing 7 TB.
#
# F1 is deliberate, not F4: the n^3.63 fit is over F1 cells, so matching its
# depth keeps the comparison about length. Depth moves peak by up to 3.5x on one
# measured cell, which would otherwise confound the exponent.
#
# No --max-rss-gb: its presence reconfigures the walk (+72.1%) and does not bound
# the peak anyway, so it would corrupt the quantity being measured.
#
# Peak is reached during table construction rather than at the end of the walk --
# P15086_F2 held 7.63 GB from 114 points through 130 to its budget stop -- so a
# short budget measures it. If construction has not finished, the figure is a
# LOWER bound, which still refutes n^2.00 if it lands well above 20.5 GB.
set -u
P=run/paper-e3
C=$P/grid/.claims
MO=$P/motifs
S=$1
FA=$P/screen/corpus-xl/PRT-05.fasta
sleep 7200 & HOLD=$!
CORE=""
for try in $(seq 1 15); do
  for c in 21 23 17; do
    f=$C/core.$c
    if [ -e "$f" ]; then p=$(cat "$f" 2>/dev/null); [ -n "$p" ] && kill -0 "$p" 2>/dev/null && continue; fi
    echo "$HOLD" > "$f"
    [ "$(cat "$f" 2>/dev/null)" = "$HOLD" ] && { CORE=$c; break; }
  done
  [ -n "$CORE" ] && break
  sleep 12
done
[ -z "$CORE" ] && { echo "no core after 15 attempts"; kill $HOLD 2>/dev/null; exit 1; }
echo "=== PRT-05 F1, 1 thread, core $CORE, holder $HOLD, started $(date -Is)"
echo "=== predictions at 1082 aa, one thread: n^2.00 -> 20.5 GB   n^3.63 -> 81.7 GB"
/usr/bin/time -f '%M %e %U' -o $S/prt05.rss \
  taskset -c $CORE timeout -k 30 2700 \
  $P/bin/derna_06605630 front --protein $FA --threads 1 \
    --motif $MO/set1-manufacturing.txt \
    --out $S/prt05_f1.jsonl --force > $S/prt05.log 2>&1
echo "=== run exit $? at $(date -Is)"
read -r kb wall user < $S/prt05.rss
python3 -I -c "
kb=$kb
gb=kb/1048576
print('peak RSS %.2f GB  (wall %.0f s, user %.0f s)' % (gb, $wall, $user))
for name, pred in (('n^2.00', 20.5), ('n^3.63', 81.7)):
    print('   vs %s prediction %5.1f GB  ->  %.2fx' % (name, pred, gb/pred))
"
echo "=== did table construction finish? last lines of the walk log ==="
tail -4 $S/prt05.log
[ "$(cat $C/core.$CORE 2>/dev/null)" = "$HOLD" ] && command rm -f -- $C/core.$CORE && echo "released core $CORE"
kill $HOLD 2>/dev/null; true
