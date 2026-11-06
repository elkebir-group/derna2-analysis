#!/usr/bin/env bash
# The peak-RSS exponent, measured on a ladder rather than argued between two fits.
#
# Two fits of peak RSS disagree -- n^2.00 over nineteen unconstrained one-thread
# cells at 511-694 aa, n^3.63 over six F1 cells at 265-478 aa -- and at 2828 aa
# they differ by 18x, which decides whether the long band's top instance runs at
# all. Neither can be extended cheaply, because both were taken from whole front
# walks and a walk at 1082 aa did not finish constructing its tables in 2700 s.
#
# `design` is the cheaper instrument. A single solve at one weight builds the same
# constrained tables, so its peak is the CONSTRUCTION peak, and construction is
# where a walk's peak is set -- P15086_F2 held 7.63 GB from 114 points through 130
# to its budget stop. One solve replaces a walk.
#
# What design does NOT measure is the walk's own storage on top of construction,
# so a design peak is a lower bound on the matching front peak. That is why the
# ladder's first rung is 478 aa, where the front peak is already known at one
# thread (4.21 GB, from 11.28 GB at four threads over the measured 2.68x): the
# rung calibrates the two instruments against each other before any rung is used
# to extrapolate. If design at 478 aa comes back far under 4.21 GB, the instrument
# measures a different quantity and the ladder is reported as the fold's curve
# rather than the walk's.
#
# Depth is held at set I (F1) across every rung, matching the n^3.63 fit, because
# depth moved peak by 3.5x on one measured cell and would otherwise confound a
# trend in length. No --max-rss-gb: its presence reconfigures the walk and does
# not bound the peak.
set -u
P=run/paper-e3
C=$P/grid/.claims
MO=$P/motifs
S=$1
BIN=$P/bin/derna_06605630
PER=1800
LAMBDA=0.5
sleep 21600 & HOLD=$!
release() {
  [ -n "${CORE:-}" ] && [ "$(cat $C/core.$CORE 2>/dev/null)" = "$HOLD" ] && \
    command rm -f -- $C/core.$CORE && echo "released core $CORE"
  kill $HOLD 2>/dev/null
}
trap release EXIT
CORE=""
for try in $(seq 1 40); do
  for c in 21 23 17; do
    f=$C/core.$c
    if [ -e "$f" ]; then p=$(cat "$f" 2>/dev/null); [ -n "$p" ] && kill -0 "$p" 2>/dev/null && continue; fi
    echo "$HOLD" > "$f"
    [ "$(cat "$f" 2>/dev/null)" = "$HOLD" ] && { CORE=$c; break; }
  done
  [ -n "$CORE" ] && break
  sleep 15
done
[ -z "$CORE" ] && { echo "no core after 40 attempts"; exit 1; }
echo "=== design ladder on core $CORE, holder $HOLD, $(date -Is)"
printf "%-10s %5s %10s %9s %8s %s\n" instance aa "peak GB" "wall s" exit note
for inst in Q9UBM8:$P/screen/corpus/Q9UBM8.fasta \
            Q6UXY8:$P/screen/corpus-xl/Q6UXY8.fasta \
            PRT-05:$P/screen/corpus-xl/PRT-05.fasta; do
  name=${inst%%:*}; fa=${inst#*:}
  [ -s "$fa" ] || { printf "%-10s %5s %10s %9s %8s %s\n" "$name" "-" "-" "-" "-" "fasta missing"; continue; }
  aa=$(grep -v '^>' "$fa" | tr -d '[:space:]' | wc -c)
  /usr/bin/time -f '%M %e' -o $S/ladder.$name.t \
    taskset -c $CORE timeout -k 30 $PER \
    $BIN design --protein "$fa" --lambda $LAMBDA --threads 1 \
      --motif $MO/set1-manufacturing.txt > $S/ladder.$name.log 2>&1
  code=$?
  read -r kb wall < $S/ladder.$name.t
  note=$([ $code -eq 0 ] && echo "solved" || echo "KILLED at ${PER}s, peak is a LOWER bound")
  printf "%-10s %5s %10.2f %9s %8s %s\n" "$name" "$aa" "$(echo "$kb" | awk '{printf "%.2f",$1/1048576}')" "$wall" "$code" "$note"
done
echo "=== done $(date -Is)"
