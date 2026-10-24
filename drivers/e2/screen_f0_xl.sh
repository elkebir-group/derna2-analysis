#!/usr/bin/env bash
# The binding screen's occurrence half for the band at 1000 aa and above.
#
# WHY THIS IS NOT screen_f0_long.sh. That script folds a whole baseline front F0
# per protein under a wall timeout, and records a protein as censored when the
# timeout kills it, because a killed front emits nothing. At these lengths EVERY
# protein would be censored: at 1082 aa a front reached 8 points and closed zero
# gaps in four hours, and a whole front is on the order of 200 core-hours. A
# screen that returns nothing for all nineteen proteins is not a screen.
#
# So `occurs` here is counted over a LAMBDA GRID of single designs rather than
# over a front. Each design is the exact optimum at its own weight, so every
# sequence scanned is a real optimum and an occurrence on one of them is a real
# occurrence. What changes is the denominator, and it changes in one direction
# only:
#
#   occurs > 0   means the same thing it means for a front -- the motif is spelled
#                by an optimum, so forbidding it must move that optimum
#   occurs = 0   is MUCH WEAKER here than on a whole front. Five weights are not
#                the hundreds of designs a front carries, so a zero is "none of
#                these five optima spells it", not "the front avoids it". Any
#                zero from this screen must be reported with the weight count
#                beside it, and must not be quoted as a front's zero.
#
# That asymmetry is the whole reason the column is worth having: a motif that
# DOES occur is settled cheaply here, and the ones that read zero are the only
# ones a front would have to settle.
#
# One file per (protein, weight), and the weights run one at a time rather than
# as a single comma list, so a kill loses the weight in flight instead of the
# protein. The parts live in a SUBDIRECTORY and are concatenated into one file
# per protein, because the occurrence counter globs `*.jsonl` and keys a protein
# by the file's basename: parts beside their merge would be counted twice and
# would not match the lattice table's accessions. Ascending in length, so the cheap majority lands first and a stop
# leaves a usable prefix of the corpus. Already finished pairs are skipped, so
# the script is re-runnable.
#
# Walls are NOT reportable: this runs on whatever cores are free beside other
# sessions' pinned work. The readouts are occurrence counts, which contention
# cannot perturb.
set -u
P=/scratch/xinyugu4/paper-e3
BIN=$P/bin/derna_06605630
CORPUS=$P/screen/corpus-xl
OUT=$P/screen/f0-xl
CORES=${CORES:?set CORES to a taskset mask this screen may use}
PER=${PER:-14400}                      # per design, seconds
LAMBDAS=${LAMBDAS:-"0.9 0.7 0.5 0.3 0.1"}
mkdir -p "$OUT/parts"
SUM=$OUT/summary.csv
[ -f "$SUM" ] || printf 'protein,aa,lambda,exit,wall_s,peak_rss_kb\n' > "$SUM"

echo "=== f0-xl screen: weights [$LAMBDAS], ${PER}s per design, cores $CORES, $(date -Is)"

for fa in $(for f in "$CORPUS"/*.fasta; do
              printf '%s %s\n' "$(grep -v '^>' "$f" | tr -d '\n\r ' | wc -c)" "$f"
            done | sort -n | cut -d' ' -f2); do
  acc=$(basename "$fa" .fasta)
  aa=$(grep -v '^>' "$fa" | tr -d '\n\r ' | wc -c)
  # ReverseCovid is a reversed synthetic sequence, kept in the corpus as a free
  # control and excluded from selection; the census over it is informative, so it
  # is screened like the rest and labelled where it is reported.
  for lam in $LAMBDAS; do
    o=$OUT/parts/$acc.l$lam.jsonl
    [ -s "$o" ] && continue
    start=$(date +%s)
    taskset -c "$CORES" timeout -k 30 "$PER" \
      "$BIN" design --protein "$fa" --lambda "$lam" --out "$o" --force \
      > "$OUT/parts/$acc.l$lam.log" 2>&1
    code=$?
    wall=$(( $(date +%s) - start ))
    if [ "$code" -ne 0 ] || [ ! -s "$o" ]; then
      rm -f "$o"
      printf '%s,%s,%s,%s,%s,\n' "$acc" "$aa" "$lam" "$code" "$wall" >> "$SUM"
      echo "=== $acc ($aa aa) lambda $lam exit $code after ${wall}s — censored" >&2
      # A protein whose FIRST weight is censored will censor the rest too, and the
      # corpus runs ascending, so stop this protein and move to the next.
      break
    fi
    kb=$(python3 -I -c "
import json
print(json.loads(open('$o').readline()).get('peak_rss_kb',''))")
    printf '%s,%s,%s,0,%s,%s\n' "$acc" "$aa" "$lam" "$wall" "$kb" >> "$SUM"
    echo "=== $acc ($aa aa) lambda $lam ${wall}s peak ${kb} kb" >&2
  done
  # One merged file per protein, named for the accession the lattice table uses.
  # Concatenation is additive and each part is a distinct weight, so this is not
  # a last-writer merge; the counter ignores the parts' header lines because it
  # keeps only records carrying a sequence.
  if ls $OUT/parts/$acc.l*.jsonl >/dev/null 2>&1; then
    cat $OUT/parts/$acc.l*.jsonl > $OUT/$acc.jsonl
    echo "=== $acc merged $(ls $OUT/parts/$acc.l*.jsonl | wc -l) weight(s) -> $acc.jsonl" >&2
  fi
done
echo "=== done $(date -Is)" >&2
