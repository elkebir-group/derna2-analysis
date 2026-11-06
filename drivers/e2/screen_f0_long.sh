#!/usr/bin/env bash
# The binding screen over the 501-700 aa band: the baseline front F0 per protein,
# whose designs are then scanned for motif occurrences. Occurrences on F0 are the
# only ones that can move the front, so they are what `occurs` means.
#
# A second corpus and a second output directory, not an extension of the first,
# so the 61-protein screen keeps its denominator. The two summaries are joined
# only where a statement is about both bands, and the band is always named.
#
# What this band answers is where the frontier sits on real proteins between the
# grid's ceiling and the shortest Cas ortholog. It does not feed the grid, which
# is fixed: see the plan under the selection rule.
#
# Runs ascending in length and writes one file per protein, so a stop at any
# point leaves a usable prefix of the corpus rather than nothing. Already
# finished proteins are skipped, so the script is re-runnable.
#
# Walls here are NOT reportable: the E-cores carry another session's build, and
# the screen's readouts are the occurrence counts, which load cannot perturb.
set -u

BIN=derna2/e2/bin/derna     # derna2 at 8603d37d
HERE=derna2/e3-plan/experiments/paper-benchmark
CORPUS=$HERE/screen/corpus-l
OUT=$HERE/screen/f0-l
CORES=24              # E2 holds core 24 alone; 16-23 and 25-27 belong to peers
MEMMAX=24G
PER=28800             # per protein, seconds; 1.8x the +1sd projection for 694 aa on the >= 400 aa fit

mkdir -p "$OUT"
SUM=$OUT/summary.csv
[ -f "$SUM" ] || printf 'protein,aa,exit,wall_s,points,hull,gaps\n' > "$SUM"

# Ascending length, so the cheap majority lands first.
for fa in $(for f in "$CORPUS"/*.fasta; do
              printf '%s %s\n' "$(grep -v '^>' "$f" | tr -d '\n\r ' | wc -c)" "$f"
            done | sort -n | cut -d' ' -f2); do
  acc=$(basename "$fa" .fasta)
  aa=$(grep -v '^>' "$fa" | tr -d '\n\r ' | wc -c)
  [ -s "$OUT/$acc.jsonl" ] && continue

  start=$(date +%s)
  systemd-run --user --scope -q -p MemoryMax=$MEMMAX -p MemorySwapMax=0 -- \
    taskset -c "$CORES" \
    timeout -k 30 "$PER" \
    "$BIN" front --protein "$fa" --threads 1 --out "$OUT/$acc.jsonl" --force \
      > "$OUT/$acc.log" 2>&1
  code=$?
  wall=$(( $(date +%s) - start ))

  if [ "$code" -ne 0 ]; then
    # A killed front emits nothing, so record the censoring and move on.
    rm -f "$OUT/$acc.jsonl"
    printf '%s,%s,%s,%s,,,\n' "$acc" "$aa" "$code" "$wall" >> "$SUM"
    echo "=== $acc ($aa aa) exit $code after ${wall}s — censored" >&2
    continue
  fi

  read -r pts hull gaps <<<"$(python3 -I -c "
import json,sys
recs=[json.loads(l) for l in open('$OUT/$acc.jsonl')]
hdr=next((r for r in recs if 'points' in r), {})
pt=[r for r in recs if 'sequence' in r]
print(hdr.get('points',len(pt)), sum(1 for r in pt if r.get('on_hull')), hdr.get('gaps',''))
")"
  printf '%s,%s,0,%s,%s,%s,%s\n' "$acc" "$aa" "$wall" "$pts" "$hull" "$gaps" >> "$SUM"
  echo "=== $acc ($aa aa) ${wall}s  $pts point(s), $hull on hull" >&2
done
echo "=== done $(date -Is)" >&2
