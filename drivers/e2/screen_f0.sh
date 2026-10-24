#!/usr/bin/env bash
# The binding screen's second half: the baseline front F0 per protein, whose
# designs are then scanned for motif occurrences. Occurrences on F0 are the only
# ones that can move the front, so they are what `occurs` means.
#
# Runs ascending in length and writes one file per protein, so a stop at any
# point leaves a usable prefix of the corpus rather than nothing. Already
# finished proteins are skipped, so the script is re-runnable.
#
# Walls here are NOT reportable: the E-cores carry another session's build, and
# the screen's readouts are the occurrence counts, which load cannot perturb.
set -u

BIN=/scratch/xinyugu4/paper-e3/bin/derna_8603d37d   # the local build; a binary built elsewhere may not link here
HERE=/scratch/xinyugu4/paper-e3
CORPUS=$HERE/screen/corpus
OUT=$HERE/screen/f0
CORES=1,3,5,7,9,11,13,15   # socket 1, NUMA node1; a mask from another host means nothing
# Single-threaded, and it stays that way: the 61 baselines in hand were walked at
# --threads 1 without --deterministic, so a re-run at a wider setting could return
# a different front and would not be comparable with them.
MEMMAX=24G
PER=21600             # per protein, seconds; the tamazirt cap at ada's measured 1.53x and a margin, not a fitted exponent

# The cap is the kernel's where a user scope can be created, and nothing where
# systemd has no user bus — then `--max-rss-gb`, read between folds, is the only
# guard and the host's headroom is what stands behind it.
if systemd-run --user --scope -q -p MemoryMax=$MEMMAX -- /bin/true 2>/dev/null; then
  GUARD=(systemd-run --user --scope -q -p "MemoryMax=$MEMMAX" -p MemorySwapMax=0 --)
else
  GUARD=()
  echo "=== no systemd user bus: running without a kernel memory cap" >&2
fi

mkdir -p "$OUT"
SUM=$OUT/summary.csv
[ -f "$SUM" ] || printf 'protein,aa,exit,wall_s,points,hull,gaps\n' > "$SUM"

# Ascending length, so the cheap majority lands first.
for fa in $(for f in "$CORPUS"/*.fasta; do
              printf '%s %s\n' "$(grep -v '^>' "$f" | tr -d '\n\r ' | wc -c)" "$f"
            done | sort -n | cut -d' ' -f2); do
  acc=$(basename "$fa" .fasta)
  aa=$(grep -v '^>' "$fa" | tr -d '\n\r ' | wc -c)
  # Output alone is not enough to skip on: a pass interrupted between the walk
  # returning and the row being written leaves a cell every later pass skips and
  # no row anywhere. Require both.
  if [ -s "$OUT/$acc.jsonl" ] && grep -q "^$acc," "$SUM"; then
    continue
  fi

  start=$(date +%s)
  # `${A[@]+"${A[@]}"}`, not `"${A[@]}"`: bash below 4.4 — ada runs 4.2 — treats
  # an empty array as an unbound variable under `set -u`.
  ${GUARD[@]+"${GUARD[@]}"} \
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
