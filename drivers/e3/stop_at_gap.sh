#!/usr/bin/env bash
# Stop a motif cell AT A GAP BOUNDARY, so the checkpoint captures the gap in
# flight instead of discarding it.
#
# A checkpoint is written only when a gap settles, so a SIGTERM mid-gap loses
# everything since the last one -- 9 to 40 minutes on the cells this is applied
# to. Waiting for the boundary keeps that work; waiting without bound defeats the
# point of stopping. So: poll for the settled-gap count to rise, signal at once
# when it does, and signal anyway at the deadline.
#
# The deadline is an argument rather than a policy, because what counts as too
# long depends on the cell: these marginal gap costs span 473 s to 3163 s and are
# climbing, so a figure that is patient for one is reckless for another.
#
# Signals the PAYLOAD, not the wrappers. timeout and /usr/bin/time then reap, the
# driver subshell writes its summary row, and the claim is released by its owner.
#
# No apostrophe appears in any ${VAR:?...} word here: one inside that word opens a
# quote that swallows the rest of the file, which bash -n reports as an EOF error
# forty lines later.
set -u
PID=${1:?pid of the derna payload}
OUT=${2:?path to the cells --out jsonl}
DEADLINE=${3:?seconds to wait for a gap boundary}
LBL=$(basename "$OUT" .jsonl)

gaps() {
  python3 -I -c "
import json
try:
    h=next(r for r in (json.loads(l) for l in open('$OUT')) if 'points' in r)
    print(h.get('gaps') or 0)
except Exception:
    print(-1)" 2>/dev/null
}

start=$(date +%s)
g0=$(gaps)
echo "=== $LBL: pid $PID at $g0 gap(s); waiting up to ${DEADLINE}s for the next to close"
while kill -0 "$PID" 2>/dev/null; do
  g=$(gaps)
  if [ "${g:-0}" -gt "${g0:-0}" ] 2>/dev/null; then
    echo "=== $LBL: gap closed ($g0 -> $g) after $(( $(date +%s) - start ))s; stopping now"
    kill "$PID" 2>/dev/null
    exit 0
  fi
  if [ $(( $(date +%s) - start )) -ge "$DEADLINE" ]; then
    echo "=== $LBL: no gap closed within ${DEADLINE}s; stopping anyway, losing the one in flight"
    kill "$PID" 2>/dev/null
    exit 0
  fi
  sleep 15
done
echo "=== $LBL: exited on its own before the deadline"
