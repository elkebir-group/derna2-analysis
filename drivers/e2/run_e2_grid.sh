#!/bin/bash
# The whole E2 grid, in the order the passes depend on each other.
#
#   run_e2_grid.sh CONFIG.json [--max-n N]
#
# Four passes, and the order is not a preference:
#
#   1. hull  --phase single   the reference walk; every replay arm's lambdas
#   2. sweep --phase single   the single-core grid, at the piloted concurrency
#   3. hull  --phase par      the multi-core walk, holding every timing core
#   4. sweep --phase par      the multi-core replay arms, one cell at a time
#
# The multi-core passes run last and alone: a cell on eight cores is every
# timing core this host has, so anything beside it is measured under
# contention the grid does not otherwise have.
#
# Each pass is resumable. A pass re-run skips the cells already recorded
# against this pin, so an interruption costs the cell it was in and nothing
# else.
set -u -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:?usage: run_e2_grid.sh CONFIG.json [--max-n N]}"
shift
EXTRA=("$@")

OUT=$(python3 -I -c "import json,sys;print(json.load(open(sys.argv[1]))['out_dir'])" "$CONFIG")
mkdir -p "$OUT"
LOG="$OUT/grid.log"

say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }

# The timing cores must be idle, and so must their hyperthread siblings: two
# threads of one core share its execution units, so a job on the sibling slows
# the cell without ever appearing on the cell's own core. A sibling is not
# visible in a cell's numbers, and every cell measured beside one is suspect.
CORES=$(python3 -I -c "import json,sys;print(','.join(map(str,json.load(open(sys.argv[1]))['cores'])))" "$CONFIG")
WATCH=$(python3 -I - "$CORES" <<'PY'
import sys, pathlib
watch = set()
for cpu in sys.argv[1].split(","):
    watch.add(int(cpu))
    siblings = pathlib.Path(
        f"/sys/devices/system/cpu/cpu{int(cpu)}/topology/thread_siblings_list")
    if siblings.exists():
        for part in siblings.read_text().strip().replace("-", ",").split(","):
            watch.add(int(part))
print(",".join(str(c) for c in sorted(watch)))
PY
)
# Sampled over a window rather than once: a job that lands on a watched core
# between two instants is still a job that will land on it during a cell, and
# one sample cannot tell an idle core from an unlucky moment.
busy=""
for _ in 1 2 3 4 5 6; do
  busy+=$(ps -eo psr,pcpu,comm --no-headers | awk -v c="$WATCH" '
    BEGIN { n = split(c, a, ","); for (i = 1; i <= n; i++) mine[a[i]] = 1 }
    $2 > 20 && ($1 in mine) { print $1 ":" $3 }')' '
  sleep 2
done
busy=$(echo "$busy" | tr ' ' '\n' | grep -v '^$' | sort -u | tr '\n' ' ')
if [ -n "$busy" ]; then
  say "cores $WATCH (timing cores and their siblings) are busy: $busy"
  say "refusing to start; quiet the host or edit the config's cores"
  exit 1
fi
say "timing cores $CORES idle, with siblings ($WATCH)"
say "config $CONFIG, output $OUT"

python3 "$HERE/drive_e2.py" plan --config "$CONFIG" "${EXTRA[@]}" | tee -a "$LOG"

for pass in "hull single" "sweep single" "hull par" "sweep par"; do
  set -- $pass
  command=$1
  phase=$2
  say "=== $command --phase $phase"
  started=$(date +%s)
  if ! python3 "$HERE/drive_e2.py" "$command" --config "$CONFIG" \
      --phase "$phase" "${EXTRA[@]}"; then
    say "$command --phase $phase failed; the passes after it need it, stopping"
    exit 1
  fi
  say "=== $command --phase $phase done in $(( ($(date +%s) - started) / 60 )) min"
done

say "=== plan after the run"
python3 "$HERE/drive_e2.py" plan --config "$CONFIG" "${EXTRA[@]}" | tee -a "$LOG"
say "grid finished"
