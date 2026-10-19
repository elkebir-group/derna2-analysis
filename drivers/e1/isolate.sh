#!/usr/bin/env bash
# isolate.sh [--apply] [RESERVED] [KEEP]
#
# Reserve a CPU set for the time pass by moving this account's own CPU-heavy
# processes off it. Without root this is the only isolation available: the
# cpuset cgroup is root-owned, and `taskset` can move only processes this
# account owns. Another account's process can still land on the reserved set,
# which is why the time pass checks those cores before and during every run.
#
# Defaults reserve socket 1 (the odd CPUs, NUMA node 1) and leave socket 0
# (the even CPUs) to everything else. A whole socket is reserved rather than a
# few cores because the last-level cache and memory bandwidth are per socket:
# a neighbour on the same socket slows a run without ever touching its cores.
#
# Dry run by default: it lists the busy processes and what else is running.
# --apply moves EVERY process this account owns, every thread (`taskset -a -p`),
# not only the busy ones: a sweep driver idles at ~0% CPU between cells, and a
# driver left on all cores would launch its next binary onto the reserved set.
# Children inherit the new mask. The benchmark driver pins its own runs with
# an explicit `taskset`, so it is unaffected by being launched from a moved shell.
# Moving a running job changes that job's own timings from that point on.
set -uo pipefail

APPLY=0
if [ "${1:-}" = "--apply" ]; then APPLY=1; shift; fi
RESERVED=${1:-1,3,5,7,9,11,13,15,17,19,21,23,25,27}
KEEP=${2:-0,2,4,6,8,10,12,14,16,18,20,22,24,26}
MIN_PCPU=${MIN_PCPU:-5}

echo "reserved for timing: $RESERVED"
echo "everything else to:  $KEEP"
echo

echo "== this account's processes above ${MIN_PCPU}% CPU"
ps -u "$(id -u)" -o pid=,pcpu=,etime=,comm= --sort=-pcpu | awk -v m="$MIN_PCPU" '$2 > m' |
while read -r pid pcpu etime comm; do
    [ "$pid" = "$$" ] && continue
    now=$(taskset -c -p "$pid" 2>/dev/null | sed 's/.*: //')
    printf '%7s %6s%% %12s %-20s affinity %s\n' "$pid" "$pcpu" "$etime" "$comm" "$now"
done

if [ "$APPLY" -eq 1 ]; then
    moved=0; failed=0
    for pid in $(ps -u "$(id -u)" -o pid=); do
        if taskset -a -c -p "$KEEP" "$pid" > /dev/null 2>&1; then moved=$((moved + 1)); else failed=$((failed + 1)); fi
    done
    echo "moved $moved processes to $KEEP ($failed had exited or could not be moved)"
fi

echo
echo "== other accounts' processes above ${MIN_PCPU}% CPU (cannot be moved; ask the owner or an admin)"
ps -e -o user=,pid=,pcpu=,psr=,comm= --sort=-pcpu | awk -v me="$(id -un)" -v m="$MIN_PCPU" '$1 != me && $3 > m'

if [ "$APPLY" -eq 0 ]; then echo; echo "dry run: re-run with --apply to move the processes listed above"; fi
