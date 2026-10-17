#!/usr/bin/env bash
for k in $(seq 0 11); do taskset -c 0-26:2 ${FIG1_WORK:-work}/bin/fold_seq --in ${FIG1_WORK:-work}/top_exh/shard$k.tsv --out ${FIG1_WORK:-work}/top_exh/shard$k.free.jsonl > /dev/null 2>&1 & done; wait; echo "free folds done"
