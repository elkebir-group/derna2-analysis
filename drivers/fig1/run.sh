#!/usr/bin/env bash
set -u
S=${FIG1_WORK:-work}
FA=/home/xinyugu4/derna2/experiments/scope-25-beam-runtime-quality/fixtures/GLPE_HUMAN_78.fasta
U5=$(cat /home/xinyugu4/derna2/experiments/scope-19-utr-support/m5/hbb5.txt); U3=$(cat /home/xinyugu4/derna2/experiments/scope-19-utr-support/m5/hbb3.txt)
run(){ arm=$1; shift; mkdir -p $S/$arm; taskset -c 0-26:2 /usr/bin/time -v $S/bin/scope35_walk_front --fasta $FA --n 78 --out $S/$arm/report.json --front $S/$arm/front.jsonl "$@" 2> $S/$arm/stderr.log; echo "$arm done"; }
run cdsonly &
run unguarded --utr5 $U5 --utr3 $U3 &
run symmetric --utr5 $U5 --utr3 $U3 --region 5-cds,cds-3 &
run full --utr5 $U5 --utr3 $U3 --region 5-cds,cds-3,5-3 &
wait; echo all-finished
