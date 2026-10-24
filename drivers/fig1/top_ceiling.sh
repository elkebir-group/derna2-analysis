#!/usr/bin/env bash
C=$1
SEQ=$(python3 -c "import json; print(json.loads(open('${FIG1_WORK:-work}/so/seed_stair.jsonl').readlines()[-1])['sequence'])")
taskset -c 0-26:2 /usr/bin/time -v ${FIG1_WORK:-work}/bin/scope37_window --fasta /home/xinyugu4/derna2/experiments/scope-25-beam-runtime-quality/fixtures/GLPE_HUMAN_78.fasta --n 78  --utr5 $(cat /home/xinyugu4/derna2/experiments/scope-19-utr-support/m5/hbb5.txt) --utr3 $(cat /home/xinyugu4/derna2/experiments/scope-19-utr-support/m5/hbb3.txt)  --table symmetric --lo-seq $SEQ --hi 0.000001 --ceiling $C --scan-threads 3 --verify-threads 3 --max-secs 3600 --max-rss-gb 30 --energy-slopes 0,0.25x,0.5x,2x,4x  --out ${FIG1_WORK:-work}/so/top/c$C.jsonl > ${FIG1_WORK:-work}/so/top/c$C.log 2>&1
echo "top ceiling $C exit $?"
