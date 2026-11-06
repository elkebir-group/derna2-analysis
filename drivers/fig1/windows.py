#!/usr/bin/env python3
"""Capped windows between adjacent known SO designs (symmetric guard), plus the bottom stretch.
Usage: windows.py STAIR.jsonl OUTDIR JOBS [k ...]   (k = window index; -1 = bottom)"""
import json, subprocess, sys, os, math, csv, collections
from concurrent.futures import ThreadPoolExecutor
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = 'derna2'
BIN = f'{S}/bin/scope37_window'
FA = f'{R}/experiments/scope-25-beam-runtime-quality/fixtures/GLPE_HUMAN_78.fasta'
U5 = open(f'{R}/experiments/scope-19-utr-support/m5/hbb5.txt').read().strip()
U3 = open(f'{R}/experiments/scope-19-utr-support/m5/hbb3.txt').read().strip()
FRONT = f'{S}/symmetric/front.jsonl'
stair = [json.loads(l) for l in open(sys.argv[1])]
out = sys.argv[2]; jobs = int(sys.argv[3]); want = [int(x) for x in sys.argv[4:]]
os.makedirs(out, exist_ok=True)
# least possible adaptation, for the bottom window's lo
by = collections.defaultdict(list)
for r in list(csv.reader(open(f'{R}/packages/derna-codon-dp/data/codon_usage_human.csv')))[1:]:
    by[r[0]].append(float(r[2]))
prot = ''.join(l.strip() for l in open(FA) if not l.startswith('>'))[:78]
LO = sum(math.log(min(by[a]) / max(by[a])) for a in prot) - 0.01
common = ['taskset', '-c', '0-26:2', BIN, '--fasta', FA, '--n', '78', '--utr5', U5, '--utr3', U3,
          '--table', 'symmetric', '--slope', 'band', '--front', FRONT, '--scan-threads', '3',
          '--verify-threads', '3', '--max-secs', '3600', '--max-rss-gb', '30'] + os.environ.get('EXTRA', '').split()
def job(k):
    tag = f'w{k}'
    if k == -1:
        args = ['--lo', str(LO), '--incumbent', stair[0]['sequence']]
    else:
        args = ['--lo-seq', stair[k]['sequence'], '--incumbent', stair[k + 1]['sequence']]
    rep = f'{out}/{tag}.jsonl'
    if os.path.exists(rep):
        return tag, 'skip'
    with open(f'{out}/{tag}.log', 'w') as fh:
        r = subprocess.run(['/usr/bin/time', '-v'] + common + args + ['--out', rep], stdout=fh, stderr=subprocess.STDOUT)
    line = json.loads(open(rep).readlines()[-1]) if os.path.exists(rep) else {}
    print(tag, 'exit', r.returncode, 'complete', line.get('complete'), 'pts', len(line.get('front', [])), flush=True)
    return tag, r.returncode
ks = want or [-1] + list(range(len(stair) - 1))
with ThreadPoolExecutor(jobs) as ex:
    list(ex.map(job, ks))
print('windows all done', flush=True)
