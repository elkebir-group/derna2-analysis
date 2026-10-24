#!/usr/bin/env python3
"""Exact SO front over every design with cai_log >= -0.5816: descending CAI, guarded fold only where F beats the best SO energy so far."""
import json, glob, subprocess, os
# FIG1_WORK: the working directory the run used (binaries in bin/, outputs beside them).
W = os.environ.get('FIG1_WORK', 'work')
D=W + '/top_exh'
cai={}; seq={}
for l in open(f'{D}/designs.tsv'):
    n,c,s=l.rstrip('\n').split('\t'); cai[n]=float(c); seq[n]=s
F={}
for p in glob.glob(f'{D}/shard*.free.jsonl'):
    for l in open(p):
        d=json.loads(l); F[d['name']]=round(d['mfe_kcal_mol'],2)
assert len(F)==len(cai), (len(F),len(cai))
order=sorted(cai, key=lambda n:(-cai[n], F[n]))
best=float('inf'); front=[]; nguard=0; i=0
def guard(names):
    with open(f'{D}/q.tsv','w') as f:
        for n in names: f.write(f'{n}\t{seq[n]}\n')
    out=subprocess.run(['taskset','-c','0-26:2',W + '/bin/guard_fold','--in',f'{D}/q.tsv','--utr5-len','50','--utr3-len','132','--threads','12'],capture_output=True,text=True,check=True).stdout
    return {l.split('\t')[0]:int(l.split('\t')[1])/100 for l in out.splitlines()}
while i<len(order):
    batch=[n for n in order[i:i+2000] if F[n]<best-1e-9][:400]
    if not batch:
        i+=2000; continue
    E=guard(batch); nguard+=len(batch)
    last=order.index(batch[-1])
    for n in order[i:last+1]:
        if F[n]<best-1e-9 and n in E and abs(E[n]-F[n])<1e-6:
            # same-CAI tie: keep only strictly better
            front.append({'cai_log':cai[n],'mfe_kcal_mol':F[n],'sequence':seq[n],'name':n}); best=F[n]
    i=last+1
print('guarded folds',nguard,'SO front points in top stretch',len(front))
json.dump({'free_min':min(F.values()),'free_max':max(F.values()),'n':len(F)},open(f'{D}/stats.json','w'))
with open(f'{D}/top_front.jsonl','w') as f:
    for p in front: f.write(json.dumps(p)+'\n')
for p in front: print(round(p['cai_log'],4),p['mfe_kcal_mol'])
