#!/usr/bin/env python3
"""Merge seeds + window finds into the SO front; report coverage of [min, max] adaptation."""
import json, glob, os, sys
S = os.path.dirname(os.path.abspath(__file__))
pts = {}
for l in open(f'{S}/seed_stair.jsonl'):
    d = json.loads(l); pts[d['sequence']] = (d['cai_log'], d['mfe_kcal_mol'], 'seed')
reports = []
for p in sorted(glob.glob(f'{S}/r*/w*.jsonl') + glob.glob(f'{S}/top/*.jsonl')):
    d = json.loads(open(p).readlines()[-1])
    reports.append((os.path.relpath(p, S), d['complete'], d['lo'], d['hi'], d.get('open')))
    for q in d['front']:
        pts.setdefault(q['sequence'], (q['cai_log'], q['energy_kcal_mol'], os.path.relpath(p, S)))
front, best = [], float('inf')
for s, (a, e, src) in sorted(pts.items(), key=lambda kv: (-kv[1][0], kv[1][1])):
    if e < best - 1e-9:
        front.append({'cai_log': a, 'mfe_kcal_mol': e, 'sequence': s, 'source': src}); best = e
front.reverse()
for r in reports: print(*r)
print('distinct SO designs', len(pts), 'front points', len(front))
with open(f'{S}/so_front.jsonl', 'w') as f:
    for p in front: f.write(json.dumps(p) + '\n')
