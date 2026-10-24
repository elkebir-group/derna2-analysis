import json,sys
def load(p,ek='mfe_kcal_mol'):
    pts={}
    for l in open(p):
        d=json.loads(l); L=round(d.get('cai_log', 0),4); E=d.get(ek, d.get('mfe'))
        pts[L]=(E,d.get('sequence') or d.get('rna_seq'))
    return pts
old=load(sys.argv[1]); new=load(sys.argv[2])
print(f'old {len(old)} pts, new {len(new)} pts')
common=set(old)&set(new)
dE=[new[k][0]-old[k][0] for k in common]
print('shared CAI levels',len(common),'energy changed',sum(abs(x)>1e-6 for x in dE), 'min dE',min(dE) if dE else None)
print('same seq at shared', sum(old[k][1]==new[k][1] for k in common))
print('only old',len(set(old)-set(new)),'only new',len(set(new)-set(old)))
