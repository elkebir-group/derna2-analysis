import csv, math, collections
T='derna2/packages/derna-codon-dp/data/codon_usage_human.csv'
rows=list(csv.reader(open(T)))
hdr=rows[0]; rows=rows[1:]
def build():
    by=collections.defaultdict(dict)
    for r in rows:
        d=dict(zip(hdr,r))
        cod=d.get('codon') or r[0]; aa=d.get('aa') or d.get('amino_acid') or r[1]; f=float(d.get('fraction') or r[2])
        by[aa][cod.replace('T','U')]=f
    w={}
    for aa,cs in by.items():
        m=max(cs.values())
        for c,f in cs.items(): w[c]=f/m
    return w
W=build()
def cai_log(cds): return sum(math.log(W[cds[i:i+3]]) for i in range(0,len(cds),3))
