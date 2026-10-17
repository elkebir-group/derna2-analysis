import csv,math,collections,json,os
# FIG1_WORK: the working directory the run used (binaries in bin/, outputs beside them).
W = os.environ.get('FIG1_WORK', 'work')
R='/home/xinyugu4/derna2'
by=collections.defaultdict(list)
for r in list(csv.reader(open(f'{R}/packages/derna-codon-dp/data/codon_usage_human.csv')))[1:]: by[r[0]].append((r[1].replace('T','U'),float(r[2])))
prot=''.join(l.strip() for l in open(f'{R}/experiments/scope-25-beam-runtime-quality/fixtures/GLPE_HUMAN_78.fasta') if not l.startswith('>'))[:78]
u5=open(f'{R}/experiments/scope-19-utr-support/m5/hbb5.txt').read().strip(); u3=open(f'{R}/experiments/scope-19-utr-support/m5/hbb3.txt').read().strip()
T=0.5815909122842027+1e-9
opts=[]
for a in prot:
    m=max(f for _,f in by[a]); opts.append(sorted(((-math.log(f/m),c) for c,f in by[a])))
# suffix min penalty is 0; DFS
out=[]
def dfs(i,pen,acc):
    if i==len(prot): out.append((pen,''.join(acc))); return
    for p,c in opts[i]:
        if pen+p<=T: acc.append(c); dfs(i+1,pen+p,acc); acc.pop()
        else: break
dfs(0,0.0,[])
print(len(out))
K=12
fs=[open(f'{W}/top_exh/shard{k}.tsv','w') for k in range(K)]
with open(W + '/top_exh/designs.tsv','w') as d:
    for i,(p,c) in enumerate(out):
        d.write(f'd{i}\t{-p:.9f}\t{u5+c+u3}\n'); fs[i%K].write(f'd{i}\t{u5+c+u3}\n')
