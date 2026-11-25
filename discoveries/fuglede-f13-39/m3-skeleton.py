"""Exact planar and small-profile controls; run with Python 3.10+."""
from itertools import combinations, product
from pathlib import Path

P = 13
INF = P
points = list(product(range(P), repeat=2))
N = len(points)
index = {pt:i for i,pt in enumerate(points)}
base = tuple(index[p] for p in [(0,0), (1,0), (0,1)])
dirs = [[None for _ in points] for _ in points]
vecs = [[None for _ in points] for _ in points]
line = [[0]*(P+1) for _ in points]
for i,(a,b) in enumerate(points):
    for j,(c,d) in enumerate(points):
        if i==j:
            continue
        x,y=(c-a)%P,(d-b)%P
        dirs[i][j]=INF if x==0 else y*pow(x,-1,P)%P
        vecs[i][j]=min((x,y),(-x%P,-y%P))
        line[i][dirs[i][j]]|=1<<j
base_dirs = {dirs[i][j] for i,j in combinations(base,2)}


def find_cap(goal, allowed=range(P+1), monochord=False):
    """Return a cap if one exists up to affine normalization of its triangle."""
    adj = [sum(line[i][s] for s in allowed) for i in range(N)]
    avail = (1<<N)-1
    for i in base:
        avail &= adj[i]
    for i,j in combinations(base,2):
        avail &= ~line[i][dirs[i][j]]
    initial = {dirs[i][j]:vecs[i][j] for i,j in combinations(base,2)}
    def rec(chosen,left,used):
        if len(chosen)==goal:
            return chosen
        while len(chosen)+left.bit_count()>=goal:
            bit=left & -left
            l=bit.bit_length()-1
            left^=bit
            insert={dirs[i][l]:vecs[i][l] for i in chosen}
            assert len(insert)==len(chosen)
            if monochord and any(k in used and used[k]!=v for k,v in insert.items()):
                continue
            nxt=left & adj[l]
            for i in chosen:
                nxt &= ~line[i][dirs[i][l]]
            found=rec((*chosen,l),nxt,used|insert if monochord else used)
            if found:
                return found
        return None
    ans=rec(base,avail,initial)
    return None if ans is None else [points[i] for i in ans]



def seven_caps():
    avail=(1<<N)-1
    for i in base:
        avail &= ~(1<<i)
    for i,j in combinations(base,2):
        avail &= ~line[i][dirs[i][j]]
    initial={dirs[i][j]:vecs[i][j] for i,j in combinations(base,2)}
    def rec(chosen,left,used):
        if len(chosen)==7:
            yield chosen
            return
        while len(chosen)+left.bit_count()>=7:
            bit=left & -left
            k=bit.bit_length()-1
            left^=bit
            insert={dirs[i][k]:vecs[i][k] for i in chosen}
            if any(d in used and used[d]!=v for d,v in insert.items()): continue
            nxt=left
            for i in chosen:
                nxt &= ~line[i][dirs[i][k]]
            yield from rec((*chosen,k),nxt,used|insert)
    return rec(base,avail,initial)

def orient(seq):
    groups = {}
    for i,j in combinations(range(7),2):
        x,y = points[seq[i]]
        a,b = points[seq[j]]
        d = dirs[seq[i]][seq[j]]
        dlt=((a-x)%13,(b-y)%13)
        group=groups.setdefault(d,[])
        group.append((i,j) if vecs[seq[i]][seq[j]]==dlt else (j,i))
    return groups

def heights(groups):
    # Affine shear makes h0=h1=h2=0: they're over a fixed affine frame.
    for rest in product(range(13),repeat=4):
        h=(0,0,0)+rest
        if all(len({(h[j]-h[i])%13 for i,j in pairs})==1 for pairs in groups.values()):
            yield h

V = list(product(range(13),repeat=3))
V = [i for i in V if next((t for t in i if t),0)==1]
vindex={v:i for i,v in enumerate(V)}
def proj(h,x,y):
    d=pow(next(t for t in (h,x,y) if t),-1,13)
    return vindex[(h*d%13,x*d%13,y*d%13)]
inc=[sum(1<<i for i,(a,b,c) in enumerate(V) if (u*a+v*b+w*c)%13==0) for u,v,w in V]
def pencilmax(S,h):
    mask=1<<proj(1,0,0)
    for i,j in combinations(range(7),2):
        x,y=points[S[j]]
        a,b=points[S[i]]
        for k in range(-2,3):
            mask|=1<<proj((h[j]-h[i]+k)%13,(x-a)%13,(y-b)%13)
    return max((mask&l).bit_count() for l in inc)

step_set={i%13 for i in range(-2,3)}
aff={(a,b) for a in range(1,13) for b in range(13) if {(a*i+b)%13 for i in step_set}==step_set}
assert aff=={(1,0),(12,0)}
n=0
hist={}
for S in seven_caps():
    g=orient(S)
    for h in heights(g):
        m=pencilmax(S,h)
        hist[m]=1+hist.get(m,0)
    n+=1
print('7-cap normalized skeletons',n,'height histogram',hist)
