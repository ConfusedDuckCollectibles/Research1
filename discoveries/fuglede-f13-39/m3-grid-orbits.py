from pathlib import Path
import sys
from collections import Counter
p=13
inv=[0]+[pow(i,-1,p) for i in range(1,p)]
def canonical(S):
    c=min(S) # centroid is not canonical: 39 = 0 in F_13.
    # The argument below uses the unique center of reflection whenever present.
    centers=[(x,y) for x,y in S if {((2*x-u)%p,(2*y-v)%p) for u,v in S}==S]
    if len(centers)!=1:
        return None
    cx,cy=centers[0]
    SS=[((x-cx)%p,(y-cy)%p) for x,y in S]
    best=None
    for a,b in SS:
        for d,e in SS:
            k=(a*e-b*d)%p
            if k:
                trial=tuple(sorted(13*((e*x-d*y)*inv[k]%p)+(-b*x+a*y)*inv[k]%p for x,y in SS))
                if best is None or trial<best:best=trial
    return best
hist=Counter()
for s in Path(sys.argv[1]).read_text().splitlines():
    if s.startswith("grid "):
        r=canonical({divmod(int(t),13) for t in s.split()[1:]})
        hist[r]+=1
for r,c in hist.items(): print(c,r)
