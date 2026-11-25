#!/usr/bin/env python3
"""Exact planar lower-direction gates over F_13, requiring a normalized triple."""
from itertools import combinations
p = 13
points = [(x,y) for x in range(p) for y in range(p)]
n = len(points)
inv = [0]+[pow(i,-1,p) for i in range(1,p)]
def slope(P,Q):
    x=(P[0]-Q[0])%p
    return (P[1]-Q[1])*inv[x]%p if x else p
ds = [[slope(P,Q) for Q in points] for P in points]
triple = [0,1,p]
base = {0,p-1,p}
def gate(goal, directions, cap):
    stats = [0]*(n+1)
    survivors = []
    for ex in combinations(sorted(set(range(p+1))-base),directions-len(base)):
        mask=base.union(ex)
        allowed = [sum(1<<j for j in range(n) if i != j and ds[i][j] in mask) for i in range(n)]
        cand=(1<<n)-1
        for i in triple: cand &= allowed[i]
        def walk(chosen, can):
            stats[len(chosen)]+=1
            if len(chosen)==goal:
                survivors.append(chosen)
                return
            while can.bit_count() >= goal-len(chosen):
                lsb=can&-can; i=lsb.bit_length()-1; can^=lsb
                if cap < goal and max(sum(ds[i][j]==s for j in chosen) for s in set(ds[i][j] for j in chosen)) >= cap:
                    continue
                walk(chosen+[i],can & allowed[i])
        walk(triple,cand)
    print((goal,directions,cap), "states",stats[3:goal+1],"solutions",len(survivors),
          "first",[[points[i] for i in s] for s in survivors[:1]],flush=True)
    return survivors
if __name__ == '__main__':
    for gate_parameters in ((8, 5, 5), (8, 7, 3), (9, 7, 4)):
        if gate(*gate_parameters):
            raise SystemExit("Unexpected candidate")
