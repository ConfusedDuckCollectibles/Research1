"""Exact finite controls for the analytic triangle-hole argument."""
from fractions import Fraction as F
from itertools import combinations, permutations, product

D = {tuple(i if k == a else j if k == b else 0 for k in range(4))
     for a, b in combinations(range(4), 2) for i, j in product((-1,1), repeat=2)}
E = tuple(tuple(int(i == j) for i in range(4)) for j in range(4))
dot = lambda u,v: sum(a*b for a,b in zip(u,v))
neg = lambda v: tuple(-x for x in v)
sup = lambda v: frozenset(i for i,x in enumerate(v) if x)
def check_code(points):
    assert len(set(points)) == len(points)
    assert all(dot(x,x)==2 for x in points)
    assert all(dot(x,y)<=1 for x,y in combinations(points,2))

frame = [tuple(v for e in E for v in (e,neg(e)))]
frame.extend(tuple(tuple(F(i,2) for i in sig) for sig in product((-1,1),repeat=4)
                   if sum(i == -1 for i in sig)%2 == par) for par in range(2))
matrix=tuple(tuple(F(i,2) for i in row) for row in
             ((1,1,1,1),(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1)))
assert {tuple(dot(r,row) for row in matrix) for r in D} == D
assert all(dot(u,v) in (-F(1,2),F(1,2))
           for a,b in combinations(frame,2) for u,v in product(a,b))

triangles = {frozenset(t) for t in combinations(D,3)
             if all(dot(x,y)==1 for x,y in combinations(t,2))}
M = frozenset(tuple(int(k in (i,j)) for k in range(4))
              for i,j in combinations(range(3),2))
Mf = frozenset(tuple(dot(r,row) for row in matrix) for r in M)
def orbit(t):
    return {frozenset(tuple(r[i]*s for i,s in zip(p,sgn)) for r in t)
            for p in permutations(range(4)) for sgn in product((-1,1),repeat=4)}
assert orbit(M) | orbit(Mf) == triangles and not orbit(M) & orbit(Mf)
assert (len(orbit(M)),len(orbit(Mf))) == (32,64)

for t in triangles:
    faces=0
    for G in frame:
        B=tuple(g for g in G if next(i for i in g if i)>0)
        assert len(B)==4 and all(dot(a,b)==int(i==j)
                                for i,a in enumerate(B) for j,b in enumerate(B))
        SS = [sup(tuple(dot(r,b) for b in B)) for r in t]
        assert len(set(SS))==3 and all(len(s)==2 for s in SS)
        faces += len(set.union(*map(set,SS))) == 3
    assert faces==1
print("96 triangles; 32 axial faces + 64 axial stars; one face frame each")

core=[r+(0,) for r in sorted(D-M)]
cap_pos=[e+(1,) for e in frame[0]]
cap_neg=[e+(-1,) for e in frame[0]]
X=(F(1),F(3,5),F(4,5),F(0),F(0))
assert {r for r in D if dot(r,X)>1}==M
check_code(core+cap_pos+cap_neg+[X])
# Symbolic vertical coordinate +/- sqrt(23)/5; all dots to roots use u only.
u=(F(3,5),)*3+(0,)
assert dot(u,u)+F(23,25)==2 and dot(u,u)-F(23,25)==F(4,25)
assert {r for r in D if dot(r,u)>1}==M
print("Exact 38-point full-cap witness and opposite-side triple-low pair")

# Nonpolar seven-high example and axial cases in each frame.
w=(F(3,5),F(4,5),F(0),F(0))
high=[e+(1,) for e in frame[0] if sup(e) <= {2,3}] + [
    neg(E[0])+(1,),neg(E[1])+(1,),w+(1,)]
check_code(core+high)
assert {r for r in D if dot(w,r)>1}=={r for r in M if sup(r)=={0,1}}
for G in frame:
    for missing in G:
        kept=set(G)-{missing}
        assert all(any(dot(r,w)==1 for w in kept) for r in M)
print("Exact nonpolar seven-high control; all 24 polar seven-high subframes")
# Height-one triple violation blocks rays strictly. A W dot has rational
# horizontal part and positive 1/sqrt(2), so a horizontal value >1 is enough.
v=(F(4,7),F(4,7),F(4,7),F(1,7))
check_code(core+[v+(1,)])
assert {r for r in D if dot(r,v)>1}==M
for r in M:
    assert all(dot(v,g)>0 for g in frame[0] if dot(r,g)==1)
    assert all(dot(v,g)>0 for G in frame[1:] for g in G if dot(r,g)==1)
    s,t=M-{r}
    ws=tuple(F(a+b,2) for a,b in zip(r,s))
    wt=tuple(F(a+b,2) for a,b in zip(r,t))
    assert dot(ws,wt)+F(1,2)==F(7,4)
    assert dot(v,ws)>1 and dot(v,wt)>1
print("Triple-violating height-one control and all W/W products 7/4")
