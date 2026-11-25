from fractions import Fraction as F
from itertools import combinations, product, permutations
from collections import Counter

def dot(a,b):
    return sum(x*y for x,y in zip(a,b))

def rank4_vertices(guards):
    """Vertices of {x in R4: guard*x<=1}. Solve each full-rank active quartet."""
    signed = [((-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4)),p)
              for p in permutations(range(4))]
    def det(m):
        return sum(s*m[0][p[0]]*m[1][p[1]]*m[2][p[2]]*m[3][p[3]]
                   for s,p in signed)
    verts=set()
    for rows in combinations(guards,4):
        d=det(rows)
        if not d:
            continue
        nums=tuple(det(tuple(tuple(1 if j==i else row[j] for j in range(4))
                             for row in rows)) for i in range(4))
        if d<0:
            d,nums=-d,tuple(-n for n in nums)
        if all(dot(nums,g)<=d for g in guards):
            verts.add(tuple(F(n,d) for n in nums))
    return verts

def roots():
    return {tuple((a if k==i else b if k==j else 0) for k in range(4))
            for i,j in combinations(range(4),2) for a,b in product((-1,1),repeat=2)}

def neg(p): return tuple(-x for x in p)

def valid(points,check_size=None):
    if check_size is not None:
        assert len(set(points))==check_size
    assert all(dot(p,p)==2 for p in points)
    assert all(dot(p,q)<=1 for p,q in combinations(points,2))
    return sum(neg(p) in points for p in points)//2

D=roots()
r=(1,1,0,0)
axis={tuple(s if i==j else 0 for i in range(4)) for j in range(4) for s in (-1,1)}
half=set(product((F(-1,2),F(1,2)),repeat=4))
U=axis|half
assert len(D)==len(U)==24
V=rank4_vertices(sorted(D-{r}))
assert V==U|{r}, (U-V,V-(U|{r}))
print("One deletion: 8,855 active systems; vertices=24 polar units + missing root")

directions=sorted({min(p,neg(p)) for p in U})
frames=[]
for frame in combinations(directions,4):
    if all(dot(p,q)==0 for p,q in combinations(frame,2)):
        frames.append(set(frame)|{neg(p) for p in frame})
even={h for h in half if sum(x<0 for x in h)%2==0}
odd=half-even
assert {frozenset(f) for f in frames}==set(map(frozenset,(axis,even,odd)))
for frame in frames:
    assert any(tuple(a+b for a,b in zip(p,q))==r for p,q in combinations(frame,2))
print("495 four-line choices: exactly axis/even/odd signed frames")
eq={p+(0,) for p in D}
caps=[(top,bottom) for top in frames for bottom in frames]
antipodes=Counter()
for top,bottom in caps:
    pts=eq|{p+(1,) for p in top}|{p+(-1,) for p in bottom}
    antipodes[valid(pts,40)]+=1
    valid(pts-{r+(0,)},39)
assert antipodes=={20:3,12:6}
print("Nine 40-point completions and nine 39-point examples omitting r;",dict(sorted(antipodes.items())))

def extra(k,sign=1):
    a,t=F(k*k-2,k*k+2),F(sign*4*k,k*k+2)
    return (a,a,0,0,t)
a,b,a_dn,hi=extra(4),extra(5),extra(4,-1),extra(3)
for x in (a,b,a_dn,hi):
    valid((eq-{r+(0,)})|{x},24)
    assert dot(x,r+(0,))>1
assert abs(a[-1])<1 and abs(b[-1])<1 and abs(hi[-1])>1
assert dot(a,b)>1 and dot(a,a_dn)==F(34,81)
for frame in frames:
    assert max(dot(a,p+(1,)) for p in frame)>1
    assert max(dot(b,p+(1,)) for p in frame)>1
print("Rational cone controls: same side conflict, opposite side dot=34/81; nonpolar high exists")

s=(1,0,1,0)
w=(1,F(1,2),F(1,2),0)
A=rank4_vertices(sorted(D-{r,s}))
print("Adjacent deletions vertex inventory",len(A),sorted(A-U),"missing polar",sorted(U-A))
# Frames admit a disjoint pair partition covered by fully retained roots at both signs.
for frame in frames:
    dirs=set(min(p,neg(p)) for p in frame)
    plist=sorted(dirs)
    matchings=0
    for j in (1,2,3):
        k,l=set(range(1,4))-{j}
        required={tuple(a*x+b*y for x,y in zip(plist[u],plist[v]))
                  for u,v in ((0,j),(k,l)) for a,b in product((-1,1),repeat=2)}
        assert len(required)==8 and required<=D
        matchings += required.isdisjoint({r,s})
    assert matchings>=1
print("All three triality frames admit a retained full-sign orthogonal matching")
assert A==U|{r,s}
# Adjacent two holes do not imply at most one low overall in each hemisphere.
x=(F(3,5),1,0,0,F(4,5))
y=(F(5,13),0,1,0,F(12,13))
extra_core={(z+(0,)) for z in D-{r,s}}
assert dot(x,y)==F(63,65)
valid(extra_core|{x,y},24)
assert dot(x,r+(0,))>1>=dot(x,s+(0,))
assert dot(y,s+(0,))>1>=dot(y,r+(0,))
print("Adjacent-hole control: two same-side lows compatible, dot=63/65")
# A low simultaneously violating both missing roots; its height is sqrt(2)/2.
seven={r+(0,),tuple(F(m) for m in (1,0,0,0,1)),tuple(F(m) for m in (0,1,0,0,1))}
seven |= {h+(1,) for h in half if h[:2]==(F(1,2),F(1,2))}
assert len(seven)==7
assert dot(w,w)+F(1,2)==2
assert all(dot(w,z)<=1 for z in D-{r,s})
assert dot(w,r)==dot(w,s)==F(3,2)
assert all(dot(x,x)==2 for x in seven)
assert all(dot(x,y)>=1 for x,y in combinations(seven,2))
for v in seven:
    # Require dot(w,v[:4])+v[4]*sqrt(2)/2 > 1, with exact sign test.
    z=1-dot(w,v[:4])
    assert z<0 or (v[4]>0 and F(v[4]**2,2)>z**2)
assert min(dot(w,z[:4]) for z in seven if z[-1])==F(1,2)
print("Mixed W: actual two-violation low, augmented ray pair dots certified over Q(sqrt(2))")
pairs=Counter(dot(p,q) for p,q in combinations(D,2))
assert pairs=={-2:12,-1:96,0:72,1:96}
print("Missing-pair inner products and counts:",dict(sorted(pairs.items())))

# Truncated polar vertex controls in the other pair types. Two signed-permutation
# templates have dot zero (same coordinate support, or disjoint).
templates=((neg(r),set()),((-1,0,1,0),set()),((1,-1,0,0),{(1,0,0,0)}),
           ((0,0,1,1),{(F(1,2),)*4}))
for k,nonvertex in templates:
    verts=rank4_vertices(sorted(D-{r,k}))
    assert verts==(U-nonvertex)|{r,k}
    assert all(dot(v,v)<2 for v in verts-{r,k})
print("Five signed-permutation pair templates: 26,26,26,25,25 vertices; norm-2 vertices are just the holes")
