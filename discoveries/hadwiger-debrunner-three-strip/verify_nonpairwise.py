#!/usr/bin/env python3
"""Independent exact verifier for a nonpairwise sharp ten-hexagon (4,3) family."""
from fractions import Fraction
from itertools import combinations

shifts=((0,0),(10,0),(0,10),(1,0),(0,1),(9,0),(9,1),(0,9),(1,9))
simple=[(a,a+10,b,b+10,a+b,a+b+10) for a,b in shifts]
# Opposite triangle x<=9, y<=10, x+y>=10, with redundant tight bounds.
simple.append((0,9,1,10,10,19))
wide=[tuple(40*x+(-1 if j%2==0 else 1)*(1+i%3) for j,x in enumerate(row))
      for i,row in enumerate(simple)]
n=len(wide)
N=(1<<n)-1

def constraints(row):
    return ((-1,0,-row[0]),(1,0,row[1]),(0,-1,-row[2]),(0,1,row[3]),(-1,-1,-row[4]),(1,1,row[5]))

def vertex(h,k):
    a,b,c=h
    d,e,f=k
    det=a*e-b*d
    if det:
        return Fraction(c*e-b*f,det),Fraction(a*f-c*d,det)

def accepts(row,p):
    x,y=p
    return row[0]<=x<=row[1] and row[2]<=y<=row[3] and row[4]<=x+y<=row[5]

def all_masks(rows):
    # Arrangement vertices suffice: each nonempty intersection of bounded
    # convex polygons has an extreme point on two independent boundary lines.
    lines=tuple(set(h for row in rows for h in constraints(row)))
    V={p for h,k in combinations(lines,2) if (p:=vertex(h,k)) is not None}
    return {sum(1<<i for i,row in enumerate(rows) if accepts(row,p)) for p in V}-{0}

def faces(max_masks):
    F=set()
    for s in max_masks:
        t=s
        while t:
            F.add(t)
            t=(t-1)&s
    return F

base,expanded=all_masks(simple),all_masks(wide)
F=faces(base)
assert F==faces(expanded)
assert [q for q in combinations(range(n),2) if sum(1<<i for i in q) not in F]==[(1,9)]
assert all(any(sum(1<<i for i in t) in F for t in combinations(q,3))
           for q in combinations(range(n),4))
assert all(a|b!=N for a in expanded for b in expanded)
assert any(a|b|c==N for a,b,c in combinations(expanded,3))
assert all(len({p for h,k in combinations(constraints(row),2)
                 if (p:=vertex(h,k)) is not None and accepts(row,p)})==6 for row in wide)
print('ten true six-edged polygons; exactly one disjoint pair; all 210 quadruples good; piercing number exactly 3')
print('nerve faces:',len(F))
for i,row in enumerate(wide):
    print(i,*row)
