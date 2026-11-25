#!/usr/bin/env python3
"""Exact small controls for the +4-off-subgroup width-28 template."""
from collections import Counter
from itertools import combinations_with_replacement, product
from math import isqrt

def orbits(h):
    u=next(x for x in range(1,37)
           if pow(x,h,37)==1 and all(pow(x,k,37)!=1 for k in range(1,h)))
    O={frozenset(x*pow(u,j,37)%37 for j in range(h)) for x in range(37)}
    assert set().union(*O)==set(range(37))
    assert sum(map(len,O))==37
    assert sorted(map(len,O))==sorted([1]+[h]*(36//h))
    return O

def orbit_sums(h):
    possible={0}
    for O in orbits(h):
        possible={x+y*len(O) for x in possible for y in (-1,1)}
    return possible

def check_small_row_norms():
    for h,P in [(4,{(3,-3):29,(5,-5):1,(13,-13):1,(21,-3):1}),
                (6,{(1,-1):20,(1,1):9,(13,-13):1,(13,-1):1,(25,-13):1})]:
        S=orbit_sums(h)
        assert sum(P.values())==32 and all(a in S and b in S for a,b in P)
        assert sum((a*a+b*b)*n for (a,b),n in P.items())==1360
        assert sum((a+b)**2*n for (a,b),n in P.items())==324
    print("Orders 4 and 6 both pass the prerequisite two pair norms (positive controls).")
    for h in (9,12,18,36):
        # Transform each pair toward nonnegative residual budgets.
        S=orbit_sums(h)
        scores={((a*a+b*b-2)//8,((a+b)//2)**2)
                for a,b in product(S,repeat=2)}
        remaining={(162,81)}
        for _ in range(32):
            remaining={(r-x,s-y) for r,s in remaining for x,y in scores if r>=x and s>=y}
        assert (0,0) not in remaining
    print("Orders 9, 12, 18, 36 fail the two pair norms (independent residual-budget check).")

def four_square_reps(n):
    to_root={x*x:x for x in range(isqrt(n)+1)}
    R=[]
    for prefix in combinations_with_replacement(range(isqrt(n)+1),3):
        rest=to_root.get(n-sum(x*x for x in prefix),-1)
        if rest>=prefix[-1]: R.append(prefix+(rest,))
    return R

def check_quotient_norms():
    R=four_square_reps(599)
    assert len(R)==14
    U=Counter(sum(max(s,4*s-24) for s in row) for row in R)
    assert U=={82:1,84:1,86:1,87:1,90:1,92:1,93:3,94:2,95:1,98:1,99:1}
    assert [r for r in R if sum(max(s,4*s-24) for s in r)==82]==[(3,5,6,23)]
    print("Order 4: 14 sorted four-square reps of 599, minimum total lower bound 82 > 81.")
    # h=6 optional quantitative row check. No light-heavy mixing at deficit <=49.
    E=sorted({(a+b)//2 for a,b in product(orbit_sums(6),repeat=2) if abs((a+b)//2)<=9})
    assert E==[-7,-6,-5,-1,0,1,5,6,7]
    possible=set()
    for L in combinations_with_replacement(E,8):
        t=sum(x*x for x in L); s=sum(L)
        if t<=81 and 8*t-s*s<=49:
            assert all(-1<=x<=1 for x in L)
            possible.add((t,s*s))
    assert possible and max(t for t,_ in possible)==8
    print("Order 6: every admissible single quotient block has t<=8, four give <81.")

def check_positives():
    # Generated sign rows on C8, N={0}, width=28, target order 60.
    X=[{0},{3,7},{5,6,7},{2,5,7}]
    R=[tuple(1-2*int(j in S) for j in range(8)) for S in X]
    assert [sum(a[x]*a[(x+s)%8] for a in R for x in range(8)) for s in range(8)]==[32]+[4]*7
    assert sum(sum(a)**2 for a in R)==60
    print("Synthetic native core: C8, target order 60, R(0)=32 and seven off-subgroup R=4.")
    # All compressed correlations for one common Singer order-7 row symmetry.
    # Coordinate 0/8 is F8=0 in the corresponding C2 slice.
    T=((-13,11,1,-1),(-13,13,-1,-1),(-13,11,1,-5),(-7,-1,1,-5))
    V=[(a,)+(c,)*7+(b,)+(d,)*7 for a,b,c,d in T]
    raw=lambda a,y: 2368 if (a,y)==(0,0) else (-28 if a in (0,8) else 4)
    target=[sum(raw(d,y) for y in range(37)) for d in range(16)]
    assert target==[1360]+[148]*7+[-1036]+[148]*7
    assert [sum(a[x]*a[x^d] for a in V for x in range(16)) for d in range(16)]==target
    for h in (1,2,3):
        S=orbit_sums(h)
        assert all(x in S for r in V for x in r)
    print("Orders 1, 2, 3 pass all 16 shifts of a Singer-symmetric compressed witness (no full lift).")

if __name__=="__main__":
    check_small_row_norms()
    check_quotient_norms()
    check_positives()
