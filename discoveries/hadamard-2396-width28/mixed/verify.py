#!/usr/bin/env python3
"""Exact quotient/Walsh, lifting, and affine checks for width-28 order 2396."""
from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement, product
from math import factorial, isqrt

F8=range(8)

def walsh(x,k):
    return (-1)**((x&k).bit_count()%2)

def coeff(row,k):
    return sum(v*walsh(x,k) for x,v in enumerate(row))

def inv_rows(trivial,amplitude):
    rows=set()
    for signs in product((-amplitude,amplitude),repeat=7):
        spectrum=(trivial,)+signs
        n=tuple(sum(c*walsh(x,k) for k,c in enumerate(spectrum)) for x in F8)
        if any(t%8 for t in n): continue
        row=tuple(t//8 for t in n)
        assert tuple(coeff(row,k) for k in F8)==spectrum
        rows.add(row)
    return rows

def support_rows():
    planes={frozenset((x,y,z,x^y^z)) for x,y,z in combinations(F8,3)}
    assert len(planes)==14 and all(len(P)==4 for P in planes)
    D=lambda n,j:tuple(n*int(x==j) for x in F8)
    odd={D(1,j) for j in F8}
    even={D(2,j) for j in F8} | {
        tuple(-1 if x==j else int(x in P) for x in F8) for P in planes for j in P}
    assert inv_rows(1,1)==odd and inv_rows(2,2)==even
    assert len(odd)==8 and len(even)==64
    print("Independent inverse Walsh: 8 odd-normalized rows and 64 even-normalized rows.")
    return odd,even,{P for P in planes if 0 in P}

def reps(N):
    # Independent two-pair join for sorted four-square representations.
    lists=defaultdict(list)
    for a,b in combinations_with_replacement(range(isqrt(N)+1),2):
        lists[a*a+b*b].append((a,b))
    return sorted((a,b,c,d) for n,P in lists.items() for a,b in P
                  for c,d in lists.get(N-n,()) if b<=c)

def restore(rows,s,t):
    # t=1 or 2 is normalized trivial coefficient, s>0 the actual one.
    sign=1 if (s-t)%8==0 else -1
    assert (s-sign*t)%8==0
    off=(s-sign*t)//8
    ans={tuple(off+sign*x for x in row) for row in rows}
    assert len(ans)==len(rows) and all(sum(row)==s for row in ans)
    return ans

def real_profiles(odd,even):
    R=reps(599)
    assert len(R)==14
    keep=[]
    nclasses=0
    for t in R:
        E=[x for x in t if x%8 in (2,6)]
        O=[x for x in t if x%8 in (1,7)]
        if len(E)==1 and len(O)==3:
            keep.append(t)
            pool=[restore(even if s%2==0 else odd,s,2 if s%2==0 else 1) for s in t]
            class_factor=1
            for P,s in zip(pool,t):
                square=1 if s%2 else 4
                for row in P:
                    assert max(map(abs,row))<=7
                    assert [coeff(row,k)**2 for k in F8]==[s*s]+[square]*7
                orbits={frozenset(tuple(row[x^j] for x in F8) for j in F8) for row in P}
                assert len(orbits)==(1 if s%2 else 8)
                class_factor*=len(orbits)
            nclasses+=class_factor
    assert keep==[(1,7,15,18),(2,9,15,17),(6,7,15,17),(7,10,15,15)]
    L=0
    for t in keep:
        places=factorial(4)
        for n in Counter(t).values(): places//=factorial(n)
        L+=places*2**4*len(odd)**3*len(even)
    assert L==44_040_192 and nclasses==32 and 4*2**4*len(odd)**3*len(even)==2_097_152
    print("Four physical profiles:",*keep)
    print("44,040,192 labeled signed quotient arrays; 32 translation/sign/permutation classes.")
    return keep

def add(mods,x,y):
    return tuple((a+b)%m for a,b,m in zip(x,y,mods))

def linear_maps(mods):
    V=list(product(*(range(n) for n in mods)))
    output=set()
    # Choose images of cyclic generators, evaluate and test surjectivity.
    for B in product(*[[t for t in V if all(n*a%m==0 for a,m in zip(t,mods))]
                       for n in mods]):
        p=tuple(tuple(sum(x[i]*B[i][j] for i in range(len(mods)))%m
                      for j,m in enumerate(mods)) for x in V)
        if len(set(p))==len(V): output.add(p)
    return V,output

def affine_types(odd,profiles):
    for mods,expect in [((8,), (4,32,1)), ((4,2),(8,64,1)), ((2,2,2),(168,1344,225))]:
        V,A=linear_maps(mods)
        idx={v:i for i,v in enumerate(V)}
        H={tuple(idx[add(mods,y,t)] for y in a) for a in A for t in V}
        H3=[p for p in H if all(p[p[p[i]]]==i for i in F8)]
        assert (len(A),len(H),len(H3))==expect
        print("Q",mods,"aut, affine, affine order-dividing-3:",expect)
        if mods==(2,2,2):
            assert Counter(sum(x==p[x] for x in F8) for p in H3)=={8:1,2:224}
            for p in H3:
                fixed=[x for x in F8 if x==p[x]]
                for profile in profiles:
                    assert any(not any(all(r[x]==r[p[x]] for x in F8)
                                       and all(abs(r[x])<=1 for x in fixed)
                                       for r in restore(odd,s,1)) for s in profile if s%2)

def lift_planes(planes):
    V,Aut=linear_maps((2,2,2))
    C=range(3)
    idx={v:i for i,v in enumerate(V)}
    T=[tuple(idx[z] for z in a) for a in Aut]
    H=frozenset(i for i,x in enumerate(V) if x[0]==0)
    U=[p for p in T if frozenset(p[i] for i in H)==H]
    assert len(T)==168 and len(U)==24
    for p in U:
        B=[V[p[idx[tuple(int(j==i) for j in C)]]] for i in C]
        # Lift to A=C4 x C2 x C2. First cyc generator has order 4; others order 2.
        lift=lambda x:tuple(sum(x[i]*B[i][j] for i in C)%m for j,m in enumerate((4,2,2)))
        A=list(product(range(4),range(2),range(2)))
        assert len({lift(x) for x in A})==16 and lift((2,0,0))==(2,0,0)
        for x in A:
            y=lift(x)
            assert idx[(y[0]%2,y[1],y[2])]==p[idx[(x[0]%2,x[1],x[2])]]
    def plane_orbits(T):
        return sorted({frozenset(frozenset(p[i] for i in P) for p in T) for P in planes},key=len)
    assert list(map(len,plane_orbits(T)))==[7]
    assert list(map(len,plane_orbits(U)))==[1,6]
    assert plane_orbits(U)[0]==frozenset((H,))
    print("Liftable common automorphisms: 8 classes on C2^4 and 12 on C4 x C2^2 (four physical profiles).")

def field_order(a):
    return next(k for k in range(1,37) if pow(a,k,37)==1)

def field_control():
    seen=0
    for s in range(1,37):
        h=field_order(s)
        if h not in (9,18,36): continue
        k=h//9
        for t in range(37):
            F=lambda y:(s*y+t)%37
            I=[]
            for y in range(37):
                z=y
                for _ in range(k): z=F(z)
                I.append(z)
            unseen=set(range(37))
            lengths=[]
            while unseen:
                x=min(unseen); y=x; l=0
                while y in unseen: unseen.remove(y); y=I[y]; l+=1
                assert y==x
                lengths.append(l)
            assert sorted(lengths)==[1,9,9,9,9]
            seen+=1
    assert seen==888
    S={c+18*x for c in (-1,1) for x in range(-2,3)}
    assert sorted({(x+y)//2 for x in S for y in S if abs((x+y)//2)<=7})==[-1,0,1]
    print("Scalar power control: all 888 field-affine maps with scalar order 9, 18 or 36 yield 1+4×9 orbits after removal of 2-part.")

def singer_control():
    P=((-13,11,1,-1),(-13,13,-1,-1),(-13,11,1,-5),(-7,-1,1,-5))
    V=[(a,)+(c,)*7+(b,)+(d,)*7 for a,b,c,d in P]
    assert [sum(r[x]*r[x^d] for r in V for x in range(16)) for d in range(16)]==[1360]+[148]*7+[-1036]+[148]*7
    Q=[tuple((r[x]+r[x+8])//2 for x in F8) for r in V]
    assert [sum(coeff(row,k)**2 for row in Q) for k in F8]==[599]+[7]*7
    assert sorted(abs(sum(r)) for r in Q)==[1,7,15,18]
    print("Prior Singer compressed control passes all +148 off-pair shifts and physical profile (1,7,15,18).")

if __name__=="__main__":
    odd,even,planes=support_rows()
    profiles=real_profiles(odd,even)
    affine_types(odd,profiles)
    lift_planes(planes)
    field_control()
    singer_control()
