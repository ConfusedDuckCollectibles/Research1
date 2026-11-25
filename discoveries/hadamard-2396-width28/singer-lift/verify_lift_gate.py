#!/usr/bin/env python3
"""Exact independent checks of two orthogonal quotient gates. Standard library."""
from itertools import product

def iterate_group(moduli):
    return tuple(product(*(range(m) for m in moduli)))

def ring_add(x,y,mods):
    return tuple((a+b)%n for a,b,n in zip(x,y,mods))

def paff(blocks, mods):
    pts = iterate_group(mods)
    ids = {x:i for i,x in enumerate(pts)}
    shf = [[ids[ring_add(x,s,mods)] for x in pts] for s in pts]
    return tuple(sum(row[i]*row[j] for row in blocks for i,j in enumerate(t)) for t in shf)

def target(v,q,zero,in_N):
    return 4*v*int(zero)-4*q*int(in_N)+4

def project_small_control():
    # Width 12: G=(C2)^2, N={0}, q=4. The first three fibres are
    # precisely (+,+,+) outside zero, and the fourth is all positive.
    rows = ((-1,1,1,1),)*3 + ((1,1,1,1),)
    R = paff(rows,(2,2))
    assert R == tuple(target(7,4,s==0,s==0) for s in range(4))
    plus = [sum(a)//2 for a in rows]
    char = (1,1,-1,-1)
    minus = [sum(x*y for x,y in zip(a,char))//2 for a in rows]
    assert plus == [1,1,1,2] and minus == [-1,-1,-1,0]
    assert sum(x*x for x in plus) == 7 and sum(x*x for x in minus) == 3
    print("positive width-12 control: direct PAF [16,4,4,4]; B4 plus/minus energies (7,3)")

def order7_profiles():
    # Enumerate all scalar sign profiles on F8 and the two possible
    # orders of B∩N, rather than assuming the bound from the note.
    norms = sorted({sum(row)**2 for row in product((-1,1),repeat=8)
                    if len(set(row[1:]))==1})
    assert norms == [36,64]
    projected_norms=[]
    for e_in_B in (False,True):
        projected = sum(target(599,8,x==0,x==0 or (x==7 and e_in_B)) for x in range(8))
        projected_norms.append(projected)
        assert 74*norms[0]>projected
    assert projected_norms == [2396,2364]
    print("order 7: one block lower 2664, full four-block upper 2396 or 2364")

def b4_native_difference():
    # G=B×C with B=(C2)^2 and C=(C2)^2×C37; N occupies 74
    # shifts of C including zero. Direct character-weighted native R.
    C = iterate_group((2,2,37))
    chi = (1,1,-1,-1)
    X=[]; Y=[]
    for e,x,y in C:
        vals=[target(599,8,b==e==x==y==0,b==x==0) for b in range(4)]
        assert sum(vals)%4 == sum(v*c for v,c in zip(vals,chi))%4 ==0
        X.append(sum(vals)//4)
        Y.append(sum(v*c for v,c in zip(vals,chi))//4)
    assert all(x-y==4 for x,y in zip(X,Y))
    assert (X[0],Y[0],sum(X),sum(Y)) == (595,591,599,7)
    btable = [(a+3*b)//2 for a,b in product((-1,1),repeat=2)]
    ctable = [(a-b)//2 for a,b in product((-1,1),repeat=2)]
    assert set(zip(btable,ctable)) == {(1,-1),(-1,1),(2,0),(-2,0)}
    admissible_neighbors={ (l,r) for l,r in product((-1,1),repeat=2) if 2*(l+r)==4}
    assert admissible_neighbors=={(1,1)}
    assert (1+len(C))**2==22201>sum(X)
    print("order 3 on transverse B4: all 148 differences +4; one defect implies 22201>599")

def c37_compressed_control():
    # Recheck the previous positive margin tuple in a different quotient;
    # importantly, it is NOT a positive original core.
    records=((-13,11,1,-1),(-13,13,-1,-1),(-13,11,1,-5),(-7,-1,1,-5))
    orbit_sums={s+3*(2*k-12) for s in (-1,1) for k in range(13)}
    assert set(sum((list(row) for row in records),[]))<=orbit_sums
    # A=C2×F2^3. Last 3 coordinates zero is the exceptional position.
    blocks=[[ ((z0 if e==0 else z1) if x==y==z==0 else (w0 if e==0 else w1))
             for e,x,y,z in iterate_group((2,2,2,2))]
            for z0,z1,w0,w1 in records]
    R=tuple(sum(target(599,8,e==u==y==0,u==0) for y in range(37))
            for e in range(2) for u in range(8))
    obs = paff(blocks,(2,2,2,2))
    assert obs==R and obs==tuple([1360]+[148]*7+[-1036]+[148]*7)
    print("old C37 compression remains positive (16/16), but every Singer lift fails the transverse norm")

if __name__=="__main__":
    project_small_control()
    order7_profiles()
    b4_native_difference()
    c37_compressed_control()
