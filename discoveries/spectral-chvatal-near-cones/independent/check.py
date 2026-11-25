#!/usr/bin/env python3
"""Independent exact congruence/inertia and adversarial controls."""
from fractions import Fraction as F
from itertools import product

def sub(s,t): return s&t==s

def down(facets):
    out=set()
    for orig in facets:
        s=orig
        while True:
            out.add(s)
            if not s: break
            s=(s-1)&orig
    return out

def inertia(entries):
    # Symmetric rational pivots. If no nonzero diagonal exists, use a
    # hyperbolic pair; Schur-complement pivots preserve exact congruence.
    M=[list(map(F,row)) for row in entries]
    result=[0,0,0]
    while M:
        n=len(M)
        i=next((i for i in range(n) if M[i][i]),None)
        if i is not None:
            q=M[i][i]
            result[0 if q>0 else 1]+=1
            I=[j for j in range(n) if j!=i]
            M=[[M[a][b]-M[a][i]*M[i][b]/q for b in I] for a in I]
            continue
        p=next(((i,j) for i in range(n) for j in range(i+1,n)
                if M[i][j]),None)
        if p is None:
            result[2]+=n
            break
        i,j=p
        q=M[i][j]
        result[0]+=1
        result[1]+=1
        I=[a for a in range(n) if a not in p]
        M=[[M[a][b]-(M[a][i]*M[j][b]+M[a][j]*M[i][b])/q
            for b in I] for a in I]
    return tuple(result)

def mob_column(s,U):
    # The transpose of the inverse zeta takes e_s to these coefficients.
    return {u:(-1)**((s^u).bit_count()) for u in U if sub(u,s)}

def dd(va,vb,C): return sum(x*y*C[a,b] for a,x in va.items() for b,y in vb.items())

def assemble(U,E,C,B,star):
    n=len(U)
    X=[row+[B[A,T] for T in E] for row,A in zip(
       [[C[A,B] for B in U] for A in U],U)]
    X += [[B[A,T] for A in U]+[0]*len(E) for T in E]
    assert len(X)==n+len(E)
    labels=list(U)+[T|star for T in E]
    for a,b in product(range(len(X)),repeat=2):
        assert X[a][b]==X[b][a]
        if labels[a]&labels[b]: assert X[a][b]==0
    return X

def check(W,i,all_inertia=True):
    star=1<<i
    U=sorted(s for s in W if not s&star)
    E=sorted(s^star for s in W if s&star)
    R=sorted(set(U)-set(E))
    B={(A,T):int(not A&T) for A in U for T in E}
    col={s:mob_column(s,U) for s in U}
    if U:
        assert all(sum(v*B[a,T] for a,v in col[s].items())
                   == (-1)**s.bit_count()*int(sub(s,T))
                   for s in U for T in E)
    chain=next(((a,b) for a,b in product(R,repeat=2)
                if sub(a,b) and (a^b).bit_count()==1),None)
    compressed=all(((A^j)|star) in W for A in U for j in
                   (1<<i2 for i2 in range(A.bit_length()) if A&(1<<i2)))
    assert compressed == (chain is None)
    if chain:
        a,b=chain
        w=col[a].copy()
        for s,v in col[b].items(): w[s]=w.get(s,0)+v
        w={s:v for s,v in w.items() if v}
        assert w and set(w)=={s for s in U if sub(s,b) and s&(a^b)}
        assert all(s&t for s in w for t in w)
        assert all(sum(v*B[s,T] for s,v in w.items())==0 for T in E)
        return False
    # Use twice the rational correction, plus doubled ordinary B.
    Rset=set(R)
    C={(A,F):0 if A&F else (
        (-1)**(1+A.bit_count())*int(A in Rset)
        +(-1)**(1+F.bit_count())*int(F in Rset))
       for A,F in product(U,repeat=2)}
    for A,F in product(R,repeat=2):
        assert dd(col[A],col[F],C)==(-2 if A==F else 0)
    if all_inertia:
        twiceB={pair:2*v for pair,v in B.items()}
        M=assemble(U,E,C,twiceB,star)
        assert inertia(M)==(len(E),len(U),0)
    return True

def all4():
    n=4
    ds=[]
    for mask in range(1<<(1<<n)):
        W={a for a in range(1<<n) if mask&(1<<a)}
        if all((a^(1<<j)) in W for a in W for j in range(n) if a&(1<<j)):
            ds.append(W)
    assert len(ds)==168
    yes=0
    for W in ds:
        for k in range(4):
            yes+=check(W,k)
    print(f"all 168 labelled downsets on four ground points, 672 cuts: "
          f"{yes} exact inertia, {672-yes} exact chain stresses")

def examples():
    W={s for s in range(128) if s.bit_count()<=3}
    assert check(W,0)
    fclq_weights=[(0,0,F(1,4),F(1,2))[s.bit_count()] for s in W]
    assert sum(fclq_weights)==F(91,4)>22
    for n2,n3 in product(range(4),repeat=2):
        if n2*2+n3*3<=7:
            assert F(n2,4)+F(n3,2)<=1
    print("EFF 7-point rank-3 skeleton: inertia (22,42,0), "
          "fractional independence 91/4 witness > 22")

    # Edges 12,13,14,15 and 23,34,45,25 induce a C4 off apex.
    W=down([3,5,9,17,6,12,24,18])
    assert len(W)==14 and check(W,0)
    print("nonshiftable C4 off universal apex: inertia (5,9,0)")

    # Interior lower diagonals must not be inadvertently counted as edges.
    # Matching is not a near-cone at i, so ordinary B is impossible.
    W=down([3,12])
    assert not check(W,0)
    U,E=[0,2,4,8,12],[0,2]
    B={(A,T):int((A,T) in {(2,0),(12,0),(0,2),(4,2),(8,2)})
       for A in U for T in E}
    C={(A,S):(-1 if A==S==0 else int({A,S} in ({2,12},{4,8})))
       for A,S in product(U,repeat=2)}
    M=assemble(U,E,C,B,1)
    assert inertia(M)==(2,5,0)
    print("disjoint 2-facets: fixed-B obstruction but alternate B inertia (2,5,0)")

if __name__=="__main__":
    all4()
    examples()
    print("PASS")
