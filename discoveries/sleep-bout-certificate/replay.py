#!/usr/bin/env python3
"""Private exact reconstruction: no empirical inputs. Python standard library only."""
from collections import Counter
from fractions import Fraction as F
from itertools import accumulate, groupby, permutations, product
from random import Random

Z, C, B, A = range(4)                 # SS, SW, WS, WW
N = 960

def tracks(j):
    return tuple(x >> 1 for x in j), tuple(x & 1 for x in j)

def words(j, K):
    return {k: Counter(zip(*(j[i:] for i in range(k)))) for k in range(1, K+1)}

def dwells(x):
    return Counter((v, sum(1 for _ in g)) for v,g in groupby(x))

def wake_dwells(x):
    return {k:v for k,v in dwells(x).items() if k[0]}

def bouts(x, minimum=10):
    cut = [0]+[i for i in range(1,len(x)) if x[i]!=x[i-1]]+[len(x)]
    return [(a,b) for a,b in zip(cut,cut[1:]) if x[a] and b-a >= minimum]

def hit(r,d,tau=F(1,2)):
    a,b=r; c,e=d
    o=max(0,min(b,e)-max(a,c))
    return o>0 and F(o,b-a+e-c-o)>=tau

def results(j, minimum=10, tau=F(1,2)):
    r,d=tracks(j)
    P,Q=bouts(r,minimum),bouts(d,minimum)
    H=[(i,k) for i,p in enumerate(P) for k,q in enumerate(Q) if hit(p,q,tau)]
    # All test thresholds and records used here give unique matched edges.
    assert len(H)==len({i for i,k in H})==len({k for i,k in H})
    return len(H),len(P),len(Q),sum(r),sum(d),Counter(j)

def show(label,j):
    h,nr,nd,wr,wd,t=results(j)
    print(f'{label}: hits={h}, bouts={nr},{nd}, awake_min={F(wr,2)},{F(wd,2)}, '
          f'confusion[SS,SW,WS,WW]={tuple(t[i] for i in range(4))}, '
          f'accuracy={F(t[Z]+t[A],len(j))}, wake_sensitivity={F(t[A],wr)}')
    return h

def unaligned_shift(x,L,gap=90):
    j=[Z]*gap
    for v in x:
        j += [B]*(L-v)+[A]*v+[C]*(L-v)+[Z]*gap
    assert len(j)<=N
    return tuple(j+[Z]*(N-len(j)))

def fixed_shift(x,L=10):
    j=[Z]*N
    for i,v in enumerate(x):
        a=72+100*i; b=a+L-v
        for k in range(a,a+L): j[k]|=B
        for k in range(b,b+L): j[k]|=C
    return tuple(j)

def formula(n,X,a,b,q):
    return {h for h in range(n+1) if h*q+(n-h)*a <= X <= h*b+(n-h)*(q-1)}

def split_pairs(m,T,d=1):
    c=lambda i: max(0,min(i+1,2*m-1-i))
    U=[(T+d*i,T+d*(i+2)) for i in range(2*m-1) for _ in range(c(i))]
    V=[(T+d*k,T+d*k) for k in range(1,2*m) for _ in range(c(k-1)-(k==m))]
    V.append((T,T+2*m*d))
    key=lambda p: (sum(p),p[0])
    U.sort(key=key); V.sort(key=key)
    assert len(U)==len(V)==m*m
    assert list(map(sum,U))==list(map(sum,V))
    assert Counter(x for p in U for x in p)==Counter(x for p in V for x in p)
    return U,V

def split_record(pairs,K,mirror=False,n=None,g=1):
    gap=max(20,K+2,g+1)
    j=[Z]*gap
    for x in ([B,C] if mirror else [B]):
        for p,q in pairs:
            j += [A]*p+[x]*g+[A]*q+[Z]*gap
    if n is not None:
        assert len(j)<=n
        j += [Z]*(n-len(j))
    return tuple(j)

def from_word_counts(w,U,minimum):
    span=4*U//3
    assert max(w)>=span+2
    total=0
    for width in range(3,span+3):
        for j,n in w[width].items():
            r,d=tracks(j)
            total += n*sum(hit(p,q) and min(p[0],q[0])==1 and max(p[1],q[1])==width-1
                           for p in bouts(r,minimum) for q in bouts(d,minimum))
    return total

def main():
    f=[unaligned_shift(x,16) for x in ((10,)*6,(11,)*5+(5,))]
    assert words(f[0],5)==words(f[1],5)
    assert list(map(wake_dwells,tracks(f[0])))==list(map(wake_dwells,tracks(f[1])))
    h={sum(x>=11 for x in v) for v in product(range(4,13),repeat=6) if sum(v)==60}
    assert h==formula(6,60,4,12,11)==set(range(6))
    assert [show('Frozen 1 '+s,j) for s,j in zip(('lo','hi'),f)]==[0,5]
    print('Frozen 1 exact support:',sorted(h))
    for n in range(1,6):
        for q in range(2,6):
            brute={}
            for v in product(range(1,6),repeat=n):
                brute.setdefault(sum(v),set()).add(sum(x>=q for x in v))
            for X,e in brute.items(): assert formula(n,X,1,5,q)==e

    x0=(3,6,5,6,6,6,6,6); x1=(3,3,4,7,7,7,7,6)
    low,high=map(fixed_shift,(x0,x1))
    assert words(low,4)==words(high,4)
    assert tracks(low)[0]==tracks(high)[0]
    assert list(map(dwells,tracks(low)))==list(map(dwells,tracks(high)))
    assert dwells(tracks(low)[1])==Counter({(0,79):1,(0,87):1,(0,89):1,(0,90):4,(0,91):1,(0,174):1,(1,10):8})
    sig=Counter((-1,0,0,0,0,1,3))
    full=[]
    for mid in product(range(3,8),repeat=6):
        xs=(3,)+mid+(6,)
        if sum(xs)==44 and Counter(y-x for x,y in zip(xs,xs[1:]))==sig: full.append(xs)
    walk=[xs for v in set(permutations(sig.elements()))
          if ((xs:=tuple(accumulate((3,)+v)))[-1]==6 and sum(xs)==44 and all(3<=x<=7 for x in xs))]
    assert set(walk)==set(full) and len(full)==11
    possible=set()
    for x in full:
        j=fixed_shift(x)
        assert words(j,4)==words(low,4)
        assert tracks(j)[0]==tracks(low)[0]
        assert list(map(dwells,tracks(j)))==list(map(dwells,tracks(low)))
        possible.add(results(j)[0])
    assert possible==set(range(5))
    assert [show('Frozen 2 '+s,j) for s,j in zip(('lo','hi'),(low,high))]==[0,4]
    print('Frozen 2 exact support:',sorted(possible),'; 11 feasible trajectories')

    for m in range(2,10):
        for K in (1,4,11,15):
            for g in (1,3):
                U,V=split_pairs(m,max(10,K-1),(g+1)//2)
                for flag in (False,True):
                    j,z=(split_record(p,K,flag,g=g) for p in (U,V))
                    assert words(j,K)==words(z,K)
                    assert list(map(dwells,tracks(j)))==list(map(dwells,tracks(z)))
                    assert flag or tracks(j)[0]==tracks(z)[0]
                    v,w=results(j),results(z)
                    assert (v[0],w[0])==(m*m*(1+flag),1+flag)
                    assert not flag or (v[1]==v[2]==3*m*m and v[3]==v[4])
    for flag,g in ((False,1),(True,1),(True,3)):
        U,V=split_pairs(2,10,(g+1)//2)
        l,r=(split_record(p,11,flag,960,g) for p in (V,U))
        assert [show(f'Split trade mirror={flag} gap={g} '+s,j)
               for s,j in zip(('lo','hi'),(l,r))]==[1+flag,4*(1+flag)]
        if flag and g==1:
            for j,want in ((l,2),(r,8)):
                assert all(results(j,tau=t)[0]==want for t in (F(241,500),F(1,2),F(14,27)))
    rng=Random(214)
    for _ in range(500):
        n=rng.randrange(8,43)
        j=(0,)+tuple(rng.randrange(4) for _ in range(n-2))+(0,)
        U=max((b-a for x in tracks(j) for a,b in bouts(x,1)),default=1)
        w=words(j,4*U//3+2)
        for minimum in (1,2,4):
            assert from_word_counts(w,U,minimum)==results(j,minimum)[0]
    r,p=tracks((Z,)+(A,)*10+(B,)+(A,)*10+(Z,))
    assert sum(hit(x,y,F(10,21)) for x in bouts(r) for y in bouts(p))==2
    assert sum(hit(x,y) for x in bouts(r) for y in bouts(p))==0
    print('Checked finite support by two independent strategies, scalable exact trades, '
          '1,500 sufficient order conditions, and collision below IoU 1/2.')

if __name__=='__main__': main()
