#!/usr/bin/env python3
"""Exact Gaussian-rational path propagation vs joint-character theorem."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import comb, prod

ZERO=(F(0),F(0)); ONE=(F(1),F(0))
PHASE=((1,0),(0,-1),(-1,0),(0,1)) # (-i)^k


def add(a,b):
    return a[0]+b[0],a[1]+b[1]


def mul(a,b):
    return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]


def abs2(a):
    return a[0]**2+a[1]**2


def bit(w,i):
    return (w>>i)&1


def syndrome(w,edges):
    return sum((bit(w,a)^bit(w,b))<<i for i,(a,b) in enumerate(edges))


def row(mask,edges):
    r=0
    for i,(a,b) in enumerate(edges):
        if bit(mask,i):
            r^=(1<<a)|(1<<b)
    return r


def amplitude_weights(slot):
    n=len(slot)
    return [prod(slot[j][bit(e,j)] for j in range(n)) for e in range(1<<n)]


def propagate(slot1,slot2,echo,mid,full,logical=(ONE,ZERO)):
    """Project separately for each measured syndrome, then square each basis entry."""
    n=len(slot1); inv=(1<<n)-1
    space=range(1<<n)
    w1,w2=amplitude_weights(slot1),amplitude_weights(slot2)
    mid_labels=[syndrome(e,mid) for e in space]
    final_labels=[syndrome(e,full) for e in space]
    w2=[w*(-1)**((e&echo).bit_count()%2) for e,w in enumerate(w2)]
    paths=defaultdict(lambda:ZERO)
    for e in space:
        for f in space:
            if w1[e] and w2[f]:
                net=e^f
                m,y=mid_labels[e],final_labels[net]
                phase=PHASE[(e.bit_count()+f.bit_count())%4]
                g=w1[e]*w2[f]
                paths[m,y,net]=add(paths[m,y,net],(phase[0]*g,phase[1]*g))
    out=defaultdict(lambda:ZERO)
    for (m,y,u),amp in paths.items():
        if amp!=ZERO:
            out[m,y,u]=add(out[m,y,u],mul(logical[0],amp))
            if logical[1]!=ZERO:
                out[m,y,u^inv]=add(out[m,y,u^inv],mul(logical[1],amp))
    law=defaultdict(F)
    for (m,y,u),amplitude in out.items():
        law[m,y]+=abs2(amplitude)
    assert sum(law.values())==1
    return dict(law), paths


def classical(slot1,slot2,mid,full):
    w1,w2=amplitude_weights(slot1),amplitude_weights(slot2)
    law=defaultdict(F)
    for e,u in enumerate(w1):
        for f,v in enumerate(w2):
            law[syndrome(e,mid),syndrome(e^f,full)]+=u*u*v*v
    return dict(law)


def moment(law,b,a):
    return sum(p*(-1)**(((m&b).bit_count()+(y&a).bit_count())%2)
               for (m,y),p in law.items())


def theorem_moment(slot1,slot2,echo,B,A,mid):
    C=A^B
    n=len(slot1)
    c1=[u*u-v*v for u,v in slot1]; s1=[2*u*v for u,v in slot1]
    c2=[u*u-v*v for u,v in slot2]; s2=[2*u*v for u,v in slot2]
    total=F(0)
    for D in range(1<<n):
        if D & ~(A&~B) or syndrome(D,mid):
            continue
        total+=((-1)**D.bit_count()
                *prod(((-1)**bit(echo,j))*s1[j]*s2[j] for j in range(n) if bit(D,j))
                *prod(c1[j] for j in range(n) if bit(C&~D,j))
                *prod(c2[j] for j in range(n) if bit(A&~D,j)))
    return total


def check(law,predict,mid,full):
    for b,a in product(range(1<<len(mid)),range(1<<len(full))):
        assert moment(law,b,a)==predict(row(b,mid),row(a,full)),(a,b)


def same(a,b):
    return all(a.get(k,0)==b.get(k,0) for k in a.keys()|b.keys())


def average(laws):
    out=defaultdict(F)
    for law in laws:
        for event,v in law.items():
            out[event]+=v/len(laws)
    return dict(out)


def chain(n):
    return [(i,i+1) for i in range(n-1)]


def pair(m,sign=1):
    return F(m*m-1,m*m+1),sign*F(2*m,m*m+1)


complex_logical=((F(3,5),F(0)),(F(12,25),F(16,25)))
assert sum(map(abs2,complex_logical))==1

# Independent heterogeneous slot inputs, signs, three midpoint partitions,
# two nonidentical choices of arbitrary logical input, complete joint moments.
n=5; full=chain(n)
slot1=[pair(m,(-1)**m) for m in (4,5,6,7,8)]
slot2=[pair(m,(-1)**m) for m in (9,4,6,11,8)]
blocks=[(full,[0]),
        ([full[i] for i in (0,2,3)],[0,2]),
        ([full[i] for i in (0,3)],[0,2,3])]
for mid,reps in blocks:
    pauli=classical(slot1,slot2,mid,full)
    laws=[]
    for bits in product((0,1),repeat=len(reps)):
        echo=sum(1<<j for j,v in zip(reps,bits) if v)
        law,_=propagate(slot1,slot2,echo,mid,full,complex_logical)
        direct_zero,_=propagate(slot1,slot2,echo,mid,full)
        assert same(law,direct_zero)
        check(law,lambda B,A:theorem_moment(slot1,slot2,echo,B,A,mid),mid,full)
        laws.append(law)
    assert same(average(laws),pauli)
    if mid==full:
        varied,_=propagate(slot1,slot2,0b10110,mid,full,complex_logical)
        assert same(varied,pauli)
    else:
        assert not same(laws[0],laws[-1])
    print('all joint moments & representative echo twirl pass;',len(mid),'midpoint checks')

# Only one physical qubit rotates: merely withholding its outcome does
# nothing (average of full projector readout); both incident checks must go.
n=5; j=2; full=chain(n)
case=[(F(1),F(0)) for _ in range(n)]
case[j]=pair(40)
e_plus=0; e_minus=1<<j
for removed,positive in (((),False),((1,),False),((2,),False),((1,2),True)):
    mid=[edge for i,edge in enumerate(full) if i not in removed]
    plus,_=propagate(case,case,e_plus,mid,full)
    minus,_=propagate(case,case,e_minus,mid,full)
    p=case[j][1]**2
    gap=(moment(minus,0,1<<1)-moment(plus,0,1<<1))/2
    assert gap==(4*p*(1-p) if positive else 0)
    if not positive:
        assert same(plus,minus)
    print('single-site omitted edges',removed,'exact contrast',gap)

# Rare whole component click at r=1,2,3: unmatched complement by odd n,
# paired path cancellation for exact D (not only leading order).
def hit(law,y):
    return sum(q for (m,f),q in law.items() if f==y)


for n,r in ((5,1),(5,2),(7,3)):
    all_edges=chain(n)
    D=(1<<2) if r==1 else (1<<r)-1
    removed=(1,2) if r==1 else (r-1,)
    mid=[edge for i,edge in enumerate(all_edges) if i not in removed]
    slots=[pair(40)]*n
    laws=[]
    for echo in (0,D & -D):
        law,raw=propagate(slots,slots,echo,mid,all_edges)
        f=syndrome(D,all_edges)
        if echo:
            assert all(v==ZERO for (m,y,u),v in raw.items() if (y,u)==(f,D))
        laws.append(hit(law,f))
    p=slots[0][1]**2
    # Count the exact leading split coefficients before multiplying
    # the common product of rotation amplitudes; minimum degree r.
    groups=defaultdict(int)
    for s in range(1<<n):
        if not s & ~D:
            groups[syndrome(s,mid)] += 1
    assert set(groups.values())=={2}
    assert sum(c*c for c in groups.values())==2**(r+1)
    assert laws[0]>laws[1]
    print(f'whole-component r={r}, no echo / echo hit probabilities:',
          f'{float(laws[0]):.10g}',f'{float(laws[1]):.10g}')
    if r==1:
        # Six all-no-echo hits by a fixed shot budget. Exact power in
        # the specified independent ideal-binomial toy control.
        N,k=2000,6
        h=(laws[0]+laws[1])/2
        label=1/(1+laws[1]/laws[0])
        fail=sum(comb(N,i)*h**i*(1-h)**(N-i) for i in range(k))
        power=(1-fail)*label**k
        assert power>F(9,10)
        print(f'ideal example: {N} shots, exact null size <=1/64, power >0.90')

# Exact local readout effect at declared synthetic p=.01, q=.02.
p,q=F(1,100),F(2,100)
delta=4*p*(1-p)
assert (1-2*q)*delta==F(2376,62500)
delta2=(1-2*q)*delta*(1-2*p)**2
assert delta2==F(1426194,39062500)
print('reported-bit click gaps: single site',float((1-2*q)*delta),
      'rotating neighbor',float(delta2))
print('All checks passed.')
