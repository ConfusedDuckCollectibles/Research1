#!/usr/bin/env python3
"""Clean-room algebra and rational determinant/Sturm audit; no optimizer."""
from fractions import Fraction as F
from itertools import permutations,product
from math import isqrt
from pathlib import Path
import json
import sys
import sympy as sp

DIRECTIONS=list(permutations(range(4)))

def parity(seq):
    return (-1)**sum(seq[i]>seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))

SIGNED=[(p,parity(p)) for p in DIRECTIONS]

def tensor(a,b):
    return [[a[i//2][j//2]*b[i%2][j%2] for j in range(4)] for i in range(4)]

def direct(a,b,u,v,d):
    # In the interior conjugate A1 by S_A=diag(1,sin(theta_A)).
    # The off-diagonals become sin(theta_A)^2 and 1. Same for B1.
    # Principal minors extend continuously to zero sin by polynomial identity.
    Z=[[1,0],[0,-1]]
    I=[[1,0],[0,1]]
    A1=[[1-u,u*(2-u)],[1,u-1]]
    B1=[[1-v,v*(2-v)],[1,v-1]]
    action=[tensor(Z,Z),tensor(Z,B1),tensor(A1,Z),tensor(A1,B1),
            tensor(Z,I),tensor(I,Z)]
    weights=[1,1,1,-1,a,b]
    return [[(2+a+b+d)*int(i==j)-sum(c*M[i][j] for c,M in zip(weights,action))
             for j in range(4)] for i in range(4)]

def det(A):
    n=len(A)
    sig=SIGNED if n==4 else [(p,parity(p)) for p in permutations(range(n))]
    work=0
    for p,k in sig:
        v=k
        for i,j in enumerate(p):
            v*=A[i][j]
        work+=v
    return work

def invariants(a,s):
    t=a+s
    b=2-t
    A=4-a*a
    B=4-b*b
    C=8-2*s+a*b
    return t,b,A,B,C

def algebra():
    a,s,u,v,d=sp.symbols('a s u v d')
    t,b,A,B,C=invariants(a,s)
    M=direct(a,b,u,v,d)
    P=d*(d+4+2*a)*(d+4+2*b)*(d+2*(a+b))
    D=P+8*u*v*(a*b*d-s*C+(A-2*v)*u+B*v)
    assert sp.expand(det(M)-D)==0
    assert sp.factor(C-(2+a)*(4-t))==0
    P0=(4+2*a)*(4+2*b)*2*(a+b)
    R0=(2+a)*(4-t)/(2*(2-a)*(2-s))
    assert sp.factor(4*(2+a)*C**2/(A*P0)-R0)==0
    assert sp.expand((R0-1)*2*(2-a)*(2-s)-(6*a+2*s-a*a-3*a*s)).cancel()==0
    assert F(19,25)*2*F(199,100)**2>6
    assert F(1,8)*(F(19,25))+1/(6*F(199,100))<F(18,100)
    assert 4*F(798,100)*F(398,100)>127
    assert F(68,8*127)<F(7,100)
    # Independent Schur constants at the far corner.
    A0,B0,Q=F(199,50),F(397,50),F(2,5)
    dd=A0*B0-Q**2
    assert dd>0
    assert B0/dd<F(3,10) and A0/dd<F(13,100) and Q/dd<F(1,50)
    cross=2*F(1,50)*F(4,5)   # sqrt(0.08*8) = 4/5
    assert F(3,10)*F(2,25)+F(13,100)*8+cross==F(137,125)<F(11,10)
    k,x=sp.symbols('k x')
    f=x*(1-x)**2/(2*(1-k*x))
    assert sp.factor(sp.diff(f,x)-(1-x)*(1-3*x+2*k*x*x)/(2*(1-k*x)**2))==0
    beta=sp.symbols('beta')
    assert sp.expand(8+beta*(2-beta)-(4-beta)*(2+beta))==0
    assert len([v for v in product([-1,1],repeat=4)
                if v[0]*v[2]+v[0]*v[3]+v[1]*v[2]-v[1]*v[3]==2])==8
    print('Independent symbolic 4x4 determinant, exact uniform constants and Schur bounds: PASS')

def trim(p):
    p=p[:]
    while p and not p[-1]:p.pop()
    return p

def neg_remainder(a,b):
    a=trim(a)
    assert b
    while len(a)>=len(b):
        k=len(a)-len(b)
        v=a[-1]/b[-1]
        for i,x in enumerate(b):
            a[k+i]-=v*x
        a=trim(a)
    return [-x for x in a]

def value(p,x):
    o=0
    for c in reversed(p):o=o*x+c
    return o

def sturm(p):
    p=list(map(F,p))  # This includes the integral leading sextic coefficient.
    q=[p,[i*p[i] for i in range(1,len(p))]]
    while q[-1]:
        q.append(neg_remainder(q[-2],q[-1]))
    return q[:-1]

def variations(seq,x):
    vals=[value(p,x) for p in seq]
    signs=[1 if a>0 else -1 for a in vals if a]
    return sum(a!=b for a,b in zip(signs,signs[1:]))

def upper_sign(a,s,d):
    t,b,A,B,C=invariants(a,s)
    assert 0<a<t<=F(1,100) and d>0
    assert (s*C+2*a*a)/B<A/2  # removes any near-opposite-corner gain
    J=s*C-a*b*d
    if J<=0:return
    P=d*(d+4+2*a)*(d+4+2*b)*(d+2*(a+b))
    # When gain is possible complete the square in u; multiply the
    # resulting lower envelope by its strictly positive denominator.
    poly=[P*A,-2*P-2*J*J,4*J*B,-2*B*B]
    interval=(F(0),J/B)
    assert all(value(poly,x)>0 for x in interval)
    seq=sturm(poly)
    assert variations(seq,interval[0])==variations(seq,interval[1])

def lower_sign(a,s,d,u,v):
    assert 0<u<2 and 0<v<2
    b=invariants(a,s)[1]
    assert det(direct(a,b,u,v,d))<0

def gigena_coefficients(a,b,corrected=True):
    # Gigena et al. (43), with Pakhunov §5's sign correction at tau5.
    return [
        a**6*(27*b*b-8)+a**4*(-54*b**4+48*b*b+32)+a*a*(27*b**6+48*b**4-400*b*b+128)
            -8*b**6+32*b**4+128*b*b-512,
        -168*a**3*b**3+60*a**5*b+160*a**3*b+60*a*b**5+160*a*b**3-576*a*b,
        -6*a**4*b*b-6*a*a*b**4-64*a*a*b*b+20*a**4+96*a*a+20*b**4+96*b*b+320,
        8*a**3*b**3-24*a**3*b-24*a*b**3+96*a*b,
        11*a*a*b*b-16*a*a-16*b*b-64,
        (-4 if corrected else 4)*a*b,
        4,
    ]

def symbolic_prior_equivalence():
    a,b,u,v,l=sp.symbols('a b u v l')
    ca,cb=1-u,1-v
    printed46=(l**4-2*l*l*(a*a+b*b+4)+8*a*b*l*(ca*(cb-1)-cb-1)
        +8*b*b*ca*(cb*cb-1)-8*cb*(a*a+(b*b-2)*cb)+(a*a-b*b)**2
        +8*ca*ca*(cb-1)*(a*a-2*cb-2))
    assert sp.expand(det(direct(a,b,u,v,l-2-a-b))-printed46)==0
    # At an interior stationary point A*u=B*v. Set x=v/A, p=A*B,
    # write J=J(d), and substitute the stationarity conditions in F.
    p,J,P,x=sp.symbols('p J P x')
    r=sp.resultant(4*p*x*x-3*p*x+J, P-8*p*p*x**3*(1-2*x),x)
    core=sp.cancel(r/(16*p**4))
    s=2-a-b
    before=core.subs({p:(4-a*a)*(4-b*b),J:s*(8-2*s+a*b)-a*b*(l-2-a-b),
        P:(l-2-a-b)*(l+2+a-b)*(l+2-a+b)*(l-2+a+b)})
    q=sum(c*l**i for i,c in enumerate(gigena_coefficients(a,b)))
    assert sp.expand(before-(a*b+2*l)**2*q)==0
    print('Symbolic equivalence: determinant = Gigena Eq. 46; stationary resultant factors as corrected sextic times stated factors: PASS')

def equivalence(a,s,lo,hi):
    b=invariants(a,s)[1]
    lam0,lam1=(2+a+b+d for d in (lo,hi))
    q=gigena_coefficients(a,b)
    assert value(q,lam0)*value(q,lam1)<0
    st=sturm(q)
    assert variations(st,lam0)-variations(st,lam1)==1
    pos_inf=[1 if p[-1]>0 else -1 for p in st]
    assert variations(st,lam1)==sum(x!=y for x,y in zip(pos_inf,pos_inf[1:]))
    uncorrected=gigena_coefficients(a,b,False)
    other=sturm(uncorrected)
    assert variations(other,lam0)==variations(other,lam1)

def original_checks(cases):
    blob=json.loads(cases.read_text())
    assert len(blob['results'])==17
    for data in blob['results']:
        a,s=map(F,(data['a'],data['s']))
        lo,hi=map(F,data['exact_delta_enclosure'])
        u,v=map(F,data['exact_witness_u_v'])
        assert 0<lo<hi and (hi-lo)/lo<F(1,10**17)
        lower_sign(a,s,lo,u,v)
        upper_sign(a,s,hi)
        equivalence(a,s,lo,hi)
    print('All 17 reported rational intervals: independent direct lower determinants and upper Sturm certificates PASS')
    print('All 17 intervals isolate the largest real root of the corrected Gigena sextic; printed-sign sextic has no root in them')

def sqrt_bracket(z, precision=90):
    p=10**precision
    n=isqrt(z.numerator*p*p//z.denominator)
    a,b=F(n,p),F(n+1,p)
    assert a*a<=z<b*b
    return a,b

def x_bracket(h):
    l,r=sqrt_bracket(9-8*h)
    return 2/(3+r),2/(3+l)

def g(x):
    return x*x*(1-x)

def extra():
    for t,h in [
        (F(1,100),F(1,10**12)),(F(1,100),1-F(1,10**12)),
        (F(1,10**12),F(1,10**30)),(F(1,10**12),1-F(1,10**30)),
        (F(1,10**40),F(1,2)),(F(1,10**40),1-F(1,10**40))
    ]:
        a,s=t*(1-h),t*h
        _,b,A,B,C=invariants(a,s)
        V=s*C/B
        k=2*V/A
        assert 0<k<1
        lo,hi=x_bracket(h)
        theta_lo,theta_hi=(g(hi)-t/5), (g(lo)+t/5)
        assert 0<theta_lo<theta_hi
        d0,d1=s*s*h*theta_lo,s*s*h*theta_hi
        v=V*sum(x_bracket(k))/2
        u=(s*C-B*v)/(2*(A-2*v))
        lower_sign(a,s,d0,u,v)
        upper_sign(a,s,d1)
        print('t =',t,', h =',h,': exact t/5 bracket PASS')

def opposite():
    uv=[2-e for e in (F(0),F(1,1000),F(1,100),F(1,50))]
    tested=0
    for a,s in [(F(1,1000),F(9,1000)),(F(1,10**30),F(1,1000))]:
        b=invariants(a,s)[1]
        for u,v in product(uv,repeat=2):
            if u*v<(2-s)**2:continue
            M=direct(a,b,u,v,F(0))
            for size in [1,2,3]:
                J=range(1,size+1)
                assert det([[M[i][j] for j in J] for i in J])>0
            tested+=1
    assert tested>10
    print(tested,'opposite-corner controls: exact compression minors all positive')

if __name__=='__main__':
    algebra()
    symbolic_prior_equivalence()
    original_checks(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('checks.json'))
    extra()
    opposite()
