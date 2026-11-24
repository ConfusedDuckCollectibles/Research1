#!/usr/bin/env python3
"""Symbolic proof identities and exact finite-s Bell-value controls."""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path

import sympy as S

HERE=Path(__file__).resolve().parent


def model_checks():
    a,b,u,v,x,y,d,c = S.symbols("a b u v x y d c")
    Ku=S.Matrix([[-u,x],[x,u]])
    Kv=S.Matrix([[-v,y],[y,v]])
    K=S.kronecker_product(Ku,Kv)
    D=S.diag(d,4+2*b+d,4+2*a+d,2*(a+b)+d)
    det=(D+K).det(method="berkowitz").expand()
    det=S.rem(S.rem(det,x*x-u*(2-u),x),y*y-v*(2-v),y)
    s=2-a-b
    formula=d*(d+4+2*a)*(d+4+2*b)*(d+2*(a+b))+8*u*v*(a*b*d-s*(8-2*s+a*b)+(4-a*a)*u+(4-b*b)*v-2*u*v)
    assert S.expand(det-formula)==0
    assert S.diff(formula,d).subs({a:0,b:2,u:0,d:0})==128
    h,x=S.symbols("h x")
    g=x*(1-x)**2/(2*(1-h*x))
    assert S.factor(S.diff(g,x)-(1-x)*(1-3*x+2*h*x*x)/(2*(1-h*x)**2))==0
    optx=2/(3+S.sqrt(9-8*h))
    simple=optx**2*(1-optx)
    assert simple.subs(h,0)==S.Rational(2,27) and simple.subs(h,1)==S.Rational(1,8)
    assert S.expand(S.series(S.sqrt(8+2*(2-c)**2)-(4-c),c,0,3).removeO())==c*c/8
    assert max(A0*B0+A0*B1+A1*B0-A1*B1 for A0,A1,B0,B1 in product([-1,1],repeat=4))==2
    assert all(C*(x+y)<D for C,x,y,D in [(Q(1),Q("0.024"),Q("1.072"),Q(11,10))])
    # Rational Schur positivity constants for t=a+s<=.01.
    R1,R2,off=Q(398,100),Q(794,100),Q(2,5)
    assert R1*R2>off**2
    assert R2/(R1*R2-off**2)<Q(3,10)
    assert R1/(R1*R2-off**2)<Q(13,100)
    assert off/(R1*R2-off**2)<Q(1,50)
    # Loss mapping, no erased trial omitted.
    EA,EB,X,Y,Z=S.symbols("EA EB X Y Z")
    Atilt,Btilt=2*(1-EB)/EB,2*(1-EA)/EA
    observed=EA*EB*Z+2*EA*(1-EB)*X+2*(1-EA)*EB*Y+2*(1-EA)*(1-EB)
    assert S.factor(observed-2-EA*EB*(Z+Atilt*X+Btilt*Y-2-Atilt-Btilt))==0
    return str(S.factor(formula.subs(d,0)))


def sqrt_interval(r, scale=10**90):
    assert r>=0
    n=isqrt(r.numerator*scale**2//r.denominator)
    low,high=Q(n,scale),Q(n+1,scale)
    assert low*low<=r<high*high
    return low,high


def certify(a,s,d):
    """Prove sign of minimum det(L+d-B) over the full angle box.
    Returns +1 if globally positive, -1 for a definite rational witness.
    No pseudo-moments, floating eigenvalues, or optimizer required.
    """
    t=a+s
    assert 0<s<t<=Q(1,100)
    assert d>=0
    b=2-t
    A, B, C=4-a*a, 4-b*b, 8-2*s+a*b
    D=s*C-a*b*d
    base=d*(d+4+2*a)*(d+4+2*b)*(d+2*(a+b))
    vcap=(s*C+2*a*a)/B
    assert vcap<A/2 and vcap<=2
    if D<=0:
        return 1, None
    # Any gain uv(D-Au-Bv+2uv)>0 forces v<=vcap. For such v
    # A-2v>0, so unconstrained u gives a valid upper envelope.
    # Its unique maximum in [0,D/B] is the first quadratic root.
    top=D/B
    p0,p1=sqrt_interval(9*A*A-16*A*D/B)
    vlo,vhi=(3*A-p1)/8,(3*A-p0)/8
    assert 0<vlo<=vhi<top and vhi<A/2
    envelope_hi=vhi*(D-B*vlo)**2/(4*(A-2*vhi))
    if base>8*envelope_hi:
        return 1, None
    v=(vlo+vhi)/2
    u=(D-B*v)/(2*(A-2*v))
    assert 0<=u<=2 and 0<=v<=2
    det=base-8*u*v*(D-A*u-B*v+2*u*v)
    if det<0:
        return -1,(u,v)
    return 0,None


def format_q(x, precision=23):
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec=precision
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def exact_case(a,s,n_bisect=94):
    lo,hi=Q(0),s*s
    _,wit=certify(a,s,lo)
    assert wit is not None
    z,_=certify(a,s,hi)
    assert z==1
    for _ in range(n_bisect):
        m=(lo+hi)/2
        z,w=certify(a,s,m)
        assert z!=0
        if z<0:
            lo,wit=m,w
        else:
            hi=m
    assert lo>0 and (hi-lo)/lo<Q(1,10**17)
    zz,_=certify(a,s,hi); zz2,_=certify(a,s,lo)
    assert zz==1 and zz2==-1
    t=a+s
    h=s/t
    q0,q1=sqrt_interval(9-8*h)
    xlo,xhi=2/(3+q1),2/(3+q0)
    # g is increasing on [1/3,1/2]; evaluation at bracket endpoints valid.
    g0=xlo*xlo*(1-xlo)
    g1=xhi*xhi*(1-xhi)
    return dict(a=str(a),s=str(s),t=str(t),h=str(h),
                exact_delta_enclosure=[str(lo),str(hi)],
                delta_display=[format_q(lo,14),format_q(hi,14)],
                scaled_delta_display=[format_q(lo/(s*s*h),15),format_q(hi/(s*s*h),15)],
                predicted_g_display=[format_q(g0,15),format_q(g1,15)],
                exact_witness_u_v=[str(x) for x in wit])


def main():
    c=model_checks()
    params=[]
    for h in (Q(9,10),Q(1,2),Q(1,3),Q(1,10),Q(1,1000)):
        for t in (Q(1,100),Q(1,10000),Q(1,10**7)):
            params.append((t*(1-h),t*h))
    # Variable ratio: s << a, both approach the apex. And a << s.
    params += [(Q(1,10000),Q(1,10**10)),(Q(1,10**12),Q(1,10**5))]
    cases=[exact_case(a,s) for a,s in params]
    check=dict(explicit_factorized_d0=c, n_exact_finite_s_cases=len(cases),
         n_successfully_global_bracketed_cases=len(cases), results=cases,
         code_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/"checks.json").write_text(json.dumps(check,indent=2)+"\n")
    for case in cases:
        print("a",case["a"],"s",case["s"],"h",case["h"],"scaled",case["scaled_delta_display"][0],"g(h)",case["predicted_g_display"][0])


if __name__=="__main__":
    main()
