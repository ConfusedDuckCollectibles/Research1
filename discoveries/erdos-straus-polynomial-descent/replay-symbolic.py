#!/usr/bin/env python3
"""Algebra sanity checks, not a search or proof of the universal lemma."""
from itertools import product
from math import gcd, isqrt
from sympy import Poly, QQ, expand, symbols

u,v,P,A,B,C,D,a,b,e,r,s,R,S,U,W,H,G,lam = symbols("u v P A B C D a b e r s R S U W H G lam")
N=u**2+v**2

def exact(p, abc):
    x,y,z=abc
    assert expand(4*x*y*z-p*(x*y+x*z+y*z))==0

# Published modular constructions substituted into the norm:
I=((3*N+1)/11,3*N+1,N*(3*N+1)/4)
II=((N+1)/3,N,N*(N+1)/3)
exact(N,I); exact(N,II)
for (i,j),trip in [((5,2),I),((1,2),II)]:
    p=int(N.subs({u:i,v:j}))
    values=tuple(t.subs({u:i,v:j}) for t in trip)
    assert all(t.is_Integer and t>0 for t in values)
    assert 4*values[0]*values[1]*values[2] == p*sum(values[i]*values[j] for i,j in [(0,1),(0,2),(1,2)])
    print("norm control", p, values)
# Check generic Type-I and Type-II algebra independently by polynomial expansion.
original=4*A*B*C-P*(A*B+B*C+C*A)
assert expand(original.subs(C,P*R) - P*((4*R-1)*A*B-P*R*(A+B)))==0
assert expand(((4*R-1)*A*B-P*R*(A+B)).subs(R,(1+P*U)/4) -
              P*(4*U*A*B-(1+P*U)*(A+B))/4)==0
assert expand((4*U*A*B-W*(A+B)).subs({A:D*a,B:D*b,W:a*b*e}) -
              D*a*b*(4*U*D-e*(a+b)))==0
assert expand((P*(a+b)-lam*(a*b*e-1)) -
              (b*(P-lam*a*e)+P*a+lam))==0
assert expand(original.subs({B:P*R,C:P*S}) -
              P**2*(A*(4*R*S-R-S)-P*R*S))==0
assert expand(((4*R-1)*(4*S-1)-1)-4*(4*R*S-R-S))==0
assert expand((A*(4*R*S-R-S)-P*R*S).subs({R:D*r,S:D*s}) -
              D*(A*(4*D*r*s-r-s)-P*D*r*s))==0
assert expand((4*D*r*s-r-s-lam*P*D).subs({r:a,s:b+e*P}) -
              (D*(4*a*b+(4*a*e-lam)*P)-(a+b+e*P)))==0
print("symbolic transformations verified")
# Countercontrol: nonconstant units make genuinely 2-variable identities possible.
M=4*u*v-u-v
J=(u*v,M*u,M*v)
exact(M,J)
ring=Poly(M,u,v,domain=QQ)
assert Poly((4*u-1)*(4*v-1)-1,u,v,domain=QQ).rem(ring).is_zero
assert not Poly(u*v,u,v,domain=QQ).rem(ring).is_ground
sq=sorted({i*i%840 for i in range(840) if gcd(i,840)==1})
assert sq==[1,121,169,289,361,529]
for i,j in product(range(1,100),repeat=2):
    p=int(M.subs({u:i,v:j}))
    if p%840 in sq and p>2 and all(p%d for d in range(2,isqrt(p)+1)):
        values=tuple(t.subs({u:i,v:j}) for t in J)
        assert all(t.is_Integer and t>0 for t in values)
        print("nonconstant-unit countercontrol on a Mordell prime",p,(i,j),values)
        break
else:
    raise AssertionError
print("no finite conjecture verification performed")
