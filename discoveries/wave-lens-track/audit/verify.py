"""Exact symbolic controls, not a replacement for the continuum proof."""
import sympy as sp
x,k,s,y,w,e,u,v = sp.symbols('x k s y w e u v', real=True)
tau=x-k*sp.log(1+x)
kap=sp.factor(1-sp.diff(tau,x)-x*sp.diff(tau,x,2))
assert kap==k/(x+1)**2
assert sp.factor(sp.diff(kap,x))==-2*k/(x+1)**3
assert sp.factor(sp.diff(kap,x,2))==6*k/(x+1)**4
assert sp.factor(sp.diff(tau,x).subs(x,k-1))==0
# The squared-radius Abel projection of (1+v/2)^(-5/2)
# has coefficient sqrt(2) B(1/2,2) and the exponent decreases to -2.
assert sp.simplify(sp.gamma(sp.Rational(1,2))*sp.gamma(2)/
                    sp.gamma(sp.Rational(5,2)))==sp.Rational(4,3)
print("Exact Plummer potential, positive-convex projection, and Abel constant")

j0=sp.besselj(0,w*y*sp.sqrt(2*x)).series(y,0,5).removeO().expand()
assert j0.coeff(y,2)==-w*w*x/2
assert j0.coeff(y,4)==w**4*x*x/16
print("Bessel offset coefficients: -w^2 x/2 and w^4 x^2/16")

# Exact common shear on synthetic local branches; s indexes delays.
lo,hi=3-s,3+s
h=e*s**5*(1-s)**5
deltaW=sp.expand((hi+h)-(lo+h)-(hi-lo))
M=lambda a,b:(b*b-a*a)/2
assert deltaW==0 and sp.expand(M(lo+h,hi+h)-M(lo,hi)-h*(hi-lo))==0
assert sp.expand(hi+h-lo-(hi-lo)-h)==0  # one-sided shift is not invisible
for sign,expected in ((1,hi+h),(-1,lo+h)):
    assert sp.simplify(M(lo+h,hi+h)/(hi-lo)+sign*(hi-lo)/2-expected)==0
print("Common shear preserves width, its shifted first moment inverts, and one-sided test fails")

# Affine-track Jacobian test on a compact fibre patch around the effective
# v-minimum c. The identity here is algebra; no physical claim for this toy.
c=sp.Rational(1,3)
sigma=(v-c)**2+u**2
H=e*(s-u**2)*(1-(s-u**2))**5
jac=sp.diff(v+H.subs(s,sigma),v)
assert sp.expand(jac-1-sp.diff(H,s).subs(s,sigma)*sp.diff(sigma,v))==0
for j in range(5):
    # Derivative of a primitive in s. v endpoints have the same sigma.
    changed=sp.Poly(sp.expand(sigma**j*(jac-1)),v)
    primitive=changed.integrate()
    assert sp.expand(primitive.eval(c+1)-primitive.eval(c-1))==0
# At exactly one delay, translating both ends by H moves the area barycentre.
a,b=c-s,c+s
h0=H.subs(s,s*s+u*u)
assert sp.simplify(M(a+h0,b+h0)-M(a,b)-(b-a)*h0)==0
print("Affine-track conditional Jacobian moments exact; first transverse moment changes")
