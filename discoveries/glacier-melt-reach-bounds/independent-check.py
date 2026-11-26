"""Independent exact replay of synthetic controls (Python standard library)."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
import json

# Sharp signed tide example: calculate primitive areas independently by
# trapezoids (including exact zero crossing of the piecewise-linear Peano kernel).
a, b_t = (-80,-50,0),(-20,10,60)
stencil = (1,-2,1)
pairs = list(zip(a, stencil)) + list(zip(b_t, (-c for c in stencil)))
atoms = defaultdict(int)
for x,c in pairs: atoms[x] += c
knots = sorted(atoms)
G = 0
W1 = 0
P1 = F(0)
H = lambda t: sum(c*max(x-t,0) for x,c in pairs)
for x,y in zip(knots,knots[1:]):
    G += atoms[x]
    W1 += abs(G)*(y-x)
    c,d = H(x),H(y)
    if c*d >= 0: P1 += (y-x)*F(abs(c)+abs(d),2)
    else: P1 += (y-x)*F(c*c+d*d,2*(abs(c)+abs(d)))
assert W1 == 160
contrast = lambda ts, fun: sum(c*fun(x) for x,c in zip(ts,stencil))
assert contrast(a, lambda x:x) == contrast(b_t, lambda x:x) == 20
hinge = lambda x: F(max(x-20,0),4)
assert contrast(a,hinge) == 0 and contrast(b_t,hinge) == 10
# Attain endpoint bounds using actual slopes in the intervals.
for sign in (-1,1):
    G = 0
    yval = {knots[0]: F(0)}
    for x,y in zip(knots,knots[1:]):
        G += atoms[x]
        yval[y] = yval[x] + (y-x)*F(int(sign*G < 0),4)
    assert sum(c*yval[x] for x,c in pairs) == sign*20

# Independent finite LP by interpolation on all 3-support distributions in a
# 1/8-km grid. Lagrange moments produce exact weights, dropping zero masses.
# Adding mean and second moment cuts gives ≤3 support at finite LP vertices.
mu, s2, ridge, cap = F(3),F(15),F(4),F(6)
def distal(l):
    return (l-ridge)**2/(2*l) if l > ridge else F(0)
lo = hi = None
ats = [F(i,8) for i in range(49)]
feasible = 0
for xs in combinations(ats,3):
    w = []
    for i,x in enumerate(xs):
        y,z = (l for j,l in enumerate(xs) if i != j)
        w.append((s2-(y+z)*mu+y*z)/((x-y)*(x-z)))
    if min(w) < 0: continue
    feasible += 1
    val = sum(p*distal(x) for p,x in zip(w,xs))
    if lo is None or lo > val: lo = val
    if hi is None or hi < val: hi = val
assert (lo,hi) == (F(3,50),F(2,15))

# High-ridge interior upper optimum: choose rational contact points,
# derive the moments, and verify the dual at the touching point exactly.
c = F(29,5)
a_tan = lambda u: u*(3*ridge-u)/(u+ridge)
a = a_tan(c)
um,us = (a+c)/2,(a*a+c*c)/2
low_atom = lambda u: (um*u-us)/(u-um)
ar = low_atom(cap)
assert a == F(899,245) and low_atom(c) == a
assert 0 <= a < ridge < c < cap and 0 <= ar < um < cap
assert ar > a_tan(cap)
assert c**3-3*ridge*c**2+(4*ridge*um-us)*c-ridge*us == 0
P = lambda t: 2*(ridge-a)*t-(t-ridge)*(t-a)
assert P(ridge)>0 and P(c)==0 and P(cap)<0
assert P(ridge)/1>0 and P(F(10))<0
inter = distal(c)/2
old = (um-ar)/(cap-ar)*distal(cap)
assert inter == F(81,580) and inter > old

# Positive-length mixture moment inversion and diagnostic.
law = [(F(2),F(1,3)),(F(5),F(2,3))]
Mj = lambda j: sum(prob*L**(j+1) for L,prob in law)/((j+1)*(j+2))
M0,M1,M2 = (Mj(j) for j in range(3))
effective_L = 3*M1/M0; effective_A = 2*M0/effective_L
assert M0*M2 > F(3,2)*M1*M1
assert effective_L != sum(p*L for L,p in law)

result = {
 "tide_W1_cm":W1, "tide_Peano_cm2":str(P1), "tide_box_endpoint_abs_cm":20,
 "moment_grid_feasible_three_supports":feasible, "moment_grid_distal":list(map(str,(lo,hi))),
 "interior_upper_rational":list(map(str,(a,c,inter,old))),
 "effective_fixed_ramp_moments":list(map(str,(M0,M1,M2,effective_L,effective_A)))
}
print(json.dumps(result, separators=(",",":")))
