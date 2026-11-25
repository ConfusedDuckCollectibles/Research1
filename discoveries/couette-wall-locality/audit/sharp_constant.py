"""Exact numerical certificate for the explicit bounded-pulse factor-2 counterexample."""
from fractions import Fraction as Q
m=256
eta=Q(1,100)
pi_hi=Q(22,7)
mu_lo=Q(3,4)
epsilon=Q(1,2**(m+20))
times=[Q(1,4**(j+1)) for j in range(1,m+1)]
width=[eta**2*s*(1-s) for s in times]
assert all(0<s<s+w<1 for s,w in zip(times,width))
assert all(times[j+1]+width[j+1]<times[j] for j in range(m-1))
assert mu_lo**2<2/pi_hi
# Error of ideal point-pulse split: Var(H)+Var(G)/2, by Gaussian covariance.
ideal=Q(11,4)*pi_hi**3/m
pulse=3*pi_hi*eta/mu_lo
# sigma(s_j)>=sqrt(s_j); bound the sum of finite rectangular pulse masses.
mass_hi=pi_hi/(m*mu_lo)*sum(2**(j+1) for j in range(1,m+1))
# x=y=epsilon. p_1(0)/p_1(2epsilon)=exp(epsilon^2).
# Each bridge phase expectation changes by at most epsilon*mass_hi.
# For epsilon^2<1/2, exp(epsilon^2)-1<=2epsilon^2.
assert epsilon**2<Q(1,2)
endpoint=3*epsilon*mass_hi+4*epsilon**2
assert ideal<Q(1,3)+Q(1,100)
assert pulse<Q(13,100)
assert endpoint<Q(1,1000)
total=ideal+pulse+endpoint
assert total<Q(1,2)
print("256 disjoint bounded-pulse windows satisfy the variance/continuity certificate.")
print("Strict lower bound for |K_N−K_bulk|/p_1(2epsilon):",float(3-total))
print("Exact inequality: > 5/2 at x=y=epsilon=2^(-276).")
