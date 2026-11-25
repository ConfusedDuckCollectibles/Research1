#!/usr/bin/env python3
"""Independent tent identity controls and certified re-evaluation of fixed witness."""
from fractions import Fraction as F
from math import factorial
from verify import I, add, mul, div, root, art, phi, log, fmt, simulate_algebra

def mt(A): return [list(x) for x in zip(*A)]
def mm(A,B): return [[sum(a*b for a,b in zip(r,c)) for c in zip(*B)] for r in A]
def ms(A,B): return [[a+b for a,b in zip(r,s)] for r,s in zip(A,B)]
def mc(a,A): return [[a*x for x in row] for row in A]
def md(d): return [[d[i] if i==j else F(0) for j in range(len(d))] for i in range(len(d))]
def part(A,I,J): return [[A[i][j] for j in J] for i in I]
def powers(A,N):
    result=[md([F(1)]*len(A))]
    for _ in range(N): result.append(mm(result[-1],A))
    return result

def ring_control():
    # Nonuniform stationary law, both hidden-internal jumps and direct V-V.
    pi=[F(k,10) for k in (1,2,3,4)]; V=[0,3]; H=[1,2]
    Q=[[F(0)]*4 for _ in range(4)]
    for i in range(4):
        Q[i][i]=-5/pi[i]; Q[i][(i+1)%4]=4/pi[i]; Q[i][(i-1)%4]=1/pi[i]
    Qr=[[pi[j]*Q[j][i]/pi[i] for j in range(4)] for i in range(4)]
    K=[[F(0)]*4 for _ in range(4)]
    for i in range(4):
        K[i][i]=Q[i][i]; K[i][(i+1)%4]=K[i][(i-1)%4]=2/pi[i]
    assert mm([pi],Q)==[[0,0,0,0]]
    assert all(sum(row)==0 for row in Qr)
    assert Q[1][2]>0 and Q[2][1]>0 and Q[0][3]>0
    assert all((K[i][j]**2==Q[i][j]*Qr[i][j]) for i in range(4) for j in range(4))
    assert Q[0][1]*Q[1][2]*Q[2][3]*pi[0]/(Q[3][2]*Q[2][1]*Q[1][0]*pi[3])==64
    Pih=md([pi[i] for i in H]); Piv=md([pi[i] for i in V])
    N=13
    # Case 1: ordinary forward observed rates (split uses Q^* on left).
    # Case 2: geometric stopped-path affinities (split uses geometric L on both halves).
    for upper,left in [(Q,Qr),(K,K)]:
        A=part(upper,V,H); B=part(upper,H,V); L=part(upper,H,H)
        Bleft=part(left,H,V); Lleft=part(left,H,H)
        assert L[0][0]!=L[1][1] and L[0][1]!=L[1][0]
        assert mm(mt(Lleft),Pih)==mm(Pih,L)
        assert mm(mt(Bleft),Pih)==mm(Piv,A)
        pl= powers(Lleft,N); pr=powers(L,N)
        for n in range(N+1):
            # Formal coefficient of tau^(n+2) from two independently truncated halves.
            split=[[F(0)]*2 for _ in range(2)]
            for k in range(n+1):
                l=mm(pl[k],Bleft); r=mm(pr[n-k],B)
                term=mm(mm(mt(l),Pih),r)
                split=ms(split,mc(F(1,factorial(k+1)*factorial(n-k+1)),term))
            # Integral of w_tau(s)*s^n/n! over s gives this coefficient.
            c=F(2**(n+2)-2,factorial(n+2))
            marked=mc(c,mm(mm(mm(Piv,A),pr[n]),B))
            assert split==marked
            if upper is K: assert split==mt(split)
    print("nonuniform hidden-jump ring: first 14 coefficients agree for ordinary and geometric tent integration")

def lower_gt(x,decimal): assert x[0]>F(decimal)
def upper_lt(x,decimal): assert x[1]<F(decimal)

def numbers():
    e,z,a=F(1,32),F(1,16),F(1,4)
    alpha,Z=simulate_algebra(e,z,a)
    t=mul(I(4),log(I(2))); p=1-z; theta=a*e/2
    # 2 * p_H = 1, since p_H=1/2; compare old upper (2*p_H F).
    new=div(mul(I(p*p*alpha),add(I(1),log(I(2/e)))),t)
    x=root(I((1-z*z)*alpha)); old=div(mul(x,art(x)),t)
    ep=add(mul(I(4*e/Z),phi(root(I(1-theta)))),
           mul(I(alpha),phi(I(e/2))))
    fraction=div(new,ep); improvement=div(new,old); asymptotic=div(I(p*p),t)
    lower_gt(new,"0.0001876917017181")
    upper_lt(old,"0.00004123477217914")
    lower_gt(fraction,"0.279420235595")
    lower_gt(improvement,"4.55178219253")
    lower_gt(asymptotic,"0.3169"); upper_lt(asymptotic,"0.3171")
    print("specialized rate lower bound:",fmt(new))
    print("optimized prior rate ceiling:",fmt(old))
    print("guaranteed improvement:",fmt(improvement))
    print("certified fraction of actual model rate:",fmt(fraction))
    print("limiting certified fraction in witness family:",fmt(asymptotic))
    print("PASS: all outward-rounded thresholds are strict for rational interval endpoints.")

if __name__=="__main__":
    ring_control()
    numbers()
