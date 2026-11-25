#!/usr/bin/env python3
"""Independent exact controls: raw transfer metric, finite-state count DP."""
from fractions import Fraction as F
from itertools import product, combinations
from math import comb

def mm(A,B):
    return tuple(tuple(sum(a*b for a,b in zip(row,col)) for col in zip(*B)) for row in A)

def tr(A):return A[0][0]+A[1][1]
def det(A):return A[0][0]*A[1][1]-A[0][1]*A[1][0]
def trans(A):return tuple(zip(*A))
def exact_inv(A):
    d=det(A)
    return ((A[1][1]/d,-A[0][1]/d),(-A[1][0]/d,A[0][0]/d))
I=((F(1),F(0)),(F(0),F(1)))

def raw(v,E=F(0)):
    return ((E-v,-F(1)),(F(1),F(0)))

def trans_prob(M,v):
    # C0.T*C0, so no radicals even for the commensurate v=1 control.
    Q=((F(1),v/2),(v/2,F(1)))
    norm=tr(mm(mm(mm(trans(M),Q),M),exact_inv(Q)))
    assert det(M)==1
    return F(4)/(norm+2)

def Tmu(mu):
    return 4*mu*mu/(1+mu*mu)**2

def k2_present_B_counts(N):
    # Count words with exactly two defects among AA, BB, AB.
    # Dynamic program tracks the stages before, within, and after defect pair.
    ns=1
    between=[0]*(N+1)
    after=[0]*(N+1)
    for k in range(N):
        # append a defect to the middle, or append either background after it
        after=[2*x+y for x,y in zip(after,between)]
        nxt=between[:]
        nxt[0]+=ns                 # first defect
        for j in range(k+1):
            nxt[j+1]+=between[j]   # append BB, distinct from AA no increment
        between=nxt
        ns*=2
    assert sum(after)==comb(N,2)*2**(N-2)
    return after

def small_word_control(N):
    tally=[0]*(N+1)
    for w in product(range(3),repeat=N):
        spots=[i for i,x in enumerate(w) if x==2]
        if len(spots)==2:tally[w[spots[0]+1:spots[1]].count(1)]+=1
    assert tally==k2_present_B_counts(N)

def benchmark(v,N):
    A,B=raw(F(0)),raw(v)
    AA,BB,D=mm(A,A),mm(B,B),mm(B,A)
    assert AA==((-1,0),(0,-1))
    assert trans_prob(BB,v)==1
    assert D==((-1,v),(0,-1))
    assert trans_prob(D,v)==F(16,25)
    assert trans_prob(mm(D,D),v)==F(4,13)
    post=k2_present_B_counts(N)
    P=I;wins=0;succ=[]
    for j in range(N+1):
       if trans_prob(mm(mm(D,P),D),v)<=F(73,100):
           wins+=post[j]
           if post[j]:succ.append(j)
       P=mm(BB,P)
    probability=F(wins,sum(post))
    assert F(786229680,10**9)<probability<F(786229682,10**9)
    p=F(1,N)
    exactly_one=N*p*(1-p)**(N-1)
    exactly_two=comb(N,2)*p*p*(1-p)**(N-2)
    per=exactly_one+exactly_two*probability
    power=1-(1-per)**4
    ceiling=1-(1-p)**(4*N)
    assert per>F(5131,10000) and power>F(9438,10000)
    print("detected BB-count entries",len(succ),"exact K2 conditional",f"{float(probability):.12f}")
    print("single lower",f"{float(per):.12f}","four lower",f"{float(power):.12f}",
          "zero-size four upper",f"{float(ceiling):.12f}")

def certs(N):
    e=F(1,10**5); kappa=(1+e)/(1-e)
    err=kappa/(1-6*N*e)-1
    assert err<F(25,1000)
    for mu,error,sense in [(1+err,-F('0.02'),1),(2-2*err,F('0.02'),-1),
                            (F('1.77')-4*err,F('0.02'),-1)]:
        assert sense*(Tmu(mu)+error-F('0.80'))>0
    assert Tmu(F('1.77'))>F('0.73')
    e=F(1,2000); kappa=(1+e)/(1-e)
    Ghi=kappa*(2+F(5,4)*e)
    Glo=(2-F(5,4)*e)/kappa
    assert Glo>F('1.99') and Ghi<F('2.01')
    # Elementary coefficients in the radial/traceless map.
    assert (Ghi**2+Ghi**(-2)-2)/4 < 1
    assert Ghi+e*Ghi**2/2 < 3
    c0=F(-3,5)
    assert max(abs(2*(c0+e/2)**2-1),abs(2*(c0-e/2)**2-1))+3*e<F('.31')
    assert 1+F(9)/F('.69')<15
    cap=15*N*e*e/(2*(1-15*N*e*e))
    roundmu=F('1.08')
    cut=F(1)/Tmu(roundmu)-1
    assert Tmu(Glo/roundmu**2)+F('.02')<F('.80')
    fail=8*cap/(1/F('.82')-1)
    fits=(1-F(1,N))**(N-1)*(1-cap/cut)
    power=1-(1-fits)**8
    assert fail<F('.028') and power>F('.95')
    print("norm perturbation factor",f"{float(err):.12f}",
          "wide false <=",f"{float(fail):.12f}","wide power >=",f"{float(power):.12f}")

def commensurate():
    v=F(1);A,B=raw(F(0)),raw(v)
    D,R=mm(B,A),mm(B,B)
    assert trans_prob(D,v)==F(3,4)
    P=I;values=[]
    for j in range(3):
        values.append(trans_prob(mm(mm(D,P),D),v))
        P=mm(P,R)
    assert P==I
    assert values==[F(3,7),F(1),F(3,7)]
    print("commensurate v=1 two-defect T values by BB count modulo 3:",values)

if __name__=="__main__":
    for n in range(2,8):small_word_control(n)
    benchmark(F(6,5),400)
    certs(400)
    commensurate()
