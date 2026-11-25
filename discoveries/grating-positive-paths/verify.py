"""Exact rational gray-Debye continuum grating and sharp cone certificate."""
from fractions import Fraction as Q
from math import factorial


def mean(a):
    # Real-basis polynomial sum a[k]*(-i*z)^k, under uniform density dz/2.
    return sum(((-1)**(k//2)*v/Q(k+1) for k,v in enumerate(a) if k%2==0),Q(0))


def mean_z(a):
    # Real number b such that expected z*f equals -i*b.
    return sum(((-1)**((k-1)//2)*v/Q(k+2)
                for k,v in enumerate(a) if k%2),Q(0))


def apply_A(a,c,N,R):
    result=[-(N+R)*v for v in a]+[Q(0)]
    for k,v in enumerate(a):
        result[k+1]+=c*v
    result[0]+=(N+R)*mean(a)
    result[1]+=3*N*mean_z(a)
    return result


def propagate(c,N,R,order=48):
    term=[Q(1)]
    exact_sum=Q(1)
    for n in range(1,order+1):
        term=[v/n for v in apply_A(term,c,N,R)]
        exact_sum+=mean(term)
    # ||e^(sA)|| <= 1 for s>=0; integral remainder, not a heuristic series tail.
    tail=(N+R+c)**(order+1)/factorial(order+1)
    assert tail<Q(1,10**12)
    return exact_sum-tail,exact_sum+tail


def trig(c,shift,order=12):
    # cos for shift=0, and sin(c)/c for shift=1.
    assert 0<=c<=Q(6,5) and c*c<Q((2+shift)*(1+shift))
    terms=[(-1)**k*c**(2*k)/factorial(2*k+shift) for k in range(order+2)]
    a,b=sum(terms[:-1]),sum(terms)
    return min(a,b),max(a,b)


def diff(a,b):
    return a[0]-b[1],a[1]-b[0]


def ratio(c1,c2):
    f1,f2=trig(c1,0),trig(c2,0)
    assert 0<f2[0]<f1[0]
    return (c2/c1)**2*f2[0]/f1[1],(c2/c1)**2*f2[1]/f1[0]


def main():
    c1,c2=Q(3,5),Q(6,5)
    N,R=Q(8),Q(7,100)
    D=[diff(propagate(c,N,R),trig(c,1)) for c in (c1,c2)]
    (l1,u1),(l2,u2)=D
    assert Q(1047,10**6)<l1<=u1<Q(1048,10**6)
    assert Q(478,10**6)<l2<=u2<Q(479,10**6)
    kl,ku=ratio(c1,c2)
    assert ku<Q(9,5)
    contrast=kl*l1-u2
    assert contrast>Q(13,10000)
    assert (ku+1)*Q(1,2500)<contrast
    print("collective model and both residuals are positive; exact D1 in (1047,1048)*1e-6")
    print("exact D2 in (478,479)*1e-6; passive two-tone inequality violation >13/10000")
    print("absolute errors of 1/2500 on each residual still leave strict separation")

    # Positive-control null: remove the normal collision projection.
    (a,b),(d,e)=[diff(propagate(c,0,R),trig(c,1)) for c in (c1,c2)]
    assert a>0 and d>=ku*b and e<=(c2/c1)**2*a
    print("resistive-only stationary-reset control satisfies both finite cone edges")

    # q^2 t^3 and q^4 t^5 collision coefficients, independently expanded.
    for c in (Q(0),Q(1,2),Q(1)):
        a=[Q(1)]
        for k in range(1,6):
            a=apply_A(a,c,N,R)
            if k==3:
                assert mean(a)==R*c*c/3
        assert mean(a)==R**3*c*c/3-(4*N+22*R)*c**4/45
    # Even degrees <=2 and <=4 make the preceding values exact polynomial checks.
    assert (4*N+22*R)/45>10*R/3
    assert (4*32+22)==150
    print("exact early-time coefficients verify strict 32-to-1 threshold for the moment test")
    print("PASS")


if __name__=="__main__":
    main()
