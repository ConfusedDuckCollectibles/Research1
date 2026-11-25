"""Exact Gaussian-free transfer and rare-event probability oracle."""
from fractions import Fraction as F
from math import comb


def multiply(A, B):
    return tuple(sum(A[2*i+k]*B[2*k+j] for k in range(2))
                 for i in range(2) for j in range(2))


I = (1, 0, 0, 1)


def power(A, n):
    X = I
    while n:
        if n%2:
            X = multiply(X, A)
        n //= 2
        A = multiply(A, A)
    return X


def transfer(v, E=0):
    return (E-v, -1, 1, 0)


def flux(M):
    return multiply(multiply((F(4,5),0,F(3,5),1),M),(F(5,4),0,F(-3,4),1))


def T(M):
    assert M[0]*M[3]-M[1]*M[2] == 1
    return F(4)/(2+sum(x*x for x in M))


def T_mu(mu):
    return 4*mu**2/(1+mu**2)**2


def finite_probability():
    N = 400
    A = flux(transfer(0))
    B = flux(transfer(F(6,5)))
    D = multiply(B, A)
    assert multiply(A,A) == (-1,0,0,-1)
    assert B == (F(-3,5),F(-4,5),F(4,5),F(-3,5))
    R = multiply(B,B)
    assert R == (F(-7,25),F(24,25),F(-24,25),F(-7,25))
    assert T(D) == F(16,25)
    assert T(I) == 1
    assert T(power(D,2)) < F(4,5)
    assert sum(x*x for x in D) == F(17,4)
    for l in range(12):
        assert T(multiply(power(R,l),I)) == 1

    score = []
    Rpower = I
    for j in range(N-1):
        score.append(T(multiply(multiply(D,Rpower),D)) <= F(73,100))
        Rpower = multiply(R,Rpower)
    # Conditional count b, exact integer weights with common denominator
    # choose(N,2)*2^(N-2). Gap g has N-1-g placements of this size.
    norm = comb(N,2) * (1 << (N-2))
    wins = 0
    mass = 0
    for g in range(N-1):
        mul = (N-1-g) << (N-2-g)
        mass += mul * (1 << g)
        wins += mul * sum(comb(g,j) for j in range(g+1) if score[j])
    assert norm == mass
    conditional = F(wins,norm)
    assert F("0.7862") < conditional < F("0.7863"), str(float(conditional))
    p = F(1,N)
    p1 = N*p*(1-p)**(N-1)
    p2 = comb(N,2)*p*p*(1-p)**(N-2)
    lower = p1+p2*conditional
    assert lower > F("0.5131"), str(float(lower))
    assert 1-(1-lower)**4 > F("0.9438")
    print(f"K=2 P(T0<=.73 | K=2): {float(conditional):.9f} (rounded, exact weights)")
    print(f"single-chain robust probability lower bound: {float(lower):.9f} (rounded)")
    print(f"four-chain robust probability lower bound: {float(1-(1-lower)**4):.9f} (rounded)")


def robustness():
    # For |E|<=e, fixed normalization block perturbation is <6e.
    e = F(1,100000)
    N = 400
    assert F(5,4)*(4+F(5,4)*e) < 6
    assert F("0.79")**2 < 1-(F(3,5)+e/2)**2
    assert F("0.61")/(2*F("0.79")) < F("0.4")
    assert F("0.4")**2+F("0.5")**2 < F("0.8")**2
    prefactor = (1+e)/((1-e)*(1-6*N*e)) - 1
    assert prefactor < F("0.025")
    # Null, one defective block, and the certified class of two-block gaps.
    assert T_mu(F("1.025"))-F(".02") > F(".8")
    assert T_mu(F("1.95"))+F(".02") < F(".8")
    assert T_mu(F("1.77")) > F(".73")
    assert F("1.77")-4*F(".025") == F("1.67")
    assert T_mu(F("1.67"))+F(".02") < F(".8")
    print("all hard calibration gates passed, including |E|<=1e-5")
    # Two in Haar: R=A+B cos, using exact extrema and one-defect values.
    R1 = F(25,16)
    A = R1**2+(R1-1)**2
    B = 2*R1*(R1-1)
    assert A == F(353,128) and B == F(225,128)
    assert 1/(A+B) == F(64,289) and 1/(A-B) == 1
    print("one- and two-defect Haar controls passed")


def wider_energy_gate():
    N = 400
    e = F(1,2000)
    c_hi,c_lo=F(3,5)+e/2, F(3,5)-e/2
    assert max(abs(2*c*c-1) for c in (c_hi,c_lo)) < F(".29")
    assert (1+e)*(2+F(5,4)*e)/(1-e) < F("2.01")
    assert (2-F(5,4)*e)*(1-e)/(1+e) > F("1.99")
    assert F("2.01")**2+F("2.01")**(-2) < F("4.4")
    assert (F("4.4")/2-1)/2 < 1
    assert F("2.01")+F("2.01")**2*e/2 < 3
    assert F(".29")+3*e < F(".31")
    assert 1+9/(1-F(".31")) < 15
    cap = (15*N*e*e)/(2*(1-15*N*e*e))
    false_eight = 8*cap/(1/F(".82")-1)
    unaffected_min = F("1.99")/F("1.08")**2
    assert unaffected_min > F("1.67")
    assert T_mu(unaffected_min)+F(".02") < F(".8")
    bad_segments = cap/(1/T_mu(F("1.08"))-1)
    p1 = (1-F(1,N))**(N-1)
    lower = p1*(1-bad_segments)
    assert lower > F(".32")
    assert false_eight < F(".028")
    assert 1-(1-lower)**8 > F(".95")
    print(f"wider-energy eight-chain bounds: false alarm <= {float(false_eight):.9f}, "
          f"power >= {float(1-(1-lower)**8):.9f} (rounded)")


def commensurate_control():
    # At B=1, k=2*pi/3. Use the rational metric to avoid sqrt(3).
    Q = (1,F(1,2),F(1,2),1)
    Qi = (F(4,3),F(-2,3),F(-2,3),F(4,3))
    D = multiply(transfer(1),transfer(0))
    R = power(transfer(1),2)
    assert power(R,3) == I
    def bare_transmission(M):
        Mt = (M[0],M[2],M[1],M[3])
        W = multiply(multiply(multiply(Q,M),Qi),Mt)
        return F(4)/(2+W[0]+W[3])
    assert bare_transmission(D) == F(3,4)
    t2 = [bare_transmission(multiply(multiply(D,power(R,j)),D)) for j in range(3)]
    assert t2 == [F(3,7),F(1),F(3,7)]
    print("commensurate control: two-defect transmission = 3/7, 1, 3/7 by gap residue")


if __name__ == "__main__":
    finite_probability()
    robustness()
    wider_energy_gate()
    commensurate_control()
