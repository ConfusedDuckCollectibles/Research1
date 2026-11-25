"""Gaussian-rational, standard-library replay of the class-D count witnesses."""
from fractions import Fraction as F
from itertools import combinations, permutations
from typing import NamedTuple
from random import Random


class C(NamedTuple):
    re: F = F(0)
    im: F = F(0)

    def __add__(self, y):
        y = as_c(y)
        return C(self.re+y.re, self.im+y.im)

    __radd__ = __add__

    def __neg__(self):
        return C(-self.re, -self.im)

    def __sub__(self, y):
        return self + (-as_c(y))

    def __rsub__(self, y):
        return as_c(y) - self

    def __mul__(self, y):
        y = as_c(y)
        a, b = self
        c, d = y
        return C(a*c-b*d, a*d+b*c)

    __rmul__ = __mul__

    def __truediv__(self, z):
        return C(self.re/F(z), self.im/F(z))

    def conj(self):
        return C(self.re, -self.im)

    def norm2(self):
        return self.re*self.re + self.im*self.im


def as_c(z):
    return z if isinstance(z, C) else C(F(z))


U = C(0, 1)
I2 = [[C(1), C(0)], [C(0), C(1)]]
O2 = [[C(0), C(0)], [C(0), C(0)]]
P = [[C(0), C(1)], [C(1), C(0)]]
SKEW = [[C(0), C(1)], [C(-1), C(0)]]
BASES = [
    ([[C(1), C(1)], [C(1), C(-1)]], 2),
    ([[C(1), -U], [C(1), U]], 2),
    (I2, 1),
]


def adj(A):
    return [list(j) for j in zip(*[[a.conj() for a in row] for row in A])]


def transpose(A):
    return [list(j) for j in zip(*A)]


def mult(A, B):
    return [[sum(a*b for a,b in zip(row,col)) for col in transpose(B)] for row in A]


def scale(A, z):
    return [[a*z for a in row] for row in A]


def diff(A, B):
    return [[a-b for a,b in zip(r,s)] for r,s in zip(A,B)]


def conjugate(A):
    return [[v.conj() for v in r] for r in A]


def join(A, B, C0, D):
    return [a+b for a,b in zip(A,B)] + [c+d for c,d in zip(C0,D)]


def nambu(A, B):
    return join(A, conjugate(B), B, conjugate(A))


def det(A):
    n = len(A)
    answer = C()
    for p in permutations(range(n)):
        sign = (-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1,n))
        product = C(sign)
        for i in range(n):
            product *= A[i][p[i]]
        answer += product
    return answer


def eye(n):
    return [[C(i == j) for j in range(n)] for i in range(n)]


def counts(A, B, basis, norm):
    assert mult(basis, adj(basis)) == scale(I2, norm)
    incident = mult(basis, A) + mult(conjugate(basis), B)
    pairs = {}
    for i,j in combinations(range(4),2):
        pairs[(i,j)] = det([incident[i], incident[j]]).norm2() / norm**2
    alone = {
        (i,): sum(a.norm2() for a in incident[i])/norm
              -sum(v for pair,v in pairs.items() if i in pair)
        for i in range(4)
    }
    full = pairs | alone
    full[()] = 1 - sum(full.values())
    assert all(0 <= w <= 1 for w in full.values()) and sum(full.values())==1
    score = sum(pairs.values())-3*(pairs[(0,3)]+pairs[(1,2)])
    # Optional neutral-pair-blind net-charge reconstruction: marginal
    # one-body returned number from the two separate single-input controls.
    charge = lambda outcome: sum((1 if i<2 else -1) for i in outcome)
    p1 = sum(mass for pattern,mass in full.items() if abs(charge(pattern))==1)
    mu = sum(a.norm2() for row in incident for a in row)/norm
    assert (mu-p1)/2 == sum(pairs.values())
    singles_by_input = [[a.norm2()/norm for a in row] for row in incident]
    return score, full, singles_by_input


def score_check(A, B):
    expected = sum(counts(A,B,*v)[0] for v in BASES)/3
    assert det(nambu(A,B)) == C(expected)
    return expected


def global_pair():
    c,s=F(4,5),F(3,5)
    K = join(scale(I2,c), scale(I2,s), scale(I2,-s), scale(I2,c))
    J0=join(scale(P,c), scale(P,s), scale(P,-s), scale(P,c))
    J1=[[(1 if i%2 else -1)*v for v in row] for i,row in enumerate(J0)]
    all_intensities=[]
    all_signs=[]
    for J in [J0,J1]:
        assert mult(K, adj(K)) == eye(4) and mult(J, adj(J)) == eye(4)
        assert det(K)*det(J) == C(1)
        assert det(diff(K,eye(4))).norm2() > 0 and det(diff(J,eye(4))).norm2() > 0
        A=scale(diff(K,scale(J,-1)), F(1,2))
        B=scale(diff(K,J), F(1,2))
        S=nambu(A,B)
        assert mult(S,adj(S)) == eye(8)
        all_intensities.append([[v.norm2() for v in row] for row in S])
        assert sum(S[i][0].norm2() for i in [2,3,6,7]) == F(9,25)
        all_signs.append(score_check([a[:2] for a in A[:2]],[b[:2] for b in B[:2]]))
    assert all_intensities[0] == all_intensities[1]
    assert all_signs == [-F(256,625),F(256,625)]
    print("full two-lead |S_ij|^2 match; transmission 9/25; det r = -256/625 and +256/625")


def lower_pair():
    a,b=F(2,5),F(7,12)
    assert a+b == F(59,60) < 1
    A=scale(I2,a)
    v_plus, v_minus = [[C(1),U]], [[C(1),-U]]
    all_counts=[]
    out=[]
    for v in [v_plus,v_minus]:
        mixed = scale(mult(transpose(v_plus),conjugate(v)), a*b/2)
        B=transpose(scale(mult(SKEW,mixed),-1/a))
        K=mult(B,adj(B))
        assert K[0][0]+K[1][1] == C(b*b) and det(K)==C()
        assert mult(scale(SKEW,a),transpose(B))==mixed
        all_counts.append([counts(A,B,*BASES[index])[1:] for index in (0,2)])
        out.append(score_check(A,B))
    assert all_counts[0] == all_counts[1]
    assert out == [F(16,625), -F(649,22500)]
    print("two-basis full local law and individual inputs match; det r = +16/625 and -649/22500")


def controls():
    assert score_check(I2,O2) == 1                # two normal reflections
    assert score_check(O2,I2) == 1                # two Andreev reflections
    A = [[C(1),C(0)],[C(0),C(0)]]
    B = diff(I2,A)
    assert score_check(A,B) == -1                # just one Andreev mode
    assert score_check(A,O2) == 0                # one fully transmitted orbital
    assert score_check(O2,O2) == 0                # all transmitted
    # Fixed exact nonspecial contractions cover phases, off-diagonal A/B and leakage.
    generator=Random(329830)
    for _ in range(12):
        A,B = ([[C(F(generator.randint(-5,5),100),F(generator.randint(-5,5),100))
                  for _ in range(2)] for _ in range(2)] for _ in range(2))
        assert sum(x.norm2() for row in nambu(A,B) for x in row)<1
        score_check(A,B)
    # Hoeffding at n=3000 and alpha=.05: ln 40 < 15/4 by lower truncated exp sum.
    def factorial(j):
        v=1
        for i in range(1,j+1):
            v*=i
        return v
    assert sum(F(15,4)**j/factorial(j) for j in range(9))>40
    assert F(9,6000)*F(15,4) == F(".075")**2
    assert 6*F(".01")+3*F(".01")+3*F(".01") == F(".12")
    assert 2*F(".195") < F(256,625)
    print("normal/Andreev/singular + 12 complex controls; neutral-pair workaround; statistical margin <0.075")


if __name__ == "__main__":
    global_pair()
    lower_pair()
    controls()
    print("all exact Gaussian-rational checks passed")
