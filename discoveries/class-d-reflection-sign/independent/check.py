#!/usr/bin/env python3
"""Independent exact class-D checks; SymPy 1.14."""
from itertools import combinations
from random import Random
import sympy as sp

R = sp.Rational
I = sp.I
E = sp.eye(2)
Z = sp.zeros(2)
SWAP = sp.Matrix([[0, 1], [1, 0]])
J = sp.Matrix([[0, 1], [-1, 0]])
X = sp.Matrix([[1, 1], [1, -1]])/sp.sqrt(2)
Y = sp.Matrix([[1, -I], [1, I]])/sp.sqrt(2)
axes = (E, X, Y)

def nambu(A, B):
    return A.row_join(B.conjugate()).col_join(B.row_join(A.conjugate()))

def two_sector(K, L):
    return nambu((K+L)/2, (K-L)/2)

def sq(x):
    return sp.simplify(x * x.conjugate())

def clean(M):
    return M.applyfunc(sp.simplify)

def wedge_law(C, loc):
    """Given orthonormal occupied free modes, count the exact local subset."""
    law = {(): 0}
    for ix in combinations(range(C.rows), 2):
        y = tuple(i for i in ix if i in loc)
        law[y] = law.get(y, 0) + sq(C[list(ix), :].det())
    law = {y: sp.simplify(p) for y,p in law.items()}
    assert sum(law.values()) == 1 and all(p>=0 for p in law.values())
    return law

def local_law(C): # algebraically integrate out arbitrary environment
    P = {ij:sq(C[list(ij),:].det()) for ij in combinations(range(4), 2)}
    marg = [sum(sq(x) for x in C[j,:]) for j in range(4)]
    f = P | {(j,):sp.simplify(marg[j] - sum(p for ij,p in P.items() if j in ij))
             for j in range(4)}
    f[()] = sp.simplify(1-sum(f.values()))
    assert all(p>=0 for p in f.values()) and sum(f.values()) == 1
    return f

def score(law, loc):
    e0,e1,h0,h1 = loc
    return sp.simplify(sum(p for y,p in law.items() if len(y)==2)
                         -3*(law[e0,h1]+law[e1,h0]))

def from_local(r, Q):
    return clean(sp.diag(Q,Q.conjugate()) * r[:,:2])

def weidenmuller_assert(S):
    K = (S-sp.eye(S.rows)).inv() * (S+sp.eye(S.rows))
    assert K+K.T == sp.zeros(S.rows)
    # H=iK, gamma=1, W=I, in the W*W convention.
    H = I*K
    assert H.H == H and H.conjugate() == -H
    S0 = sp.eye(S.rows) - 2*I*(-H + I*sp.eye(S.rows)).inv()
    assert clean(S0) == S

def fixed_probability_ambiguity():
    a,b = R(4,5), R(3,5)
    # Bottom rows phase-flipped relative to involutive construction, to
    # realize by full-rank *ordinary* charge-preserving lead-dot coupling.
    Q = sp.diag(-1,1)*SWAP
    T = sp.Matrix([[a,b],[-b,a]])
    K = sp.kronecker_product(T,E)
    Ns = [sp.kronecker_product(T,P) for P in (SWAP,Q)]
    Smaj = [sp.diag(K,N) for N in Ns]
    S = [two_sector(K,N) for N in Ns]
    assert S[0].applyfunc(sq) == S[1].applyfunc(sq)
    loc = (0,1,4,5)
    for k,d in enumerate((-a**4,a**4)):
        assert Smaj[k].T*Smaj[k] == sp.eye(8) and Smaj[k].det()==1
        weidenmuller_assert(Smaj[k])
        q = S[k]
        assert q.H*q == sp.eye(8) and q.det()==1
        assert sum(sq(q[i,0]) for i in (2,3,6,7)) == b*b
        r = q[list(loc),list(loc)]
        assert r.det() == d
        all_laws = []
        for U in axes:
            V = sp.diag(U,E,U.conjugate(),E)
            all_laws.append(wedge_law(V*q[:,:2],loc))
        scores = [score(L,loc) for L in all_laws]
        assert sum(scores)/3 == d
        assert all(sum(v for y,v in L.items() if len(y)==2) == a**4 for L in all_laws)
        print("same fixed Nambu |S_ij|² and leakage", b*b, "| det", d, "| x/y/z scores",
              [scores[1],scores[2],scores[0]])

def two_outputs_fail():
    # Any two axes lie in one plane; choose the orthogonal direction z
    # and take x/y as two exact representative in-plane examples.
    a,b = R(2,5),R(7,12)
    assert a+b == R(59,60)
    Ms = [sp.Matrix([[a*b,0],[0,0]]),sp.Matrix([[0,a*b],[0,0]])]
    rs = [nambu(a*E,(-J*M/a).T) for M in Ms]
    ds = tuple(r.det() for r in rs)
    assert ds == (R(16,625),-R(649,22500))
    for V in (X,Y):
        assert local_law(from_local(rs[0],V)) == local_law(from_local(rs[1],V))
        for pol in (sp.Matrix([1,0]),sp.Matrix([0,1]),sp.Matrix([3,4*I])/5):
            for inp in [pol.col_join(sp.zeros(2,1)),sp.zeros(2,1).col_join(pol)]:
                outputs = [(sp.diag(V,V.conjugate())*r*inp).applyfunc(sq) for r in rs]
                assert outputs[0] == outputs[1]
    assert local_law(from_local(rs[0],E)) != local_law(from_local(rs[1],E))
    for r in rs:
        assert sum(score(local_law(from_local(r,V)),(0,1,2,3)) for V in axes)/3 == r.det()
    print("two exact local laws (both mixed bases) equal for opposite determinants",*ds)

def quasi_controls():
    rng = Random(23858)
    S = [nambu(E,Z),nambu(sp.diag(1,0),sp.diag(0,1)),
         nambu(sp.diag(R(1,2),R(1,3)),sp.diag(R(1,2),R(1,4)))]
    for k in range(15):
        matrices = [sp.Matrix(2,2,[R(rng.randint(-6,6),40) + I*R(rng.randint(-6,6),40)
                                   for _ in range(4)]) for j in range(2)]
        assert sum(sq(x) for M in matrices for x in M) < R(1,2)
        S.append(nambu(*matrices))
    assert [r.det() for r in S[:3]] == [1,-1,0]
    for r in S:
        laws = [local_law(from_local(r,U)) for U in axes]
        assert sum(score(L,(0,1,2,3)) for L in laws)/3 == r.det()
        assert abs(r.det()) <= sum(v for y,v in laws[0].items() if len(y)==2)
    print("isolated +/- unit norm, nontrivial singular, and 15 complex rational contraction controls: PASS")

if __name__ == "__main__":
    fixed_probability_ambiguity()
    two_outputs_fail()
    quasi_controls()
