"""Pure first-jet and actual three-register parity controls; no numerical optimizer."""
import sympy as S

R = S.Rational
I = S.I
eye = S.eye
col0 = S.Matrix([1, 0])
col1 = S.Matrix([0, 1])
simp = S.simplify


def eq(a, b):
    if isinstance(a, S.MatrixBase):
        assert (a-b).applyfunc(simp) == S.zeros(*a.shape), (a, b)
    else:
        assert simp(a-b) == 0, (a, b)


def norm(v):
    return simp((v.H*v)[0])


def symbolic():
    s, x, y, z, w, d = S.symbols("s x y z w d", real=True)
    ell, q = 1-s, d*d
    W, L = 1-ell*x, ell*x
    j0, j1 = I*(s*x-ell*d*(z+I*w)), I*ell*d*(z+I*w)
    A0, A1 = s*x+2*s*d*z+q*(1-ell*y), ell*q*y
    q0 = A0/W-j0*S.conjugate(j0)/W**2
    q1 = A1/L-j1*S.conjugate(j1)/L**2
    T = q+(s*x*(1-x)+2*s*d*z-ell*q*z*z/x)/W
    classical = (ell*d*w)**2*(1/W+1/L)
    eq(2*S.re(j0), 2*ell*d*w)
    eq(2*S.re(j1), -2*ell*d*w)
    eq(T, W*q0+L*q1+classical)
    eq(q0.subs(w, 0)-q0, ell**2*q*w**2/W**2)
    eq(q1.subs(w, 0)-q1, q*w**2/x**2)
    eq(q1.subs({w: 0, y: z*z/x}), 0)
    eq(q1.subs({w: 0, y: 1-z*z/(1-x)}),
       q*(x*(1-x)-z*z)/(x*x*(1-x)))
    eq(T-q, s*x/ell-ell*q*(z-s*x/(ell*d))**2/(x*W))
    e = S.sqrt(x*(1-x))
    eq(q0.subs({w: 0, z: e, y: 1-x}), s*(d+e)**2/W**2)
    eq(q1.subs({w: 0, z: e, y: 1-x}), 0)
    print("symbolic: Bellman mean, CFI, real reduction, capacity, square, saturation")


def isometry(x, y, c):
    """Exact Cholesky isometry S-M-A, including the coherent basis tag."""
    c2 = simp(c*S.conjugate(c))
    assert 0 <= x <= 1 and 0 <= y <= 1
    assert c2 <= x*y and c2 <= (1-x)*(1-y)
    if x in (0, 1):
        assert c == 0
        L1 = S.diag(S.sqrt(x), S.sqrt(y))
        L0 = S.diag(S.sqrt(1-x), S.sqrt(1-y))
    else:
        L1 = S.Matrix([[S.sqrt(x), c/S.sqrt(x)], [0, S.sqrt(y-c2/x)]])
        L0 = S.Matrix([[S.sqrt(1-x), -c/S.sqrt(1-x)],
                       [0, S.sqrt(1-y-c2/(1-x))]])
    M = S.Matrix([[x, c], [S.conjugate(c), y]])
    eq(L1.H*L1, M)
    eq(L0.H*L0, eye(2)-M)
    V = S.zeros(8, 2)
    for a, A in enumerate((L0,L1)):
        for b in range(2):
            V[4*a+2*b+a,:] = A[b,:]
    eq(V.H*V, eye(2))
    return V


def branches(s, d, x, y, c):
    V = isometry(x, y, c)
    B = [S.diag(1, S.sqrt(s)), S.Matrix([[0,S.sqrt(1-s)],[0,0]])]
    P = S.kronecker_product(S.diag(0,1), eye(4))
    out = []
    for j, Bj in enumerate(B):
        K = S.kronecker_product(Bj, eye(4))
        v = (K*V*col0).applyfunc(simp)
        dv = (I*(P*v+d*K*V*col1)).applyfunc(simp)
        # Tag syndrome is observable without access to a Kraus environment.
        for r in range(8):
            if (((r >> 2) & 1) ^ (r & 1)) != j:
                eq(v[r], 0)
                eq(dv[r], 0)
        p, inner = norm(v), simp((v.H*dv)[0])
        dp = simp(2*S.re(inner))
        if p:
            jet = (dv-v*inner/p)/S.sqrt(p)
            eq((v.H*jet)[0], 0)
            out.append((p, dp, norm(jet), v/S.sqrt(p), jet))
        else:
            eq(p, 0)
            eq(dv*v.H+v*dv.H, S.zeros(8))
            out.append((p, dp, None, v, dv))
    eq(out[0][0]+out[1][0], 1)
    eq(out[0][1]+out[1][1], 0)
    if 0 < x < 1 and 0 < s < 1:
        W = 1-(1-s)*x
        T = d*d+(s*x*(1-x)+2*s*d*S.re(c)-(1-s)*d*d*S.re(c)**2/x)/W
        quantum = sum(p*q for p, dp, q, _, _ in out)
        classical = sum(dp*dp/(4*p) for p, dp, _, _, _ in out)
        eq(quantum+classical, T)
    return out


def aligned_recovery(out):
    """Construct the recovery Kraus restrictions on occupied nominal tangent spans."""
    rho = S.zeros(2)
    drho = S.zeros(2)
    for p, dp, q, psi, jet in out:
        if not p:
            continue
        K = col0*psi.H
        if q:
            K += I*col1*jet.H/S.sqrt(q)
        eq(K*psi, col0)
        eq(K*jet, I*S.sqrt(q)*col1)
        eq(K*K.H, eye(2) if q else col0*col0.H)
        v, dv = K*psi, K*jet
        rho += p*v*v.H
        drho += dp*v*v.H+p*(dv*v.H+v*dv.H)
    eq(rho, col0*col0.H)
    return drho.applyfunc(simp)


def pure_tangent(d):
    return I*d*(col1*col0.H-col0*col1.H)


def witness():
    s = R(1,4)
    x, y, c = R(3,5), R(27,32), R(1,4)
    weights, weighted_roots, T = [], [], 0
    target = [(R(11,20), R(8407,15488)), (R(9,20), R(355,1152))]
    lower = [R(736,1000), R(555,1000)]
    for j, (p, dp, q, _, _) in enumerate(branches(s, R(1,2), x, y, c)):
        eq(dp, 0)
        assert (p,q) == target[j]
        assert lower[j]**2 < q
        T += p*q
        weights.append(p)
        weighted_roots.append(p*S.sqrt(q))
    eq(T, R(1539,3520))
    value = 4*(T+sum(weighted_roots)+R(1,4))
    lb = 4*(T+sum(a*b for a,b in zip(weights,lower))+R(1,4))
    ub = 4*(R(1,4)+R(1,3)+R(1,2)+R(1,4))
    eq(lb, R(590377,110000))
    eq(ub, R(16,3))
    eq(lb-ub, R(11131,330000))
    assert lb-ub > R(337,10000)
    print("physical witness (1,1/4,1): saturated <= %s; actual > %s; gap > %s"
          % (ub, lb, lb-ub))
    print("witness value (display only):", value.evalf(14))


def controls():
    s, d = R(1,4), R(1,2)
    out = branches(s, d, R(1,2), R(1,2), R(1,4)+I/8)
    real = branches(s, d, R(1,2), R(1,2), R(1,4))
    assert out[0][1] > 0 and out[1][1] < 0
    assert all(a[2] < b[2] for a,b in zip(out, real))
    missing = branches(s, d, 0, 1, 0)
    assert missing[1][0] == 0 and norm(missing[1][4]) > 0
    eq(missing[0][2], s*d*d)
    eq(branches(s, d, 0, 0, 0)[0][2], d*d)
    for d in (0,R(1,2),S.sqrt(R(355,1152))):
        eq(branches(1, d, R(1,2), R(1,2), R(1,2))[0][2],
           (d+R(1,2))**2)
    print("physical controls: complex score, null-probability jet, parking, unitary end")


def schedule():
    s = R(1,4)
    seed = branches(s, 0, R(1,2), R(1,2), R(1,2))
    eq(aligned_recovery(seed), pure_tangent(S.sqrt(s)/2))
    # A warmup slot using the seed. No measurement history between corrections.
    q, d, x, z = s/4, S.sqrt(s)/2, R(1,4), 0
    W = 1-(1-s)*x
    T = q+s*x*(1-x)/W
    y = x*T/q
    warm = branches(s, d, x, y, z)
    eq(warm[0][2], T)
    eq(warm[1][2], T)
    eq(aligned_recovery(warm), pure_tangent(S.sqrt(T)))
    # Rational high-information example: saturates completed square but not z bound.
    kappa = s/(1-s)
    q, d = 100, 10
    x = 1-2*kappa/d
    z = kappa*x/d
    T = q+kappa*x
    y = x*T/q+z*z/x
    high = branches(s, d, x, y, z)
    eq(high[0][2], T)
    eq(high[1][2], T)
    eq(aligned_recovery(high), pure_tangent(S.sqrt(T)))
    print("no-record controls: seed q=%s, next warmup q=%s; high-q increment=%s"
          % (s/4, warm[0][2], kappa*x))


if __name__ == "__main__":
    symbolic()
    witness()
    controls()
    schedule()
    print("PASS")
