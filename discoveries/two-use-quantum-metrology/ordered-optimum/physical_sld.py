#!/usr/bin/env python3
"""Direct exact qubit/tag Kraus, phase derivative, and SLD checks."""
import sympy as S
Q=S.Rational
I=S.I
tup = [(0,0),(0,1),(1,0),(1,1)]
h = S.diag(0,0,1,1)  # sensor excitation
parity = S.diag(0,1,1,0)
col = lambda a,b:S.Matrix([a,0,0,b])
def simp(M): return M.applyfunc(S.simplify)
def az(M):return simp(M)==S.zeros(*M.shape)

def kraus(s):
    return [S.kronecker_product(K,S.eye(2)) for K in (
      S.diag(1,S.sqrt(s)), S.Matrix([[0,S.sqrt(1-s)],[0,0]]))]

def sectors(s,psi,gamma):
    """Physical tagged branch outputs at theta=0; v′=iPv+Bγ."""
    V=[]
    for j,K in enumerate(kraus(s)):
        v=K*psi; w=I*h*v+K*gamma
        assert az(parity*v-j*v) and az(parity*w-j*w)
        V.append((simp(v),simp(w)))
    return V

def sld(inputs):
    rho=S.zeros(4);dr=S.zeros(4)
    for v,w in inputs:
        rho+=v*v.H
        dr+=w*v.H+v*w.H
    rho=simp(rho);dr=simp(dr)
    L=S.zeros(4);cdf=0
    for inds in ([0,3],[1,2]):
        R=rho.extract(inds,inds);D=dr.extract(inds,inds)
        mass=S.trace(R);mass_d=S.trace(D)
        if mass==0:
            assert az(D)
            continue
        assert mass.is_positive
        assert az(R*R-mass*R)
        # On a pure subnormalized block, radial term requires a correction
        # to the common constant-weight formula L=2D/mass.
        X=2*D/mass-mass_d*R/mass**2
        cdf+=mass_d**2/mass
        for i,k in enumerate(inds):
            for j,l in enumerate(inds): L[k,l]=X[i,j]
    assert az((rho*L+L*rho)/2-dr)
    return S.simplify(S.trace(rho*L*L)),S.simplify(cdf)

def flag_checks():
    s=Q(81,100);d=Q(3,10);X=Q(1,2)
    psi=col(1/S.sqrt(2),1/S.sqrt(2))
    phi=col(-I/S.sqrt(2),I/S.sqrt(2))
    value,cfi=sld(sectors(s,psi,I*d*phi))
    W=1-(1-s)*X
    assert value==4*(d*d+s*X*(1-X)/W)
    assert cfi>0 # demonstrates nonzero classical flag score in admissible controls
    print(f"Direct complex-overlap control: F={value}, classical flag term={cfi}>0.")
    d=Q(1,3)
    v0=sectors(s,col(1,0),I*d*col(0,1))
    assert az(v0[1][0]) and not az(v0[1][1])
    assert sld(v0)==(4*s*d*d,0)
    assert sld(sectors(s,col(0,1),I*d*col(1,0)))==(4*d*d,0)
    print(f"Zero-mass jump at x=0 verified: SLD F={4*s*d*d}; x=1 gives {4*d*d}.")

def physical_witness():
    a=Q(16,25);s=Q(81,100)
    x=Q(1,2);t=Q(9,10);y=Q(10,19)
    first=sectors(a,col(S.sqrt(1-x),S.sqrt(x)),S.zeros(4,1))
    qfi=0;prob=0
    for j,(v,w) in enumerate(first):
        W=S.simplify((v.H*v)[0]);prob+=W
        if j==1:
           assert az(w)
           v2=S.sqrt(W)*col(S.sqrt(1-y),S.sqrt(y))
           d2=S.zeros(4,1)
        else:
           alpha=v/S.sqrt(W)
           beta=col(-alpha[3],alpha[0])
           psi=col(t,1)/S.sqrt(1+t*t)
           phi=col(-1,t)/S.sqrt(1+t*t)
           # Isometry acts on the accessible no-jump code, and may be
           # extended arbitrarily on its complement.
           C=psi*alpha.H+phi*beta.H
           assert az(C.H*C-S.diag(1,0,0,1))
           v2=C*v; d2=C*w
        v2,d2=simp(v2),simp(d2)
        f,c=sld(sectors(s,v2,d2))
        assert c==0
        qfi+=f
    assert prob==1 and qfi==Q(203724139,66974525)
    print(f"Physical two-call parity-conditioned circuit: F={qfi}, total probability exactly one.")

if __name__=="__main__":
    flag_checks();physical_witness()
