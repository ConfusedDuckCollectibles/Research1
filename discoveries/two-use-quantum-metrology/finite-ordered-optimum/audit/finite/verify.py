"""Independent exact controls for local tangent compression; Python + SymPy."""
from sympy import I, Matrix, MatrixBase, Rational as R, conjugate, eye, im, kronecker_product
from sympy import re, simplify as S, sqrt, trace, zeros

def eq(a,b):
    if isinstance(a, MatrixBase):
        return all(S(t)==0 for t in a-b)
    return S(a-b)==0

P=kronecker_product(Matrix.diag(0,1),eye(4))
P_even=Matrix.diag(1,1,0,0,0,0,1,1)
P_odd=eye(8)-P_even
e0,e1=Matrix([1,0]),Matrix([0,1])

def physical_jet(R0,R1,s,d):
    """Explicit sensor⊗tag⊗logical, then accessible parity dephasing."""
    W=zeros(8,2)
    W[0:2,:]=R0
    W[6:8,:]=R1
    assert eq(W.H*W,eye(2))
    B0=kronecker_product(Matrix.diag(1,sqrt(s)),eye(4))
    B1=kronecker_product(Matrix([[0,sqrt(1-s)],[0,0]]),eye(4))
    v=[B*W*e0 for B in (B0,B1)]
    t=[I*(P*b+d*B*W*e1) for B,b in zip((B0,B1),v)]
    rho=sum((b*b.H for b in v),zeros(8))
    drho=sum((c*b.H+b*c.H for b,c in zip(v,t)),zeros(8))
    for j,tag in enumerate((P_even,P_odd)):
        assert eq(tag*rho*tag,v[j]*v[j].H)
        assert eq(tag*drho*tag,t[j]*v[j].H+v[j]*t[j].H)
    # The physical state is already parity-block diagonal at the local point,
    # including its derivative (not assumed across all nearby phases).
    assert eq(P_even*drho*P_odd,zeros(8))
    return list(zip(v,t)),R1.H*R1

def mass_jet(v,t):
    q=S((v.H*v)[0])
    qprime=S(2*re((v.H*t)[0]))
    if q==0:
        return q,qprime,None
    n=S((t.H*t)[0])
    b=(v.H*t)[0]
    return q,qprime,S(n/q-b*conjugate(b)/q**2)

def exact_sld(v,t):
    # Strip zero padding from the individually known parity block.
    indices=[j for j in range(len(v)) if v[j]!=0 or t[j]!=0]
    if not indices:
        return 0
    v,t=v.extract(indices,[0]),t.extract(indices,[0])
    rho=v*v.H
    dot=t*v.H+v*t.H
    q=S(trace(rho))
    if not q:
        assert eq(dot,zeros(len(v)))
        return 0
    proj=rho/q
    L=(2*dot-proj*dot*proj)/q
    assert eq(L,L.H)
    assert eq((rho*L+L*rho)/2,dot)
    return S(trace(rho*L**2))

def compression_jet(M,s,d):
    x,y,m=M[0,0],M[1,1],M[0,1]
    l=1-s
    q=(1-l*x,l*x)
    b=(I*(s*x-l*d*m),I*l*d*m)
    n=(s*x+d*d*(1-l*y)+2*s*d*re(m),l*d*d*y)
    return [(S(qj),S(2*re(bj)),S(nj/qj-bj*conjugate(bj)/qj**2) if qj else None)
            for qj,bj,nj in zip(q,b,n)]

def local_terminal(tab):
    return S(sum(qp**2/q+4*q*d2 for q,qp,d2 in tab if q))

def feasible(M):
    return M[0,0]>=0 and M[1,1]>=0 and M[0,0]<=1 and M[1,1]<=1 and S(M.det())>=0 and S((eye(2)-M).det())>=0

# Arbitrary strict effect with both real and imaginary off-diagonal parts.
rot=Matrix([[R(3,5),-R(4,5)],[R(4,5),R(3,5)]])
U=Matrix.diag(1,(3+4*I)/5)*rot
A=U*Matrix.diag(R(4,5),R(3,5))*U.H
B=U*Matrix.diag(R(3,5),R(4,5))*U.H
A,B=A.applyfunc(S),B.applyfunc(S)
s,d=R(16,25),R(3,5)
first,M=physical_jet(A,B,s,d)
assert eq(M[0,1],R(-252,3125)+I*R(336,3125))
tab=[mass_jet(*j) for j in first]
assert tab==compression_jet(M,s,d)
assert feasible(M) and all(0<q<1 for q,_,_ in tab)
F1=S(sum(exact_sld(*j) for j in first))
CFI=S(sum(qp**2/q for q,qp,_ in tab))
assert CFI>0 and F1==local_terminal(tab)
Mr=M.applyfunc(re)
tr=compression_jet(Mr,s,d)
assert feasible(Mr)
assert local_terminal(tab)==local_terminal(tr)
assert all(r[2]>t[2] for t,r in zip(tab,tr))
Mp=Mr.copy()
Mp[0,1]=Mp[1,0]=abs(Mp[0,1])
tp=compression_jet(Mp,s,d)
assert feasible(Mp) and tp[0][2]>tr[0][2] and tp[1][2]==tr[1][2]
print("Complex effect direct SLD:",F1,"positive flag CFI:",CFI)
print("Realification preserves classical + base QFI and raises both tangents; sign flip raises surviving tangent")

# Two uses, with different physical effects after the two observed outcomes.
# Multiplying explicitly the branch amplitude and its derivative tests the flag chain rule.
Tplus=Matrix([[R(1,2),R(1,2)],[R(1,2),R(1,2)]])
Tminus=eye(2)-Tplus
total=0
rhs=CFI
for j,(q,qp,d2) in enumerate(tab):
    # Nominal isometry sends the horizontal pair to the next code. For j=0
    # choose a scalar effect, and for j=1 choose a rank-one tilted effect.
    r0,r1=(R(4,5)*eye(2),R(3,5)*eye(2)) if j==0 else (Tminus,Tplus)
    u=R(144,169)
    leaves,Mn=physical_jet(r0,r1,u,sqrt(d2))
    sub=compression_jet(Mn,u,sqrt(d2))
    assert [mass_jet(*leaf) for leaf in leaves]==sub
    rhs+=q*local_terminal(sub)
    for v,t in leaves:
        total+=exact_sld(sqrt(q)*v,sqrt(q)*t+qp*v/(2*sqrt(q)))
assert eq(total,rhs)
print("Two-call physical block SLD agrees with recursive classical chain rule; total =",S(total))

# Rank-changing boundary: present density derivative in no-jump block;
# strictly absent density derivative in jump block, which has nonzero ket derivative.
edge, E=physical_jet(Matrix.diag(1,0),Matrix.diag(0,1),s,d)
v1,t1=edge[1]
assert eq(v1,zeros(8,1)) and not eq(t1,zeros(8,1))
assert exact_sld(v1,t1)==0
qfi_edge=sum(exact_sld(*leaf) for leaf in edge)
assert eq(qfi_edge,4*s*d**2)
# Positive x recovers prior tangent norm through the small jump branch in the limit.
X=R(1,100)
near=physical_jet(Matrix.diag(sqrt(1-X),0),Matrix.diag(sqrt(X),1),s,d)[0]
q,qp,d2=mass_jet(*near[1])
assert qp==0 and q*d2==(1-s)*d*d
assert qfi_edge<4*d*d
print("Rank-change control: zero-mass pointwise SLD",qfi_edge,"while positive-mass limit",4*d*d)
# Zero and unit survival share the same projector dilation and omit zero mass without dividing.
for survival in (0,1):
    ends, Z=physical_jet(A,B,survival,d)
    assert [mass_jet(*j) for j in ends]==compression_jet(Z,survival,d)
    assert S(sum(exact_sld(*j) for j in ends))==local_terminal(compression_jet(Z,survival,d))
print("Endpoint survivals 0 and 1 agree with physical direct SLD")
# Multi-slot optimizing control for all-noiseless schedule including a dummy zero slot:
# the centered-generator bound is reached by a half-excited aligned projection
# on every unit-survival call and by M=0 on dummy slots.
D=d
slots=(1,1,0,1,0,1)
for survival in slots:
    R0,R1=(Tminus,Tplus) if survival else (eye(2),zeros(2))
    jet,E=physical_jet(R0,R1,survival,D)
    q,qp,D2=mass_jet(*jet[0])
    assert q==1 and qp==0 and mass_jet(*jet[1])==(0,0,None)
    assert D2==(D+R(survival,2))**2
    assert exact_sld(*jet[0])==4*D2
    D+=R(survival,2)
assert 4*D*D==(2*d+sum(slots))**2
print("Six-slot noiseless/dummy saturation of the centered-generator bound:",4*D*D)
