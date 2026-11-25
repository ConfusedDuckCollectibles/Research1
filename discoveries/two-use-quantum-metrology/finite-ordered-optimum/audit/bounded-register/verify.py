"""Exact finite controls, independent of the Bellman optimiser."""
from fractions import Fraction as F
from sympy import Matrix as Mat, Rational as R, sqrt, simplify, expand, eye, zeros
from sympy import kronecker_product as kron

def S(x):
    return simplify(expand(x))
def clean(m): return m.applyfunc(S)
def eq(a,b): return all(S(x)==0 for x in a-b)

Q=kron(Mat.diag(0,1),eye(4))
pe=Mat.diag(1,1,0,0,0,0,1,1)
po=eye(8)-pe

def square_factor(x,y,z):
    if x==0:
        assert z==0
        return Mat.diag(0,sqrt(y))
    return Mat([[sqrt(x),z/sqrt(x)],[0,sqrt(y-z*z/x)]])

def raw_branch_vectors(W,d,s):
    # W prepackages state/tangent, with derivative i d W[:,1].
    B0=kron(Mat.diag(1,sqrt(s)),eye(4))
    B1=kron(Mat([[0,sqrt(1-s)],[0,0]]),eye(4))
    B=(B0,B1)
    v=[clean(b*W[:,0]) for b in B]
    h=[clean(Q*a + d*b*W[:,1]) for a,b in zip(v,B)]
    return list(zip(v,h))

def map_branch_to_carrier(v,h,P):
    assert eq(P*v,v) and eq(P*h,h)
    prob=S(v.dot(v))
    assert prob>0
    state=clean(v/sqrt(prob))
    hor=clean(h/sqrt(prob)-state*(state.dot(h)/sqrt(prob)))
    norm=S(hor.dot(hor))
    if norm==0:
        return Mat.vstack(state.T,zeros(1,8)),prob,0
    normroot=sqrt(norm)
    basis=clean(hor/normroot)
    assert S(state.dot(basis))==0 and S(basis.dot(basis))==1
    return Mat.vstack(state.T,basis.T),prob,norm

def physical_erase(jets):
    # 2x8 compression Kraus operators are completed by reset of the orthogonal
    # remainder. This computes the actual density jet after discarding the flag.
    comp=[map_branch_to_carrier(v,h,P) for (v,h),P in zip(jets,(pe,po))]
    C=[a[0] for a in comp]
    T=clean(sum((a.T*a for a in C),zeros(8)))
    # T is an orthogonal range projector, even for the zero-tangent seed branch.
    assert eq(T*T,T)
    residual=eye(8)-T
    E=[Mat.vstack(residual[j,:],zeros(1,8)) for j in range(8)]
    Kraus=C+E
    assert eq(sum((a.T*a for a in Kraus),zeros(8)),eye(8))
    rho=sum((v*v.T for v,h in jets),zeros(8))
    antisym=sum((h*v.T-v*h.T for v,h in jets),zeros(8))
    rr=clean(sum((k*rho*k.T for k in Kraus),zeros(2)))
    dd=clean(sum((k*antisym*k.T for k in Kraus),zeros(2)))
    assert eq(rr,Mat.diag(1,0))
    dnew=S(dd[1,0]); assert eq(dd,Mat([[0,-dnew],[dnew,0]]))
    # Actual SLD: L=2i dd. Trace(rho L^2)=4 dnew².
    assert eq(rr*dd+dd*rr,dd)
    return [a[1:] for a in comp],dnew

def sample(s,d,x,z):
    q=d*d; a=s/(1-s)
    p0=1-(1-s)*x
    T=q+(s*x*(1-x)+2*s*d*z-(1-s)*q*z*z/x)/p0
    y=z*z/x+x*T/q
    y_upper=1-z*z/(1-x)
    assert z*z/x<=y<=y_upper
    R0=square_factor(1-x,1-y,-z); R1=square_factor(x,y,z)
    assert eq(R1.T*R1,Mat([[x,z],[z,y]]))
    assert eq(R0.T*R0+R1.T*R1,eye(2))
    W=zeros(8,2); W[:2,:]=R0; W[6:,:]=R1
    jets=raw_branch_vectors(W,d,s)
    comp,dnew=physical_erase(jets)
    assert comp==[(p0,T),(1-p0,T)]
    assert S(dnew**2)==T and dnew>0
    return T,y,4*T

for s in (R(1,4),R(1,2),R(16,25)):
    W=zeros(8,2)
    W[0,0]=W[6,0]=sqrt(R(1,2))
    W[0,1]=sqrt(R(1,2)); W[6,1]=-sqrt(R(1,2))
    comp,dnew=physical_erase(raw_branch_vectors(W,0,s))
    assert S(dnew-sqrt(s)/2)==0
    print("Seed physical erasure",s,"carrier q",S(dnew*dnew))

out=sample(R(1,4),R(1,4),R(1,4),R(0))
assert out==(R(25,208),R(25,52),R(25,52))
print("Startup physical erasure q, y, SLD:",out)
out=sample(R(1,2),R(4),R(1,2),R(1,8))
assert out==(R(33,2),R(35,64),R(66))
print("Large-mode physical erasure q, y, SLD:",out)

def verify_arithmetic(s):
    a=s/(1-s); W=1-(1-s)/4; seed=s/4
    threshold=max(4*a,F(1))**2
    gain=3*s/(4*(3+s))
    k=max(0,-((seed-threshold)//gain))  # exact ceil((threshold-seed)/gain)
    Qbase=seed+k*gain
    assert seed>=s/(16*W) and Qbase>=threshold
    if k:
        assert Qbase-gain<threshold
    r=max(4*a,F(1)); x=1-2*a/r; z=a*x/r
    assert 0<x<1 and z*z<=x*(1-x)
    q=r*r; W=1-(1-s)*x
    T=q+(s*x*(1-x)+2*s*r*z-(1-s)*q*z*z/x)/W
    ymax=1-z*z/(1-x)
    lossmax=q*(ymax/x-z*z/x**2)
    assert T==q+a*x and lossmax-T==a*(r*(2/x-F(1,2))-x)>0
    return k
for s in (F(1,100),F(1,4),F(1,2),F(9,10),F(99,100)):
    print("Sufficient exact startup steps after seed",str(s)+":",verify_arithmetic(s))
