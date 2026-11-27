from itertools import combinations, permutations, product
from collections import defaultdict
from fractions import Fraction
from math import comb
import sympy as s
from sympy import QQ
from sympy.polys.matrices import DomainMatrix as DM
Z=[34,43,-22,103,-14,-28,4,-8,-10,-94,37,8,13,-17,-18,-38,-14,-16,-6,-13,0,-35,-92,-91,-90,13,93,0,5,53,-6,6,58,-6,49,73,-11,-10,10,-6,56,74,-18,8,55,-4,60,88,-16,32,38,-29,12,-7,17,35,38,-15,-10,-24,-4,-72,29,69,3,15,-37,48,74,24,2,13,20,17,-1,-91,8,24,-11,-47,-7,-3,76,115,2,-87,-69]
m=s.Symbol('m')
def bn(x,k):
  if k<0: return 0
  return s.expand(s.prod(x-j for j in range(k))/s.factorial(k))
op=[4,5,3]
per=list(permutations(range(3)))
def act(T,p):
  a,z,k=T
  mp=list(p)+[op[p[op.index(j)]] for j in range(3,6)]
  return (sum(1<<mp[j] for j in range(6) if a>>j&1),z,k)
D={11,22,37}
S=sorted([(u,0,3-u.bit_count()) for u in range(64) if u.bit_count()<=3 and u not in D]+[(1<<i,1,0) for i in range(3)])
N=len(S)
gact=[{x:act(x,p) for x in S} for p in per]
def no(x,y):return not(x[0]&y[0] or x[1]*y[1])
def ky(x,y):return min(tuple(sorted([g[x],g[y]])) for g in gact)
ks=sorted({ky(x,y) for x,y in product(S,repeat=2) if no(x,y)})
K=len(ks);ki={x:i for i,x in enumerate(ks)}
assert (N,K)==(42,75)
at={(x,y):ki[ky(x,y)] for x,y in product(S,repeat=2) if no(x,y)}
eq=[]
for anchor in [True,False]:
  for y in S:
    if anchor and y[0]&1:continue
    R=[0]*K
    for x in S:
      if (x[0]&1 if anchor else x[2]>=1) and no(x,y):
        R[at[x,y]]+=(1+x[1])*bn(m-y[2]-int(not anchor),x[2]-int(not anchor))
    eq.append(R)
A=s.Matrix(eq);As=A[list(A.subs(m,10).T.rref()[1]),:]
target=(m+5)*(m+4)/2
base=s.Matrix([(m+5)/(m+3) if x[1] or y[1] else (m+5)*(m+4)/((m+2)*(m+1)) for x,y in ks])
Dall=s.prod(m+i for i in [1,2,3,4,6])
vm=[(i,j) for i,(x,y) in enumerate(ks) for j in range(3-max(0,x[2]+y[2]-4))]
R=[];rhs=[]
for i in range(As.rows):
  ps=[s.Poly(x,m) for x in As.row(i)]
  rp=s.Poly(s.cancel((target-(As.row(i)*base)[0])*Dall/(m+5)),m)
  for k in range(5):
    row=[ps[a].nth(k-j) if k>=j else 0 for a,j in vm]
    if any(row) or rp.nth(k):R.append(row);rhs.append(rp.nth(k))
u,piv=s.Matrix(R).row_join(s.Matrix(rhs)).rref()
n=len(vm);free=[i for i in range(n) if i not in piv]
assert len(free)==len(Z)
vals=[s.Integer(0)]*n
for f,z in zip(free,Z):vals[f]=s.Integer(z)
for i,p in enumerate(piv):vals[p]=u[i,n]-sum(u[i,f]*vals[f] for f in free)
cor=[0]*K
for (a,j),v in zip(vm,vals):cor[a]+=v*m**j
q=[s.Poly(s.cancel(Dall*(base[i]-1)+(m+5)*cor[i]),m) for i in range(K)]
di=s.Poly(Dall,m)
for R in eq:
  assert s.expand(sum(a*(x.as_expr()+Dall) for a,x in zip(R,q))-target*Dall)==0
print('orbit and star identities',len(eq),flush=True)
B=defaultdict(list); seen=set(); tt=gact[1]
for x in S:
  if x in seen:continue
  orb=sorted({g[x] for g in gact});seen.update(orb);k=x[2]
  B['triv'].append((k,dict.fromkeys(orb,1)))
  if len(orb)==3:
    B['std'].append((k,{a:2 if tt[a]==a else -1 for a in orb}))
  if len(orb)==6:
    B['sgn'].append((k,{g[x]:(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3)) for p,g in zip(per,gact)}))
    pairs=sorted({tuple(sorted([y,tt[y]])) for y in orb})
    B['std'].append((k,{a:c for P,c in zip(pairs,[2,-1,-1]) for a in P}))
    B['std'].append((k,{a:c for P,c in zip(pairs,[0,1,-1]) for a in P}))
from functools import lru_cache
@lru_cache(None)
def gam(d,k,l,t):
  # Norm of harmonic lift times eigenvalue of exact-intersection operator.
  return s.Poly(bn(m-2*d,k-d)*sum((-1)**(d-p)*comb(d,p)*comb(k-d,t-p)*bn(m-d-k,l-d-t+p) for p in range(d+1) if 0<=t-p<=k-d and l-d-t+p>=0),m)
def bm(d,bs):
  ret=[]
  for k,V in bs:
    row=[]
    for l,W in bs:
      v=s.Poly(0,m)
      for x,c in V.items():
        for y,e in W.items():
          if not c*e:continue
          for t in range(min(k,l)+1):
            a=gam(d,k,l,t)
            if x==y and t==k:
              h=di*s.Poly(target/(1+x[1])-1,m)
            elif t or not no(x,y):h=-di
            else:h=q[at[x,y]]
            v+=a*h*(c*e)
      row.append(v)
    ret.append(row)
  return ret
def mk(M,x,ids):
  return s.Matrix([[M[i][j].eval(x) for j in ids] for i in ids])
def chk(M,i):
  p=list(M.rref()[1]);assert len(p)==M.rows-i,(M.rows,i,len(p));L,D=M[p,p].LDLdecomposition();assert all(D[j,j]>0 for j in range(len(p)))
def toq(v):
  n,d=v.as_numer_denom();return QQ(int(n),int(d))
def det(A):return DM(A,(len(A),len(A)),QQ).det()
pc=0;bc=[]
for d in range(4):
  for sec in ['triv','sgn','std']:
    bs=[x for x in B[sec] if x[0]>=d]
    if not bs:continue
    M=bm(d,bs);D0=len(bs);nz=2 if (d,sec)==(0,'triv') else int((d,sec) in [(0,'std'),(1,'triv')])
    for i in range(1,6):
      ids=[j for j,(k,_) in enumerate(bs) if d<=k<=i-d]
      if len(ids)>nz: chk(mk(M,i,ids),nz)
    mat=mk(M,6,range(D0));ids=list(mat.rref()[1]);assert len(ids)==D0-nz
    # A Newton positivity check is independent of the submitted power-basis conversion.
    degs=[sum(max(M[i][j].degree() for j in ids[:k]) for i in ids[:k]) for k in range(1,len(ids)+1)]
    ev=[[[toq(M[i][j].eval(x+6)) for j in ids] for i in ids] for x in range(max(degs)+1)]
    for k,deg in enumerate(degs,1):
      y=[det([row[:k] for row in a[:k]]) for a in ev[:deg+1]]
      assert y[0]>0
      while len(y)>1:
        y=[b-a for a,b in zip(y,y[1:])]
        assert y[0]>=0,(d,sec,k)
    pc+=len(ids);bc.append((d,sec,len(ids)))
    print('PSD block',bc[-1],flush=True)
print('finite bases 7–11 and all-size minors',pc,bc,flush=True)
