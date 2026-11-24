"""Exact structured random samples; no numerical results claimed as theorem."""
from collections import deque
import itertools
import random
import sys
from bound import parity_collision, theoretical, profile
N=int(sys.argv[1]) if len(sys.argv)>1 else 9
B=int(sys.argv[2]) if len(sys.argv)>2 else 500
rng=random.Random(428929)
case_cfg=[((2,0,0,2),(0,1)),((2,2,0,0),(0,1)),((2,0,0,2),(1,1))]
def build(n,kind):
    x,indeg=case_cfg[kind]
    L=[j for j in range(3,n) for _ in range(2)]+[0]*indeg[0]+[1]*indeg[1]+[2]*(2-sum(indeg))
    rng.shuffle(L)
    return (x[0],x[2],*L[::2]),(x[1],x[3],*L[1::2])
def is_sc(T):
    n=len(T[0])
    for s in (0,):
        R={s}; todo=[s]
        for x in todo:
            for f in T:
                y=f[x]
                if y not in R: R.add(y);todo.append(y)
        if len(R)<n:return False
    R={0};todo=[0]
    for x in todo:
        for f in T:
            for y,z in enumerate(f):
                if z==x and y not in R:R.add(y);todo.append(y)
    return len(R)==n
def tab(T):
    n=len(T[0]);times=1<<n
    out=[]
    for f in T:
        t=[0]*times
        for i,j in enumerate(f): t[1<<j]|=1<<i
        for mask in range(1,times):
            if mask & (mask-1):
                u=mask&-mask
                t[mask]=t[mask^u]|t[u]
        out.append(t)
    return out
def control_forward(T): # Distinct representation, sets with forward actions.
    ini=frozenset(range(len(T[0])))
    dist={ini:0};Q=[ini]
    for K in Q:
        if len(K)==1:return dist[K]
        for f in T:
            J=frozenset(f[i] for i in K)
            if J not in dist:dist[J]=dist[K]+1;Q.append(J)
    return None
def backwards_exts(P,weight):
    S=len(P[0]);n=(S-1).bit_count()
    BC=[s.bit_count() for s in range(S)]
    W=[sum(weight[i] for i in range(n) if s>>i&1) for s in range(S)]
    R=[[] for _ in range(S)]
    for f in P:
        for i,y in enumerate(f):R[y].append(i)
    def BFS(src):
        D=[-1]*S;Q=list(src)
        for x in Q:D[x]=0
        for x in Q:
            for y in R[x]:
                if D[y]<0:D[y]=D[x]+1;Q.append(y)
        return D
    # From mask K, shortest composition of preimage operators to bigger (mass or card)
    # is via multi-source distance in the reverse preimage graph.
    DS=BFS([S-1])
    rt=min((DS[1<<i] for i in range(n) if DS[1<<i]>=0),default=-1)
    if rt<0:return rt,{},{}
    E={}
    for k in range(1,n):
        D=BFS([s for s in range(S) if BC[s]>k])
        for s in range(S):
            if BC[s]==k:E[s]=D[s]
    M={}
    for v in range(1,sum(weight)):
        D=BFS([s for s in range(S) if W[s]>v])
        for s in range(S):
            if W[s]==v:M[s]=D[s]
    return rt,E,M
if __name__=='__main__':
    for n in range(3,min(10,N)+1):
        T=(tuple((i+1)%n for i in range(n)), tuple(i if i<n-1 else 0 for i in range(n)))
        P=tab(T)
        assert all(all({i for i,y in enumerate(a) if s>>y&1}==
                           {i for i in range(n) if p[s]>>i&1}
                       for s in range(1<<n)) for a,p in zip(T,P))
        assert backwards_exts(P,[2]*n)[0]==control_forward(T)==(n-1)**2
    print("Cerny controls passed")
    for c in range(3):
        best_m=[0]*8; worst_ext=0; survivors=0; bad=None; hi=0; maxrt=0
        prefix_certified=0; tau_max=0; profile_certified=0
        seen=set()
        for _ in range(B):
            T=build(N,c)
            if T in seen:continue
            seen.add(T)
            if not is_sc(T):continue
            w=[1 if c<2 else 3,1]+[2]*(N-2)
            assert all(sum(w[x] for f in T for x in range(N) if f[x]==y)==2*w[y] for y in range(N))
            rt,E,M=backwards_exts(tab(T),w)
            if rt<0:continue
            assert rt==control_forward(T)
            kind='ABC'[c]
            tau=parity_collision(T,kind)
            assert tau is not None
            tau_max=max(tau_max,tau)
            bound,Eb,Hb=theoretical(N,kind,tau)
            level={}
            cutoff=3 if kind=='C' else 2
            for K,m in M.items():
                if m>N-1:raise AssertionError((T,K,m))
                v=sum(w[i] for i in range(N) if (K>>i)&1)
                if v<cutoff:continue
                level[v]=max(level.get(v,0),m)
            assert all(m<=Eb[v] for v,m in level.items() if not v&1)
            assert all(m<=max(Hb) for v,m in level.items() if v&1)
            assert all(sum(v&1 and m>=h for v,m in level.items())<=b for h,b in Hb.items())
            assert rt<=sum(level.values())<=bound
            refined,pe,po=profile(T,kind)
            assert all({v for v,m in level.items() if v%2==i and m>=h}<=V
                       for i,P in enumerate((pe,po)) for h,V in enumerate(P,1))
            assert rt<=sum(level.values())<=refined<=bound
            profile_certified+=refined<=(N-1)**2
            prefix_certified+=bound<=(N-1)**2
            maxrt=max(maxrt,rt)
            survivors+=1
            for K,m in M.items():
                best_m[K&7]=max(best_m[K&7],m)
            for K,d in E.items():
                if d>worst_ext:worst_ext=d;bad=(K,T,rt)
                if d>N:hi+=1
        print(N,c,'attempts',B,'distinct',len(seen),'SC sync',survivors,'max tau',tau_max,'bound certified count',prefix_certified,'profile certified',profile_certified,'maxrt',maxrt,'max ext',worst_ext, 'violations',hi,'m per pattern',best_m,'max case',bad)
