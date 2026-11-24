"""Exact standalone sanity checks for orientation, strictness, and C4 types."""
from itertools import product
def verify(n,A):
    assert all(u!=v and (v,u) not in A for u,v in A)
    N=[{v for u,v in A if u==i} for i in range(n)]
    I=[{u for u,v in A if v==i} for i in range(n)]
    NN=[set().union(*(N[j] for j in N[i]))-N[i]-{i} for i in range(n)]
    return N,I,NN

for r in range(1,31):
    # Published sharp limiting construction; the 3r+1 endpoint matters.
    n=3*r+1
    A={(i,(i+a)%n) for i in range(n) for a in range(1,r+1)}
    N,I,NN=verify(n,A)
    assert all(len(N[i])==len(I[i])==len(NN[i])==r and not I[i]&NN[i] for i in range(n))
    assert all(len(N[i])+len(I[i])+len(NN[i])+1==n for i in range(n))
    # Lower the order by one while preserving all steps: the r-step triangle appears.
    n-=1
    A={(i,(i+a)%n) for i in range(n) for a in range(1,r+1)}
    N,I,NN=verify(n,A)
    assert (0,r) in A and (r,2*r) in A and (2*r,0) in A
    assert any(I[i]&NN[i] for i in range(n))
print('30 exact cyclic controls passed (sharpness and triangle at changed order)')

C={(i,(i+1)%4) for i in range(4)}
T=[]
for z in product(range(3),repeat=4): # 0 sends to cycle, 1 nonadjacent, 2 receives
    if z.count(1)!=z.count(2):continue
    if any(z[i]==0 and z[j]==2 for i,j in C):continue
    T.append(z)
assert sorted(z.count(2) for z in T)==[0]+[1]*4+[2]*6
for i in range(4):
    z=next(z for z in T if z.count(2)==1 and z[i]==2)
    assert z[(i-1)%4]==1 and z[(i+1)%4]==z[(i+2)%4]==0
    for v in T:
        if v.count(2)!=2:continue
        admissible=not any(v[j]==2 and z[j]==0 for j in range(4))
        assert admissible==(v[i]==v[(i-1)%4]==2)
print('all 11 C4 local attachment types and exact compatibility checked')
