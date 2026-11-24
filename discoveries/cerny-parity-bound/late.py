"""An exact limitation: an unbounded near-parity family rejected by the bound."""
from collections import deque
from bound import theoretical, profile, parity_collision
from examples import graph_cycle_lengths, reset_word, strong
from pilot import tab, backwards_exts, control_forward

def build(n):
    assert n>=9 and n%2
    a=[2,0,*range(3,n),1]
    b=list(range(n));b[1]=2
    u,v=n-5,n-3
    a[u],a[v]=a[v],a[u]
    b[u],b[v]=b[v],b[u]
    return a,b

def independently_distinguish(T):
    n=len(T[0]);init=frozenset(((0,1),(1,0),(2,0)))
    now=init
    for i in range(1,2*n):
        nxt=now|{(f[q],p^eps) for q,p in now for f,eps in zip(T,(1,0))}
        if len({q for q,p in nxt})<len(nxt):return i
        if nxt==now:return None
        now=nxt
    raise AssertionError('joint parity reachability did not terminate')

for n in range(9,62,2):
    T=build(n)
    assert strong(T) and graph_cycle_lengths(T[0])==[2,n-2]
    assert parity_collision(T,'A')==independently_distinguish(T)==n-4
    rw=reset_word(T)
    assert rw
    image=set(range(n))
    for c in map(int,rw):image={T[c][q] for q in image}
    assert len(image)==1
    p,_,_=profile(T,'A')
    assert p==5*((n-1)//2)**2-16<=theoretical(n,'A',n-4)[0]
    if n in (9,13,21,31,61):print(n,n-4,p,(n-1)**2)
T=build(13)
rt,E,M=backwards_exts(tab(T),[1,1]+[2]*11)
assert rt==control_forward(T)<=144
K,d=max(E.items(),key=lambda x:x[1])
print('n=13 exact reset threshold',rt,'max exact cardinal extension',d,'set mask',K)
