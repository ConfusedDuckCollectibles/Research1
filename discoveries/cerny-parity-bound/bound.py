"""Exact integer bounds from PROOF.md, no subset search."""
from collections import deque

KINDS = ('A','B','C')

def parity_collision(T,kind):
    """Shortest depth by which a state receives both parity colors."""
    n=len(T[0])
    sigma=(1,1 if kind=='B' else 0)
    dist=[[None]*2 for _ in range(n)]
    Q=deque()
    for q,b in enumerate((1,0,0)):
        Q.append((q,b));dist[q][b]=0
    while Q:
        q,b=Q.popleft()
        d=dist[q][b]
        for f,v in zip(T,sigma):
            p=f[q]; c=b^v
            if dist[p][c] is None:
                dist[p][c]=d+1;Q.append((p,c))
    both=[max(x,y) for x,y in dist if x is not None and y is not None]
    return min(both,default=None)

def profile_limit(n,kind,h):
    if h<=2:return n-1 if kind=='C' else n-2
    if kind=='C':return n-2*(h//2)
    return n+1-2*((h+1)//2)

def theoretical(n,kind,tau=None):
    if n<3 or kind not in KINDS:raise ValueError('out of scope')
    u=min(n-1,tau+1) if tau is not None else n-1
    if kind=='C':
        E={2*j:max(j-1,n-j-1,1) for j in range(2,n)}
    else:
        E={2*j:max(j,n-j-1,1) for j in range(1,n-1)}
    H={h:max(0,min(profile_limit(n,kind,1),profile_limit(n,kind,h))) for h in range(1,u+1)}
    O=n*u-(u*u)//2 if kind=='C' else n*u-(u*u+1)//2-2
    assert sum(H.values())==O
    return sum(E.values())+O,E,H

def profile(T,kind):
    """Sharper instance bound from actual parity/ordinary ball sizes."""
    n=len(T[0]); sigma=(1,1 if kind=='B' else 0)
    w=(3,1)+(2,)*(n-2) if kind=='C' else (1,1)+(2,)*(n-2)
    W=sum(w); low=3 if kind=='C' else 2
    P={(0,1),(1,0),(2,0)}
    even=[set(range(2*((low+1)//2),W,2))]
    odd=[set(range(low|1,W,2))]
    for h in range(2,n):
        R={q for q,b in P}
        direct=sum(w[q] for q in R)
        free=n-len(R)
        levels=range(0,2*free+1,2)
        even.append(({z for j in levels for z in (j,direct+j) if low<=z<W})&even[-1])
        if len(P)==len(R):
            w0=sum(w[q] for q,b in P if not b)
            w1=direct-w0
            possible={z for j in levels for z in (w0+j,w1+j) if low<=z<W}
        else:
            possible=set()
        odd.append(possible&odd[-1])
        P=P|{(f[q],b^eps) for q,b in P for f,eps in zip(T,sigma)}
    return sum(map(len,even))+sum(map(len,odd)),even,odd

if __name__=='__main__':
    print('n | uniform AB,C | sufficient tau AB,C (assuming weak early test)')
    for n in (4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,24,26,30,32,50,100,1000):
        out=[]
        for k in ('A','C'):
            u,*_=theoretical(n,k)
            t=max((t for t in range(1,n) if theoretical(n,k,t)[0]<=(n-1)**2),default=0)
            out.append((u,t))
        print(n,out)
