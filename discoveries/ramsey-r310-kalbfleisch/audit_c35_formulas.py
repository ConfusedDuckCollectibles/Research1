#!/usr/bin/env python3
"""Independent CNF semantic and input audit for the induced C35 branch."""
from collections import Counter,defaultdict
from itertools import combinations,product
from pathlib import Path
import argparse,json
H=Path(__file__).parent

def U(n,E):
    return [set(v if i==u else u for u,v in E if i==u or i==v) for i in range(n)]

def data(n):
    return set(combinations(range(n),2))

def g6(path):
    s=path.read_text().strip().removeprefix('>>graph6<<')
    n=ord(s[0])-63
    x=list(int(b) for c in s[1:] for b in f'{ord(c)-63:06b}')
    return n,{(i,j) for b,(i,j) in zip(x,((i,j) for j in range(n) for i in range(j))) if b}

def core_data():
    C={(a,b) for a,b in data(35) if min((b-a)%35,(a-b)%35) in (1,7,11,16)}
    assert len(C)==140
    assert all(not (U(35,C)[i]&U(35,C)[j]) for i,j in C)
    assert all(U(35,C)[i]&U(35,C)[j] for i,j in data(35)-C)
    return C

def iso(C,E,V):
    V=sorted(V)
    E={tuple(sorted((i,j))) for i,j in E if i in V and j in V}
    W={v:set() for v in V}
    A=U(35,C)
    for i,j in E:
        W[i].add(j);W[j].add(i)
    left={i:{j:len(A[i]&A[j]) for j in range(35)} for i in range(35)}
    right={i:{j:len(W[i]&W[j]) for j in V} for i in V}
    P={0:V[0]}
    free=set(V)-{V[0]}
    def go():
        if len(P)==35:return True
        xs=(i for i in range(35) if i not in P)
        i,Q=min(((i,[r for r in sorted(free) if all(
                        (j in A[i])==(s in W[r]) and left[i][j]==right[r][s] for j,s in P.items())])
                   for i in xs),key=lambda r:len(r[1]))
        for r in Q:
            P[i]=r;free.remove(r)
            if go():return True
            free.add(r);del P[i]
        return False
    assert go()
    assert E=={tuple(sorted((P[i],P[j]))) for i,j in C}
    return [P[i] for i in range(35)]

def count_independent(n,E,maximum):
    # Different recursion orientation and representation: sets of surviving
    # nonneighbors, maintaining the unique largest-to-smallest ordering.
    A=U(n,E)
    out=defaultdict(set)
    def go(T,P):
        out[len(T)].add(tuple(sorted(T)))
        if len(T)==maximum:
            return
        while P:
            x=max(P)
            P=P-{x}
            go(T+(x,),P-A[x])
    go((),set(range(n)))
    return out

EXT={'star':({(0,4),(1,4),(2,4),(3,4)},(1,1,1,1,4)),
     'bipartite':({(0,2),(0,3),(0,4),(1,2),(1,3),(1,4)},(3,3,2,2,2)),
     'pentagon':({(0,1),(0,4),(1,2),(2,3),(3,4)},(2,)*5)}

def read_dimacs(path):
    L=[]
    for line in path.read_text().splitlines():
        if not line or line[0]=='c':continue
        if line[0]=='p':
            h=line.split()
            assert h[:2]==['p','cnf']
            var,m=map(int,h[2:])
            continue
        vals=list(map(int,line.split()))
        assert vals[-1]==0 and all(0<abs(a)<=var for a in vals[:-1])
        assert len(set(vals[:-1]))==len(vals[:-1]) and not set(vals[:-1])&{-a for a in vals[:-1]}
        L.append(tuple(vals[:-1]))
    assert len(L)==m
    return var,L

def k(i,j):return 5*i+j+1
def graph_only(C,F,I):
    A=U(35,C);B=U(5,F)
    expected=Counter()
    def push(*cl):expected[tuple(sorted(cl))]+=1
    for i in range(35):
        for a,b in data(5):push(-k(i,a),-k(i,b))
    for i,j in C:
        for a in range(5):push(-k(i,a),-k(j,a))
    for size in range(2,6):
        for T in combinations(range(5),size):
            if any((a,b) in F for a,b in combinations(T,2)):continue
            for S in I.get(10-size,()):
                push(*(k(i,a) for i,a in product(S,T)))
    for i,a in product(range(35),range(5)):
        push(k(i,a),*(k(j,a) for j in A[i]),*(k(i,b) for b in B[a]))
    push(1)
    return expected

def audit_arithmetic(n,L,degree):
    gates=defaultdict(list)
    terminal=Counter()
    G=[]
    for cl in L:
        peak=max(map(abs,cl),default=0)
        if peak<=175:G.append(cl)
        elif len(cl)==1:terminal[cl[0]]+=1
        else:gates[peak].append(cl)
    last=175
    expect_terminal=Counter()
    for a,d in enumerate(degree):
        cutoff=10-d
        ids={(i,j):None for i,j in product(range(36),range(cutoff+1))}
        def term(i,j):
            return True if j==0 else ids[i,j] if j<=i else False
        for i in range(1,36):
            x=k(i-1,a)
            for j in range(1,min(i,cutoff)+1):
                last+=1
                s=last
                ids[i,j]=s
                A=term(i-1,j);B=term(i-1,j-1)
                cl=gates.pop(s)
                support={v for v in (A,B,s,x) if type(v) is int}
                assert all(abs(v) in support for a in cl for v in a)
                for vals in product((False,True),repeat=len(support)):
                    M=dict(zip(sorted(support),vals))
                    ok=all(any(M[abs(z)]==(z>0) for z in q) for q in cl)
                    aver=A if type(A) is bool else M[A]
                    bver=B if type(B) is bool else M[B]
                    assert ok==(M[s]==(aver or (bver and M[x])))
        if d<4:expect_terminal[ids[35,4-d]]+=1
        expect_terminal[-ids[35,cutoff]]+=1
    assert not gates and n==last and terminal==expect_terminal
    return G

def main():
    C=core_data()
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--formulas-only',action='store_true',help='No optional third-party graph inputs.')
    args=p.parse_args()
    if not args.formulas_only:
        nn,M=g6(H/'mckay-kalbfleisch.g6')
        assert nn==35
        P=iso(C,M,range(35))
        e=(H/'exoo39.txt').read_text().split()
        E={(a,b) for a,b in data(39) if e[a][b]==e[b][a]=='0'}
        assert len(E)==167
        I39=count_independent(39,E,10)
        assert I39[9] and not I39.get(10) and all(not (U(39,E)[i]&U(39,E)[j]) for i,j in E)
        Q=iso(C,E,set(range(39))-{22,35,37,38})
        (H/'baseline-isomorphisms.json').write_text(json.dumps({'core_to_mckay':P,'core_to_exoo':Q},separators=(',',':'))+'\n')
    I=count_independent(35,C,9)
    assert [len(I[i]) for i in range(9)]==[1,35,455,2905,10150,20265,22995,13760,3360] and not I.get(9)
    for n,(F,d) in EXT.items():
        v,L=read_dimacs(H/f'{n}.cnf')
        G=audit_arithmetic(v,L,d)
        assert Counter(tuple(sorted(a)) for a in G)==graph_only(C,F,I)
        print(n,'checked every graph clause and every degree-counter truth table:',v,len(L))
    if not args.formulas_only:
        print('Exoo 39/167 positive and both C35 graph isomorphisms checked.')
if __name__=='__main__':
    main()
