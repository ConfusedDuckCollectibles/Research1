#!/usr/bin/env python3
"""Find all relative links and exclude the exceptional case before any gluing."""
import argparse
from collections import defaultdict,Counter
from itertools import combinations,permutations,product
from pathlib import Path

if not __debug__:
    raise RuntimeError("Python assertions must be enabled")

def as_edges(bits):
    return tuple(frozenset(i for i in range(9) if m & (1<<i)) for m in bits)

def degree(E,n):
    a=[0]*n
    for e in E:
        for i in e:
            a[i]+=1
    return a

def signature(points,degrees):
    return tuple(sorted(degrees[i] for i in points))

def getmaps(C,X,G,Y):
    """All bijections X->Y carrying hyperedge list C to G."""
    src=defaultdict(list)
    dst=defaultdict(list)
    dc,dg=degree(C,9),degree(G,9)
    for i in X:
        src[dc[i]].append(i)
    for i in Y:
        dst[dg[i]].append(i)
    if any(len(src[i])!=len(dst[i]) for i in set(src)|set(dst)):
        return
    keys=sorted(src)
    s=tuple(i for k in keys for i in src[k]);target=set(G)
    for blocks in product(*(permutations(dst[k]) for k in keys)):
        m=dict(zip(s,(i for b in blocks for i in b)))
        if {frozenset(m[i] for i in e) for e in C}==target:
            yield m

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('cards',type=Path);ap.add_argument('catalog',type=Path)
    ap.add_argument('links',type=Path)
    p=ap.parse_args()
    rows=[tuple(map(int,line.split())) for line in p.cards.read_text().splitlines()]
    assert len(rows)==242
    assert all(len(a)==len(set(a))==15 and all(0<e<512 and e.bit_count()==4 for e in a) for a in rows)
    raw=[as_edges(a) for a in rows]
    checked=0; rejected_edge=0; rejected_cdeg=0
    for A in raw:
        d=degree(A,9)
        if max(d)>10:
            continue
        for t in combinations([i for i in range(9) if d[i]==7],3):
            checked+=1
            if any(not (1<=len(e.intersection(t))<=2) for e in A):
                rejected_edge+=1
                continue
            if any(sum(i in e and j in e for e in A)!=2 for i,j in combinations(t,2)):
                rejected_cdeg+=1
                continue
            raise AssertionError(f'Exceptional case candidate: {A}, {t}')
    assert (checked,rejected_edge,rejected_cdeg)==(2728,2423,305)
    print('exceptional marks tested',checked,'rejected by edge sizes',rejected_edge,
          'rejected by codegree',rejected_cdeg,'survivors 0')
    catalog=defaultdict(list)
    for line in p.catalog.read_text().splitlines():
        m,*ms=map(int,line.split()); G=as_edges(ms)
        assert len(G)==m
        catalog[m,tuple(sorted(degree(G,8)))].append((tuple(ms),G))
    links=defaultdict(set)
    y=tuple(range(8))
    for A in raw:
        d=degree(A,9)
        for i,di in enumerate(d):
            if di not in (4,5,6):
                continue
            C=tuple(e for e in A if i not in e)
            L=tuple(e-{i} for e in A if i in e)
            x=tuple(j for j in range(9) if j!=i)
            for ms,G in catalog[len(C),signature(x,degree(C,9))]:
                for m in getmaps(C,x,G,y):
                    links[ms].add(tuple(sorted(sum(1<<m[j] for j in e) for e in L)))
    assert len(links)==24
    assert sum(len(ls)*(len(ls)+1)//2 for ls in links.values())==12018
    print('templates:',len(links),'link counts by overlap:',
          [(z,sum(len(ls) for ms,ls in links.items() if len(ms)==z)) for z in (9,10,11)])
    print('total unordered link pairs', sum(len(ls)*(len(ls)+1)//2 for ls in links.values()))
    with p.links.open('w') as f:
        for ms,ls in sorted(links.items()):
            print(len(ms),len(ls),*ms,file=f)
            for l in sorted(ls):
                print(len(l),*l,file=f)

if __name__=='__main__':
    main()
