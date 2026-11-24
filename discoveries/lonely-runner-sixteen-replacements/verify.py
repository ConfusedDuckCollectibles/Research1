#!/usr/bin/env python3
"""Exact streamed verification for the independent C++ witness generator.

Does not compute or intersect good regions: each prefix only verifies one
proposed interval analytically, and each leaf verifies a rational witness.
"""
import argparse
from fractions import Fraction
import sys

def check(lines, partial=False):
    it = iter(lines)
    def row(tag, num):
        r = next(it).split()
        if len(r)!=num+1 or r[0]!=tag:
            raise ValueError(f"expected {tag} record, got {r!r}")
        return tuple(map(int, r[1:]))

    n,K=map(int,next(it).split())
    if n!=16 or not 2<=K<=5:
        raise ValueError("wrong task")
    caps=[0]*K
    visited=[0]*(K+1)
    counts=[0,0]
    minimum=Fraction(1,2)
    largest=0

    def walk(S,lo,r):
        nonlocal minimum,largest
        depth=K-r
        if r==0:
            a,b=row('W',2)
            if not 0<a<b: raise ValueError("bad rational time")
            gap=b
            for v in S:
                z=(a*v)%b
                d=min(z,b-z)
                if n*d<=b: raise ValueError("leaf witness not strict")
                gap=min(gap,d)
            if Fraction(gap,b)<minimum: minimum=Fraction(gap,b)
            largest=max(largest,S[-1])
            visited[depth]+=1
            return
        a,b,c,d=row('I',4)
        if not (b>0 and d>0 and 0<=a*d<c*b<=b*d):
            raise ValueError("bad rational endpoints")
        p,q=c*b-a*d,b*d
        abmid=a*d+c*b
        for w in S:
            j=w*abmid//(2*q)
            if n*(w*a-j*b)<b or n*(w*c-j*d)>(n-1)*d:
                raise ValueError("interval crosses bad band")
        bound=2*r*(n-2)*q//(n*(n-2*r)*p)
        visited[depth]+=1
        caps[depth]=max(caps[depth],bound)
        for w in range(lo,bound+1): walk(S+(w,),w+1,r-1)

    try:
        for mask in range(1<<(n-1)):
            if mask.bit_count()!=K:
                continue
            counts[0]+=1
            if row('R',1)!=(mask,):
                raise ValueError("missing or out-of-order removal set")
            S=tuple(v for v in range(1,n) if not ((mask>>(v-1))&1))
            walk(S,n,K)
            counts[1]+=1
    except StopIteration:
        if not partial: raise ValueError("truncated certificate") from None
        status="VERIFIED PREFIX ONLY"
    else:
        if next(it,None) is not None:
            raise ValueError("unexpected records after final root")
        status="VERIFIED COMPLETE"
    print(status,"k",K,"verified finished roots",counts[1],"roots expected through",counts[0],
          "records checked",visited,"caps",caps,"witness minimum",minimum,"max terminal speed",largest)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("input",nargs="?",default="-")
    p.add_argument("--partial",action="store_true",help="report only a verified prefix if input terminates early")
    a=p.parse_args()
    if a.input=="-": check(sys.stdin,a.partial)
    else:
        with open(a.input,encoding="ascii") as f: check(f,a.partial)
