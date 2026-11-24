#!/usr/bin/env python3
"""Exact arithmetic checks and graph controls; no packages or inputs."""
from fractions import Fraction
from itertools import combinations


def c2(x):
    # Signed polynomial extension, intentional for x < 0.
    return x*(x-1)//2


def g(q, points, lines, marks, s):
    return s*(points+1-s)-(marks-(lines-1))+c2(q+1-s)


def budget(q, points, lines, marks):
    return c2(points)-lines*c2(q+1)+q*((q+1)*lines-marks)


def arithmetic():
    shape_checks = 0
    for q in range(9,401):
        N = q*(q+1)
        V = N+1
        j0 = q//2
        j = j0+1
        B = (q+1)*N
        E = B-j
        t = q-j
        K = (q-1)*(2*j-q)//2
        orig = budget(q,V,N,E)
        assert orig == N//2+q*j
        assert g(q,V,N,E,q+2)==q*q+j
        assert orig-g(q,V,N,E,q+2)==K
        assert 0<K<=q-1 and K//(q-1)<=1
        assert g(q,V,N,E,q+3)-g(q,V,N,E,q+2)==q*q-q-1>K
        assert g(q,V,N,E,V)>g(q,V,N,E,q+3)>orig
        side = budget(q,N,V,E)
        assert side==q*j
        assert g(q,N,V,E,q+2)-side==(q-1)*(q-j)-1>0
        assert g(q,N,V,E,N)>g(q,N,V,E,q+2)
        if q<=80:
            for Q,M,low,limit in [(V,N,q+3,orig),(N,V,q+2,side)]:
                for s in range(low,Q+1):
                    assert g(q,Q,M,E,s)>limit
                    shape_checks+=1
        assert V-q-2==q*q-1
        assert q*(q-1)<q*q-2
        assert (q-1)**2+q-1<q*q-2
        assert (q+1)*V-E==q+1+j
        assert (q+2)*(q+1)-(q+2+N-2)==q+2
        assert (q+2)*(q+1)-(q+2+N-1)==q+1
        for h in range(q):
            gap=(q-3)*j-2*(q-1-h)-c2(h)
            assert gap>0
        assert (q-3)*j-2*(q-1)>=Fraction(q*q-6*q+1,2)>0
        assert (q-3)*j-c2(q-1)>=Fraction(q-5,2)>0
        assert N-j>=q*q+4  # Dow threshold is below q^2 at q>=9.
        assert t>=2
        assert t*t+t+2>q+1+j
        assert 3*(q+1-j//2)-1>q+1+j
        base_gap=budget(q,V,N,B-j0)-g(q,V,N,B-j0,q+2)
        assert base_gap == (0 if q%2==0 else -(q-1)//2)
        assert 2*(N-j0)>=2*q*q+q # Stinson's bound is strict.
    print(f"arithmetic: q=9..400, {shape_checks} explicit high-degree checks")
    for q,value in [(10,1203),(14,3141),(30,28813)]:
        assert (q+1)*q*(q+1)-q//2-2==value
        print(f"q={q}: {value}")


def graph_identity(q, points, lines):
    marks=sum(map(len,lines))
    xs=sum(c2(q+1-len(line)) for line in lines)
    pairs=set()
    for line in lines:
        for p in combinations(sorted(line),2):
            assert p not in pairs
            pairs.add(p)
    assert c2(points)-len(pairs)+xs==budget(q,points,len(lines),marks)
    for i,A in enumerate(lines):
        s=len(A)
        others=[line for k,line in enumerate(lines) if i!=k and not A&line]
        observed=s*(points-s)-sum((a in A)!=(b in A) for a,b in pairs)
        right=s*(points+1-s)-(marks-(len(lines)-1))+sum(len(x)-1 for x in others)
        assert observed==right
        assert observed+xs>=g(q,points,len(lines),marks,s)+(q-1)*len(others)
    # Cross-check transposition without assuming degrees are bounded.
    dual=[{i for i,line in enumerate(lines) if p in line} for p in range(points)]
    assert sum(map(len,dual))==marks
    assert c2(len(lines))-sum(c2(len(line)) for line in dual)+sum(c2(q+1-len(line)) for line in dual)==budget(q,len(lines),points,marks)


def controls():
    q=3
    # Affine lines with their distinct slope points at infinity, no infinity
    # line: the incidence graph of the prime-order plane minus one line.
    pts=q*q+q+1
    finite=[{x*q+(m*x+b)%q for x in range(q)}|{q*q+m}
           for m in range(q) for b in range(q)]
    vertical=[{x*q+y for y in range(q)}|{q*q+q} for x in range(q)]
    assert len(finite+vertical)==pts-1
    graph_identity(q,pts,finite+vertical)
    # A linear hypergraph with negative line defect and an empty line tests
    # the signed extension and intersections with truly disjoint lines.
    q=9
    pts=q*q+q+1
    A=set(range(q+2))
    B=set(range(q+2,q+2+q+1))
    other=[{a,b} for a in sorted(A) for b in range(q+2,pts)]
    graph_identity(q,pts,[A,B,set()]+other[:pts-4])
    print("crossing controls: plane-minus-line and sparse signed defects, both orientations")


if __name__=="__main__":
    arithmetic()
    controls()
    print("PASS")
