#!/usr/bin/env python3
"""Independent all-state and all-cut LP checks. Synthetic instances only."""
from fractions import Fraction as F
from itertools import product,permutations
import numpy as np
from scipy.optimize import linprog

def pal(V,edges,steps):
    seq=[]
    for p in permutations(V,steps):
        arc=list(zip(p[:-1],p[1:]))
        if all(e in edges for e in arc):seq.append((p,[edges.index(e) for e in arc]))
    return seq
def receive(paths,X,k):
    a=np.ones(len(paths))
    for i,(p,e) in enumerate(paths):
        if k not in p: a[i]=0
        else:
            for h in range(p.index(k)):
                if not X[e[h]]: a[i]=0;break
    return a
def allstate(V,E,p,X,w,lab):
    paths=pal(V,E,3)
    assert np.isclose(sum(w),1)
    A=np.array([tuple([1]+list(x)) for x in X],dtype=float).T
    Y=np.array([sum(receive(paths,x,k) for k in V) for x in X])
    sol=linprog(Y@np.array(w),A_eq=A,b_eq=[1]+[float(p)]*len(E),bounds=(0,None),method='highs')
    assert sol.success
    return F(sol.fun).limit_denominator(100000)
def allcut(V,E):
    paths=pal(V,E,3)
    if not paths:return None
    K=len(paths)
    C=[]
    for x in product([0,1],repeat=len(E)):
        # A state and its complement at .5/.5 are the extreme cut representations.
        X=np.array([x,tuple(1-v for v in x)])
        a=np.zeros(K)
        for y in X:
            a+=sum(receive(paths,y,k) for k in V)/2
        C.append(list(-a)+[1])
    sol=linprog([0]*K+[-1],A_ub=C,b_ub=[0]*len(C),A_eq=[[1]*K+[0]],b_eq=[1],
                bounds=[(0,None)]*K+[(None,None)],method='highs')
    assert sol.success
    return F(sol.x[-1]).limit_denominator(100000)
def check():
    small=(
      (list(range(3)),[(0,1),(1,2),(2,0)],F(5,3)),
      (list(range(5)),[(0,1),(1,2),(2,3),(3,4),(1,3)],F(8,5)),
      (list(range(4)),[(0,1),(1,2),(2,3),(3,0)],F(3,2)),
      (list(range(4)),[(0,1),(1,2),(2,3)],F(3,2)))
    P=[]
    for V,E,z in small:
        paths=pal(V,E,3);X=list(product([0,1],repeat=len(E)))
        a=allstate(V,E,F(1,2),X,[F(1,len(paths))]*len(paths),'expected')
        b=allcut(V,E)
        assert a==b==z
        print('separate all-state and all-cut LP certificate',len(paths),len(X),a)
    V,E,_=small[0]; X=list(product([0,1],repeat=len(E)))
    checks=0
    for i in range(61):
        p=F(i,60)
        expected=1+p+max(0,p-F(1,3),2*p-1)
        assert allstate(V,E,p,X,[F(1,3)]*3,'triangle sweep')==expected
        checks+=1
    print('fresh p grid independently checked',checks)
    paths=pal(V,E,3)
    allplan_count=len(paths)+1
    # Last column is an optional atomic 3-cycle. Shorter histories only reduce
    # the realized receipt set relative to extending the same start.
    Y=np.array([list(sum(receive(paths,x,k) for k in V))+[3*all(x)] for x in X])
    I={k:np.array([list(receive(paths,x,k))+[all(x)] for x in X]) for k in V}
    cs=[]; cb=[]
    for j,y in enumerate(X):
        if y[0]:continue
        k=X.index(tuple(1-v for v in y))
        cs.append(list(-(Y[j]+Y[k])/2)+[1]);cb.append(0)
        for i,b in enumerate([F(3,5),F(7,20),F(7,20)]):
            cs.append(list(-(I[i][j]+I[i][k])/2)+[0]);cb.append(-float(b))
    sol=linprog([0]*allplan_count+[-1],A_ub=cs,b_ub=cb,
                 A_eq=[[1]*allplan_count+[0]],b_eq=[1],
                 bounds=[(0,None)]*allplan_count+[(None,None)],method='highs')
    assert sol.success
    assert [F(x).limit_denominator(10000) for x in sol.x]==[F(7,15),F(4,15),F(4,15),F(0),F(49,30)]
    print('priority full-plan LP including optional atomic cycle:',[F(x).limit_denominator(10000) for x in sol.x])

if __name__=='__main__':check()
