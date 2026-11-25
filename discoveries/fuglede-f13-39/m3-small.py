from itertools import combinations, product
from collections import Counter
p=13
T=list(combinations(range(p),3))
def apsteps(A):
    return {d for d in range(1,7) if any({c,(c+d)%p,(c-d)%p}==set(A) for c in A)}
steps={A:apsteps(A) for A in T}
ap=0
for A in T:
    for B in T:
        n=len({(a-b)%p for a in A for b in B})
        assert n>=5
        if n==5:
           assert steps[A]&steps[B]
           ap+=1
print('triple differences: exact tight pairs',ap)
for d in (1,2,3):
    hist=Counter()
    extreme=[]
    # For d >= 2, translation uniquely eliminates the next coefficient.
    # For d = 1, every function is a permutation and this normalization is harmless too.
    for a in range(1,p):
        for b in product(range(p),repeat=d-1):
            c=b+(0,a)
            v=[sum(y*x**j for j,y in enumerate(c))%p for x in range(p)]
            hist[sum(v)]+=1
            if sum(v)<=39:extreme.append(c)
    print('degree',d,'mass histogram',sorted(hist.items()),'extreme',extreme)
