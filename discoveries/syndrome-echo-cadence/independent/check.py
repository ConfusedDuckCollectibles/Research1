#!/usr/bin/env python3
"""Exact independent state-vector replay, Python 3.10+; no dependencies."""
from collections import Counter
from fractions import Fraction
from math import prod

P = (3, 4, 5), (5, 12, 13), (15, -8, 17), (24, 7, 25), (20, 21, 29)

def masks(n):
    return range(1 << n)

def subsets(x):
    i = x
    while True:
        yield i
        if i == 0:
            return
        i = (i-1) & x

def sign(x):
    return (-1) ** x.bit_count()

def edge(i,j):
    return (1<<i)|(1<<j)

def syndrome(x, edges):
    return sum(((x & e).bit_count() % 2) << a for a,e in enumerate(edges))

def sqrt_norm(v):
    return v[0]**2 + v[1]**2

def local_rotation(v, site, c, s):
    # Integer numerator of (c I - i s X) |psi>. All sites normalize separately.
    w = []
    for j, (x, y) in enumerate(v):
        u, t = v[j^(1<<site)]
        w.append((c*x+s*t, c*y-s*u))
    return w

def pauli_Z(v, z):
    return [((sign(j&z)*x), (sign(j&z)*y)) for j,(x,y) in enumerate(v)]

def rotate(v, slots):
    for j, (c,s,d) in enumerate(slots):
        assert c*c+s*s == d*d
        v = local_rotation(v,j,c,s)
    return v

def initial_state(n,init):
    u,v,d = init
    a = [(0,0)] * (1<<n)
    a[0],a[-1] = u,v
    assert sqrt_norm(u) + sqrt_norm(v) == d*d
    return a

def exact_joint(n, mid, last, s1, s2, z, init=((1,0),(0,0),1)):
    first = rotate(initial_state(n, init), s1)
    counts = Counter()
    for m in {syndrome(j,mid) for j in masks(n)}:
        prj = [v if syndrome(j,mid)==m else (0,0) for j,v in enumerate(first)]
        # Conjugate by actual diagonal Z before and after the second slot.
        end = pauli_Z(rotate(pauli_Z(prj,z),s2),z)
        for j,v in enumerate(end):
            counts[m,syndrome(j,last)] += sqrt_norm(v)
    d2 = (init[-1] * prod(d for c,s,d in (*s1,*s2))) ** 2
    assert sum(counts.values()) == d2
    return {k: Fraction(v,d2) for k,v in counts.items()}

def twirl(n,mid,last,s1,s2):
    def weight(word,slot):
        return prod(s**2 if (word>>j)&1 else c**2 for j,(c,s,d) in enumerate(slot))
    counts = Counter()
    for u in masks(n):
        for v in masks(n):
            counts[syndrome(u,mid),syndrome(u^v,last)] += weight(u,s1)*weight(v,s2)
    d2 = prod(d*d for c,s,d in (*s1,*s2))
    return {k:Fraction(v,d2) for k,v in counts.items()}

def char_product(mask, edges):
    a = 0
    for j,e in enumerate(edges):
        if (mask>>j)&1:
            a ^= e
    return a

def supports(edges):
    return {char_product(i,edges):i for i in masks(len(edges))}

def opp(x,y):
    return (sign(x&y))

def expected_char(mid, A, B, s1, s2, z):
    trig = lambda c,s,d: (Fraction(c*c-s*s,d*d), Fraction(2*c*s,d*d))
    l1 = [trig(*t) for t in s1]
    l2 = [trig(*t) for t in s2]
    C = A ^ B
    ans = 0
    for D in subsets(A&~B):
        if syndrome(D,mid) != 0:
            continue
        ans += sign(D)*sign(D&z)*prod(l1[j][bool(D&(1<<j))] for j in range(len(s1)) if C&(1<<j))*prod(l2[j][bool(D&(1<<j))] for j in range(len(s2)) if A&(1<<j))
    return ans

def exact_moment(law, bmask, amask):
    return sum(v*opp(m,bmask)*opp(a,amask) for (m,a),v in law.items())

def push_readout(law, edge_list, num=1, den=7):
    # Actual independent flips on nonredundant final chain digit labels.
    k = len(edge_list)
    out=Counter()
    for (m,a), v in law.items():
        for e in masks(k):
            out[m,a^e] += v * Fraction(num**e.bit_count()*(den-num)**(k-e.bit_count()),den**k)
    return out

def marginal(law,edges,err):
    out = syndrome(err,edges)
    return sum(v for (m,a),v in law.items() if a==out)

def all_checked_replay():
    trials = chars = 0
    for n,last in ((3,(edge(0,1),edge(1,2),edge(0,2))), (5,(edge(0,1),edge(1,2),edge(2,3),edge(3,4),edge(4,0)))):
        ss1 = P[:n]
        ss2 = (P[1:]+P[:1])[:n]
        last_support = supports(last)
        for f in masks(len(last)):
            mid = tuple(e for j,e in enumerate(last) if (f>>j)&1)
            mid_support = supports(mid)
            tw = twirl(n,mid,last,ss1,ss2)
            avg = Counter()
            for z in masks(n):
                l = exact_joint(n,mid,last,ss1,ss2,z)
                h = exact_joint(n,mid,last,ss1,ss2,z,((3,0),(0,4),5))
                g = exact_joint(n,mid,last,ss1,ss2,z,((3,0),(4,0),5))
                assert l==h==g
                if len(mid_support)==len(last_support):
                    assert l==tw
                for A,a in last_support.items():
                    for B,b in mid_support.items():
                        assert exact_moment(l,b,a) == expected_char(mid,A,B,ss1,ss2,z)
                        chars += 1
                avg.update(l)
                trials += 1
            assert {k:v/(1<<n) for k,v in avg.items()} ==tw
    return trials, chars

def local_and_rare():
    l = (edge(0,1),edge(1,2))
    p=Fraction(16,25)
    slot=(P[0],(1,0,1),(1,0,1))
    a=exact_joint(3,(l[1],),l,slot,slot,0)
    b=exact_joint(3,(l[1],),l,slot,slot,1)
    w=marginal(a,l,1)-marginal(b,l,1)
    assert w==4*p*(1-p)==Fraction(576,625)
    assert marginal(a,l,1)==Fraction(576,625)
    assert marginal(b,l,1)==0
    assert exact_joint(3,l,l,slot,slot,0)==exact_joint(3,l,l,slot,slot,1)
    ar=push_readout(a,l)
    br=push_readout(b,l)
    parity=lambda law: sum(v for (m,e),v in law.items() if e&1)
    assert parity(ar)-parity(br)==Fraction(5,7)*w
    print(f"single endpoint: true no-echo minus echo gap={w}; after readout q=1/7, gap={parity(ar)-parity(br)}")
    for r in (2,3,4):
        n=2*r+1
        l=tuple(edge(i,i+1) for i in range(n-1))
        mid=l[:r-1]+l[r:]
        slot=P[:r]+((1,0,1),)*(n-r)
        err=(1<<r)-1
        a=exact_joint(n,mid,l,slot,slot,0)
        b=exact_joint(n,mid,l,slot,slot,1)
        exact=2**(r+1)*prod(Fraction(c*c*s*s,d**4) for c,s,d in slot[:r])
        assert marginal(a,l,err)==exact
        assert marginal(b,l,err)==0
        print(f"matched component r={r}: no-echo hit={exact}; echoed hit=0")

def negative_controls():
    n=4
    l=tuple(edge(i,i+1) for i in range(n-1))
    slots=(P[0],)*n
    even_a=exact_joint(n,l,l,slots,slots,0)
    even_b=exact_joint(n,l,l,slots,slots,1)
    assert even_a!=even_b
    tv=sum(abs(even_a[k]-even_b[k]) for k in even_a)/2
    print(f"even n=4, fully checked between slots: exact echo total variation={tv}; not an odd-block theorem")
    n=3
    l=tuple(edge(i,i+1) for i in range(n-1))
    slots=(P[0],(1,0,1),(1,0,1))
    a=exact_joint(n,(l[1],),l,slots,slots,0)
    b=exact_joint(n,(l[1],),l,slots,slots,1)
    erase=lambda law: sum(law.values())
    assert a!=b and erase(a)==erase(b)==1
    print("singular erasing readout removes real signal, as expected")

if __name__ == '__main__':
    trials,chars=all_checked_replay()
    print(f"Independent computational-vector exact replay: {trials} instruments, {chars} characters; coherent initial real/imaginary tests and twirl checks pass")
    local_and_rare()
    negative_controls()
