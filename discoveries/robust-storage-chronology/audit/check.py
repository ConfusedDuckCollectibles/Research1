#!/usr/bin/env python3
"""Independent Fraction-only microchronology coverage at integer breakpoints.

With rho=ln(2), output weighted deficits multiplied by rho so every
microstep integration has integer weight 2**k. The accompanying proof,
not sampling, reduces the continuous model to finite cuts.
"""
from fractions import Fraction as F
from itertools import accumulate, product

LO, HI, PC = -1, 8, 6
ETA_C, ETA_D = F(3,4), F(1,2)
LEVELS = (-1, 0, 2, 3, 5, 6, 8)

def flow(x):
    return ETA_C * min(max(x, 0), PC) + min(x, 0) / ETA_D

def brute(widths, heights, q):
    bybin = []
    for width, h in zip(widths, heights):
        total = width*LO+(HI-LO)*h
        bybin.append([p for p in product(LEVELS, repeat=width)
                      if sum(p) == total])
    n = sum(widths)
    best = {(i, k): -10**9 for i in range(n) for k in range(n+1)}
    witness = {}
    for ps in product(*bybin):
        p = sum(ps, ())
        f = tuple(map(flow, p))
        for i in range(n):
            d = F(0)
            for k in range(n+1):
                if d > best[i, k]:
                    best[i, k], witness[i, k] = d, p
                d -= q**k * f[(i+k) % n]
    return best, witness, tuple(map(len, bybin))

def formula(widths, heights, i, k, q):
    # Independently allocate a bin's deficit mass on its largest weights;
    # each slot represents one exact integrated unit-interval weight.
    n = sum(widths)
    weight = { (i+off) % n: q**off for off in range(k) }
    end = list(accumulate(widths, initial=0))
    ans = F(0)
    for j, h in enumerate(heights):
        w = sorted((weight.get(t, 0) for t in range(end[j], end[j+1])),
                   reverse=True)
        nneg = widths[j] - h
        ans += -flow(LO) * sum(w[:nneg]) - flow(HI) * sum(w[nneg:])
    return ans

def vertices(widths, heights):
    n = sum(widths)
    bs = list(accumulate(widths, initial=0))
    starts = set(bs) | {s+h for s,h in zip(bs,heights)}
    ends = set(bs) | {t-h for t,h in zip(bs[1:],heights)}
    ends |= {t+n for t in list(ends)}
    lengths = {0,n}
    lengths |= {d-h for d,h in zip(widths,heights)}
    lengths |= {n-h for h in heights}
    pts = {(a,b-a) for a,b in product(starts,ends)}
    pts |= set(product(starts,lengths))
    pts |= {(b-k,k) for b,k in product(ends,lengths)}
    return {(a%n,k) for a,k in pts if 0<=a<=n and 0<=k<=n}

def periodic_greedy(p, cap, q):
    # Piecewise-constant continuous flow makes motion monotone within each
    # unit segment until cap. A period map is min(x/q**n+J, M).
    # Compose the cap-limited one-step maps to obtain exact scaled x0.
    n = len(p)
    J, M = F(0), cap
    for f in map(flow,p):
        A = F(1,q)
        w = (1-A)*f if q != 1 else f
        J, M = A*J+w, min(cap,A*M+w)
    x = min(J/(1-F(1,q)**n), M) if q != 1 else M
    assert min(x,cap) == x
    init = x
    for f in map(flow,p):
        assert 0 <= x <= cap
        A = F(1,q)
        w = (1-A)*f if q != 1 else f
        x = min(cap,A*x+w)
    assert x == init

def check(widths, heights, q):
    best, wit, counts = brute(widths, heights, q)
    for (i,k), d in best.items():
        assert d == formula(widths, heights, i, k, q)
    cand = vertices(widths, heights)
    cap = max(best.values())
    n = sum(widths)
    cy = max(best[i,n] for i in range(n))
    assert cap == max(best[pt] for pt in cand)
    assert cy == max(best[pt] for pt in cand if pt[1]==n)
    assert cap >= 0
    if cy <= 0:
        for p in set(wit.values()):
            periodic_greedy(p, cap, q)
    edge = {i % n for i in accumulate(widths,initial=0)}
    edge_cap = max(d for (i,k),d in best.items() if i in edge and (i+k)%n in edge)
    max_pts = [pt for pt,d in best.items() if d==cap]
    print(f"bins {widths} heights {heights}, base {q}: {counts} microtrace factors, "
          f"cap={cap}, max cycle deficit={cy}, edges-only cap={edge_cap}, "
          f"max arcs={max_pts}")
    # Here cap is merely the one-period arc deficit if a cycle fails.
    return cap, cy, edge_cap

if __name__ == "__main__":
    assert check((4,), (3,), 2) == (2, F(-31,2), 0)
    assert check((4,), (2,), 2)[:2] == (F(21,2), F(21,2))
    assert check((4,), (2,), 1)[:2] == (4, -5)
    # Even the all-start feasibility test can fail at an *unobserved*
    # interior phase although each bin-edge phase passes.
    assert check((2,3), (1,2), 2) == (F(33,2), F(33,2), 0)
    b,w,_ = brute((2,3),(1,2),2)
    assert b[(0,5)] == -F(45,2) and b[(2,5)] == -F(19,2)
    assert w[(1,5)] == (-1,8,8,8,-1)
    # Independent replay of the author's symmetric-rate controls.
    q = F(11,10)
    traces = [p for p in product((-1,1),repeat=6)
              if sum(p[:3]) == sum(p[3:]) == 1]
    vals = {(a,k): max(-(q-1)*sum(q**i*p[(a+i)%6] for i in range(k))
                        for p in traces)
            for a in range(6) for k in range(7)}
    assert max(vals.values()) == F(21,100)
    assert max(vals[a,6] for a in range(6)) == -F(156639,10**6)
    assert vals[0,6] == vals[3,6] == -F(207459,10**6)
    r=F(1001,1000)
    assert 2*r**101-r**200-1 < 0
    # Repeating the same bin interior is a real hypothesis: when it
    # varies by cycle, a last and then a first deficit can touch.
    assert -flow(-1)*(1+2) == 6  # versus repeated-trace reserve 2.
    print("PASS: author's symmetric controls and repetition countercontrol")
    for widths, heights in [((4,2),(3,2)),
                            ((2,2,2),(1,2,1)), ((1,3,1),(1,2,1)),
                            ((2,3,1),(0,3,1)), ((2,3,1),(2,3,1))]:
        check(widths, heights, 2)
    print("PASS: exhaustive finite profile/arc maxima, candidate cuts, and feasible greedy witnesses")
