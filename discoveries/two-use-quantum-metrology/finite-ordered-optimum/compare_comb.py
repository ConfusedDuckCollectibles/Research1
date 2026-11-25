"""Independent author 3-use comb data smoke test, NOT a certified optimizer."""
from argparse import ArgumentParser
from hashlib import sha256
from math import sqrt
from pathlib import Path
from urllib.request import Request, urlopen
from scipy.optimize import brentq, minimize, minimize_scalar

# The authors' plot.py uses uniform decay p in [0,1], 201 points, and
# the actual AD model adds exp(-i theta Z/2), differing from exp(i theta P)
# by only a global phase. Phase and axial AD commute at the channel level.
COMMIT = "07748a4e94d7b29dc888d2cbf41d0b0f421c58dc"
URL = ("https://raw.githubusercontent.com/qiushi-liu/strategies_in_metrology/"
       +COMMIT+"/data/QFI_seq_rangep_0_1_N_3.txt")
SHA = "756927f516d70321fbb454538b76d1a1ac95147f8964424768591c134590b4e0"


def final(q, s):
    d = sqrt(q)
    t = brentq(lambda t: s-t**4-d*(1-s)*t*(1+t*t), 0, 1, xtol=3e-15)
    x = 1/(1+t*t)
    return s*(d+sqrt(x*(1-x)))**2/(1-(1-s)*x)


def middle_value(x, a, q, s):
    if x <= 0 or x >= 1:
        return final(q, s)
    W = 1-(1-s)*x
    z = a*sqrt(x*(1-x))
    T = q+(s*x*(1-x)+2*s*sqrt(q)*z-(1-s)*q*z*z/x)/W
    R = q*(1-a*a)/x
    b = min(T,max(0,R))
    return W*final((T-(1-s)*x*b)/W, s)+(1-W)*final(b, s)


def middle(q, s):
    # Very small deterministic multistart, independent of reference values.
    starts = sorted(((middle_value(x,a,q,s),x,a) for x in (0.25,0.55,0.80,0.95)
                     for a in (0.20,0.55,0.88,0.97)), reverse=True)
    sols = [minimize(lambda xa: -middle_value(*xa,q,s), [x,a],
                     method="Nelder-Mead", bounds=((1e-8,1-1e-8),(0,1)),
                     options=dict(xatol=1e-10,fatol=1e-13,maxiter=200))
            for _,x,a in starts[:3]]
    assert any(sol.success for sol in sols)
    return max(-sol.fun for sol in sols)


def first(s):
    loss_tail = middle(0,s)
    def val(x):
        W = 1-(1-s)*x
        return W*middle(s*x*(1-x)/(W*W),s)+(1-W)*loss_tail
    sol = minimize_scalar(lambda x: -val(x), method="bounded",
                          bounds=(0,1), options=dict(xatol=1e-10))
    assert sol.success
    return 4*max(-sol.fun, val(0), val(1))


def main():
    arg = ArgumentParser(description=__doc__)
    group = arg.add_mutually_exclusive_group(required=True)
    group.add_argument("--fetch", action="store_true",
                       help="Fetch pinned original anonymously into memory, no local writes")
    group.add_argument("--source", type=Path,
                       help="Read an independently downloaded original data file")
    cli = arg.parse_args()
    if cli.fetch:
        req = Request(URL, headers={"User-Agent": "Read-only academic source replay"})
        with urlopen(req, timeout=20) as resp:
            raw = resp.read(100001)
            assert resp.status == 200 and 0 < len(raw) <= 100000
    else:
        raw = cli.source.read_bytes()
    assert sha256(raw).hexdigest() == SHA
    arr = [float(line) for line in raw.splitlines()]
    assert len(arr) == 201
    print("p=decay | three Bellman slots | author comb | absolute difference")
    for i in (20,80,160):
        p=i/200
        ans = first(1-p)
        err = abs(ans-arr[i])
        assert err < 5e-6,(p,ans,arr[i])
        print(f"{p:.2f}   {ans:.10f}   {arr[i]:.10f}   {err:.3g}")
    print("PASS (numerical comparison only)")

if __name__ == "__main__":
    main()
