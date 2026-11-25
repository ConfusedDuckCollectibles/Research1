#!/usr/bin/env python3
"""Independent exact source and benchmark replay; no third-party dependencies."""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import argparse

RSR_SHA="b0a86a7e86b07713cecfe9c853292d56f19da6c385de45f6265b676cc1330f2a"

def native(rsr,gator):
    data=Path(rsr).read_bytes()
    assert sha256(data).hexdigest()==RSR_SHA
    vals=[tuple(map(F,s.split())) for s in data.decode("ascii").splitlines()
          if s.strip() and not s.startswith("#")]
    assert len(vals)==548
    assert [v[0] for v in vals] == [F(j,100) for j in range(253,801)]
    lo=min(r*(1-4*sd/1000) for x,r,sd in vals if F("4.10")<=x<=F("4.50"))
    h=[(x,r*(1+4*sd/1000)) for x,r,sd in vals]
    # Fractional trapezoid for linearly interpolated upper response,
    # with rigorous 1/x replacement at the left endpoint.
    ii=sum((b-a)*(r+s)/(2*a) for (a,r),(b,s) in zip(h,h[1:]))
    assert lo > F(".45654") > F(".45")
    assert ii < F(".179") < F(".185")
    print("native RSR: 548 ratios, lower", float(lo), "rigorous integral upper",float(ii))
    raw=Path(gator).read_text()
    row=[s.split() for s in raw.splitlines() if s.lstrip().startswith("J230630.02-050234.1")]
    assert len(row)==1
    d,ra,dc,clon,clat,w2,dw2,snr,art,q,m,n,*_=row[0]
    assert (d,F(ra),F(dc),F(w2),F(dw2),F(snr),art,q,m,n)==(
        "J230630.02-050234.1",F("346.625105"),F("-5.0428184"),
        F("9.807"),F(".020"),F("54.6"),"0000","AAAU","21","21")
    assert "AllWISE Source Catalog (allwise_p3as_psd)" in raw
    print("independent IRSA Gator source and 21 covered/detected W2 frames match exactly")

def ln_lo(x, n=128):
    t=(x-1)/(x+1)
    p=2*sum(t**(2*j+1)/(2*j+1) for j in range(n))
    e=2*abs(t)**(2*n+1)/((2*n+1)*(1-t*t))
    return p-e

def bi_lo(beta,size):
    return beta*ln_lo(beta/size)+(1-beta)*ln_lo((1-beta)/(1-size))

def atan(x,n):
    s=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(n))
    other=s+(-1)**n*x**(2*n+1)/(2*n+1)
    return min(s,other),max(s,other)

def benchmark():
    # 9.7 Vega mag means 10^(-97/25); the comparisons are integer exact.
    assert 7500**25 < 10**97
    h=F("6.62607015e-34")
    quiet=26*200*F("1e-26")*F(".185")/(7500*F(".45")*h)
    assert quiet < 4_302_000
    phi=2*4_302_000
    p,H,u,d,L=map(F,[".068",".01",".00045",".0001","2"])
    pp=16*atan(F(1,5),14)[0]-4*atan(F(1,239),14)[1]
    assert pp>F("3.14159265358")
    alpha=F(".05")
    div80=bi_lo(F(".8"),alpha)
    assert div80>F("1.90644")
    def kl(eps, H_effective):
        return phi*L*L*d*d*H_effective*H_effective/(pp*pp*p*(1+p)*u*(1-eps))
    base=kl(F(".041"),H)
    assert 10*base<F("1.113")<bi_lo(F(".591"),alpha)
    assert 17*base<F("1.892")<div80
    # Uniform sharper epsilon: proof ratio increases up to p=1.
    eps=F(".0013")
    assert p>H/(2*(1-H))
    assert (L*(2*H+H*H))**2*H < (eps*pp)**2
    # Area-coordinate shell width improvement, valid for weighted mark RMS.
    assert d <= p*H and d<2*p*H+H*H
    out=kl(eps,H-d/(4*(p+H)))
    assert 19*out<F("1.903")<div80
    assert 10*out<F("1.0011")<bi_lo(F(".556"),alpha)
    print("conditional ideal incident-rate cap: <",float(2*quiet))
    print("baseline: 10 KL <",float(10*base), "17 KL <",float(17*base), 
          "power <.591, >=18 visits necessary")
    print("optional refined: 10 KL <",float(10*out),"19 KL <",float(19*out),
          "power <.556, >=20 visits necessary")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--native-files", nargs=2, metavar=("RSR","GATOR"))
    args=p.parse_args()
    if args.native_files: native(*args.native_files)
    benchmark()
    print("PASS")
