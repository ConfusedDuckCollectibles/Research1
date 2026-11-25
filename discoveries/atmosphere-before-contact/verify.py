"""Rational certificate; --native also refetches two public IRSA calibrations."""
import csv
from fractions import Fraction as F
from hashlib import sha256
from io import StringIO
import sys
from urllib.request import urlopen

SRC_URL = ("https://irsa.ipac.caltech.edu/TAP/sync?"
           "QUERY=SELECT+TOP+8+designation,ra,dec,w2mpro,w2sigmpro,w2snr,cc_flags,ph_qual,nb,ext_flg"
           "+FROM+allwise_p3as_psd+WHERE+CONTAINS%28POINT%28%27ICRS%27,ra,dec%29,"
           "CIRCLE%28%27ICRS%27,346.6265,-5.0435,0.006%29%29=1&FORMAT=csv")
RSR_URL = "https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/figures/RSR-W2.txt"
SRC_SHA = "6a98146f7c226d20de6aa8c639459d60ab8f4ea811da97372c49490e4cc7a470"
RSR_SHA = "b0a86a7e86b07713cecfe9c853292d56f19da6c385de45f6265b676cc1330f2a"


def read_url(url, digest):
    with urlopen(url, timeout=20) as response:
        data = response.read(250000)
    assert sha256(data).hexdigest() == digest
    return data.decode("ascii")


def native(source, response):
    rows = list(csv.DictReader(StringIO(source)))
    assert len(rows) == 1 and rows[0]["designation"] == "J230630.02-050234.1"
    row = rows[0]
    assert F(row["w2mpro"]) == F("9.807") and F(row["w2sigmpro"]) == F(".020")
    assert row["cc_flags"] == "0000" and row["ph_qual"] == "AAAU"
    assert row["nb"] == "1" and row["ext_flg"] == "0"
    assert F(row["w2mpro"]) - 5 * F(row["w2sigmpro"]) > F("9.70")
    points = [tuple(map(F, line.split())) for line in response.splitlines()
              if line.strip() and not line.startswith("#")]
    assert len(points) == 548 and (points[0][0], points[-1][0]) == (F("2.53"), F("8.00"))
    assert all(len(p) == 3 and 0 <= p[1] <= 1 and 0 <= p[2] for p in points)
    assert all(b[0] - a[0] == F(".01") for a, b in zip(points, points[1:]))
    lo_pb = min(v * (1 - 4*e/1000) for x, v, e in points if F("4.10") <= x <= F("4.50"))
    def upper(point):
        _, v, e = point
        return v * (1 + 4*e/1000)
    hi_int = sum(max(upper(a), upper(b)) * (b[0] - a[0]) / a[0]
                 for a, b in zip(points, points[1:]))
    assert lo_pb > F(".45654") > F(".45")
    assert hi_int < F(".179187") < F(".185")
    print("native IRSA: source identity and magnitudes; 548 response knots; exact band/integral envelopes verified")


def atan_interval(x, n):
    partial = sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(n))
    following = (-1)**n*x**(2*n+1)/(2*n+1)
    return (min(partial, partial+following), max(partial, partial+following))


def log_interval(x, n=128):
    z = (x - 1) / (x + 1)
    partial = 2*sum(z**(2*j+1)/(2*j+1) for j in range(n))
    error = 2 * abs(z)**(2*n+1) / ((2*n+1)*(1-z*z))
    return partial-error, partial+error


def d_bern_lo(p, q):
    return p*log_interval(p/q)[0] + (1-p)*log_interval((1-p)/(1-q))[0]


def main():
    if sys.argv[1:] == ["--native"]:
        native(read_url(SRC_URL, SRC_SHA), read_url(RSR_URL, RSR_SHA))
    elif len(sys.argv) == 4 and sys.argv[1] == "--files":
        source, response = (open(path, "rb").read() for path in sys.argv[2:])
        assert sha256(source).hexdigest() == SRC_SHA
        assert sha256(response).hexdigest() == RSR_SHA
        native(source.decode("ascii"), response.decode("ascii"))
    else:
        assert not sys.argv[1:], "usage: verify.py [--native | --files source.csv W2.txt]"
    assert 10**97 > 7500**25
    pi_lo = 16*atan_interval(F(1,5), 8)[0] - 4*atan_interval(F(1,239), 8)[1]  # Machin
    assert pi_lo > F("3.1415")
    pi_simple = F("3.1415")
    phi_quiet = 26 * 200 * F("1e-26") * F(".185")/(7500*F(".45")*F("6.62607015e-34"))
    assert phi_quiet < 4_302_000
    phi = 2 * phi_quiet
    assert phi < 8_604_000
    assert 2 * (2 * F(".01") + F(".01")**2) < F(".041")  # p<=1, L<=2
    k = phi * 2**2 * F(".0001")**2 * F(".01")**2 / (pi_simple**2 * F(".068") * F("1.068") * F(".00045") * F(".959"))
    assert 10*k < F("1.113")
    assert 17*k < F("1.892")
    assert d_bern_lo(F(".8"), F(".05")) > F("1.9064")
    assert d_bern_lo(F(".591"), F(".05")) > F("1.1149")
    assert 17*k < d_bern_lo(F(".8"), F(".05"))
    assert 10*k < d_bern_lo(F(".591"), F(".05"))
    print("10-visit information ceiling < 1.113 nats; 17-visit ceiling < 1.892 nats")
    print("size 0.05: power < 0.591 at 10 visits; >=18 visits necessary for power 0.8")
    print("all rational certificate checks passed")


if __name__ == "__main__":
    main()
