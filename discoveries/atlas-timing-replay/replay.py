#!/usr/bin/env python3
"""Replay a versioned 3I/ATLAS timing/photometry audit. Publishes no input data."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from hashlib import sha256
from html import unescape
import json
import re
from urllib.request import Request, urlopen

import numpy as np
from scipy.optimize import minimize_scalar

PAPER1 = "https://arxiv.org/html/2601.10860v1"
PAPER2 = "https://arxiv.org/html/2602.18512v2"
JD_UNIX = 2440587.5


def fetch(url, limit=2_000_000, partial=False):
    headers = {"User-Agent": "AtlasMethodsReplay/1.0 (versioned scientific reproducibility)"}
    if partial:
        headers["Range"] = f"bytes=0-{limit-1}"
    with urlopen(Request(url, headers=headers), timeout=25) as r:
        if partial and r.status != 206:
            raise ValueError(f"Expected HTTP 206 for FITS header range: {url}")
        return r.read(limit)


def rows(doc, table_id, columns):
    table = re.search(fr'<table id="{re.escape(table_id)}"[^>]*>([\s\S]*?)</table>', doc)
    if not table:
        raise ValueError(f"Missing stable arXiv HTML table ID: {table_id}")
    rr = re.findall(r'<tr[^>]*>([\s\S]*?)</tr>', table.group(1))
    vals = [re.findall(r'<td[^>]*>([\s\S]*?)</td>', r) for r in rr]
    return [[unescape(re.sub(r'<[^>]*>', '', c)).strip() for c in row] for row in vals if len(row) == columns]


def fit(period, t, y, sd, order, groups=None):
    """Weighted linear least squares at fixed period; each night gets a free mean."""
    B = [np.ones(len(t))] if groups is None else [groups == g for g in np.unique(groups)]
    ph = 2 * np.pi * (t-t[0]) / period
    X = np.column_stack(B + [f(k*ph) for k in range(1, order+1) for f in (np.cos, np.sin)])
    Xw, yw = X / sd[:, None], y / sd
    b = np.linalg.lstsq(Xw, yw, rcond=None)[0]
    err = yw - Xw @ b
    return float(err @ err), list(b[-2*order:]) if order else []


def profile(t, y, sd, order, groups=None, lo=6.8, hi=7.4, points=6001):
    score = lambda P: fit(P, t, y, sd, order, groups)[0]
    p = np.linspace(lo, hi, points)
    es = np.array([score(P) for P in p])
    ids = np.flatnonzero((es[1:-1] < es[:-2]) & (es[1:-1] <= es[2:]))+1
    peaks = []
    for i in ids:
        q = minimize_scalar(score, bounds=(p[i-1], p[i+1]), method="bounded",
                            options={"xatol": 2e-11})
        peaks.append((float(q.x), float(q.fun)))
    return sorted(peaks, key=lambda q: q[1])


def fit_local(t, y, sd, order, groups=None, target=7.136):
    peaks = profile(t, y, sd, order, groups)
    nearest = min(peaks, key=lambda v: abs(v[0]-target))
    return dict(best=peaks[0], nearest_target=nearest, top5=peaks[:5],
                chi2_at_target=fit(target, t, y, sd, order, groups)[0])


def iso(s):
    return datetime.fromisoformat(s).replace(tzinfo=timezone.utc)


def fractional_date(s):
    yr, mo, dd = map(float, s.split())
    return datetime(int(yr), int(mo), int(dd), tzinfo=timezone.utc) + timedelta(days=dd-int(dd))


def at_jd(dt):
    return dt.timestamp()/86400 + JD_UNIX


def hst_header(root):
    key = root.lower()
    url = f"https://stpubdata.s3.amazonaws.com/hst/public/{key[:4]}/{key}/{key}_flt.fits"
    hdr = fetch(url, limit=14_400, partial=True)
    if not hdr.startswith(b"SIMPLE"):
        raise ValueError(f"Invalid FITS header for {root}")
    c = {}
    for i in range(0, len(hdr), 80):
        card = hdr[i:i+80].decode("ascii", errors="replace")
        fld = card[:8].strip()
        if fld == "END":
            break
        if card[8:10] == "= ":
            c[fld] = card[10:].split(" / ", 1)[0].strip().strip("'").strip()
    return (root, iso(c["DATE-OBS"]+" "+c["TIME-OBS"]), float(c["EXPSTART"]), float(c["EXPEND"]))


def frame_set(doc, table_id, first_four):
    v = np.array([[float(r[i]) for i in [0, 1, 3, 4]] for r in rows(doc, table_id, 7)])
    jd, mag, snr, zp = v.T
    t = (jd - jd[0])*24
    night = np.cumsum(np.r_[False, np.diff(t)>6])
    keep = night < 4 if first_four else np.ones(len(night), dtype=bool)
    return jd[keep], mag[keep], snr[keep], zp[keep], night[keep]


def rebin(jd, mag, snr, zp, night):
    t = (jd-jd[0])*24
    y = mag + np.median(zp)-zp
    w = (snr/1.0857)**2
    columns = []
    for n in np.unique(night):
        ti = t[night == n]
        ix = np.floor(2*(ti-ti[0])+1e-8)
        for m in np.unique(ix):
            full = np.flatnonzero(night == n)[ix == m]
            # Average relative JDs for accuracy, add origin after.
            columns.append((np.average(t[full], weights=w[full]), np.average(y[full], weights=w[full]),
                            float(np.sum(w[full])**-0.5), n))
    return np.array(columns).T


def main():
    with ThreadPoolExecutor(2) as ex:
        p1, p2 = list(ex.map(fetch, [PAPER1, PAPER2]))
    d1, d2 = p1.decode("utf-8"), p2.decode("utf-8")
    r1, r3 = rows(d1, "S2.T1.2", 2), rows(d1, "S2.T3.2", 5)
    newer, expo = rows(d2, "S2.T1.6", 7), rows(d2, "S2.T2.2", 9)
    diffr = rows(d2, "S2.T4.2", 6)
    assert len(r1) == len(newer) == len(expo) == 21 and len(r3) == len(diffr) == 14
    assert all(float(a[1]) == float(b[2]) and fractional_date(a[0]) == fractional_date(b[0])
               for a, b in zip(r1, newer))
    g = np.array([a[0][:10] for a in r3])
    jd, b, s = np.array([[float(a[i]) for i in (1,3,4)] for a in r3]).T
    t = (jd-jd[0])*24
    j0, m0, sr0, z0, n0 = frame_set(d2, "S6.T9.2", True)
    assert len(j0) == 286
    rt, rb, rs, rn = rebin(j0,m0,sr0,z0,n0)
    dts = float(np.max(np.abs((rt + (j0[0]-jd[0])*24)-t))*3600)
    dms = float(np.max(np.abs(rb-b)))
    dss = float(np.max(np.abs(rs-s)))
    assert len(rt) == 14 and dts < .0043 and dms < .00043 and dss < .00049
    p = {
        "source_versions": [PAPER1,PAPER2],
        "source_html_sha256": [sha256(p1).hexdigest(),sha256(p2).hexdigest()],
        "recovered_input_max_rounding_difference": {"seconds":dts,"magnitudes":dms,"formal_sigma_mag":dss},
    }
    p["original_14_bin_fit"]={}
    for k in (1,2):
        for name,err in (("jitter",np.hypot(s,.089)),("formal_only",s)):
            p["original_14_bin_fit"][f"K{k}_{name}"]=fit_local(t,b,err,k,g)
    # The later five-night product is a different data reduction. Dec 28 UT at the tail
    # is the continuation of the observing visit that began Dec 27, not a sixth visit.
    tx,yx,sx=np.array([[float(a[i]) for i in (1,3,4)] for a in diffr]).T
    nx=np.array([a[0][:10].replace("2025-12-28","2025-12-27") for a in diffr])
    assert len(np.unique(nx))==5
    p["followup_v2_different_14_bin_5_night_fit"] = {
        "formal_sem":fit_local((tx-tx[0])*24,yx,sx,1,nx),
        "sem_plus_jitter":fit_local((tx-tx[0])*24,yx,np.hypot(sx,.089),1,nx),
    }
    p["followup_v2_large_aperture_control"]={}
    for four,n_expected in ((True,285),(False,398)):
        jd1,m1,snr1,z1,n1=frame_set(d2,"S6.T10.2",four)
        assert len(jd1)==n_expected
        tj=(jd1-jd1[0])*24; ej=1.0857/snr1; zj=np.median(z1)-z1
        best=profile(tj,m1+zj,ej,1,n1)[0]
        fitm=lambda y:fit(best[0],tj,y,ej,1,n1)[1]
        p["followup_v2_large_aperture_control"][str(n_expected)]={
            "best_local_hours_and_chi2":best,
            "amplitude_tabulated":float(np.linalg.norm(fitm(m1))),
            "amplitude_zeropoint_term":float(np.linalg.norm(fitm(zj))),
            "amplitude_corrected":float(np.linalg.norm(fitm(m1+zj))),
        }

    with ThreadPoolExecutor(5) as ex:
        hst=list(ex.map(hst_header,[r[0] for r in expo]))
    assert all(h[1] == iso(a[6]) and h[0]==a[0] for h,a in zip(hst,expo))
    tt0=[fractional_date(a[0]) for a in r1]
    tf=np.array([(q-tt0[0]).total_seconds()/3600 for q in tt0])
    te=np.array([float(a[-1]) for a in newer])
    tm=np.array([.5*(a[2]+a[3]) for a in hst]); tm=(tm-tm[0])*24
    ya=np.array([float(a[1]) for a in r1]); sd=np.full(len(ya),3.)
    p["jet_fit"]={
        "published_fractional_UT":fit_local(tf,ya,sd,2,target=7.2),
        "followup_elapsed_from_its_own_origin":fit_local(te,ya,sd,2,target=7.2),
        "native_midtimes":fit_local(tm,ya,sd,2,target=7.2),
        "native_drop_first_angle":fit_local(tm[1:],ya[1:],sd[1:],2,target=7.2),
        "exposure_id":hst[0][0],
        "published_first_UTC":tt0[0].isoformat(),
        "native_first_start_UTC":hst[0][1].isoformat(),
        "start_minus_fractional_seconds":(hst[0][1]-tt0[0]).total_seconds(),
        "followup_first_elapsed_interval_hours":float(te[1]-te[0]),
        "native_first_mid_to_mid_interval_hours":float(tm[1]-tm[0]),
    }
    # Runtime cross-checks prevent a changed parse or input from silently earning a pass.
    assert abs(p["original_14_bin_fit"]["K1_jitter"]["chi2_at_target"]-51.322)<.01
    assert abs(p["original_14_bin_fit"]["K2_jitter"]["best"][0]-6.88416)<1e-4
    assert abs(p["followup_v2_different_14_bin_5_night_fit"]["formal_sem"]["nearest_target"][0]-7.135967)<1e-4
    assert abs(p["jet_fit"]["start_minus_fractional_seconds"]-4316.808)<.001
    assert abs(p["jet_fit"]["native_midtimes"]["best"][0]-7.371635)<1e-4
    print(json.dumps(p,indent=2))


if __name__ == "__main__":
    main()
