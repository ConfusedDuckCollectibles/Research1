#!/usr/bin/env python3
"""Exact rational outer certificates for counts of recorded KM deaths.

Only published aggregate at-risk counts, rounded survival, and a total-count
cap are inputs. No figure coordinates, individual records or reconstructed data.
"""
from fractions import Fraction as Q
import re


def percent_band(printed):
    """Closed nearest-rounding band, clipped to [0,1]. Keep trailing zeros."""
    if not isinstance(printed, str) or not re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", printed):
        raise ValueError("pass a printed decimal percent string, without %")
    level = Q(printed) / 100
    if not 0 <= level <= 1:
        raise ValueError("printed percent outside [0,100]")
    dps = len(printed.partition(".")[2])
    half_unit = Q(1, 200 * 10**dps)
    return max(Q(0), level - half_unit), min(Q(1), level + half_unit)


def factor_bounds(start, end, deaths):
    """Sharp extrema at exact start/end risks; outer limits for risk bounds."""
    if not 1 <= end <= start or not 0 <= deaths <= start - end:
        raise ValueError("require a positive ending risk and feasible removals")
    return Q(end, end + deaths), Q(start - deaths, start)


def candidate_sets(n, landmarks, at_risk, total_cap=None, *, link=True):
    """Return cumulative-death outer sets at sorted *positive* landmarks.

    landmarks: ascending (time, rounded percent string) pairs.
    at_risk: {time: count}, time unit/origin identical to landmarks.
    total_cap: same-lock upper bound, e.g. all reported deaths; defaults to n.

    A risk count printed at the current landmark can describe either side of
    the step. The first strictly later count supplies a lower risk bound.
    """
    cap = n if total_cap is None else total_cap
    if type(n) is not int or n < 1 or type(cap) is not int or not 0 <= cap <= n:
        raise ValueError("invalid cohort size or death cap")
    anchors = sorted(at_risk.items())
    if any(t < 0 or type(k) is not int or not 0 <= k <= n for t, k in anchors):
        raise ValueError("invalid risk-table time or count")
    if any(x[1] < y[1] for x, y in zip(anchors, anchors[1:])):
        raise ValueError("risk table is not nonincreasing")
    prepared = []
    previous_time = None
    previous_band = (Q(1), Q(1))  # S(0-), not a printed S(0)
    for t, printed in landmarks:
        if t <= 0 or (previous_time is not None and t <= previous_time):
            raise ValueError("landmarks must be strictly increasing and positive")
        band = percent_band(printed)
        lower = next((k for v, k in anchors if v > t), 0)
        if lower <= 0:
            raise ValueError("each landmark needs a positive strictly later risk anchor")
        if link and previous_band[0] <= 0:
            raise ValueError("each previous survival band must have positive lower limit")
        # The first interval starts from N even if risk is printed after t=0.
        upper = n if previous_time is None else min(
            (k for v, k in anchors if v <= previous_time), default=n)
        prepared.append((band, previous_band, upper, lower))
        previous_time, previous_band = t, band

    result = []
    prev = frozenset({0})
    for (lo, hi), (prev_lo, prev_hi), upper, lower in prepared:
        possible = set()
        ratio_lo, ratio_hi = (lo / prev_hi, min(Q(1), hi / prev_lo)) if link else (None, None)
        for d in range(min(cap, n - lower) + 1):
            a, b = factor_bounds(n, lower, d)
            if a > hi or b < lo:
                continue
            if not link:
                possible.add(d)
                continue
            for p in prev:
                if 0 <= d - p <= upper - lower:
                    a, b = factor_bounds(upper, lower, d - p)
                    if a <= ratio_hi and b >= ratio_lo:
                        possible.add(d)
                        break
        prev = frozenset(possible)
        result.append(prev)
    return tuple(result)


def condensed(values):
    """Print contiguous candidate blocks without implying point precision."""
    values = sorted(values)
    if not values:
        return "none"
    blocks = [[values[0], values[0]]]
    for x in values[1:]:
        if x == blocks[-1][1] + 1:
            blocks[-1][1] = x
        else:
            blocks.append([x, x])
    return ",".join(str(x) if x == y else f"{x}-{y}" for x, y in blocks)


def run_examples():
    # PERSEUS: EMA Table 12 and Figure 6 (ITT, 1 August 2023).
    # Canada's Drug Agency Clinical Review Table 18 and Figure 4 corroborate
    # total counts and risk row. All inputs below remain aggregate numbers.
    specs = (
        ("PERSEUS", "D-VRd", 355, 34, 3, (12,24,36,48),
         ("97.7","96.3","94.2","89.4"),
         (355,347,343,341,338,335,331,329,329,326,325,323,321,316,312,284,135,21,1,0),
         ("8","13","20","22-34")),
        ("PERSEUS", "VRd", 354, 44, 3, (12,24,36,48),
         ("95.1","92.5","89.0","87.5"),
         (354,343,337,334,328,322,322,319,317,315,310,307,303,298,296,263,127,27,1,0),
         ("17","26","38","39-43")),
        # DESTINY-Breast09: EMA Table 30 and Figure 18 (FAS, 26 Feb 2025).
        ("DESTINY-Breast09", "T-DXd + pertuzumab", 383, 58, 3, (6,12,18,24,30),
         ("96.6","94.7","91.5","87.9","84.3"),
         (383,373,365,361,358,352,342,302,265,239,176,96,47,23,9,0),
         ("13","20","31-32","41-45","45-56")),
        ("DESTINY-Breast09", "THP", 387, 68, 3, (6,12,18,24,30),
         ("98.7","95.8","91.3","87.1","82.3"),
         (387,381,375,369,362,353,337,297,255,220,158,88,43,15,4,0),
         ("5","16","31-33","42-48","48-62")),
    )
    for trial, arm, n, total, risk_step, times, percents, risks, expected in specs:
        points = tuple(zip(times, percents))
        at_risk = {risk_step * i: r for i, r in enumerate(risks)}
        result = candidate_sets(n, points, at_risk, total)
        labels = tuple(map(condensed, result))
        if labels != expected:
            raise AssertionError((trial, arm, labels, expected))
        print(f"{trial} | {arm} | N={n}, OS deaths at data cutoff={total}")
        print("  elapsed months → cumulative recorded-death outer sets:",
              ", ".join(f"{t}: {x}" for t, x in zip(times, labels)))
        if trial == "PERSEUS":
            only = candidate_sets(n, points, at_risk, total, link=False)
            print(f"  36mo using baseline-to-36mo alone: {condensed(only[2])}")
            print(f"  recorded deaths with elapsed time >36mo by cutoff: {total - min(result[2])}")


if __name__ == "__main__":
    run_examples()
