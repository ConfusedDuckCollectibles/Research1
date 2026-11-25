#!/usr/bin/env python3
"""Independent synthetic exhaustions; no clinical data used here."""
from fractions import Fraction as Q
from itertools import product
from unittest import TestCase, main

from replay import candidate_sets, factor_bounds, percent_band


def km_schedule(n, schedule):
    """Block KM oracle. The schedule is tied (deaths, censors) per clock tick."""
    surv = Q(1)
    cumulative = 0
    for d, c in schedule:
        before = n
        if not 0 <= d + c < n:
            raise ValueError("toy test schedules leave at least one participant")
        surv *= Q(n - d, n)
        n -= d + c
        cumulative += d
        yield (before, n, surv, cumulative)


def render_percent(q, decimals, half_even):
    """Independent exact integer rounding oracle (not the implementation)."""
    scale = 10**decimals
    x = q * 100 * scale
    m, r = divmod(x.numerator, x.denominator)
    if 2 * r > x.denominator or (2 * r == x.denominator and (not half_even or m % 2)):
        m += 1
    digits = str(m).zfill(decimals + 1)
    return digits if decimals == 0 else digits[:-decimals] + "." + digits[-decimals:]


class TestCountCertificate(TestCase):
    def test_exact_extrema_all_orders(self):
        checked = 0
        for start in range(1, 10):
            for total_removals in range(start):  # leave positive risk
                grouped = {}
                for order in product((0,1), repeat=total_removals):
                    factor = Q(1)
                    for j, death in enumerate(order):
                        if death:
                            factor *= Q(start - j - 1, start - j)
                    grouped.setdefault(sum(order), []).append(factor)
                    checked += 1
                for d, factors in grouped.items():
                    self.assertEqual((min(factors), max(factors)),
                                     factor_bounds(start, start - total_removals, d))
        self.assertEqual(checked, 1013)

    def test_coverage_tied_grouped_schedules(self):
        # Four time ticks, including day zero and events exactly at landmarks 1
        # and 2. t=3 is the strictly-later risk anchor for the latter landmark.
        blocks = ((0,0),(1,0),(0,1),(1,1),(2,0),(0,2),(2,1),(1,2))
        counts = roundings = 0
        for n in (3, 4, 5, 6):
            for schedule in product(blocks, repeat=4):
                if sum(d+c for d,c in schedule) >= n:
                    continue
                rows = tuple(km_schedule(n, schedule))
                counts += 1
                # Both print conventions, including a day-zero risk that is
                # already depleted. Mixed conventions between tick labels
                # aren't required; testing all of them verifies a stronger case.
                for sides in product((0,1), repeat=4):
                    risks = {t: row[sides[t]] for t, row in enumerate(rows)}
                    for dec in (0, 1):
                        for half_even in (False, True):
                            points = [(t, render_percent(rows[t][2],dec,half_even))
                                      for t in (1,2)]
                            out = candidate_sets(n, points, risks, rows[3][3])
                            for where, opts in zip((1,2), out):
                                self.assertIn(rows[where][3], opts,
                                              (n, schedule, sides, points, out))
                            roundings += 1
        self.assertEqual(counts, 1520)
        self.assertEqual(roundings, 97280)

    def test_rounding_boundaries_and_clipping(self):
        self.assertEqual(percent_band("0.0"), (Q(0), Q(1,2000)))
        self.assertEqual(percent_band("100.00"), (Q(19999,20000), Q(1)))
        # Test every exact half-unit tie for 0-, 1- and 2-decimal percents.
        for decimals in range(3):
            den = 100 * 10**decimals
            for k in range(den):
                q = Q(2*k+1, 2*den)
                for half_even in (False, True):
                    lo, hi = percent_band(render_percent(q,decimals,half_even))
                    self.assertLessEqual(lo, q)
                    self.assertGreaterEqual(hi, q)

    def test_day_zero_included_with_after_zero_risk(self):
        # N=5, death then censor both at day zero. The printed post-day-zero
        # risk of 3 must NOT replace N as the start of the first interval.
        self.assertEqual(candidate_sets(5, [(1,"80.0")], {0:3, 1:3, 2:3}, 1),
                         (frozenset({1}),))
        self.assertEqual(candidate_sets(5, [(1,"80.0")], {1:3, 2:3}, 1),
                         (frozenset({1}),))

    def test_reject_unsupported_domain(self):
        with self.assertRaisesRegex(ValueError, "positive strictly later"):
            candidate_sets(5, [(1,"80")], {1:3, 2:0})
        with self.assertRaisesRegex(ValueError, "positive lower limit"):
            candidate_sets(5, [(1,"0.0"),(2,"0.0")], {1:3,2:2,3:1})
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            candidate_sets(5, [(0,"80")], {1:3,2:2})
        with self.assertRaisesRegex(ValueError, "nonincreasing"):
            candidate_sets(5, [(1,"80")], {1:3,2:4})


if __name__ == "__main__":
    main()
