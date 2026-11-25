"""Exact rational controls for the fixed-time converse; standard library only."""
from fractions import Fraction as F


def occupation_lengths(start, interval_lengths):
    """Integrate each mode's occupation over a unit period at rational start."""
    ends = [sum(interval_lengths[:j], F(0)) for j in range(len(interval_lengths) + 1)]
    start %= 1
    out = []
    for lo, hi in zip(ends, ends[1:]):
        occupied = sum(
            max(F(0), min(start + 1, epoch + hi) - max(start, epoch + lo))
            for epoch in (0, 1)
        )
        out.append(occupied)
    return out


def triangular_period(prob_low):
    # Q(p(u)) is low near both ends of the phase interval and high in between.
    return (prob_low / 2, 1 - prob_low, prob_low / 2)


def main():
    # Martingale law: Z=+-1/2 with equal weights, conditional A=+-1.
    for z in (F(-1, 2), F(1, 2)):
        p = (1 - z) / 2
        for raw_phase in (F(-7, 3), F(0), F(1, 11), F(1, 4), F(17, 27), F(87, 10)):
            low_a, high, low_b = occupation_lengths(raw_phase, triangular_period(p))
            assert low_a + low_b == p and high == 1 - p
            assert high - low_a - low_b == z
    assert ((1 - F(-1, 2)) / 2 + (1 - F(1, 2)) / 2) / 2 == F(1, 2)
    print("periodic test: same conditional occupations and averages at all test shifts")

    # c1=arccos(3/4) lies in (0, pi/4) since (3/4)^2 > 1/2.
    a, d1, d2 = F(3, 4), F(1, 8), F(3, 8)
    assert F(1, 2) < a * a < 1
    k = 4 * (2 * a * a - 1) / a
    assert 0 < k * d1 < d2 < 4 * d1
    assert d1 < 1 - a and d2 < 2 * (1 - a * a)
    assert d2 < 4 * a * d1 + 2 * d1 * d1 == F(13, 32)
    assert 2 * (1 + a) * d1 == F(7, 16)
    print("two-speed test: passes the short cone and each marginal range, but is impossible jointly")
    print("fixed-time converse exact controls passed")


if __name__ == "__main__":
    main()
