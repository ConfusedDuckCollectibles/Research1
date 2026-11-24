"""Check numerical inputs to the counting proofs; no external dependencies."""
from math import comb


def g(n, e, s):
    return s * (n + 1 - s) - (e - n + 1) + (11 - s) * (10 - s) // 2


def budget(n, e):
    return comb(n, 2) - 55 * n + 10 * (11 * n - e)


for n, e, leave, ends in (
    (111, 1212, 90, (99, 4059)),
    (110, 1199, 55, (99, 3970)),
    (111, 1215, 60, (96, 4056)),
):
    assert budget(n, e) == leave
    assert (g(n, e, 12), g(n, e, n)) == ends
    assert min(g(n, e, s) for s in range(12, n + 1)) == min(ends) > leave
    print(f"{n}×{n} at E={e}: budget {leave}, endpoint bounds {ends}")

# Dow: b > 106 - 2 sqrt(13); 99 qualifies without floating point.
assert (106 - 99) ** 2 < 4 * 13 < (106 - 98) ** 2
assert 111 - (11 * 111 - 1212) >= 99
assert 111 - (11 * 111 - 1215) == 105  # prior Stinson obstruction

# Force unit defects on side 111; audit the hub argument's integer budgets.
assert 2 * 11 - 9 == 13 and 2 * 13 > 2 * 9
assert all(r * r + 1 > 9 for r in range(3, 12))
for lows in range(2, 7):
    for arms in range(2, 9 - lows):
        leaves_one = 13 - lows - arms
        assert leaves_one >= 5
        if leaves_one <= 1 + lows:
            assert lows >= 4
assert 4 + 1 + 3 * 2 > 9

# Clean-pencil count on side 110 and the classical equality benchmark.
assert 2 * (110 - 11) > 10 * 11
assert 10 + 10 * 11 == 120
assert 99 * 11 - (1199 - 120) == 10
assert 99 - 9 * 10 == 9
assert 10 - (10 - 9) == 9  # complete large groups
assert 9 - 2 == 7  # Latin squares
assert 10**2 * 11 + 10 * 10 == 1200
print("PASS: numerical inputs")
