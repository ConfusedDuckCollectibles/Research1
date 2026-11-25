"""Independent exact 4 by 4 Fourier-minor gate in Q(zeta_13)."""
from itertools import combinations, permutations

P = 13
def orbit_representative(s):
    return min(tuple(sorted((a*t+b) % P for t in s))
               for a in range(1, P) for b in range(P))

def sign(s):
    return (-1)**sum(s[i] > s[j] for i in range(4) for j in range(i+1, 4))

def determinant(a, b):
    counts = [0]*P
    for s in permutations(range(4)):
        counts[sum(a[i]*b[s[i]] for i in range(4)) % P] += sign(s)
    return counts

reps = sorted({orbit_representative(s) for s in combinations(range(P),4)})
for a in reps:
    for b in reps:
        counts = determinant(a, b)
        # Phi_13 has degree 12 and all coefficients equal to one.
        assert len(set(counts)) > 1, (a, b, counts)
print("affine representatives:", reps)
print("nonzero exact minor representatives:", len(reps)**2)
