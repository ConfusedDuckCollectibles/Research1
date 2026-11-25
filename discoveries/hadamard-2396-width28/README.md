# Ruling out high-symmetry routes to Hadamard order 2396

A Hadamard matrix is a square grid of +1s and −1s whose rows are perfectly orthogonal. In the research literature we checked, 2,396 is the smallest size for which no such matrix is known.

A construction that worked at other sizes starts with four sign patterns and adds a border of 28 rows and columns. For this route at size 2,396, we proved that a specific kind of repeating symmetry can have order only **1, 2, or 3**. This rules out all six larger permitted orders of that symmetry, for every underlying abelian group allowed by this construction.

The proof does not rule out other symmetries, a core with no such symmetry, or other ways of building the matrix. **Order 2,396 remains open.**

- [Proof and primary sources](NOTE.md)
- [Exact replay](verify.py): run `python3 verify.py` with the standard library. The order-2396 positive control passes only after compression; it is not a matrix of the unsolved size.

### Mixed symmetries and search patterns

[The follow-up](mixed/ADDENDUM.md) rules out certain ninth-order scalar symmetries even when they also move slices or swap sign patterns. It works out the 32 normalized search patterns for one of the possible small quotients, respecting the actual row sums. Those are compressed search inputs, not complete matrices.
