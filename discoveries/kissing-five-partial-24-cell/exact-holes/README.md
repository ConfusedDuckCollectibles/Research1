# What missing contacts force in dimension five

Fix one regular 24-cell as a template for contacts in four of the five dimensions. **If exactly 22 template vertices are present, the maximum number of contacts overall is 38.** The only ways to attain 38, when the two holes are fixed, are to delete those vertices from the known aligned D5 or Leech L5 configurations. With 39 or more contacts and at least 22 template vertices, at least one missing vertex must be present after all. With 40, both must be present.

If exactly 21 template vertices are present and the three holes do not form a triangle of adjacent vertices, the overall maximum is exactly 37. This covers 1,928 of 2,024 choices of holes. Every triangular choice admits a 38-contact example; we do not give a matching upper bound for those cases.

These are fixed-template statements. They do not settle the unrestricted five-dimensional kissing number.

- [Proof and prior work](COMPANION.md)
- [Independent audit](audit/AUDIT.md)
- Exact replays: `python3 verify.py` and `python3 audit/verify.py` (standard library only).
