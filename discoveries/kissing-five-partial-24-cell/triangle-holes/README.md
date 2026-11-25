# Three missing vertices still force completion at 40

The earlier two-hole and nontriangular three-hole bounds now cover every deletion pattern once combined with this proof. **If at least 21 contacts lie on one named regular 24-cell, the configuration has at most 40 contacts. At 40, all 24 template vertices must be present and the full configuration is one of the known aligned D5 or Leech L5 arrangements.**

With exactly 21 vertices present and a triangular hole, the maximum is 38 or 39. All 96 triangular holes admit a 38-contact construction. The proof excludes 39 if both high layers are full, or if a full high layer faces a six-point high layer. The remaining triangular cases and the unrestricted five-dimensional kissing number are open here.

- [Proof and attribution](TRIANGLE.md)
- [Independent audit](audit/AUDIT.md)
- Exact checks: `python3 verify.py` and `python3 audit/verify.py` (standard library only).
