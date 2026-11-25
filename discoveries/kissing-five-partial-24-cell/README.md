# How 22 contacts can cap a five-dimensional kissing configuration

How many equal balls can touch one central ball without overlapping? The answer in five dimensions is still unknown. This note solves a restricted version: **if 22 of the contacts lie at vertices of one regular 24-cell, the total is at most 40.** Known configurations show that 40 is attainable. If 23 of those vertices are present and the total reaches 40, the missing vertex must be there too, and the full arrangement is one of the known \(D_5\) or Leech \(L_5\) configurations, up to isometry.

This extends the previously studied case with all 24 vertices present. It does not improve the unrestricted upper bound. The [exact-hole companion](exact-holes/) settles the remaining comparison when exactly 22 named vertices are present: the sharp maximum is 38, with equality only in the two known aligned completions with those vertices deleted. It also gives a sharp 37 for 1,928 of the 2,024 choices of three missing vertices.

- [Continuous proof and prior work](CHECKPOINT.md)
- [Independent audit and edge cases](audit/AUDIT.md)
- Exact controls: `python3 verify.py` and `python3 audit/verify.py` (standard library only).
