# The storage problem hidden inside an average

Suppose a meter tells you the net energy over each time bin, but not exactly when surpluses and deficits happened. Can you tell whether one lossy store will cover every repeating power profile that could have produced those readings?

In the stated scalar model, this proof gives an **exact feasibility and capacity test** using those totals and known power bounds. It reduces the continuous timing problem to a finite set of cuts. In one exact example, checking every bin boundary says feasible, but checking a point inside a bin proves no finite store covers every allowed profile. In another example, keeping separate surplus and deficit totals instead of only their difference changes the worst-case requirement from infinite to finite.

The models repeat each candidate profile and impose no ramp constraint. This is not an assessment of a real electricity grid or an operating recommendation.

- [Proof and prior work](PROOF.md)
- [Independent audit](audit/AUDIT.md)
- [Exact replays](REPLAY.md): `python3 verify.py` and `python3 audit/check.py` (Python standard library only).
