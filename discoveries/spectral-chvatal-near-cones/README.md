# A sharp inertia bound for one-coordinate-compressed downsets

A *downset* contains every subset of each of its members. In its
disjointness graph, the sets containing one fixed element form an
independent family called a *star*. Chvátal's original conjecture says
that a largest intersecting family can be a star. It was proved in
2026. The stronger spectral question asks whether an explicit
supported symmetric matrix can always certify that number by
counting its nonnegative eigenvalues.

**What this note proves.** For every downset compressed toward one
coordinate, we construct such a real symmetric matrix exactly. It has
no zero eigenvalues, and its number of positive eigenvalues equals the
largest-star size. For a fixed natural cross block, we also prove a
converse: the certificate works exactly when the residual family is an
antichain. A two-set inclusion in the residual family produces an
explicit rank-one obstruction to that cross block. It is not an
obstruction to every possible matrix.

This covers every shifted downset, including those with maximal sets
of different sizes. The ordinary largest-star result on this same
class is due to [Hunter Snevily (1992)](https://doi.org/10.1016/0097-3165(92)90059-4).
The new scope here is the explicit spectral certificate and exact
limit of its fixed cross block, not a new class for the ordinary
conjecture. [Ellis, Filmus and Friedgut, arXiv:2609.28404v3, §5](https://arxiv.org/html/2609.28404v3)
formulate the spectral questions. Their general Conjectures I and H
are not settled by this note.

- [Proof and precise scope](PROOF.md)
- [Optional graph-join extension](ADDENDUM.md)
- [Standard-library exact replay](replay.py)

Run the replay with Python 3.10 or later:

    python3 replay.py

It checks the matrix construction and inertia on 268 coordinate
choices across the four-coordinate downsets, the complementary 324
rank-one obstructions, and controls that distinguish the hypotheses.
The proof is analytic and does not rely on this finite enumeration.

Research note, AI-assisted and not peer reviewed. No priority claim.
