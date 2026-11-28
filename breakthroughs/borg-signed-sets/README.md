# The equality cases in Borg’s signed-set conjecture

**Result:** We prove exactly when a maximum intersecting family of signed sets must fix a common labelled element. The proof also gives the exact weighted Lovász theta value for arbitrary complete label cells.

## The question, in plain language

Start with any list of teams. Give each person on a team one of several badge colors, and allow **every** color assignment on every listed team. Now pick as many colored teams as possible, with one rule: every two chosen teams must share at least one person wearing the **same** color. Sharing a person in different colors does not count.

A simple strategy is to pick one person and one color and keep every team containing that exact combination. This is called a **star**. The proof shows you can never beat the best star, no matter which teams are on the original list. With three or more colors, every best solution is a star. With two colors, there is an exact exception: at least three people must always appear together in the listed teams, and their stars must be largest. For those three people, keeping teams where at least two wear red ties the best star without requiring any one person to wear red.

## What is proved

**Relation to prior work.** The ordinary bound already follows from [Chang, Liu and Liu (2026), Proposition 5.3](https://arxiv.org/html/2609.19123v1#S5), applied via [Borg’s 2011 Lemma 4.3](https://staff.um.edu.mt/peter.borg/publications/conjectures%20gen.pdf#page=8). We do not claim priority for that bound. The contributions here are a sharp semidefinite certificate and the classification of all maximizing families, including the strictness clause of Borg’s original conjecture.

- **Weighted upper bound:** Even if different teams have different weights and different people have different numbers of available colors, the weighted Lovász theta of the disjointness graph equals the largest star weight. Every person must have at least two colors, and all colorings of any one team must be present at the same weight.
- **Exact equality cases:** For a common number of colors, the proof classifies all maximum intersecting families. This gives both parts of Peter Borg’s [2007 Conjecture 1.4](https://doi.org/10.37236/959), including the binary exception. The star-existence part was restated as [2011 Conjecture 3.1](https://staff.um.edu.mt/peter.borg/publications/conjectures%20gen.pdf#page=5).
- **Not claimed:** This is not a proof for arbitrary missing color assignments, intersection in two or more points, general spectral Chvátal H, or the more general weighted hereditary conjecture from Borg’s 2011 paper.

[Read the self-contained proof](proof.md) · [Read the equality and strictness proof](equality.md)

## Check it yourself

The argument is deductive. These scripts check exact arithmetic and examples independently of the general proof; they do not replace it.

```sh
node dual-replay.cjs > /tmp/borg-dual.json
diff -u dual-expected.json /tmp/borg-dual.json
node equality-replay.cjs > /tmp/borg-equality.json
diff -u equality-expected.json /tmp/borg-equality.json
```

The first script checks complete rational dual certificates, including examples with negative entries in the support matrix. The second enumerates maximal intersecting families on all 127 nonempty support families on three binary groups, plus four additional binary and ternary cases. It checks the largest star bound and tests the predicted structure of every optimum.

**Review status:** AI-assisted mathematical research, independently re-derived and exact-checked, but not yet externally peer reviewed. A bounded literature search found earlier special cases but did not find this general equality classification or exact theta result. This is not a priority claim.
