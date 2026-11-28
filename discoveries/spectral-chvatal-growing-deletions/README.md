# Weighted Spectral Chvátal with growing triple deletions

Start with every possible team of up to three old people. Cross out some three-person teams, with overlap allowed. Then add one newcomer to the outer list, paired with any chosen old people. Under the two-level weighting below, picking teams that all share one person is optimal in the covered range. A spectral certificate proves the upper bound, so this is more than a search for good teams.

**The result.** Let $n$ be the number of old coordinates. Assume $n\ge100$, each old coordinate is in at most $n/8$ deleted triples, and at least one old coordinate is untouched. The inner family is all nonempty old sets of sizes 1–3 except the deleted triples. The outer family contains all nonempty old sets of sizes 1–3, the new marker singleton, and its bridge pairs with any selected old anchors. Every inner face has weight $u+v$ and every outer-only face weight $v$, with $u,v\ge0$. Weighted Lovász theta of the outer disjointness graph equals the maximum weighted star. With $S_0=n+\binom{n-1}{2}$ and $k$ the smallest deletion degree among selected anchors, the common value is $S_0(u+v)+\max(0,v-ku)$. If there are no anchors, the correction is zero.

**Why it extends the earlier result.** The pattern need not be a cycle, the deletions can share tips or have unequal degrees, and it need not stay fixed as the ground set grows. For example, using three disjoint blocks of size $s=\lfloor n/8\rfloor$ and deleting triples indexed by $(i,j,i+j\bmod s)$ gives $s^2$ deletions and maximum degree $s$. The theorem permits arbitrary patterns within the degree and untouched-point conditions, not just this example.

The proof gives the analytic positive semidefinite certificate. Its independent replay checks the exact inequalities and algebraic identities, including an irregular 112-point instance with 515 deletions and exactly one untouched point. It does not run a numerical semidefinite solver or use sampled cases to establish the theorem.

- [Proof and independent checks](proof.md)
- [Standalone exact replay](replay.cjs)
- [Expected output](expected.json)
- [Earlier Research1 cycle family](../../breakthroughs/spectral-chvatal-cycles/)

```sh
node replay.cjs > actual.json
diff -u expected.json actual.json
```

**Scope.** This is a restricted two-level weighted-theta theorem. It does not settle [general spectral Conjecture H](https://arxiv.org/html/2609.28404v3#S5), arbitrary downsets, or arbitrary weights. The cutoff is sufficient, not sharp. The ordinary non-spectral Chvátal conjecture was already [proved separately](https://arxiv.org/abs/2609.19123).

AI-assisted research note with independent computational and deductive rechecks. No external peer review or priority claim.
