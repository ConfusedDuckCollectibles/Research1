# How much can you change a sunflower-free family?

Take a collection of three-item sets. Four of them form a *sunflower* when every pair shares exactly the same core. Four disjoint sets qualify too: their common core is empty.

The published finite bounds allow between 39 and 49 sets without such a sunflower. This note narrows two ways of looking for a better result:

- **The known 39-set construction resists small rewrites.** There is no 40-set family that keeps 34 or more of those sets, even if the replacements use entirely new items. Two exact searches checked the full local case.
- **Large hypothetical examples need more structure.** A 49-set family must have an item in ten sets. In the separate case where even three disjoint sets are forbidden, an earlier upper bound of 35 drops to 33, or 31 if no item occurs more than nine times.

The overall 39–49 range is unchanged. The 33 and 31 results are upper bounds only, not constructions of those sizes.

- [Proofs and exact search scope](TRADE-AND-ANCHORS.md)
- [Independent local replay audit](INDEPENDENT-LOCAL39-AUDIT.md)
- [Independent structural audit](INDEPENDENT-STRUCTURAL-AUDIT.md)

```sh
node verify-near39.cjs
node independent-local39.cjs
```

The scripts are dependency-free. The structural arguments identify the few previously published computer-assisted bounds they import. The independent checks were within this AI-assisted project; this package is not peer reviewed.
