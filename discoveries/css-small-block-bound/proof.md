# Finite asymmetric distance bound for CSS codes with pair-generated X regions

## Setup and certificate

Use a binary CSS **stabilizer** code specified by X- and Z-generator matrices H_X, H_Z satisfying H_X H_Z^T=0. Gauge-generator weights of subsystem codes, propagation and measurement schedules are outside scope. All distances below are physical Hamming weights:
```
k   = n - rank(H_X) - rank(H_Z);
d_X = min wt(x) for x in ker H_Z outside row(H_X);
d_Z = min wt(z) for z in ker H_X outside row(H_Z).
```
Use a proposed positive target D for d_Z, and assume k>0 and d_Z≥D. W_X,W_Z denote the maximum **physical row weights** in the stated generator presentation; no column constraint or 2D geometry is assumed. Redundant generators may be listed, and dense changes of presentation do not preserve W_X,W_Z.

Construct a graph of physical positions with an edge ij precisely when e_i+e_j belongs to row(H_X). Use low-weight stabilizers intrinsically, not only weight-two rows visibly listed. Contract its connected components into blocks, with physical sizes w_i. First remove blocks pinned by an X-singleton stabilizer. Use the quotient check matrices A and B defined below, then remove blocks fixed by the unit vector e_i in row(B); these are encoding-zero tensor factors, with no effect on distances. Let S collect the remaining “small” blocks w_i<D and L the remaining “large” ones w_i≥D. Let n_S and n_L be their total physical sizes; n_0 is the number of removed physical positions. Define V_S as the coordinate subspace supported on S and
```
ell = dim((ker A ∩ V_S)/(row(B) ∩ V_S))
    = |S| - rank(A restricted to S) - rank(B) + rank(B restricted to L).
```
Thus ell counts independent Z logical classes representable entirely inside the small blocks. It can be computed from ranks without optimizing the physical distance.

**Theorem.** If each Z-generator row touches at most two blocks in L, then
```
n_L ≥ D d_X (k-ell), and consequently n ≥ n_0 + n_S + D d_X (k-ell).
```
In particular the row-touch condition is automatic if W_Z<3D, since Z rows must cover each block they touch in full. The actual projected row-touch property can be used even when some original rows weigh ≥3D. More generally, the same bound holds if the projected Z row space is spanned by its weight-one and weight-two vectors, irrespective of the displayed basis. See [the intrinsic rank criterion and a 24-qubit control](audit.md#3-presentation-free-projected-rank-test). Define r_S=rank(A restricted to S). The rank-only bound without needing the B intersection is
```
0≤ell≤|S|-r_S≤(W_X-1)r_S ,
n ≥ n_0+n_S+D d_X max(0,k-(W_X-1)r_S).
```
When n_S<D, ell=0. If W_Z<2D and d_X≥2, no X logical can be represented entirely on L, so ell=k; in particular k≤(W_X−1)r_S. If the X stabilizers are generated entirely by weight at most two, then A=0, no remaining small block is possible, and the main bound reduces to the sharp obstruction n≥k d_X D for W_Z<3D. In the pair-only case d_X≥2 is impossible when W_Z<2D. All statements are at fixed target D; choosing D=d_Z is optional, not required.

More specifically, under the row-touch condition one can compute the **exact minimum X logical weight constrained to L**: form a graph on L with an edge for every two-block Z projection, and pin any block with a singleton Z projection. For each unpinned graph component C, use its compressed all-ones vector 1_C. The answer is the smallest number of blocks |C| with 1_C not in row(A), or infinity if there are none. This is an upper witness for the global X distance when finite, not a general formula for the global distance when S is nonempty.

## Proof

1. **Collapse the pair-generated physical subgroup.** Products of pair stabilizers along paths give all pairs within a connected block, hence all even binary X vectors in the block. Nothing else is needed to change a physical X representative to one chosen position per odd-parity block. Any Z vector commuting with these pairs is constant across each block. A singleton X stabilizer pins its whole component: putative Z vectors vanish there and all its X vectors lie in the stabilizer space. Remove such components. The quotient X checks A are blockwise parities of H_X's rows. The compressed Z checks B encode whether each H_Z row is a full block or zero. We have AB^T=0. These maps send unrestricted logical X representatives to ker B modulo row(A) with minimum lift weight the number of nonzero coordinates, and unrestricted logical Z representatives to ker A modulo row(B) with weight the sum of the sizes at its nonzero coordinates.

2. **Remove Z-fixed blocks.** If e_i is in row(B), every row of A has ith coordinate zero and X commuting with B has zero block parity there. All Z logical representatives can subtract that stabilizer. Together with the intrablock pair stabilizers, precisely that whole physical block is a separate encoding-zero stabilizer factor: X-even and the whole-block Z. Deleting its column from A,B and qubits from the lift preserves logical dimension and minimum weights. All stated physical row and column upper bounds can only decrease under restricting the given physical checks. For the residual blocks, rank accounting gives k=m-rank A-rank B with m=|S|+|L|. In fact row(A) has no nonzero weight-one or weight-two vector: using intrablock pairs such a representative would produce a previously identified pin or intrinsic pair joining blocks.

3. **Count logicals on the two sides of the block cut.** Let V_L be the complementary coordinate subspace. Rank-nullity gives `dim(ker A ∩ V_S) = |S| - rank(A|S)` and `dim(row B ∩ V_S) = rank(B) - rank(B|L)`. These yield ell. Similarly the dimension of X classes representable inside L is
```
dim(ker B ∩ V_L) - dim(row A ∩ V_L)
= |L| - rank(B|L) - rank(A) + rank(A|S)
= k-ell.
```
This is the CSS cleaning dimension identity in this contraction; it makes no geometric separation assumption. In particular A and B can couple the two regions. If n_S<D, every Z representative on S has physical weight below the distance, so ell=0.

4. **Use the complementary graph, not only the rank identity.** A vector supported on L is in ker B exactly when its coordinates are equal on each two-block projection edge and zero on every component touched by a singleton projection. Its compressed kernel basis is the list of disjoint unpinned component indicators 1_C. Their logical classes span the k−ell dimensional quotient above, so choose k−ell independent nonzero component classes. Each gives a nontrivial physical X logical representative with weight exactly |C|, one per block. Thus |C|≥d_X. Each chosen block costs ≥D physical qubits, and the chosen supports do not overlap. This proves the lower bound on n_L. The smallest nontrivial component also gives the exact support-restricted distance: any nontrivial sum must include at least one nontrivial component, and the constituent support weights add.

5. **Show the exceptional rank cost.** After the fixed-block deletions, each small block i has e_i not in row(B). If A's ith column were zero, e_i would give a physical Z logical of weight w_i<D, impossible. Every small column is therefore hit by some row of A|S. A basis of that row space can be chosen from the presented rows, each supported on at most W_X small blocks because parity cannot inflate support. The union of the basis supports must cover S, so |S|≤W_X r_S. This proves the relaxed bound. If W_Z<2D every projected Z row has weight at most one and every nonzero unpinned graph component is a singleton, so d_X≥2 forces k−ell=0. If A=0, all small would have zero column, hence S is empty.

## Exact controls

Run `python verify.py`; no dependencies beyond standard-library Python. For each constructed check matrix, verify commutation, exact ranks, row and column maxima, and both distances. Full physical kernels are enumerated when small; otherwise exhaustive physical support enumeration goes through the first witness weight, ruling out all lighter candidates. A separate compressed computation verifies both distances under their two different metric conventions, ell by enumerating supported Z classes, and the graph-constrained distance by quotient enumeration. The code stores the actual generator supports and exact witnesses in `checks.json`.

| Control | Physical code | max X/Z row | max X/Z col | n_S, ell (D) | Bound on n |
|---|---|---|---|---|---|
| Rectangular repetition/Shor | [[18,1,(3,6)]] | 2/12 | 2/2 | 0, 0 (6) | 18 |
| Doubled Z-distance Steane | [[14,1,(3,6)]] | 4/8 | 4/3 | 14, 1 (6) | 14 |
| Their sum, exact exceptional correction | [[32,2,(3,6)]] | 4/12 | 4/3 | 14, 1 (6) | 32 |
| Add one X check crossing small/large regions | [[32,1,(3,12)]] | 6/12 | 4/3 | 14, 0 (6) | 32 |
| Same sum after an X basis change | [[32,2,(3,6)]] | 6/12 | 5/3 | 14, 1 (6) | 32 |
| One repetition block | [[6,1,(1,6)]] | 2/0 | 2/0 | 0,0 (6) | 6 |
| Product factors attached to the 18-qubit code | [[21,1,(3,6)]] | 2/12 | 2/2 | 0,0 (6); n_0=3 | 21 |
| Three-block outer even parity | [[9,2,(2,3)]] | 2/9 | 2/1 | 0,0 (3) | **Condition fails** |

The final control has three blocks of size D=3 with one Z row on all three. Here W_Z=3D exactly, d_X=2,d_Z=3 and k=2, so n=9 violates a hypothetical extension of n≥k d_X D=12. This gives necessity of the strict universal W_Z threshold for this conclusion. The case with a crossing X check takes the union of one shortest X logical from each component as a new stabilizer; k drops by one while a Z class confined to the small part disappears. This stress-checks the complementary rank term rather than assuming the code decomposes along the cut. None of these standard-derived control matrices are offered as code records.

For an exact sum of displayed per-sector column restrictions, use the actual matrix incidences: per-sector column maxima are not generally attained at the same physical qubit. No hardware, circuit-level errors, threshold or decoder accuracy conclusions follow from these distance statements.
