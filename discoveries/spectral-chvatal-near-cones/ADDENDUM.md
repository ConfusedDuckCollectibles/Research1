# Optional extensions and controls

## Disjoint-ground unions via the graph join

The following standard linear-algebra closure observation is included to
illustrate why the fixed-cross obstruction is not a counterexample to the
full inertia conjecture. Suppose each of two graphs has a real symmetric
supported matrix Mi with exactly si nonnegative eigenvalues and independence
number si. Diagonalize the two Mi separately. On the graph join, all
cross-block entries are allowed, so in these independent eigenbases the
cross block can be chosen arbitrarily. Pair min(s1,s2) nonnegative
eigenvalues from the two blocks. For each pair (a,b), set just its
cross coefficient t greater than sqrt(a*b). That pair then has one
positive and one negative eigenvalue, because its 2-by-2 determinant
is negative. All unpaired blocks are unchanged. The resulting supported
matrix has exactly max(s1,s2) nonnegative eigenvalues, as desired.
This includes zero eigenvalues in the count.

For downsets on disjoint ground sets, first restrict each certificate
to the nonempty sets. Eigenvalue interlacing cannot increase the number
of nonnegative eigenvalues, and a largest star of nonempty sets still
forces the lower bound. The Kneser graph on the combined nonempty sets
is a graph join. The common empty set can then receive its own -1
diagonal (and zero cross couplings), since its disjointness loop is
allowed. This proves sharp inertia for disjoint-ground unions of
downsets with sharp inertia certificates. For example, the union of
the powersets of two disjoint two-sets passes even though the ordinary
cross-block recipe fails for every choice of the star coordinate.

## Alternative certificate in the single-parity special case

On any downset the ordinary square disjointness matrix factors as
Z diag((-1)^|A|) Z^T. Its inertia is therefore the number of even and
odd faces. If the residual antichain R at a coordinate consists of
faces of a single parity, pair every linked face A with A union {i}.
Each pair has one face of each parity and the remaining faces all
have one parity, so one of the two signs of the ordinary disjointness
matrix already gives sharp inertia. This accounts for homogeneous
truncations and pure shifted downsets. The main certificate is stronger
because R may contain both parities.
