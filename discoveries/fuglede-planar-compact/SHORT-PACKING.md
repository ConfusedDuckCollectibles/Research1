# Short certificate for the 4,320-square spectral non-tile

This is an alternate verification of the [published construction here](PROOF.md), not a different shape. It replaces the packing clique search with a parity argument and verifies the same geometric window condition with 2,460 integer corner checks.

## At most ten centers in the finite group

In the positive 60-point seed of [Tao Zhang](https://arxiv.org/html/2607.15632v1), the finite difference set covers all nonzero CRT layers, because \(V_0-V_1=H\) and the ordered differences between their layers cover \(2,3,4\) and their negatives modulo 5. Therefore packing centers lie in one layer \(H=\mathbb{Z}_{12}^2\). In this layer the forbidden differences are the union of \(V_i-V_i\). The three forbidden forms for a difference \((u,v)\) are:

1. \(u=0\);
2. \(u-v=0\);
3. \(v\) even, \(u+v=0\pmod3\), and \(u\ne2\pmod4\).

The first two follow from Zhang's vertical and diagonal subgroups. For the third, Zhang's \(V_2\) is a six-point subgroup \(\langle(4,2)\rangle\) and its translate by \((1,2)\). Its differences are the subgroup and its two translates by \(\pm(1,2)\). Those three cosets give exactly the third form.

Suppose there are eleven centers. The labels \(u\) and \(u-v\) are each injective. Split the centers into six classes by \((v\bmod2,\ u+v\bmod3)\). In a class, each pair must differ by 2 modulo 4 in its first coordinate, so a class contains at most two points and any pair has equal \(u\)-parity. Thus five classes are pairs and one is a singleton. Translate all the centers so the singleton has even \(v\). Each of the two injective labelings uses either five or six odd residues. Each paired class contributes an even number to each odd count. The singleton contributes the same parity to both, so the two total odd counts must be equal. Yet each of the three odd-\(v\) paired classes changes the count by \(+2\) or \(-2\) when switching from \(u\) to \(u-v\). Their sum cannot be zero. This contradiction excludes eleven and any larger packing. The ten-center witness is checked in the script below.

## Check only the extreme pairs in the 437 × 99 window

The 4,320-square board repeats its 60-cell motif with periods \((60,12)\) eight by nine times. On a fixed modular conflict class, each seed difference in one coordinate with period \(n\) is represented by \(r\) or \(r-n\), where \(0\le r<n\). Write a proposed window difference as \(r+hn\). With \(k\) repeats, the two seed representatives accept respectively \(h\in[-(k-1),k-1]\) and \(h\in[-k,k-2]\).

For each coordinate use the lowest and highest allowable quotient indices in the window. Every intervening index accepts every choice available at **at least one** of its two endpoints: the two acceptance intervals differ by one step and overlap. Choose such an endpoint independently in each coordinate. If the resulting corner is realizable by an actual pair of seed differences, the original interior index pair is realizable by the same pair. It is therefore enough to check four corners for each of 615 modular conflicts. The executable checks all 2,460, including the zero residue.

```sh
node verify-window-corners.js
```

Thus every lattice packing again has at most ten origins per \(437\cdot99=43{,}263\) window. A tiling would require origin density \(1/4{,}320>10/43{,}263\). The original proof supplies the spectral construction and transfers the density obstruction to arbitrary real translations. The conclusion and seed attribution are unchanged; this shorter certificate is independently replayed, not externally peer reviewed.
