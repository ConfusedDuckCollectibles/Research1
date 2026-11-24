# The least period at a power-of-two torus width

**Theorem.** Let τ(m,n) be the maximum size of a set in `(ℤ/m) × (ℤ/n)` with no three distinct points on a translate of a cyclic subgroup. For every integer a ≥ 1, the least positive period of the sequence `n ↦ τ(2^a,n)`, n ≥ 1, is **2^(2a−1)**.

Skotnica’s [2016 thesis](https://dspace.cuni.cz/handle/20.500.11956/83781), Theorem 3.20 and Corollary 4.4, respectively give attainment of 2^(a+1) points at this height and identify it as a period. The Introduction, printed p. 4, leaves its minimality open. We prove the necessary lower bound on height. The known attainment construction is given credit below.

A translate of a cyclic subgroup can be used in place of an image of a primitive integer line: an integral lift of any nonzero generator, after dividing by the gcd of its two coordinates, gives a containing primitive line. Conversely, the image of a primitive integer line is a cyclic coset. Thus the two conventions give the same caps (sets with no collinear triple).

## 1. A maximum-size cap contains a graph with a valuation gap

Set M = 2^a and N = 2^b. Assume b ≥ a and that a cap S of size 2M exists on the M-by-N torus. For each integer t modulo M, the M fibers of x − ty modulo M are cyclic cosets: one has direction (t,1), and M divides N. Each fiber has at most two selected points and so has exactly two. In particular, each column x holds two.

For 1 ≤ r ≤ a, let F_r(x,u) count points of S in column x with ordinate y ≡ u (mod 2^r). A Fourier coefficient at frequencies (A,B) is the sum over S of

```
exp(2πi (A x/M + B y/2^r)).
```

If A is nonzero modulo M and v₂(A) ≤ a−r, choose t satisfying At + B·2^(a−r) ≡ 0 (mod M). Such t exists because gcd(A,M) divides B·2^(a−r). The coefficient is then a sum of a nontrivial character of x − ty over fibers of constant size two, so it vanishes. All remaining x frequencies, including zero, are divisible by 2^(a−r+1). By Fourier inversion, F_r(x,u) depends on x only through x modulo 2^(r−1).

At r = 1 choose one common ordinate residue modulo 2 across all columns. At every subsequent level r, choose a lift of the previously chosen residue in each x class modulo 2^(r−1), independently for each class. A lift is available in one column of the class because the preceding residue was represented there. The constancy of F_r ensures it is then available in every column in that class. Finally choose a selected point with the prescribed ordinate residue modulo M in each column. Write this cap graph as (x,f(x)). If x and x′ have difference of valuation k < a, their choices at level k+1 agree. Therefore

```
v₂(f(x′)−f(x)) ≥ v₂(x′−x)+1                  (1)
```

for distinct x and x′; the valuation of zero is infinity. Only divisibilities modulo the given powers of two are used.

## 2. Count the lines through a fixed point of that graph

A line with odd y-direction is a fiber of x − ty modulo M. By (1), it cannot contain a pair of distinct graph points: the horizontal difference has smaller valuation than t times the vertical difference. All other primitive directions have odd x-component and even y-component, and can be represented by (1,s) with s even, by multiplying their parameters by an odd unit. Such a line through an anchor (x₀,y₀) has equation

```
(y−y₀) − s(x−x₀) ≡ 0 (mod q),      q = gcd(N,Ms).         (2)
```

Indeed the line parameters for x−x₀ are that residue plus multiples of M. Two slope labels describe the same line through the anchor exactly when both q and s modulo q agree, seen by comparing the anchor column and the next column.

Put d = b−a. If d ≥ 1, the labels fall into the following disjoint families:

- For 1 ≤ e < d, take s = 2^e u with u odd modulo M and q = 2^(a+e), giving 2^(a−1) labels for each e.
- All higher valuations, including s = 0, form one family with s = 2^d u, u arbitrary modulo M, and q = N, giving M labels.

Thus there are (d+1)·2^(a−1) eligible lines through the anchor. When d = 0, q = M and there are 2^(a−1) labels, just the even s modulo M.

Now fix a graph partner with horizontal difference of valuation k < a and vertical difference of valuation T ≥ k+1. It lies on exactly 2^k eligible lines. To see this, for d ≥ 1 and T < k+d, only the low family with e = T−k can work, since other families give unequal valuations in (2) below the modulus. Dividing (2) by 2^(e+k) gives a unique class of u modulo 2^(a−k); this class is odd and has 2^k lifts modulo M. If T ≥ k+d, only the high family can work, and division at e = d gives the same count with no parity constraint. This includes the case of zero vertical difference. If d = 0, divide (2) by 2^k instead: the unique class of s modulo 2^(a−k) is even by (1), and again has 2^k lifts.

There are 2^(a−1−k) partners of valuation k at each k = 0,…,a−1. Different partners cannot share an eligible line through the anchor, or that line would contain three graph points. Thus the total number of lines required is a·2^(a−1). For d ≥ 1 this gives a ≤ d+1, hence **b ≥ 2a−1**. For d = 0, the available 2^(a−1) lines imply a ≤ 1, also consistent with b ≥ 2a−1. For b < a, the row bound τ(M,N) ≤ 2N already excludes size 2M.

## 3. Attainment and least period

Skotnica’s Theorem 3.20 attains 2M at height H = 2^(2a−1). Its recursive construction, in bit form for x = Σ_i x_i 2^i, is

```
(x, ε + Σ_{i=0}^{a−2} x_i 2^(2i+1)),      ε = 0, 1.
```

The sum is empty for a = 1. The thesis starts with integer lifts (i,2i+ε) at width two and adds an optional shift (2^(a−1),2^(2a−1)) at each step. Unwinding gives the displayed formula because the last term vanishes modulo H. This is the thesis construction, not a new witness.

For completeness, the reduction from saturation to the least period is short. For n = 2^b h with h odd, split the torus by the Chinese remainder theorem into its 2-primary part `(ℤ/M) × (ℤ/2^b)` and the cyclic factor `ℤ/h`. A triple is collinear in the product if and only if it is collinear in the 2-primary part, allowing repeated points in the projection: the product of cyclic groups of coprime orders is cyclic. A cap of size at least three therefore projects injectively to a cap, and every cap in the 2-primary part lifts in one h-fiber. Since that part has at least two points, the optimal sizes are equal.

When b increases, the 2-primary torus embeds by multiplying the second coordinate by the corresponding power of two. Such an embedding preserves the cyclic or noncyclic character of any subgroup generated by two differences, so τ(M,2^b) is nondecreasing in b. It never exceeds 2M, and hence is constant at 2M for b ≥ 2a−1. It follows that τ(M,n) depends on n only through min(v₂(n),2a−1), so H is a period. Any smaller least positive period would divide H, by subtraction of periods, and would itself have to attain the same value 2M as H. The necessary height proved above rules that out. Thus H is the least period.

## References and status

- Michael Skotnica, *Maximální množiny bodů na diskrétní torické mřížce bez trojic bodů ležících na stejné přímce* (Charles University bachelor’s thesis, 2016), [institutional record](https://dspace.cuni.cz/handle/20.500.11956/83781), [full PDF](https://dspace.cuni.cz/bitstream/handle/20.500.11956/83781/BPTX_2015_1_11320_0_410101_0_174164.pdf?isAllowed=y&sequence=1); printed pp. 4 and 28–29 (PDF pp. 8 and 32–33).
- Michael Skotnica, “No-three-in-line problem on a torus: periodicity,” *Discrete Mathematics* 342 (2019), 111611, [DOI](https://doi.org/10.1016/j.disc.2019.111611), [arXiv v2](https://arxiv.org/abs/1901.09012v2), Theorems 2–4. This proves periodicity and connects the period to first attainment.
- Jim Fowler, Andrew Groot, Deven Pandya, and Bart Snapp, “The no-three-in-line problem on a torus,” [arXiv:1203.6604](https://arxiv.org/abs/1203.6604) (2012), and Aleksander Misiak, Zofia Stępień, Alicja Szymaszkiewicz, Lucjan Szymaszkiewicz, and Maciej Zwierzchowski, “A note on the no-three-in-line problem on a torus,” *Discrete Mathematics* 339 (2016), 217–221, [DOI](https://doi.org/10.1016/j.disc.2015.08.006), develop the toric problem and determine several prior special cases.

This AI-assisted note is not peer reviewed. A targeted search did not find a matching proof of binary minimality; no priority claim is made. The result concerns toric cyclic cosets, not lines on ordinary unwrapped square grids.
