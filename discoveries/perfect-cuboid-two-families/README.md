# Neither internal-ratio parametric edge-cuboid family yields a perfect cuboid

AI-assisted research note, not peer reviewed. No priority claim. [De Grey, Gibbs and Helm, arXiv:2401.06784v1](https://arxiv.org/html/2401.06784v1#S7.SS4)
give two internal-ratio edge-cuboid parametrizations in §7.4 (the two
constructions in the proof of Theorem 5). In
[§10.2](https://arxiv.org/html/2401.06784v1#S10.SS2) they ask whether either
parametric family ever contains a perfect cuboid. This note excludes
both families **for all parameters**, without a search bound.
It does not exclude other cuboids with the same internal aspect ratios,
and does not answer the unrestricted perfect-cuboid problem. Their separate
face-ratio family is already excluded by their Theorem 3.

## Definitions directly from the published tuples

Write
`P=p²+q², D=p²−q², U=pq, F=P⁴−4U²D²`
with p>q>0. The common factor in the integer edges and space diagonal
is the product `(P²+2UD)(P²−2UD)=F` in the two printed tuples.
Use a,b for the two rational edges, c for the possibly irrational edge
and f for the diagonal of the face with edges b,c. The following comes
from the lengths in §7.4, not from supposing a factorization of c²:

| Internal family, in source order | a | b | f | body diagonal |
|---|---|---|---|---|
| I (second overall) | DF | 4UD²P² | 2UF | PF |
| II (third overall) | 2UF | 8U²DP² | DF | PF |

In either family c²=f²−b². We initially take p,q to be coprime of
opposite parity, as in the paper. If an integer is the square of a
rational number then it is an integer square (reduce the denominator).
This makes it enough to rule out the following normalized radicands.

### I. A residue class excludes the missing edge

$$
R_1=\frac{c^2}{(2U)^2}=F^2-4D^4P^4.
$$

Here P and D are odd and U is even, so F ≡ 1 (mod 16) and
`R_1 ≡ 1−4 ≡ 13 (mod 16)`, not a square.

### II. A coprime factor forces Adrain's impossible quartic

$$
R_2=\frac{c^2}{D^2}=(F-8U^2P^2)(F+8U^2P^2).
$$
For a real candidate with rational missing edge, this radicand would be a positive integer square; if it is negative, the tuple does not give a real box. Both factors are odd. Any odd prime dividing both divides their
difference `16U²P²`. If it divides U, then `F ≡ P⁴` and P is
a unit at that prime: gcd(P,U)=1 follows from coprimality of p,q.
If it divides P, use `D²=P²−4U²`, giving `F ≡ 16U⁴`,
again nonzero at that prime. Thus the two factors are coprime.

A positive square product forces the positive factor to be square, but
$$
F+8U^2P^2=P^4+(2U)^2P^2+(2U)^4.
$$
For nonzero rational A,B, (A^4+A^2B^2+B^4) cannot be a rational square.
This is the (d=1) case in [Ezra Brown's account,
§1](https://doi.org/10.1017/S0017089500007862), which attributes it to
**Robert Adrain**, not Euler. For a short independent proof, see below.
Hence the second missing edge cannot be rational.

### Addendum: parameterization with arbitrary positive rational inputs

Each normalized radicand is homogeneous of degree 16 in p,q. Rational
rescaling of p,q therefore changes it by a rational square factor, so
we can first take coprime positive integers. If both p,q are odd, set
(r=(p+q)/2), (s=(p-q)/2); these are coprime integers of opposite
parity, since r+s=p is odd and any common divisor divides p and q.
If (P',D',U',F') refer to r,s, then (P=2P',D=4U',U=D',F=16F').
Consequently `R_1(p,q)=256R_2(r,s)` and `R_2(p,q)=256R_1(r,s)`.
Thus both results also hold for every rational (p>q>0), if the same
homogeneous formulas are used, and on any nonzero rational scaling
of the cuboid itself. We exclude only nondegenerate configurations.

The same factorization gives a unified proof: take A=D for I or A=2U
for II and B as the other leg of the Pythagorean triangle. Then
P²=A²+B², F=P⁴−A²B² and the normalized radicand in either case is
(F−2A²P²)(F+2A²P²). For a primitive integer right triangle these
are odd factors with gcd 1 by the same argument: a common odd prime
would divide A or P, while F is nonzero modulo either such prime.
The positive factor is P⁴+A²P²+A⁴, so Adrain's result applies directly.
Any rational right triangle can be scaled to primitive integers and
the radicand is multiplied by a rational square. This also gives the
rational-input extension without using the parity swap.

## Independent proof of the Adrain step

If (W^2=A^4+A^2B^2+B^4) with A,B nonzero rationals, put (t=A/B);
then ((x,y)=(t^2,tW/B^2)) is a point with x>0 on
$$
E: y^2=x^3+x^2+x.
$$
The standard 2-isogeny descent (see
[Dixit–Murty–Pathak, §2](https://www.imsc.res.in/~anupdixit/rankproblem.pdf))
takes (E': y^2=x^3-2x^2-3x=x(x-3)(x+1)).
For (y^2=x(x^2+a x+b)), the square-class map is x away from the
origin and infinity, b at (0,0), and 1 at infinity. Its possible
square classes are signed squarefree divisors of b. On E, b=1 and
(x^2+x+1>0) on the real line, so all classes equal 1. On E', b=-3
and the images of its three nonzero rational 2-torsion points have
classes -3,3,-1. Its image plus 1 thus has size 4, the maximum
possible. The rank formula
$2^{r+2}=|\operatorname{im}\alpha_E|\,|\operatorname{im}\alpha_{E'}|$
yields r=0.

The discriminant of E is -48; reduction at 5 and 7 is good and each
reduced group has order 8 (enumerated in the exact replay).
Prime-to-p torsion injects at good p, so these two independent primes
force the rational torsion order to divide 8, including for the
5- and 7-primary parts. E has just one rational point of order 2, (0,0).
A rational point of order 4 would double to (0,0). Its doubling formula
has (x(2Q)=(x(Q)^2-1)^2/(4y(Q)^2)), so x(Q)=1 or -1. Negative x is
impossible, while x=1 gives y²=3. There is therefore no order-4 point
and no higher 2-power torsion. Hence $E(\mathbb Q)=\{O,(0,0)\}$, a
contradiction to x>0. This proof uses no rank database, GRH, or BSD.

## Reproduction and scope

`python3 verify.py` uses the standard library. It checks the algebra as
polynomial identities, independently starting from the expanded formulas
printed in §7.4; all four squared length relations in both tuples are
recovered. It also checks the parity-swap substitution, the 128
opposite-parity mod-16 pairs, the two finite-field point counts, and the
authors' example pairs:
- I: p,q=2,1 gives 1443, 1800, sqrt(461776), body diagonal 2405.
- II: p,q=4,1 gives 552968, 554880, sqrt(767093529825), body diagonal 1175057.

The examples and finite-field counts are controls, not sample-based
exclusions. The proofs above close both displayed infinite families.
No matching prior application was located in the sources checked; the
quartic obstruction itself is classical. Other recent reductions, such as
[Peschmann’s](https://arxiv.org/abs/2604.28072) and
[Cipollini’s](https://arxiv.org/abs/2609.14526) fixed-parameter fiber
exclusions, address other restrictions and are not used in this proof.
