# When maximal orthogonality need not be complete

This result is separate from the Fuglede equivalence. A family of frequencies can be **maximal** (no frequency can be added while preserving orthogonality) without being a **spectrum** (its exponentials need not span all functions on the set).

For every integer $m\ge4$ there is an $m$-point tile in a finite cyclic group with a maximal incomplete orthogonal family. For $m=1,2,3$, every maximal orthogonal family of a cyclic tile is complete. This gives the exact cardinalities for the universal completeness statement, and addresses the prime and prime-power indications of [Kolountzakis, Lev and Matolcsi, §7 (2025)](https://doi.org/10.1016/j.exmath.2025.125682), [accessible manuscript](https://arxiv.org/html/2410.12387v2#S7). These are not Fuglede counterexamples: the tiles are spectral.

## Prime sizes

For a prime $p\ge5$, take $A=\{0,2,3,\ldots,p-1,p+1\}$ in the cyclic group of order $6p$. The set is a complete residue system modulo $p$, so it tiles with the subgroup of multiples of $p$ and has the spectrum $(1/p)\mathbb Z/\mathbb Z$. Its mask is

$$
(1+X+\cdots+X^{p-1})(1-X+X^2)=\Phi_p(X)\Phi_6(X).
$$

Thus the nonzero frequency differences allowed for orthogonality are all nonzero points of the $p$-grid and $\pm1/6$. The pair $\{0,1/6\}$ is maximal: subtracting $1/6$ from a nonzero point of the $p$-grid gives order $6p$, and subtracting it from $-1/6$ gives order three. Neither is an allowed difference. The pair is incomplete because it has two rather than $p$ frequencies.

## Composite sizes

Write $m=pq$ with $p$ prime and $q>1$, and choose an integer $b>1$ coprime to $m$. In the cyclic group of order $bm$, take

$$
A=b\{0,\ldots,p-1\}+p\{0,\ldots,q-1\}.
$$

The first summand chooses each class modulo $p$ exactly once; the second then covers that class modulo $m$. Thus $A$ is a complete residue system modulo $m$. It tiles with the subgroup of multiples of $m$ and has the spectrum $(1/m)\mathbb Z/\mathbb Z$. The displayed representatives lie in $[0,bm)$ because $bm-b(p-1)-p(q-1)=(b-1)p(q-1)+b>0$.

Write $G_d=(1/d)\mathbb Z/\mathbb Z$. The two geometric-series factors in the mask show that its full circle zero set is

$$
(G_{bp}\setminus G_b)\cup(G_m\setminus G_p).
$$

The $p$ frequencies $j/(bp)$, $0\le j<p$, form an incomplete orthogonal family. For any candidate $t=r/(bp)$ in the first zero set outside the family, choose $j\equiv r\pmod p$ in that range. Then $t-j/(bp)$ is a nonzero member of $G_b$, in neither zero set since $\gcd(b,m)=1$. For a candidate $t$ in the second zero set, $t-1/(bp)\notin G_m$ since $1/(bp)\notin G_m$. Nor can it belong to $G_{bp}$, since $G_m\cap G_{bp}=G_p$ and $t\notin G_p$. This proves maximality. For example, $m=4,p=2,b=3$ gives $A=\{0,2,3,5\}$ in the cyclic group of order 12 with maximal orthogonal family $\{0,1/6\}$.

## The small cases and a real-line lift

A one-point tile has a one-element basis. At size two or three, a singleton cannot be maximal because the tile has a spectrum, which can be translated to include that frequency. At size two every orthogonal pair is complete. At size three, translate an orthogonal pair to $\{0,t\}$. Three unit complex numbers summing to zero form a rotated set of cube roots; their squares also sum to zero. Hence $2t$ is an orthogonal difference and $\{0,t,2t\}$ is a spectrum. These frequencies are distinct: $t=0$ was excluded, and frequency of order two cannot give a vanishing sum of three signs.

In either construction, $\Omega=A+[0,1)$ tiles the real line and is spectral by the complete-residue construction. The lift $\Lambda+\mathbb Z$ of the exhibited finite family is orthogonal and still maximal: a frequency in a new class modulo one has no zero from the interval factor in the Fourier transform, so that class would have to be addable to the finite family against the same mask. None is. The lift is incomplete because a nonzero finite vector on $A$ orthogonal to the finite family, viewed as a piecewise constant function on the intervals, is orthogonal to each lifted exponential. This is a classification of the universal completeness property by cardinality, not a classification of all incomplete families or a claim of priority.
