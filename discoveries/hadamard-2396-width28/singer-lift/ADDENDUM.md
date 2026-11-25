# Order 2396: a complete obstruction to the transverse Singer lift

AI-assisted Research1 verification note, 25 November 2026. Not peer reviewed. No full order-2396 matrix, and no nonexistence claim outside these symmetries.

## Scope and published convention

[Egan–Kotsireas–London–Ó Catháin, arXiv:2610.27175v2](https://arxiv.org/html/2610.27175v2), Theorem 7 and Table 5, give the bordered template and list 2396 as the first missing order. [Munemasa–Núñez Ponasso, arXiv:2611.27591](https://arxiv.org/pdf/2611.27591), submitted November 20, likewise describes 2396 as open. Use the width-28 case (quotient index q=8) over an abelian group G of order 592. Write G=A×F37 (additive field), with |A|=16 and N=⟨e⟩×F37 for an involution e∈A. A core consists of four signs a_f:G→{−1,+1}; if v=|G|+q−1, its necessary summed periodic autocorrelation is

    R(s) = Σ_f Σ_x a_f(x)a_f(x+s) = 4v[s=0] − 4q[s∈N] + 4.

The value off N really is +4, not zero. In this case v=599 and R(0)=2368. The [2026 construction preprint](https://hexagonmath.org/pdf/2611.00637v1), equation (4) and Theorem 5, confirms the same convention. Its Section 5.8 uses order-five multipliers on F11 for the **physical** generators at order 1436, including finite coset-sum classifications and a projection congruence. Appendix C's Singer cycle organizes the **frame** on a quotient, not a Singer invariance imposed upon the four full physical generators. No equality of those symmetries is assumed here. The projection identity below is standard; its consequence is specific to the stated profiles and template. No matching obstruction was found in the bounded source check.

## Complete no-lift: even one physical Singer block is impossible

For a subgroup B≤G, compress a_f by summing it over each coset of B to obtain z_f on G/B. Expanding the square of each coset sum, or collecting pairs of indices with difference s∈B, gives the exact diagonal identity

    Σ_f Σ_{u∈G/B} z_f(u)^2 = Σ_{s∈B} R(s)
                           = 4(|G|+|B|−1−q(|B∩N|−1)).

Suppose even one sign block has, in each B-coset for an order-8 subgroup B, a singleton and seven signs all equal among themselves. Its coset sum is one of ±6,±8. At order 2396 that single block alone contributes at least 74·36=2664 to the nonnegative left hand side. The right hand side is at most 4(592+7)=2396 (exactly 2396 if B∩N=0, or 2364 if B∩N has order 2). This is impossible. More generally, this one-block order-8 test excludes the same profile whenever |G|>56 in the stated template. Different singleton positions on different cosets are harmless, provided this profile actually holds.

For the frozen positive C37 quotient at G=C2×F8×F37, independent order-seven scalar action on F8 makes **each** physical block have exactly this profile with B=F8. All 296 coset entries would contribute at least 10,656. Hence no genuine C7×C3 lift exists, for any of the frozen C37 margins or any other margins. The C3 restriction on F37 is not needed. The previous h=3 projected table indeed satisfies all 16 shifts after summing along F37; it is the distinct, transverse compression along F8 that now decides the lift.

### Single-block automorphisms and affine scope

A genuine linear automorphism of G with order seven on A and identity on F37 forces such a profile if it fixes even one physical block. To see this, let τ∈Aut(A) have order seven. More generally, for any prime ℓ odd and automorphism τ of order ℓ on the 2-group A, the averaging endomorphism P=(1+τ+⋯+τ^{ℓ−1})/ℓ is defined (multiplication by ℓ is invertible). It splits A as C×B, with C the fixed subgroup and B=ker P. All nonzero τ-orbits on B have size ℓ. For ℓ=7, |B| divides 16 and is congruent to 1 modulo 7; the only nontrivial choice is 8. In fact B is elementary abelian, since transitivity forces its nonzero elements to have the same order and B contains an involution. Thus C has order two and each B-coset has the forced profile. N need not be τ-invariant for this single-block conclusion.

The same holds for one block invariant up to a **constant** sign under a genuine affine product map (x,y)↦(τx+t, λy+u) with τ of order seven, λ≠0, and constant t,u. Take the 2·37·36-th power, coprime to seven: this removes the affine F37 component (its group order is 37·36), constant sign, and translation on C of order two. On B the remaining affine action has linear part a generator of ⟨τ⟩, hence has a unique fixed point (1−τ^j invertible) and one orbit of seven. It is the same physical profile. This is not a conjugacy claim for arbitrary slice-dependent maps whose centers or relative shifts evolve between different fibers.

## Additional reusable prune: all four cannot share a transverse Klein-four rotation

Suppose G=B×K, N≤K, with B elementary abelian of order four. Impose a common order-three rotation on B, identity on K, on **all four** physical generators. Denote the sign at 0∈B by c_f(u) and the common other sign by d_f(u). Write

    U_f=(c_f+3d_f)/2,   W_f=(c_f−d_f)/2 .

These are half the B-trivial and any nontrivial B-character sums. For each shift s∈K, let P_U(s), P_W(s) be their summed periodic autocorrelations over K. The product rule for a real B-character χ gives P_χ(s)=¼ Σ_{b∈B} χ(b)R(b,s): this expands the two B-index sums and groups them by difference b. The delta term and N indicator are supported on b=0 and cancel in the difference. The +4 term gives

    P_U(s)−P_W(s) = ¼ Σ_{b∈B}(1−χ(b))·4 = 4

for every s, including zero. Each fiber is either ordinary, U=±1 and W=−U, or a defect, U=±2 and W=0. The shift-zero difference makes exactly one defect across all four blocks. At nonzero s only the products using this defect differ between P_U and P_W; in its block at u0, the difference is U(u0)[U(u0+s)+U(u0−s)]. For it to be 4, both neighboring signs must agree with the defect sign. This holds for every s≠0, including involutions (then the two products are counted twice). The defect block has row sum ±(|K|+1). But by total compression, Σ_f(Σ_u U_f(u))²=¼ Σ_s R(s)=v. Therefore such a rotation is impossible whenever (|K|+1)²>v=4|K|+q−1. Here |K|=148 and 149²=22,201>599. Unlike the order-seven test, the argument uses all four blocks sharing the pure rotation.

### Consequence for common linear multipliers

Every common genuine linear multiplier acting only on A preserves N (because R singles out its nonzero elements by the value −28), so fixes e. If a prime-order ℓ-odd such automorphism is nontrivial, its splitting A=C×B above has |C|≥2, thus |B|≤8. The only possibilities are |B|=4 and ℓ=3, or |B|=8 and ℓ=7: ℓ must divide |B|−1. The transitive B is elementary in both cases. The complement K=C×F37 contains N, so the two respective tests exclude both. By Cauchy's theorem, the group of common pure linear A multipliers is a **2-group**, for every abelian A of order 16. This does not eliminate coupled order-three actions on A and F37. If one also uses the [separate mixed-companion exclusion of an order-nine projected multiplier](https://github.com/ConfusedDuckCollectibles/Research1/blob/4bd4c6fbf791bba6afec2d904afec51ffb4349b6/discoveries/hadamard-2396-width28/mixed/ADDENDUM.md), the projection of any group M of common unsigned linear multipliers to F37^× has 2-group kernel and image with odd part at most 3; hence the odd part of |M| is at most 3. The companion hypothesis is separate from the two new arguments here. No claim is made about arbitrary signed automorphisms of a full Hadamard matrix.

## Reproducibility and stopping

    python3 verify_lift_gate.py

This standard-library verifier uses direct group autocorrelation and replays the earlier 16/16 C37-projected positive. It exhausts the four single-block order-seven sign profiles, both possible values of the native projected norm, and the direct B4 character differences at all 148 shifts. A small genuine width-12 **core** at target order 28 (not a new matrix construction) is the positive control, with R=[16,4,4,4] and the two Klein-four energies 7 and 3. No solver was launched: the full frozen family fails by the exact norm identity even before ordering its 208 orbit positions. No broader order-2396 search was done.
