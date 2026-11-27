# Exact PSD for shared-marker cycle anchors when n ≥ 32

Research1 proof note, November 2026. Not peer reviewed. This proves the weighted residual theta bound for every simple 2-regular anchor graph with r≥5 and n≥max(32,2r+1), connected or disconnected. In particular it covers the minimum n=2r+1 for every r≥16. The n<32 cases are not established by this proof. No general Chvátal H result is claimed. The [arbitrary-cycle range theorem](range.md) solved affine feasibility already; the new step here is positive semidefiniteness of the complete coupled Gram. Compare [Ellis–Filmus–Friedgut, §5](https://arxiv.org/html/2609.28404v3#S5).

## 1. Model, normalized Gram, and the coupled rank part

There are n old points, r anchors, and r tips (one distinct tip for each edge of the 2-regular graph). Delete from the full old triple layer each triple consisting of an anchor edge and its tip. Let S be the remaining triples, all of weight 1. Add r weight-2 bridges: each consists of its anchor and the same new marker. Every anchor and untouched point is tight, with weighted star size K=C(n−1,2); each tip star has K−1. The marker has weight 2r≤n−1<K. Let X denote the incidence matrix of S with ALL old points, x_i its columns, and Y the bridge-by-tight-point weighted incidence (one entry 2 for each bridge's own anchor). Write a=C(n−2,2), b=C(n−3,2), q=C(n−4,2), and L=K−2r>0. J is an all-ones matrix of the indicated dimensions.

We seek a PSD block matrix with diagonal reciprocal weights and zeros on distinct intersecting faces:
    H = [ B+P   C ; C^T   I_r/2 ].
The old/bridge block C vanishes when the old triple contains the bridge anchor. Prescribe for each anchor i:
    C_i^T x_i=0,
    C_i^T x_j=1 for other tight j,
    C_i^T x_j=λ=(K−1)/K for every tip j.
These conditions can be met by taking C_i=1/a on surviving triples avoiding i, and adding an incidence-linear function on those triples. Indeed their incidence Gram on the n−1 retained coordinates is the complete triple Gram minus the deleted triple incidences. For m coordinates the full Gram is C(m−2,2) I +(m−2)J. Its inverse row absolute-sum norm is ≤2/C(m−2,2). The deletion perturbation has row absolute sums ≤6, since each old coordinate is in at most two deletions, and so does its spectral norm. Thus a surviving incidence Gram has inverse row norm ≤2/[C(m−2,2)−12], by a Neumann bound. The incidence residual to restore is ≤2/a in magnitude (each coordinate loses at most two counts; tips lose at most one and we target λ). Adding the triple sum of the solved coefficients gives the claimed exact C.

For each tip t choose a column V_t with value 1 on triples containing t and such that X^T V_t=(K−1)e_t. On triples avoiding t solve the analogous incidence-linear problem, with demands minus the original t,j codegrees, absolute size at most n−2. This determines a valid column; in particular |V_t(U)|≤v=6(n−2)/(b−12) off incidence. Let V collect the tip columns, and put
    g=1−2C1_r,
    F=VV^T/[K(K−1)],
    B=2CC^T+gg^T/L+F.
For a tight scalar coordinate i, g^T x_i=L, and Bx_i=1−2C_i for an anchor, or 1 for untouched. For a tip t, g^T x_t=λL, hence Bx_t=λ1+V_t/K. In particular (Bx_i)_U=1 whenever i belongs to U, for every old i including tips. These facts ensure the old fixed entries have the correct row support.

If P is PSD, P X=0, and the old locked entries of B+P equal I on the diagonal and zero on distinct intersections, then the complete H above works. Its Schur complement on eliminating bridges is P+gg^T/L+F⪰0 (this controls all leaf/nonstar coupling, not just isolated blocks). The tight weighted star rows give H [X_tight;Y]=J by construction. For a weighted tight star s, s^T H s=K and Hs=1, hence KH−J⪰0 by Cauchy–Schwarz. This is a weighted theta dual with value K: its diagonal at a face is K/w−1 and its entries on distinct intersecting faces are −1. A tight star supplies the matching lower bound.

## 2. Projected Kneser core and exact masked inverse

Let P* be the full old triple Kneser residual
    P*=I+A_disj/q−J/K.
It is zero on the span of the old incidence columns, and has spectral gap at least μ=1−1/q on their orthogonal complement. For a short self-contained derivation, f_R(U)=1_{R⊆U} for |R|=j≤3 spans the successive inclusion layers; they are nested by summing f_(R+w). A_disj f_R(U)=C(n−3−j,3−j) Σ_(T⊆R)(−1)^|T| f_T(U). Self-adjointness gives the scalar (−1)^j C(n−3−j,3−j) on the part orthogonal to lower inclusion layers. The j=0 contribution cancels with J/K; j=1 also vanishes. The j=2,3 eigenvalues are respectively 1+(n−5)/q and 1−1/q.

Write P_s for its survivor restriction, Π=X(X^TX)^−1 X^T, E=I−Π. Pad any survivor vector in ker X^T by zeros on deleted triples; its full old incidence inner product is still zero. Thus P0=E P_s E is PSD with gap at least μ on ker X^T and annihilates X.

Let M mask a symmetric matrix, preserving precisely positions (U,W) with nonempty intersection, including its diagonal. Call a masked symmetric matrix D compatible if (DX)_(U,j)=0 for every j in U. Let
    D0=M(I−B−P0),   T(D)=M(EDE).
D0 is compatible: for any incident j every entry in its star sum lies on the mask, and B supplies 1 there while P0 kills X. The map T takes compatible matrices to compatible matrices because EDE kills X. On this subspace it is self-adjoint for Frobenius inner product, with
    <D,T(D)> = ||EDE||_F² ≥ ||D||_F²−2||ΠD||_F².
For the inequality, split rows and columns orthogonally using Π,E. The two off-diagonal blocks have equal squared Frobenius norm and ΠDΠ is counted twice on the right, yielding the inequality directly.

To bound ΠD, take a column d supported on triples meeting its indexed triple U, and pad it by zero on deleted triples. Let X_in be the three complete incidence columns for U restricted to these local triple rows, X_out the other columns. The vector d is orthogonal to X_in by compatibility. The local matrices satisfy
    X_in^T X_in = a I_3+(n−2)J_3,
    X_in^T X_out = (n−2) J_(3,n−3),
    X_out^T X_out = 3(n−4) I_(n−3)+3 J_(n−3).
The last diagonal is 3(n−3) (all triples containing an outside point and meeting U); offdiagonal is 3. Projecting X_out off X_in gives Gram
    3(n−4) I + [3−3(n−2)²/(a+3(n−2))]J.
The bracket equals 3(7−n)/(n+3)≤0. Therefore ||X^T d||²≤3(n−4)||d||²; padding by zero preserves this argument. The smallest eigenvalue of X^TX is at least a−6. Summing the column projection estimates gives
    T ≽ κ I in Frobenius inner product,   κ=1−6(n−4)/(a−6).
For n≥32, κ>0; T is invertible on the compatible subspace and the compatible solution D=T^−1(D0) obeys ||D||_F≤||D0||_F/κ. Everything is rational: compatibility is a rational linear subspace and T a rational finite invertible linear system. Set P=P0+EDE. It has the exact locked entries, and is PSD whenever ||D0||_F< μ κ, since ||EDE||₂≤||D||_F and both terms kill X. No numerical solver is used.

## 3. Uniform quantitative bound

Use the following nonnegative rational bounds, with v defined in §1:
    c = (n−2)/(K a)+12/[a(b−12)];
    B_bd = 2(3/K²+r c²)+4(3/K+r c)²/L+(3+r v²)/[K(K−1)];
    P_bd = 2[(6/K+9(1/q−1/K))/a +120/(a(a−12)q)]
             +16(1+r/q)/a².
They give
    ||D0||_F ≤ (B_bd+P_bd) sqrt(M_n),
    M_n = C(n,3)[C(n,3)−C(n−3,3)].
Here M_n bounds the number of ordered locked entries even before deletions. Proof of both bounds follows.

**Rank-block cancellation.** Let W=C−J/K, h=W1. Since g=(L/K)1−2h, exact expansion gives
    2CC^T+gg^T/L = J/K+2WW^T+4hh^T/L.
Each triple row has at most 3 incident anchor coordinates where W=−1/K; the others have magnitude at most c, since 1/a−1/K=(n−2)/(Ka) and the C adjustment is at most 12/[a(b−12)]. Each V row has at most 3 incident tip entries equal to 1 and the remaining magnitudes ≤v. Applying Cauchy–Schwarz rowwise now proves |B−J/K|≤B_bd entrywise. The removal of all terms linear in W is essential for the uniform threshold.

**Projected Kneser perturbation.** Let X_D be the deleted-triple incidence, G=X^TX and G_f=aI+(n−2)J the complete old-incidence Gram. The deleted-interaction row
    u_U = X_U G^−1 X_D^T
has full-Gram approximation whose D-entry is (o_D−ω)/a with o_D=|U∩D| and ω=6/(n−1). The row X_U G_f^−1 has ℓ1 norm at most 5/a: explicitly subtract 2/(n−1) from every coordinate of the 3-indicator before dividing by a. Since G=G_f−X_D^T X_D, the inverse identity gives
    ||X_U(G^−1−G_f^−1)||_1 ≤ 60/[a(a−12)].
The mixed block identity P_s X=−P*_(S,D) X_D follows from P* killing full incidence. Thus ΠP_s=−XG^−1 X_D^T P*_(D,S). For any fixed W in S the magnitude of P*_(D,W) is 1/K when D meets W, 1/q−1/K otherwise. There are at most six meeting D. Also o_D≤2 (survivors and deletions differ), Σ_D o_D≤6 and rω≤3. For any set of at most six D, the sum of |o_D−ω| is ≤6: if S' is the sum of overlaps on those D, at least S'/2 terms are positive overlaps, giving sum ≤S'+ω(6−S')≤6. Across all deletions the sum is ≤9 by triangle inequality. Each coordinate belongs to at most two D and |P*_(D,W)|≤1/q. It follows that
    |ΠP_s| ≤ (6/K+9(1/q−1/K))/a +120/[a(a−12)q].
Moreover, ||u_U||₂≤4/a: the full-Gram approximation norm is ≤sqrt(12+18/(n−1))/a≤(18/5)/a when n≥25; the error term is at most sqrt(2)·60/[a(a−12)]≤(2/5)/a, using n≥25 and sqrt(2)≤3/2. The identity X^T P_s X = X_D^T P*_(D,D) X_D again follows from P* killing full incidence. The deleted principal block has spectral norm ≤1+r/q by its absolute row sums. Therefore |ΠP_s Π| ≤16(1+r/q)/a². Combining the two first-order terms and this quadratic term proves |P0−P_s|≤P_bd. On locked positions P_s=I−J/K, so B_bd+P_bd bounds every entry of D0 as claimed.

**Threshold and tail.** All bounds increase in r; replace it by r_max=(n−1)/2, which can be half-integral and is only used in the inequality. Then n⁴(B_bd+P_bd) is nonincreasing for real n≥32. Here are explicit checks:
* K/n², (K−1)/n², a/n², (a−12)/n², (b−12)/n² and q/n² all increase by their displayed binomial forms, as does unscaled L=K−n+1. Rewrite (n−2)/(Ka)=4/[(n−1)(n−2)(n−3)]; hence n³c decreases, the second term of n³c being 12/[n(a/n²)((b−12)/n²)]. The factor r_max/n²=(n−1)/(2n²) decreases. These facts give decay for n⁴ times 3/K²+r_max c².
* n²r_max c decreases: its first term simplifies to 2n²/[(n−2)(n−3)], and its other term is 6(n−1)/[n²(a/n²)((b−12)/n²)]. Thus n²(3/K+r_max c) decreases, handling the h-term. The function nv=12n(n−2)/(n²−7n−12) decreases, so r_max v² decreases and so does the scaled tip term.
* For the linear projection part use K−q=3(n−3), so the scaled difference term is a positive constant times (n−3)/n² times the reciprocal increasing normalized denominators K/n², q/n², a/n². The other terms decrease via their normalized denominators and (where appropriate) negative powers of n. In the quadratic projection part r_max/q also decreases.
* sqrt(M_n)/n⁴ decreases: C(n,3)/n⁴ and [C(n,3)−C(n−3,3)]/n⁴ separately decrease (differentiate the degree-three and degree-two polynomials divided by n⁴). Both μ and κ increase; for κ use (n−4)/(a−6)=2(n−4)/[(n−6)(n+1)], whose derivative is negative.

At n=32 the exact squared rational inequality
    (B_bd+P_bd)^2 M_n < μ² κ²
holds at r_max=31/2. Numerically only for orientation, its unsquared two sides are about 0.58614 and 0.60678. The accompanying rational replay checks the inequality without square roots or rounding, along with small exact construction controls. The monotonic comparisons above prove the strict inequality at every larger n and every allowed r. Consequently P and the full coupled H are PSD throughout the claimed scope. The proof permits arbitrarily many disconnected anchor cycles. It does not require a numerical/inexact Gram or a finite family interpolation.

## 4. Full two-layer downset crossing, if required

In the usual two-layer model let the inner family be the old rank-at-most-three downset minus the deleted triples, and let the outer be the full old layer plus r bridges and their shared marker singleton. Weights are a_0 on membership in the inner plus b_0 on membership in the outer. The largest star is
    (a_0+b_0)(n+K)+max(0,b_0−2a_0).
Untouched/anchor stars attain the stated sides of the crossing; tip and marker stars cannot exceed them under n≥2r+1. At a_0=1,b_0=2 split the weights into the old full triple certificate of weight 2, the weight-1 survivor/weight-2 bridge certificate just proved, and the full old singleton and pair layers of weight 3. These combine at s=3(n+K) using a common tight anchor. Explicitly, if a face has weights w_v^(k) from the portions with total w_v, take its normalized Gram vector to be the direct sum of (w_v^(k)/w_v) u_v^(k), omitting zero portions. Its squared norm is 1/w_v and a common weighted star is mapped to all ones, as its dot product with this vector is Σ_k w_v^(k)/w_v=1. The marker singleton intersects only the bridges; in the equivalent dual Gram Q=sG−J their r×r principal metric is (d+1)I−J with d=s/2−1≥r. Its projection onto their span is −Σ_i q_i/(d+1−r), with squared norm r/(d+1−r)≤d, and an independent residual gives exactly its required squared norm d.

The pure inner endpoint is restriction of full old rank-layer Kneser certificates. At the pure outer endpoint the rank-three part costs K. The remaining pair graph (old complete plus r anchor-marker pairs) and singleton faces can be colored with n+1 disjoint classes: choose odd M in {n,n+1}, label old points as distinct residues modulo M, color old pairs by i+j mod M and anchor-marker pairs by 2i mod M, unused on old pairs incident to i and distinct at the marker. Each singleton takes an unused color at its point, possible since the underlying pair-degree is at most n. Each color is a disjoint clique in Kneser terminology. The per-rank Kneser or clique certificates give total n+K+1, attained by an anchor. Weighted theta is homogeneous and subadditive in weights from its handle maximum formulation; interpolate the crossing with the appropriate endpoint with a shared maximizing center. This optional consequence proves all nonnegative two-layer ratios for n≥32, not all downsets.

Replay files: [replay.cjs](replay.cjs) and [expected.json](expected.json). The computation is a rational threshold check and finite identity controls; it does not generate the full Gram at n=32. Its existence follows from the analytic masked inverse. For n<32 the resulting positivity problem remains open here.
