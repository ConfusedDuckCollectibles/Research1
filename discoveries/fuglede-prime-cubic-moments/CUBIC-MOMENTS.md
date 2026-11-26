# A sharp cubic rigidity lemma for difference directions

Mathematical proof and primary-source comparison, 26 November 2026. This does not settle the full conjecture over F17 or any new exact intermediate cardinality. The point is a uniform constraint available for **every prime p>3** without enumerating pencils. The singular-locus formula used below is elementary and compatible with standard cubic surface theory; we do not claim that cubic surface taxonomy is new.

## The geometric theorem

Let E be a finite set in an affine 3-space over a field of characteristic different from 2 and 3, with full affine span, and let c be its greatest size of intersection with any affine plane. Suppose there is a nonzero homogeneous cubic F in three variables with F(x−y)=0 for every x,y in E. Then

    |E| <= max(4,c+1,floor(3c/2)).

Because full affine span implies c>=3, the right side is floor(3c/2). This bound is best possible on F_p for each c in [4,p], using the axis controls below. The quantifiers refer to a formal homogeneous polynomial and actual affine difference vectors. A polynomial vanishing on a projective direction set vanishes on all representatives by homogeneity.

### Step 1: Promote differences to a cubic surface with singularities at E

Translate E so that 0 belongs to it, and choose a basis b1,b2,b3 among its points. Polarize F to a symmetric trilinear form T with F(x)=T(x,x,x). This uses characteristic different from 2 and 3. Each x in E has F(x)=0. Expanding F(x−y)=0 gives T(x,x,y)=T(y,y,x) for x,y in E. Define a symmetric bilinear form B by

    B(b_i,b_j) = T(b_i,b_i,b_j).

Symmetry follows from the previous equality on basis pairs, and linearity in the third slot shows for any x in E and every i that T(x,x,b_i)=T(b_i,b_i,x)=B(b_i,x). Thus T(x,x,•)=B(x,•), and B(x,x)=F(x)=0, on E. The nonzero homogeneous projective cubic

    P(X,w) = F(X) − (3/2) w B(X,X)

has value and all partial derivatives zero at every (x:1), x in E: its derivative in X is 3T(x,x,•)−3B(x,•) and its derivative in w is −(3/2)B(x,x). The restriction w=0 equals nonzero F, so the surface polynomial is not identically zero.

### Step 2: Four singular axes reduce to four square-free monomials

The affine basis points 0,b1,b2,b3 are projectively independent. Send them to the four coordinate axes of projective 3-space. This is an invertible projective linear transformation; plane intersections with E still contain at most c points, since a projective plane pulls back to an affine plane or to the original infinity plane, the latter disjoint from E. Singularity at coordinate axis i excludes the coefficients of z_i^3 and z_i²*z_j for j≠i. Hence

    P(z) = sum_i alpha_i * product_{j != i} z_j.

There are three exhaustive classes:

1. All four alpha_i nonzero. On the torus all four z_i nonzero, set P=(product z_i)*sum(alpha_i/z_i). At P=0, the derivative w.r.t. z_i is the nonzero product times −alpha_i/z_i², so there is no singular point there. Exactly one z_i zero already makes P nonzero. Exactly two z_i zero gives a nonzero derivative w.r.t. either zero coordinate, from the term indexed by the other zero coordinate. The singular points are only the four axes, giving |E|<=4.
2. Exactly three alpha_i nonzero. Say alpha_0=0; then P=z_0*Q(z_1,z_2,z_3), where Q=az_1z_2+bz_1z_3+cz_2z_3 and abc nonzero. The 3-by-3 derivative matrix of Q has determinant 2abc, nonzero. If z_0≠0 and all partials of P vanish, z_1=z_2=z_3=0. Otherwise the point lies in the plane z_0=0 (indeed on its conic Q=0). There are at most c+1 points of E.
3. One or two alpha_i nonzero. Then P factors into three independent linear forms (for two terms, two common coordinates times a linear combination in the other two). The singular locus consists of the three projective lines given by pairs of these linear forms equal to zero; their common intersection is a point. Put s=1 if this point is in E, else 0, and let r_i count the E points on line i except the common point. Every pair of lines lies in a plane, hence r_i+r_j+s<=c. Summing yields 2Σr_i+3s<=3c, so |E|=Σr_i+s<=floor((3c−s)/2).

### Sharpness control

Choose approximately half c points on each coordinate axis, not counting the origin. If c is even, take c/2 on each axis and omit the origin; if c is odd, take (c−1)/2 on each and include the origin. Two axes lie in one plane with exactly c selected points. Any affine plane containing only one axis intersects the other axes at most in the origin; any containing none meets each in at most one point. For c>=4, E has two distinct nonzero points on each axis, so its affine span is full. It has floor(3c/2) points, and every difference kills F=X1*X2*X3. This realizes the bound whenever the field has at least ceil(c/2) nonzero elements; in particular in the stated prime-field control range. The replay also checks the exact projective singular locus sizes and point-by-point derivative identities for representative alpha patterns in F5 and F17.

## Spectral consequence

For p>3 and a dimension-three spectral pair (E,A) of any intermediate size mp, 2<=m<p, the general spectral plane cap is p (Fallon–Mayeli–Villano; [Ferguson–Iosevich, Lemma 8.2](https://arxiv.org/html/2610.30712v1)). Thus both sets affine-span (mp>p), and satisfy mp>=2p>floor(3p/2). No nonzero cubic vanishes on D(A) or D(E), by the theorem. The same holds for degrees one and two, since multiplying a nonzero polynomial of lower degree by an appropriate power of a linear form would give a nonzero cubic vanishing there.

For j=1,2,3 define the moment polynomial M_j(v)=Σ_{x in E}(v·x)^j. Orthogonality of the spectrum implies p equal parallel-slice populations for normals in D(A), so there M_j=mΣ_{t in F_p}t^j=0 because j<p−1. Therefore M_j is identically zero. By polarization or expansion in coordinates, every mixed moment of total degree j vanishes, as the multinomial coefficients are nonzero. Degree zero vanishes because |E|=mp. Interchange E and A. This argument is uniform and needs no p17-specific or F13-specific casework.

The conclusion is not "uniform marginals": in a nonuniform pencil an occupancy g can still have degree up to p−2−3 as a reduced polynomial. The conclusion also does not make weighted incidence integral by itself: q(v)=Σ_t binom(g(t),2) depends on squares of the profile counts, not on the point second moment Σ_t t²g(t).

## Primary comparison and scope

- [Ferguson–Iosevich arXiv:2610.30712v1](https://arxiv.org/html/2610.30712v1), §8.1 and individual prime-eleven cases, establishes the moments it needs via separate incidence bounds or other conforming case-specific arguments; §4's cubic moment tensors are instead in characteristic two. We did not find this prime-uniform cubic vanishing conclusion in that paper.
- [Fallon–Mayeli–Villano arXiv:1902.02936](https://arxiv.org/abs/1902.02936) supplies basic line/plane and layering geometry for p=5,7; [Malikiosis 2024 §§2–6](https://doi.org/10.1007/s43670-023-00078-2) supplies projective blocking LP. This note uses these as prior foundations.
- Searches for a cubic on every difference direction, singular-locus promotion of an affine set, and universal prime-field cubic moment vanish did not locate matching earlier primary treatment. This is a bounded search, not a verified priority claim.
