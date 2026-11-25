# The final linked-layer contradiction at size 78

Step 4 of the size-78 exclusion. Independently reproduced; in combination with the preceding nonlayer, mixed and layered-structure proofs, this excludes size 78 in F13³. Inputs and final scope are separate from the implication that the other spectral cardinalities are already decided.

The [earlier checked reductions](../03-linked-layers/PROOF.md) put a hypothetical spectral pair in linked dual coordinates
```
E = { (x,f_s(x),s) : x in F13, s in S },
A = { (u,t,g_t(u)) : u in F13, t in T },
|S|=|T|=6,
```
with line cap five in both. In the nonhorizontal difference chart, let I_r be the set of vertical slopes (f_s(x)−f_{s'}(x'))/(s−s') across all pairs with (x−x')/(s−s')=r. A slope r is active if some complex column C_t(r)=Σ_u ζ^(ru+g_t(u)) is nonzero. The earlier Chebotarev argument gives |I_r|=5 on active r and at least ten active r; two suffice here. For every x₀ and active r, the plane x−rs=x₀ has one E point from each layer, so gives a six-point graph with abscissa set S whose slopes belong to I_r. That graph has at least five directions by the audited planar size-six/cap-five gate. Its finite direction set is therefore exactly I_r.

Each such planar graph is injective in s. The audited six-point abscissa atlas implies that S is affinely equivalent to one of the five supports with masks 111, 187, 215, 311, 427 (a 13-bit mask lists the indicated field elements). Put S in that representative once and for all, by an affine change to E's layer coordinate and the dual change on A. These changes preserve spectrality and the linked setup; the active chart is then defined in the new coordinates. There is no relabelling per r or per slice.

## One same pair for every feasible five-slope palette

For a listed representative S and each five-element I⊂F13, enumerate all functions h:S→F13 with h(0)=0 and exactly slope set I. Once h(0)=0, each h(s) must belong to sI for s≠0, by its difference with the origin. The exact recursion appends each allowed height and checks every pair slope against I. It accepts if all five slopes have occurred. This covers all functions, all 1287 palettes per support, and does not assume a symmetry of I. The result is:
```
S mask   number of feasible I   heights per feasible I   fixed layer pair
 111             156                       1                (0,1)
 187              78                       2                (3,5)
 215             156                       1                (0,1)
 311             156                       1                (0,1)
 427              78                       2                (3,5)
```
In both two-shape cases, h(5)−h(3) is the same on the two shapes for *each* fixed I, checked exhaustively. Thus for every S there is a single pair of different labels s,s′ whose height difference is determined by I. The pair is fixed independently of I; the value of the difference can vary with I. A common height translation does not change it.

## Constant nonzero-step derivative contradicts the cap

Fix two distinct active r,r′. For each fixed r and each x₀, subtract the slice height at label zero to put its six points into the checked class of functions for palette I_r. For the fixed same label pair s,s′ this gives, for some C_r,
```
f_s(x+rs) − f_{s′}(x+rs′) = C_r    for every x in F13.
```
Put y=x+rs: f_s(y)=f_{s′}(y+r(s′−s))+C_r. Compare with r′. The finite derivative of f_{s′} at displacement (r−r′)(s′−s) is a constant. The step is nonzero, so its multiples traverse F13. Iterating the equation makes f_{s′} affine, even if the constant is zero. Its graph then contains 13 collinear points, contradicting line cap five. Therefore no linked pair exists, closing the final size-78 case.

The adjacent exact JavaScript checks every finite palette. The independent checker uses affine height normalization followed by direct residual-height enumeration, then restores the omitted shear and scale.
