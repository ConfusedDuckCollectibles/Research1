# Reproduce the finite checks

Use Node.js 18 or newer. There are no npm packages or outside input datasets. From this directory run:

```sh
node replay.cjs
```

The runner executes all twelve scripts and checks their complete deterministic JSON records against the bundled exact transcripts; it prints only a status line and elapsed time per gate. Give the runner one or more gate IDs to replay a subset, for example `node replay.cjs planar-hexagon marked-quintic`. All phases of the default scripts are required for their stated exhaustive gates. In particular the seven-uniform script runs `g,6,7,8,9`, and the exactly-three-rich script defaults to the entire half-open interval `[0,19)` with all branches. The all-gate run takes several minutes on a typical single core; runtime is not used in any inference.

The finite scripts prove their stated negative results by complete enumeration within the mathematical reductions, not by an elapsed-time cutoff. The last two horizon obstructions are analytic, and no nonsplit-orbit or height census enters them. See the [proof](PROOF.md) for where each arithmetic gate is needed.

<a name="a1"></a>

## Large pencils · [A1](A-AFFINE.md#a1)

- **initial-pencils**: Initial pencil gates through K = 9. [Source](verify/initial-pencils.cjs), [expected JSON](expected/initial-pencils.json), 4 records.
- **eight-uniform**: K = 8, degrees 6, 7, 8. [Source](verify/eight-uniform.cjs), [expected JSON](expected/eight-uniform.json), 1 record.

<a name="a2"></a>

## Seven uniform normals · [A2](A-AFFINE.md#a2)

- **seven-uniform**: Seven-point planar arc gate and K = 7. [Source](verify/seven-uniform.cjs), [expected JSON](expected/seven-uniform.json), 5 records.

<a name="a3"></a>

## Planar hexagon palette · [A3](A-AFFINE.md#a3)

- **planar-hexagon**: Five- and six-point planar caps. [Source](verify/planar-hexagon.cjs), [expected JSON](expected/planar-hexagon.json), 1 record.

<a name="b1"></a>

## No rich family · [B1](B-MOMENTS.md#b1)

- **no-rich-family**: Boolean quotient with no rich line family. [Source](verify/no-rich-family.cjs), [expected JSON](expected/no-rich-family.json), 5 records.

<a name="b2"></a>

## At most two rich families · [B2](B-MOMENTS.md#b2)

- **at-most-two-rich**: Boolean quotient with at most two rich line families. [Source](verify/at-most-two-rich.cjs), [expected JSON](expected/at-most-two-rich.json), 38 records.

<a name="b3"></a>

## Exactly three rich families · [B3](B-MOMENTS.md#b3)

- **exactly-three-rich**: Boolean quotient with exactly three rich line families. [Source](verify/exactly-three-rich.cjs), [expected JSON](expected/exactly-three-rich.json), 20 records.

<a name="b4"></a>

## Marked quintic · [B4](B-MOMENTS.md#b4)

- **marked-quintic**: Marked five-plane, exterior pair, and quintic gate. [Source](verify/marked-quintic.cjs), [expected JSON](expected/marked-quintic.json), 1 record.

<a name="c2"></a>

## Quadrics · [C2](C-CURVES.md#c2)

- **quadric-controls**: Exact finite binary-form inputs and additional quadric controls. [Source](verify/quadric-controls.cjs), [expected JSON](expected/quadric-controls.json), 1 record. Includes an exact hexagon binary-form table used in C2; remaining outputs are cross-checks.

<a name="c3"></a>

## Cubic linkage · [C3](C-CURVES.md#c3)

- **quartic-twisted-cubic**: Quartic-plus-twisted-cubic branch. [Source](verify/quartic-twisted-cubic.cjs), [expected JSON](expected/quartic-twisted-cubic.json), 1 record.

<a name="c4"></a>

## Catalecticants and genus two · [C4](C-CURVES.md#c4)

- **genus-two**: Genus-two point bound over the prime field. [Source](verify/genus-two.cjs), [expected JSON](expected/genus-two.json), 1 record.
- **catalecticant-controls**: Supplementary rank and binary convention checks. [Source](verify/catalecticant-controls.cjs), [expected JSON](expected/catalecticant-controls.json), 1 record. Diagnostic only: checks representative local algebra; the rank floor in C4 is proved symbolically.

The full transcripts include all fields emitted by the scripts. The initial pencil script omits an unrelated cardinality overview; its four mathematical records are unchanged. Timing shown by the runner is supplementary. The seventh-moment, residue pairing, and reduced/nonreduced horizon case analysis in [D](D-HORIZONS.md) do not require computation.
