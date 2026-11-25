# When two different models of battery aging give the same readings

Can a partial voltage scan and impedance at one charge level tell *how* a battery aged? In the idealized model here, not always. Two mathematically feasible states give exactly the same measured voltage throughout the window and the same linear impedance spectrum at its center. Their inferred loss of available lithium differs by **14.08 percentage points**. The inferred material losses at the two electrodes differ by **16 points each**.

The proof explains why this can happen even when both fits look well-conditioned locally. In the model, an independently calibrated electrode reference can resolve it. For quadratic electrode curves, a calibrated split measurement at a second charge level also resolves it. The model specifies what calibration and electrode information those statements require.

No battery was tested: the values are exact synthetic examples, not measurements or charging advice. The linked public lab archive documents relevant channels, but its native raw rows were not fitted or used to validate chemical losses.

- [Proof, assumptions, and primary sources](PROOF.md)
- [Independent audit](independent/AUDIT.md)
- Replays: `node check.js` and `node independent/check.js`, with expected JSON outputs.

This note has not been peer reviewed.
