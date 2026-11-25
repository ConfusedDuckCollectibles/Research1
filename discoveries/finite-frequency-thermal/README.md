# What repeated heating can tell you about steady warming

Suppose you heat a surface in known cycles and measure both the size and timing of its temperature response. Can you infer how strongly it loses heat at equilibrium? Even if the measurements are exact and the total hidden heat capacity is known, different passive models can match every observed cycle and disagree about that equilibrium coefficient.

This proof gives the **largest possible disagreement** from any finite set of distinct frequencies in the stated reciprocal linear model. It constructs pairs of physically valid thermal networks that match all the measurements and approach that maximum. If the hidden heat-transfer rates are known to lie in a strictly positive finite band, the sharp maximum is attained. An explicit estimator reaches the best worst-case absolute error. The proof also shows why tightly clustered frequencies cannot beat the single-frequency limit under a stated bounded calibration-error model.

This is an identifiability result for an idealized heat-transfer model, not an estimate of Earth's response or a claim about all climates.

- [Proof and primary prior work](PROOF.md)
- [Independent proof audit](audit/AUDIT.md)
- [Replays](REPLAY.md): `python3 verify.py` and `python3 audit/verify.py` (standard library only).
