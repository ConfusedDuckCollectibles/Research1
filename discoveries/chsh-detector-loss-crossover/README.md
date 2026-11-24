# How much Bell-test margin does a nearly perfect detector cost?

An exact threshold says whether a violation is possible. It does not
say how far above the classical limit an experiment can get. Near an
asymmetric threshold, treating the better detector as perfect can
substantially overestimate that margin.

We prove a uniform formula for the largest attainable CHSH excess in
this corner. The model counts every trial, assigns +1 to each missed
detection, and assumes independent losses that do not depend on the
measurement setting, outcome or state. It is a model of attainable
quantum scores, not a premise required to interpret a valid violation
in recorded data.

| Residual better-detector loss in tilt units, divided by the distance above threshold | Leading margin relative to the perfect-detector approximation at the same distance |
| ---: | ---: |
| 0 | 100% |
| 0.25 | 68.2% |
| 1 | 36.1% |
| 4 | 12.7% |
| 9 | 6.13% |

The result is the closed-form transition between the known quadratic
and cubic regimes, with a uniform error. Earlier work established
the threshold and the exact optimum, and already identified the
qualitative crossover. This note is not a new detection threshold,
a hardware test, or a guarantee of the amount of entanglement.
See the [proof](proof.md), [uniform error proof](uniform-proof.md), [independent audit](audit.md),
and [finite rational verification](checks.json).

### Reproduce

Checked with Python 3.12 and SymPy 1.14.0.

```sh
python3 verify.py      # Python 3 with SymPy
python3 exact_audit.py # independent replay; Python 3 with SymPy
```

The paper proof establishes the result for all parameter paths in its
domain. The finite tests check the algebra and independently bracket
representative cases; they are not a substitute for the proof.

### Closest prior work

[Gigena et al.](https://doi.org/10.1038/s41534-025-01029-6)
give the exact doubly tilted optimization and derive the loss
mapping. [Pakhunov, §2](https://arxiv.org/html/2607.13774v1#S2)
gives the fixed-interior cubic coefficient and discusses the
transition to the singly tilted case. The determinant and exact
optimization agree with the prior characterization. The contribution
here is the explicit two-scale profile and quantitative, uniform
error across the transition. This note has been independently
replayed but is not peer reviewed.
