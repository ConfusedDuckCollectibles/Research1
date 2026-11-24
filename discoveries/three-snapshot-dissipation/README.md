# Detecting hidden dissipation from three snapshots

A system can dissipate energy even when the part you observe looks the
same forwards and backwards. This note proves that three equally
spaced readings of a binary signal can sometimes place a sharp lower
bound on that hidden dissipation. It uses the population probabilities
of the readings, not a count of everything the system did in between.

In a four-state mathematical test, the new bound accounted for **94.4%
of the true dissipation rate**. The earlier complete-positivity bound
accounted for **74.8%** using the same observations. This is a constructed
test case, not a measurement or evidence that real experiments will
achieve those percentages. A zero bound does not prove equilibrium.

The result assumes a stationary finite-state Markov jump process and
ordinary time reversal. The contribution is the equal-marginal closed form of the Bernoulli
optimization for this three-snapshot problem. The underlying
variational program, balanced special case and related persistence tests are prior work.
See the [proof and prior comparisons](proof.md) and [independent audit](audit.md).
Not peer reviewed.

### Reproduce

```sh
python3 verify.py
python3 independent.py
```

Both replays are mathematical checks. Any inference from a finite
record would need an appropriate uncertainty and dependence analysis.
