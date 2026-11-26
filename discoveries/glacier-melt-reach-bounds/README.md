# The same total glacier melt can reach different places

Two melt patterns can remove the same amount of ice overall but deliver very different melt beyond an inland ridge. One can reach it repeatedly; the other might never reach it at all.

**What is proved.** For a specific one-dimensional model, in which the melt rate tapers linearly inland, we derive the exact minimum and maximum melt that can occur beyond a fixed inland distance. Knowing total melt and its near-edge intensity is not enough to locate that melt. An extra measurement, the melt-weighted second moment of inland reach, gives an exact threshold for whether any distal melt is unavoidable. It also gives the sharp range when the ridge is in the outer third of the allowed reach.

For a **synthetic** example, allow reaches up to 6 km and fix a melt-weighted mean of 3 km. An unchanging 3 km reach and an intermittent 0-or-6 km reach can have the same total and peak budgets, yet the first melts nothing past 4 km and the second does. These distances are demonstration inputs, not measurements of a real ridge or intrusion.

A separate result checks a common radar comparison: two sets of measurements can have the same reported change in tide and still produce different ice motion under the same nonlinear response. The proof identifies what has to be matched or independently bounded for that comparison to be meaningful.

- [Proof and model assumptions](PROOF.md)
- [Measurement and published model scope](SOURCES.md)
- [Independent mathematics audit](INDEPENDENT-AUDIT.md)

```sh
node replay.cjs
python3 independent-check.py
```

Both run offline with no external inputs. They replay the synthetic examples, not a glacier model. No new Thwaites satellite data were fitted, no sea-level forecast was made, and the result has not been externally peer reviewed.
