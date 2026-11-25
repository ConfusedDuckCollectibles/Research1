# Two gravitational lenses, exactly the same signal along one track

If a source moves behind a gravitational lens, can its changing wave signal tell us how the lens's mass is distributed? **Not always.** In the thin-lens model, this proof constructs different positive, finite-mass lenses that give exactly the same complex amplification at every frequency and at every point along a straight source track. The example preserves geometric image data even in a narrow strip around that track.

The geometry of the measurements matters. For lenses that differ only in a bounded region, a calibrated two-dimensional patch of source positions can distinguish them at just one frequency. For a circular lens in the stated class, the change in the signal near perfect alignment supplies the extra information needed to recover the radial lens profile.

These are exact statements about the model and the data specified in the proof. They are not claims about a particular astronomical observation or how reliably noisy, finite measurements can distinguish real lenses.

- [Proof and primary prior work](PROOF.md)
- [Independent audit](audit/AUDIT.md)
- [Exact replay](REPLAY.md): `python3 verify.py` uses only the standard library; `python3 audit/verify.py` uses SymPy.
