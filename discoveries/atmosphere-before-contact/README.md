# The cost of a clean planet-atmosphere test

A mottled star complicates the interpretation of a dip during an exoplanet transit. Just before the opaque planet overlaps its star, however, stellar spots cannot cause an occultation by that core. An atmosphere extending past the core could already absorb starlight. That is a cleaner geometric comparison, but how much signal can it carry?

In a frozen TRAPPIST-1 e-scale benchmark, ten perfectly observed transits still cannot attain 80% power at a 5% false-positive rate by using only the two brief exterior-contact windows. The certified power ceiling is 59.1%. The result uses the finite star and planet geometry, even gives the analysis the correct baseline and contacts, allows a conservative independent WISE-derived photon count, and assumes perfect throughput with no additional measurement noise. It gives at least 18 perfect visits as a necessary, not sufficient, condition under the stated parameters.

This is a limit for a specific test, not a result about whether TRAPPIST-1 e has an atmosphere. Brightness at the edge of the star, whether the true core contact is known, how far the absorber extends, and the photon-weighted strength at the observation's spectral resolution all matter. The broad in-transit spectrum and the intervals of partial core overlap could contain additional information under independent stellar constraints.

- [Exact geometry and statistical bound](PROOF.md)
- [Benchmark frozen before looking up photon calibration](FROZEN.md)
- [Primary work and native calibration provenance](SOURCES.md)

Offline: `python3 -B verify.py | diff -u replay.txt -`. To refetch and hash-check the two public calibration files: `python3 -B verify.py --native | diff -u native_replay.txt -`. The replay does not read or fit a science transit.
