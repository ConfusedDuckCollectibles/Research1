# When a nearly perfect conductor still has a finite attenuation length

A one-dimensional lattice with disordered **pairs** of sites can transmit perfectly at one special energy. A rare break in the pairing destroys that perfect-conduction behavior in a long enough sample. Yet a finite sample can still look perfect because it contains no breaks, or because scattering from two breaks interferes constructively.

The result proves when this lattice produces a previously modeled random-phase transmission law. The expected number of breaks stays fixed as the chain grows. The number of breaks tends to a Poisson law. Crucially, the relative phase between them becomes random in a precisely derived limit: for two breaks, the same defect strength can give transmission anywhere from about 22% to nearly 100%. This is not described by multiplying a fixed transmission loss per break.

There are two exact model benchmarks at 400 independently disordered two-site blocks per chain and an expected one break per chain. Under the tightly calibrated energy window, four independently prepared chains give **over 94.38% detection power** with a hard zero false-alarm bound. With a 50-fold wider energy tolerance, eight independently prepared chains give **over 95% power** and **under 2.8% false alarms**, if the probe energy is fixed independently of the disorder. Both benchmarks allow an independent hard calibration bound of two percentage points on measured transmission. These are predictions for the declared model, not measured instrument tolerances.

The calculations cover the specified coherent one-particle chain and matched contacts. They do not establish localization in any material or explain finite-time arrest in interacting systems.

[Proof](PROOF.md) · [Frozen model](FROZEN.md) · [Prior work](SOURCES.md)

Offline exact replay (Python standard library): `python3 -B verify.py | diff -u replay.txt -`.
