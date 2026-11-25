# Offline check

Run python verify.py using Python 3. No data, dependency install, optimization or external connection is required. Compare stdout with replay.txt.

The script certifies rational signs of the relevant Airy series, with exact Taylor and Gamma-integral tail bounds. It then certifies the finite-channel wall excess without relying on a floating-point PDE discretization. The exact operator and asymptotic results are the proof in PROOF.md, not outputs of the script.

Physical scope: linear passive-scalar Fourier amplitude with prescribed time-dependent Couette shear, homogeneous absorbing/reflecting scalar walls, and the L² worst-case initial condition. The maintained claim is the sharp finite-width estimate and the finite-time wall singular-mode characterization, not general hydrodynamic stability or fixed-preparation echo measurements.

Run sha256sum -c manifest.sha256 to check the local core.
