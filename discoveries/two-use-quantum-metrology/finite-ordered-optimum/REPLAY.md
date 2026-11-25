# Replay and source boundary

The proof and all included code are original local calculations. The two controls are separate:

1. Exact (offline): run python verify.py, with SymPy. It verifies symbolic identities and constructs actual three-register sensor/tag states for the three-call witness, phase-dependent flag probabilities, a nominal zero-mass branch, and no-record jet alignment. Compare stdout with replay.txt.
2. Optional numerical comparison: run python compare_comb.py --fetch, with SciPy. This anonymously reads one public author data file into memory, authenticates its SHA-256, and compares three fixed control points with the Bellman equations. The file is **not included** in this package. To use an independently downloaded copy instead, pass --source PATH. Compare stdout with comb-replay.txt. This is not a mathematical certificate or a check of the mixed-survival witness.

The exact data pointer is the original authors' repo, revision 07748a4e94d7b29dc888d2cbf41d0b0f421c58dc:

https://raw.githubusercontent.com/qiushi-liu/strategies_in_metrology/07748a4e94d7b29dc888d2cbf41d0b0f421c58dc/data/QFI_seq_rangep_0_1_N_3.txt

Expected SHA-256: 756927f516d70321fbb454538b76d1a1ac95147f8964424768591c134590b4e0. Author plotting metadata in plot.py at the same revision specifies 201 equally spaced decay probabilities from zero through one. Reference paper: https://doi.org/10.1103/PhysRevLett.130.070803.

Run sha256sum -c manifest.sha256 to check the included original artifacts. Tool versions used here appear in dependencies.txt. Tiny machine-level rounding differences in the optional comparator can alter the last displayed digit; the script independently checks an absolute tolerance of 5e-6.
