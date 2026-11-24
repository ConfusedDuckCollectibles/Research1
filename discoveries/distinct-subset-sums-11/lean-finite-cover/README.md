# Lean verification of the finite exact covers in the one-fifth subset-sum proof

This package checks the finite exact-cover step in the [one-fifth subset-sum proof](https://github.com/ConfusedDuckCollectibles/Research1/tree/96875f5b71501a82cf69fe181194f7b4def0354f/discoveries/distinct-subset-sums-11). It does not formalize the full real inequality or add a new mathematical result.

## What Lean proves, and what it does not

1. **General checker soundness.** `Checker.lean` proves that `check board certificate = true` implies that the given site-indexed board has no exact cover. A residual cover selects mutually disjoint rows avoiding occupied sites, with a selected, registered row through each uncovered required site. The format is a preorder list of sites. At a valid uncovered site, the checker recomputes all compatible rows through that site and consumes one child for each. At a leaf there are no such rows. Trailing indices, missing children, attempts to choose an occupied site, and out-of-range sites are rejected. The proof removes the row chosen by a hypothetical cover and recurses on its corresponding checked child.
2. **Concrete kernel replay.** Each generated theorem `Batch…accepted` is proved by `rfl`, not by `native_decide` or a compiler result. It reduces the Boolean checker on a literal certificate and the full literal candidate list at each site named by that certificate. Rows elsewhere are supplied by `Template.board` in `Rows.lean`, but the certificates never query them. The theorem `Batch…refuted` instantiates soundness: the concrete site-indexed board with those cached rows admits no exact cover. Every named site is checked for presence in the cache before replay. The cache includes all candidate rows through that site, not just those compatible with some specific branch.
3. **Row correspondence.** `Rows.lean` independently implements floor-mask sums, shifted tags, retained windows, and the mask exclusions. `CachedEntry.board_eq_of_agrees` and `cached_batch_transfer` prove, in the kernel, that equality of the cached buckets with this reconstruction suffices to transfer the obstruction to the reconstructed board. The production build checks *all* bucket equalities using compiled Lean `#eval`; it does not use those outputs as proof terms. The optional `bench_agreement.sh` experiment attempts three worst-case equalities by raw kernel reduction. That run did not finish within its 70-second cap, so it does not certify those equality proofs. The compressed ban tables and complete normalized tail enumeration are externally reconstructed by `audit_data.py`, not proved in Lean. The Python audit independently regenerates the cube order, mirrored rank intervals, exact bans, signed tags, complete rows at every cited site, and the complete list of distinct-floor tails, then runs its own bitset verifier. The C++ producer is not trusted to refute an exact cover of the declared board; a false row list, by contrast, would change the board the theorem speaks about.
4. **Analytic bridge.** The formalized layer does not establish that a real separated-sums counterexample gives one of these covers. That deduction still uses the ordinary proofs of the saturated cuts, Harper/Räty neighborhood rigidity, prefix profiles, dimension transfer, bounded rounding, and duplicate-floor dismissal. Notik's small-dimension results are also outside this package. This is not a Lean proof of the full bound or any of those bridges.

For general soundness, the cached-to-reconstructed transfer theorem, and the concrete refutation batches, the actual Lean axiom report is `[propext, Quot.sound]`. No `sorryAx`, `Classical.choice`, `native_decide` axiom, or compiled checker result occurs in these proofs.

## Reproduce

Use **Lean 4.19.0**, release commit `6caaee842e94`, available from the [official Lean release](https://github.com/leanprover/lean4/releases/tag/v4.19.0), and a C++17 compiler, Bash and Python 3.10+ (standard library only). No Mathlib or other Lean package is required. The Linux x86-64 archive used here had byte count 343,842,845 and observed SHA-256 `6fe3ce97a58f44e2b3567d455b994eacec5bfe9ae7774f2a573444480ba813fe`. The official release API response for this old asset did not provide a digest; the value identifies the local build input, not an upstream signed checksum. The replay checks the compiler's reported version and release commit. No network is needed for replay after installation.

```sh
LEAN=/path/to/lean JOBS=4 bash replay.sh
# Optional bounded direct Lean kernel reductions, without cached rows:
LEAN=/path/to/lean bash bench_agreement.sh
```

`replay.sh` compiles `generate.cpp`, regenerates `certificates.json`, verifies its frozen SHA-256 in `corpus.sha256`, runs the independent Python audit, and passes the corpus through `emit_lean.py`. The generated corpus is not checked into the repository. The emitted modules import `Templates` (literal ban arrays) and `Rows` (which imports `Checker`, which imports core Lean). The script compiles the controls, all 81 theorem modules, and a manifest importing them; it checks each compiled row-agreement result and each axiom report. It uses a temporary build directory and keeps no precompiled theorem objects as source inputs. `maxRecDepth` and `maxHeartbeats` only set computation limits, not proof assumptions.

| Template | Distinct-floor certificates | Repeated-floor cases outside Lean cover check | Preorder site indices | Largest certificate |
|---|---:|---:|---:|---:|
| 9 | 62 | 4 | 156 | 10 |
| 10 (first three low sites omitted) | 3,755 | 455 | 6,329 | 12 |
| 11 (reference for larger dimensions) | 3,755 | 455 | 7,406 | 26 |
| Total | 7,572 | 914 | 13,891 | 26 |

For each template, an untouched certificate and truncated, extended, and invalid-first-site variants are tested using **kernel reduction on the non-cached board**. The C++ and independent Python replays both reject all 22,716 such corruptions across the corpus. The all-batch check finished in 4m44.227s with four workers (GCC 14.2.0, Python 3.12.14, Linux x86-64); the independent Python replay completed in 2.216s. The adjacent logs contain both timings and the capped direct kernel row-equality benchmark.

A shorter handwritten argument using only one-step Harper is another candidate for the analytic dependency; it is not verified in Lean here.

A concrete next step for removing the compiled row-agreement boundary is a faster kernel representation of the ban table (for example, a balanced decision tree proved equal to the literal table once for each template), followed by per-query equality proofs and application of the existing transfer theorem. The current slow spot is direct definitional reduction of subset-mask arithmetic and large array lookups, not the certificate logic.

For review, report this as **a sound Lean checker with kernel-checked concrete cached exact-cover refutations and independent row reconstruction checks**, not as an end-to-end Lean proof of the one-fifth inequality or a kernel proof of the full row-data bridge.
