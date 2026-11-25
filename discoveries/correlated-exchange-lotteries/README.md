# A cap on the value of a private lottery under correlated exchange failures

**Synthetic theorem, not a clinical proposal.** For unit value per performed exchange, certain starter links, half-success marginal risks on all non-loop paired arcs, simple chains of at most three recipients, and no fallback, the best lottery's worst expected count can be **at most 10/9** of the best pure plan's worst expected count. The bound is exact. It uses one fixed but otherwise arbitrarily correlated edge law, chosen knowing the lottery but before its independent draw. It is not a claim about actual success rates.

The mechanism is structural. An adversary can kill the late step for any single chosen chain by placing its two fragile links on opposite outcomes of one fair bit. It cannot always do that for all alternative chains with one common law. In the one-chain/no-cycle class, project each full chain onto a pair of fragile arcs. A lottery gives a positive bonus exactly when these pairs form a nonbipartite graph. Acyclicity of the original exchange graph is not enough to rule out a bonus. The exact size of the bonus is determined by a standard fractional cut-cover parameter of this derived graph. PROOF.md gives the general theorem, topology certificate, exact probability-floor calculation and failure modes.

| Synthetic cap-3 benchmark | Best pure worst expected count | Best mixed worst expected count |
|---|---:|---:|
| Directed triangle (optional atomic cycle allowed) | 3/2 | 5/3 (+11.11%) |
| Five-arc acyclic shortcut graph | 3/2 | 8/5 (+6.67%) |
| Four-edge even circle or simple line, chains only | 3/2 | 3/2 (zero) |
| Triangle, with separate receipt probability floors | No feasible pure plan at those floors | 49/30, vs **unconstrained** pure 3/2 |

The synthetic floor example reaches minimum individual receipt probabilities .60, .50, .40 while counting every performed exchange equally. It does not assume protected classes, weighted clinical priority, truthful participation or value compensation. The triangle's unprioritized bonus also survives when marginal success is varied across [.45,.55]. The universal ratio does not extend by assertion to these different marginals.

**Longer synthetic falsifier.** We checked all 4,096 edge masks for the complete bidirected cap-4 model and 1,048,576 for cap-5 (single chain only in both). Their exact worst-law lottery counts are 5/3 and 17/10, below the corresponding rank-law upper bounds 41/24 and 103/60. The cap-three exact extremum cannot simply be extrapolated by asserting complete longer graphs attain the factorial bounds. This is a synthetic holdout, not external empirical validation.

## Nearest primary work and priority limits

- [PICEF](https://arxiv.org/abs/1606.01623), [failure-aware exchange](https://doi.org/10.1287/mnsc.2018.3026), [inhomogeneous risks and CVaR](https://proceedings.mlr.press/v124/bidkhori20a.html), and [recent improved failure-aware formulation preprint](https://ssrn.com/abstract=6397118) establish the framing around chain length and cancellations. This product does not improve a deployment solver or estimate risks.
- [Anderson et al., 2015 PNAS](https://doi.org/10.1073/pnas.1421853112), [SI Stochastic Optimization section](https://rma350.github.io/assets/publications/PNAS-2015-FindingLongChainsInKidneyExchangeTSP.pdf), already permits *arbitrary* dependence when one can jointly sample edge states, and includes testing/reoptimization. [Fallback alternatives](https://pmc.ncbi.nlm.nih.gov/articles/PMC6136670/) also preexist. Our restriction instead has only marginals, a fixed independent draw and no recourse.
- [St-Arnaud, Carvalho and Farnadi](https://doi.org/10.1016/j.ejor.2025.02.014) already optimize lotteries over plans and ex ante fairness. [Zhang, Dai and Li](https://doi.org/10.1007/s12561-026-09514-y) use lottery-based calibration and independent-edge random graphs in their theoretical price result; [2026 robust participation contracts](https://doi.org/10.1016/j.econlet.2026.112934) model ambiguity-aware compatible-pair participation in a different dynamic setting. “Randomization for fairness” is not new.
- [McElfresh, Dickerson, Ren and Bidkhori, *Distributionally Robust Cycle and Chain Packing With Application to Organ Exchange*](https://doi.org/10.1109/WSC52266.2021.9715358) already applies distributionally robust optimization to uncertain match quality using CVaR and moment or data-driven models. [McElfresh, Bidkhori and Dickerson, *Scalable Robust Kidney Exchange*](https://cdn.aaai.org/ojs/3899/3899-13-6958-1-10-20190702.pdf) instead models discrete hard failure budgets. Neither is the common marginal-only binary law with a private independent plan draw used here.
- [Delage and Saif](https://doi.org/10.1287/ijoc.2020.1042) already establish broader gains and bounds from independent randomized decisions in DRO, and [Šámal](https://doi.org/10.46298/dmtcs.2134) already studies fractional cut covers and equivalent weighting/dual mathematics. The narrowly distinct contribution checked here is the simple-path rank-law ceiling, exactly attaining directed/acyclic constructions, and the resulting reporting/priority interpretation for prefix-only exchange. Priority is not established by this bounded search.

## Reproduce

Exact fractional proofs and the exact combinatorial outputs use no scientific Python dependencies. The independent small-case check uses NumPy/SciPy to replay the all-state and all-cut linear programs and check the closed-form rational answers. Two independent C++ enumeration approaches reproduce the complete-graph minima and explicit adverse masks. No real matching, private or restricted data, or deployment. This is an AI-assisted research note; it has not been peer reviewed.

The [independent audit](independent/AUDIT.md) verifies the proof and replays its own exact small-state and orbit enumerations.

    set -o pipefail
    sha256sum -c LOCK.sha256
    sha256sum -c SHA256SUMS
    g++ -O3 -std=c++17 complete.cpp -o /tmp/research1-kpe-complete
    g++ -O3 -std=c++17 brute.cpp -o /tmp/research1-kpe-brute
    /tmp/research1-kpe-complete | diff -u complete.txt -
    /tmp/research1-kpe-brute | diff -u complete.txt -
    python -B replay.py | diff -u replay.txt -
    python -B independent_lp.py | diff -u independent.txt -
