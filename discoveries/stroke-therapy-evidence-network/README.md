# A reproducibility check on a stroke-treatment evidence network

**Two reported details in one 2026 review do not match the randomized trial network described in the paper.** The INSPIRES trial is assigned to the 90-day combination-treatment group even though the trial stopped one of the two active drugs after day 21. Separately, the review reports an independent indirect comparison that cannot be formed from its listed trials.

Scope: published reporting in [Zou et al., 2026, article 108710](https://doi.org/10.1016/j.jstrokecerebrovasdis.2026.108710). This is a source and graph audit, not a new clinical analysis or a treatment recommendation.

## 1. Treatment duration and follow-up are different

The [review’s full text](https://www.strokejournal.org/article/S1052-3057%2826%2900173-4/fulltext) lists INSPIRES as a 21-day combination regimen in Table 1, but the Network geometry section groups it with POINT in the 90-day combination node. The [original INSPIRES trial](https://doi.org/10.1056/NEJMoa2309137), Trial Procedures, states that the combination arm continued clopidogrel to day 90 but replaced aspirin with placebo after day 21. POINT randomized the 90-day combination strategy. An outcome measured at day 90 does not change either trial’s exposure schedule.

## 2. The claimed indirect path uses the very comparison it is meant to check

Write A for aspirin; C21 and C90 for 21- and 90-day clopidogrel plus aspirin strategies; and T21 and T30 for the two distinct ticagrelor plus aspirin strategies. TM is ticagrelor monotherapy.

| Edge in the published primary network | Trial reports |
| --- | --- |
| A to C21 | [CHANCE](https://doi.org/10.1056/NEJMoa1215340) |
| A to C90 | [POINT](https://doi.org/10.1056/NEJMoa1800410), [INSPIRES](https://doi.org/10.1056/NEJMoa2309137) as classified in the review |
| A to T30 | [THALES](https://doi.org/10.1056/NEJMoa1916870) |
| A to TM | [SOCRATES](https://doi.org/10.1056/NEJMoa1603060) |
| C21 to T21 | [CHANCE-2](https://doi.org/10.1056/NEJMoa2111749), [PRINCE](https://doi.org/10.1136/bmj.l2211) |

The review describes [TARDIS](https://pubmed.ncbi.nlm.nih.gov/29274727/) as a separate subnetwork, excluded from this primary analysis. The five distinct edges connect six nodes in a tree. The only neighbor of T21 is C21. A [standard node-splitting check](https://doi.org/10.1002/jrsm.1167) holds out the direct studies and asks whether the remaining randomized comparisons connect the same pair. Remove CHANCE-2 and PRINCE and T21 is disconnected. **There is no remaining randomized path from C21 to T21.** Table 3’s note describes reaching T21 from aspirin via CHANCE and CHANCE-2, already using one of the direct studies. Correcting where INSPIRES sits leaves this conclusion unchanged.

Merging T21 and T30 would create a connected indirect comparison through THALES. That would be a different, coarser question, no longer comparing the two durations as distinct treatment strategies. Trials with a first-21-day analysis also need explicit outcome and follow-up alignment, not merely relabeling.

## Reproduce the check

Run `node verify.cjs` with no dependencies. It builds the published graph, the graph with INSPIRES reclassified, and a positive control that merges the two ticagrelor combination nodes. See [expected output](expected.json).

The audit uses the journal’s publicly indexed article text and original trial reports. The supplement was listed but not independently retrieved; the fitted model and inputs were not available to this audit. No linked correction was found as of November 27, 2026. The statements above identify inconsistencies in the published description. They do not establish what the model actually used, provide corrected effect estimates, assess author intent, or invalidate the underlying trials. Not externally peer reviewed.

Research direction and public archive: James Bazan / ConfusedDuckCollectibles. AI-assisted source audit and independent computational verification.
