# Three endpoint-label and denominator checks in reports on refractory out-of-hospital cardiac arrest

**Source verification: 28 November 2026.** This note reconciles published aggregate descriptions. It does not reanalyze patient data, compare treatment effects across trials, or test the validity of the 2023 pooled outcome. Throughout, arms are written as **standard/conventional, then assigned intervention**, regardless of the order in the original publication.

## 1. Prague OHCA baseline rhythm: 39%, not 49%

The [Morrison et al. NHLBI workshop report, JAHA 2025](https://www.ahajournals.org/doi/pdf/10.1161/JAHA.124.036108), printed PDF page 5, in the text continuing the description of the Prague trial from page 4, reports 49% with a nonshockable initial rhythm. The [full-text copy in PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC12132745/) contains the same discussion and table. Its Table 1 lists shockable proportions of 64% (standard) and 58% (intervention).

The [Prague OHCA primary randomized report](https://jamanetwork.com/journals/jama/fullarticle/2789313), Table 1, “Initial rhythm” (as determined by EMS), gives these integer cells:

| Assigned arm | All analyzed | Ventricular fibrillation | Asystole | Pulseless electrical activity | Nonshockable total |
|---|---:|---:|---:|---:|---:|
| Standard | 132 | 84 | 24 | 24 | 48 |
| Invasive strategy | 124 | 72 | 31 | 21 | 52 |
| Both | 256 | 156 | 55 | 45 | 100 |

Thus 100/256 = 39.0625%. The trial’s [same-cohort rhythm analysis](https://pubmed.ncbi.nlm.nih.gov/36243225/) explicitly gives 156 shockable and 100 nonshockable, and defines “initial” as the first documented EMS rhythm in its Methods. This is not a comparison of clinical outcomes. Nor can 49% be recovered by switching to the 264 recorded randomizations shown in the primary trial’s Results/Figure 1: even treating all eight who were not in the analyzed cohort as initially nonshockable gives at most 108/264 = 40.91%. That ceiling is an arithmetic control, not an imputation of the eight rhythms. The workshop paragraph points to the enrolled trial population, not a separately defined screened or postrandomization subgroup. What generated the 49% figure is not established.

## 2. INCEPTION follow-up: two identical headings, different native visit rows

The [workshop report’s Table 1, printed page 5](https://www.ahajournals.org/doi/pdf/10.1161/JAHA.124.036108#page=5) contains two consecutive rows both headed “Favorable Neuro at 180 d (CPC 1, 2), n (%).” This table gives counts plus rounded percentages, not a denominator for each measured visit. Compare them with [Suverein et al., INCEPTION, NEJM 2023, updated January 26, p. 307, Table 4](https://www.nejm.org/doi/pdf/10.1056/NEJMoa2204511?articleTools=true#page=9). See [the same published Table 4 hosted by Reanimate](https://reanimateconference.com/wp-content/uploads/2023/04/Prehosp-ECMO-RCT.pdf#page=9) if the publisher PDF is inaccessible.

| Source and time label | INCEPTION conventional | INCEPTION extracorporeal CPR |
|---|---:|---:|
| Workshop first “180 d” row | 9 (14%) | 12 (18%) |
| Workshop second “180 d” row | 10 (16%) | 14 (20%) |
| Native Table 4: 30 days | 10/62 (16%) | 14/70 (20%) |
| Native Table 4: 3 months | 9/63 (14%) | 12/68 (18%) |
| Native Table 4: 6 months | 10/63 (16%) | 14/70 (20%) |

The first workshop pair matches only the **3-month** row among these native data, not the published 6-month row. The second workshop pair is consistent with the **6-month** row, but the printed information is identical to the **30-day** numerators and rounded percentages; the workshop also prints this pair under its 30-day heading. Its precise source visit cannot be recovered from those cells alone. The native Methods explicitly put the 6-month assessment after the arrest; it should be cited at its stated visit, not assumed to fall on a common fixed randomization day in every trial.

The first Prague pair printed alongside is the full cohort (29/132, 39/124; [primary report Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC8864504/)). The second Prague pair is the initially shockable subgroup (28/84, 35/72; [same-trial rhythm analysis](https://pubmed.ncbi.nlm.nih.gov/36243225/)). That establishes the Prague cells’ provenance; it does not make the adjacent INCEPTION cells proven shockable-only subgroup estimates. We did not infer a common across-trial population from either heading.

## 3. ARREST: survival to discharge is not survival to discharge with CPC 1–2

The [2023 Prague–ARREST pooled report, Introduction](https://pmc.ncbi.nlm.nih.gov/articles/PMC10184044/) describes the original ARREST results as 7% versus 43% with **neurologically favorable** survival to discharge (in the arm order used here). [Yannopoulos et al., original Lancet ARREST Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC7856571/) ([published-layout public copy](https://broomedocs.com/wp-content/uploads/2021/01/10.1016%40S0140-67362032338-2.pdf#page=4)) instead labels 1/15 and 6/14 as the primary outcome, **survival to hospital discharge**. Its Results and Figure 2C report the lone control survivor's CPC as 4 at discharge. Table 2 reports the six intervention survivors’ discharge CPC mean (SD) as 2.5 (0.5). The [2024 ILCOR CoSTR summary, Table 12, Yannopoulos row](https://doi.org/10.1161/CIR.0000000000001288) ([visible published-layout copy](https://www.cardiacresponse.ca/uploads/b/4806fab0-ec83-11ec-a69a-172c716404c3/53c6b060-abe7-11ef-9844-155b324cf894.pdf#page=25)) lists **0/15 and 3/14 (21%)** for favorable functional outcome separately from the survival counts. The table uses a combined discharge/30-day column across trials; the [official ILCOR topic’s Consensus on Science](https://costr.ilcor.org/document/extracorporeal-cardiopulmonary-resuscitation-ecpr-for-cardiac-arrest-als-sr) explicitly identifies the ARREST observation as **at discharge**. The original aggregate CPC row also independently gives the same count if its one-decimal summary uses ordinary rounding: the six integer survivor scores sum to 15. If a score were 1 or 4, the sum of squared deviations from 2.5 would be at least 7/2, making either sample or population SD too large to round to 0.5. Hence all scores are 2 or 3, three each.

This check identifies a mismatch in the reported *historical discharge description*, not a finding of erroneous pooled 180-day data. ARREST’s original [Figure 2 caption, Methods and Table 2](https://pubmed.ncbi.nlm.nih.gov/33197396/#figures) distinguish the survival curve (whose horizontal axis is time since randomization) from CPC/mRS measured at discharge and at 3 and 6 months **after discharge**. The pooled study specifies a 180-day CPC 1–2 endpoint and says it obtained individual outcome data and definitions from both trials. The main publication does not document the exact visit alignment in enough detail to audit it here. This published-source comparison cannot exclude valid alignment in the analysis dataset and gives no corrected estimate of that endpoint.

## Reproducibility and limits

The replay [`replay.cjs`](replay.cjs), run with `node replay.cjs` and compared against [`expected.json`](expected.json), checks the arm sums, worst-case Prague denominator control, all three INCEPTION visit possibilities at the displayed precision, and the original rounded ARREST CPC moment consistency. It uses only hand-read public aggregate integers and does not fetch documents or inspect any individual record. No source PDFs, figures or publisher tables are redistributed here; links identify the specific source objects. The packet records observed versions, including the visibly updated January 26 NEJM version. Searches by journal DOI, exact title, “correction,” and “erratum,” and checks of indexed versions and available PubMed records did not locate a notice correcting these specific three items by the verification date. That search is not proof that no correction exists.

These three comparisons are useful for citation and endpoint extraction. They are not a clinical methods validation, new treatment-effect discovery, recommendation about resuscitation, or reanalysis of pooled data.

AI-assisted source reconciliation, independently checked. Not externally peer reviewed.
