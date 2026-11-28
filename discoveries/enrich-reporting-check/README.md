# Two source-label discrepancies in reporting ENRICH

This reproducibility note documents an arm-label discrepancy in an ENRICH trial table and a follow-up timing discrepancy in one guideline sentence. No treatment or pooled-analysis conclusions are drawn.

## 1. Observed subgroup denominators in the trial report

Pradilla et al. (NEJM 2024;390:1277–1289, [DOI 10.1056/NEJMoa2308440](https://www.nejm.org/doi/full/10.1056/NEJMoa2308440)) report 147 surgery and 139 control participants with an observed modified Rankin scale (mRS) score at day 180. The two location rows in Table 2 instead put the observed subgroup counts under the opposite arm headings:

| Location | Randomized, S / C (Table 1) | Observed at 180 days, S / C (Results) | Observed as printed, S / C (Table 2) |
| --- | ---: | ---: | ---: |
| Lobar | 102 / 106 | 100 / 99 | 99 / 100 |
| Anterior basal ganglia | 48 / 44 | 47 / 40 | 40 / 47 |
| Total | 150 / 150 | 147 / 139 | 139 / 147 (sum of subgroup rows) |

The Table 2 overall row itself is 147 / 139. The printed location rows reverse those arm totals, and 47 observed anterior basal ganglia controls would exceed the 44 randomized controls in that location. Exchanging the four subgroup count cells reconciles both checks. The investigators’ [AANS abstract](https://2023aansannualmeeting.eventscribe.net/fsPopup.asp?PresentationID=1237542&efp=QkJBUkFWVkQxNTIwMA&mode=presInfo) matches the Results counts in another format from the same trial. This identifies a discrepancy in displayed counts, not in the source data, randomization or Bayesian efficacy analysis.

**Where to verify:** Results, p. 1283; Table 1, p. 1284; Table 2 count cells, p. 1286 (footnotes continue on p. 1287). Both subgroup rows and the Results are visible in the publicly indexed NEJM article, which records an update on May 10, 2024. They also appear in the [May 10 updated journal-layout PDF hosted by the American Heart Association](https://www.heart.org/en/-/media/Files/Professional/Quality-Improvement/ICH-Stroke_Stryker/ENRICH-Trial.pdf?sc_lang=en) and the [earlier journal-layout copy hosted by Emory](https://med.aws.emory.edu/assets/enrich-primary-paper-nejm-2024.pdf). Direct retrieval from the live NEJM site was blocked in the audit environment; the current indexed excerpt and the May 10 PDF were checked separately.

## 2. Timing of the guideline’s ENRICH mortality comparison

The ESO/EANS [guideline](https://journals.sagepub.com/doi/pdf/10.1177/23969873251340815), PICO 4.1.3, p. 1037, identifies ENRICH counts of 30/147 versus 35/139 as deaths by day 30. The [open journal record](https://pmc.ncbi.nlm.nih.gov/articles/PMC12098356/#fig37-23969873251340815) has the same text. ENRICH Table 3, p. 1288, instead reports 14/150 versus 27/150 as the 30-day safety endpoint, measured from enrollment (Methods, p. 1281). The 30 and 35 deaths are reported at day 180 (Results, p. 1284); the denominators 147 and 139 are those with observed day-180 mRS scores (Figure 2, p. 1285). ENRICH’s separate day-30 mRS row has 147 and 146 observed scores, not 147 and 139. Restricting the sample cannot account for the discrepancy at the same day-30 horizon: 30 and 35 each exceed the corresponding full assigned-arm day-30 event count.

The guideline’s single-trial crude odds ratio is exactly \(16/21=0.761905\), with log-Wald 95% interval 0.4376–1.3267, from the 180-day observed counts. Table 3’s day-30 counts instead give \(287/612=0.468954\), interval 0.2352–0.9350, under the same simple arithmetic. These calculations identify the inputs; they are not new inferential tests and do not replace the original analysis.

**Scope of the pooled figure:** The guideline describes Figure 37 as pooling mortality at **1–12 months** (p. 1036); the figure caption does not call it a day-30 analysis. The [ESO supplementary PICO 4.1.3 mortality forest table and GRADE profile](https://eso-stroke.org/wp-content/uploads/2025-ESO-ICH-guideline-ESOC.pdf) explicitly use the 1–12-month scope and show the ENRICH row with 30/147 and 35/139. The verified discrepancy is the time label in the ENRICH-only sentence. It does not establish an error in pooled input selection, the pooled estimate, or the recommendation.

Exact DOI/title and correction/erratum searches on November 28, 2026 located no notice resolving either discrepancy. This is a bounded search, not a claim that no correction exists. The [replay](replay.cjs) uses only public printed counts and does not use participant-level data.

```sh
node replay.cjs > actual.json
diff -u expected.json actual.json
```

Compare the [expected output](expected.json). AI-assisted, independently rechecked internally, and not externally peer reviewed. This is not clinical guidance.
