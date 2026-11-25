# Counting recorded deaths from published Kaplan–Meier tables

A report can publish the Kaplan–Meier survival probability and the number still at risk without saying how many people died by that time. Multiplying the survival percentage by the randomized cohort size generally will not give the count: censoring changes the denominator at each death. But the published numbers can sometimes **force a single integer count**. This note gives a table-only audit for those cases. It does not use curve pixels, confidence intervals or reconstructed participants.

## Two worked examples

**PERSEUS, first interim cutoff on 1 August 2023.** For the overall-survival (OS) analysis, the published values force these numbers of deaths recorded by that cutoff:

| Elapsed time on the trial’s OS clock | D-VRd (N=355) | VRd (N=354) |
|---|---:|---:|
| Through month 12 | 8 | 17 |
| Through month 24 | 13 | 26 |
| Through month 36 | **20** | **38** |
| All recorded deaths at this cutoff (published) | 34 | 44 |
| Recorded deaths with elapsed time over 36 months (subtraction) | 14 | 6 |

The last row is an accounting split of deaths **recorded by the historical cutoff**. Follow-up lengths and numbers under observation differ between arms. These are not late risks or hazards and do not imply a reversal or change in treatment effect. The [authors’ 2026 ASH abstract](https://onlineprogram.hematology.org/paper-details.aspx?pid=688050) reports a later analysis at 81.3 months of median follow-up in the **same** randomized cohort. The table above addresses only the 2023 report.

Sources: [EMA PERSEUS report, Table 12 (KM levels and totals) and Figure 6 (printed risk counts), pp. 41–42](https://www.ema.europa.eu/en/documents/variation-report/darzalex-h-c-004077-ii-0072-epar-assessment-report-variation_en.pdf#page=41); [Canada’s Drug Agency clinical review, Figure 4 and Table 18, pp. 75–76](https://canjhealthtechnol.ca/index.php/cjht/article/download/PC0388r/PC0388r#page=75) independently reproduces the risk row, confirms the lock and verifies the source count totals. Both reports draw on the same underlying trial, not independent ascertainment.

**DESTINY-Breast09, first interim cutoff on 26 February 2025.** The unchanged rule gives singletons at months 6 and 12 but not at 18 or later. The ranges below are outer limits, not a claim that every integer is jointly feasible.

| Elapsed month | T-DXd + pertuzumab (N=383) | THP (N=387) |
|---|---:|---:|
| 6 | **13** | **5** |
| 12 | **20** | **16** |
| 18 | 31–32 | 31–33 |
| 24 | 41–45 | 42–48 |
| 30 | 45–56 | 48–62 |

Source: [EMA DESTINY-Breast09 report, Table 30 and Figure 18, pp. 72–73](https://www.ema.europa.eu/en/documents/variation-report/enhertu-vr-0000322236-epar-assessment-report-variation_en.pdf#page=72), the same full analysis set and cutoff. This is a second trial used to test the reporting format, not a replication of the first treatment comparison or validation against individual records.

## Proof and assumptions

**1. Bound one KM interval.** Take an ordinary unweighted product-limit estimator, right-continuous, with no delayed entry or re-entry. Suppose risk immediately before the first possible removal in an interval is `R`, risk after the interval is `r > 0`, and `d` counted deaths occur. The exact endpoints of the range with only those counts fixed are:

```
r / (r + d)  ≤  interval survival factor  ≤  (R − d) / R.
```

A tied block of deaths expands into consecutive deaths with exactly the same factor by telescoping. At risk `n`, interchanging a censor followed by a death to a death followed by a censor increases the KM factor from `(n−2)/(n−1)` to `(n−1)/n`; the risk remaining for subsequent steps is unchanged. All censors first minimizes survival and all deaths first maximizes it. This proves the extremes; it does not imply that every value between them is attainable.

**2. Allow either risk-table boundary convention.** A printed count at the interval start landmark `s` is at least the risk after `s`, whether it reports before or after removals tied at `s`. Earlier risk anchors can also upper-bound the start. Use a count at a strictly **later** time than the end landmark `t` as a lower bound for the risk after `t`; that remains safe under either convention because risk only decreases. A larger assumed start risk can only raise the upper factor bound, and a smaller assumed end risk can only lower the lower bound. The first interval always starts with the analyzed cohort `N` *before* any time-zero removals, from survival `S(0−)=1`; it counts deaths at time zero if the source clock includes them. All later intervals are `(s,t]`, with any deaths tied at `t` included. Numeric estimates and risk tables must use the same arm, endpoint, analysis set, cutoff and time scale. The code does not convert published months into an independently assumed day count.

**3. Link rounded estimates.** A percent printed to one decimal has a closed band of ±0.05 percentage points, clipped to [0,1]. Closed bands allow either rounding decision at ties. For neighboring survival bands `[Ls,Us]` and `[Lt,Ut]`, the factor lies within `[Lt/Us, Ut/Ls]` (also bounded above by 1). Keep a possible cumulative death count at `t` only if (a) it passes the baseline-to-`t` bound and (b) its difference from at least one retained predecessor passes the interval bound. Both checks are exact rational intersections. Inductively the true recorded count cannot be removed. The linked code requires a positive later risk anchor at every tested landmark and a strictly positive lower band endpoint for every ratio denominator. It does not extrapolate past exhaustion or accept truncation or unexplained rounding as nearest rounding.

**By hand: why PERSEUS VRd month 36 is 38.** The baseline check permits 37, 38 or 39, and the 24-month count is 26. With 317 at risk at month 24 and 298 at month 39, the rounded OS percentages (92.5, 89.0) require a 24-to-36-month factor between about 0.961102 and 0.963223. For 11 deaths the factor is at least `298/309 ≈ 0.964401`; for 13 it is at most `304/317 ≈ 0.958991`. Both fail. Only 12 more deaths survives the necessary tests, forcing 38 at month 36.

The procedure is a finite **outer sieve**. A singleton certifies the only count if the source satisfies the inputs; a non-singleton gives coverage, not a sharp identified set. Zero candidates suggests an inconsistent input or assumption, not a corrected value. There is no independent-censoring assumption for the bookkeeping result. Independent censoring can still matter for interpreting a survival estimate as a population risk. Never apply these formulas unchanged to weighted, delayed-entry or covariate-adjusted estimators, cumulative-incidence curves, unknown or misaligned time conventions, or unobserved deaths.

## Closest prior work

[Fojo and Labadie (2024)](https://doi.org/10.1093/oncolo/oyae070) already used the all-censors-before and all-censors-after extremes with published KM and risk counts, and [their follow-up](https://doi.org/10.1093/oncolo/oyaf003) applies interval accounting to investigate censoring in a trial. The core interval bound is not new. The practical contribution here is an exact-rational implementation that chains rounded table entries, allows for unknown risk-count conventions at ties, and reports an exact count only for a singleton. It does not require a plotted curve trace. [Guyot et al. (2012)](https://doi.org/10.1186/1471-2288-12-9), [Titman (2026)](https://doi.org/10.1002/sim.70474) and [kmextract exact mode](https://github.com/open-computational-biology/kmextract#exact-mode-km_extract_exact) use digitized or measured curve information for reconstruction or exact enumeration conditional on figure-derived inputs. This note makes no claim to be the first inverse-KM method.

## Replay

```sh
python3 -B test_replay.py
python3 -B replay.py > /tmp/os-count-replay.txt
diff -u replay.txt /tmp/os-count-replay.txt
sha256sum -c SHA256SUMS
```

Python standard library only. The synthetic tests exhaust 1,013 death/censor orders, 1,520 grouped schedules with day-zero and landmark ties (97,280 checks crossing risk conventions and rounding), and 22,200 exact rounding-edge checks. No individual records, reconstructed participant histories, or actual event dates are used; the demonstrations use only published aggregate counts and percentages. The trials are demonstrations of the accounting check, not treatment recommendations.

### Appendix: a different counting shortcut that fails

In the same DESTINY-Breast09 analysis, subtracting the deaths counted as the first centrally assessed PFS event from total OS deaths does **not** give deaths following documented central progression. The regulator’s [Table 28, pp. 69–70](https://www.ema.europa.eu/en/documents/variation-report/enhertu-vr-0000322236-epar-assessment-report-variation_en.pdf#page=69) and [Table 30, p. 72](https://www.ema.europa.eu/en/documents/variation-report/enhertu-vr-0000322236-epar-assessment-report-variation_en.pdf#page=72) share the same groups and cutoff. There are 58 and 68 total OS deaths, of which 18 and 5 are counted as deaths without progression in the primary central PFS analysis. The PFS table explicitly lists another 11 and 17 deaths with no documented central progression, censored for PFS after at least two missed visits. Thus at least 11 of the 40 and 17 of the 63 in the subtraction remainders lack documented earlier central progression. Disease development during missed assessments cannot be decided from these counts.
