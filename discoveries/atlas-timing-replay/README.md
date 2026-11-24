# Timing and period-fit checks for 3I/ATLAS

**Reproducibility note | 24 November 2026**  
A reanalysis of published numerical inputs and Hubble exposure metadata, not a new measurement of the comet’s rotation. The scope is [arXiv:2601.10860v1](https://arxiv.org/abs/2601.10860v1) (15 January 2026) and [arXiv:2602.18512v2](https://arxiv.org/abs/2602.18512v2) (24 June 2026).

## Hubble exposure time

The first jet-angle entry in the January preprint’s Table 1 gives a fractional date corresponding to **19:25:00.192 UTC on 30 November 2025**. The later preprint supplies the exposure ID, **IFLE01FBQ**, and lists **20:36:57 UTC** in its clock-time column and its separate exposure table. All 21 exposure start times in that separate table agree with the public Hubble FITS headers. The first exposure started at 20:36:57 UTC, about **72 minutes** after the time implied by the fractional date. [Hubble program 18152’s visit log](https://www.stsci.edu/hst-program-info/visits/?program=18152) and the independent exposure summary in [Hui et al., Table 1](https://doi.org/10.3847/2041-8213/ae471c) are consistent with the native time.

The later preprint's elapsed-time column is not a separate check: its first interval is 1.317 hours, versus about 0.121 hours between the archived exposure midpoints. In a 6.8–7.4-hour comparison window, a two-harmonic fit to the printed fractional dates favors 7.200 hours over 7.374 hours (χ² 2.938 versus 8.965). With archived midpoints and the same angles, **7.372 and 7.200 hours fit almost equally well** (χ² 2.783 versus 2.990). Omitting the first angle also yields a near tie. This identifies a sensitivity of the published table replay, not the time convention used in the authors’ original code. **7.20 hours is not ruled out.**

## Small-aperture light curve

We rebuilt the January preprint’s four-night, 14-bin table from 286 frame measurements in the June preprint’s Table 9. Per-frame zeropoint adjustment, weighting by the stated SNR-derived errors, and 30-minute bins starting at the first frame of each visit reproduce all 14 binned times, magnitudes, and formal uncertainties to the printed precision.

The January paper reports **7.136 ± 0.001 hours** and describes both a one-harmonic and a two-harmonic fit. In a weighted least-squares replay with a free mean for each night and the stated 0.089-mag jitter added in quadrature:

| Model | χ² at 7.136 h | Better fit within 6.8–7.4 h |
| --- | ---: | ---: |
| One harmonic | 51.322 | 6.89294 h, χ² 16.301 |
| Two harmonics | 24.506 | 6.88416 h, χ² 12.544 |

Using only formal errors does not make 7.136 hours the best fit either, although the two-harmonic fit has a higher-χ² local minimum at 7.13764 hours. These comparisons do not establish any of the other trial periods as physical signals.

The June paper also supplies a **different five-night dataset** (Table 4, groups of 30 consecutive frames rather than 30-minute bins). With one harmonic and its tabulated formal standard errors, it does produce a **secondary** minimum at **7.13597 hours**, which rounds to 7.136, but 7.37040 hours fits better in the same window (χ² 42.956 versus 62.346). Adding the stated jitter moves the secondary minimum near 7.136 to 7.13324 hours. This nearby minimum comes from different inputs and does not reproduce the optimum or precision for the January four-night dataset. As a separate check on the replay method, applying the stated correction to the June revision’s larger aperture reproduces its 7.1423- and 7.1583-hour local fits for the four- and five-night frame subsets, respectively.

## Scope and relationship to other work

This AI-assisted note has not been peer reviewed and makes no priority claim. This is a check of two specified preprint versions, not a correction of an independently measured rotation, a claim about the true period, or a statement about any undocumented analysis. A later [dedicated TESS study](https://doi.org/10.3847/PSJ/aeac2c) found no persistent brightness signal at its tested periods in January, a different observing interval. [Ahuja et al.](https://doi.org/10.1093/mnras/stag1800) found multiple aliases in a separate set of late-December observations. Neither result is presented here as a replication of the timing and input checks above.

## Replay from the original sources

With Python 3.10+, run `python -m pip install -r requirements.txt` and then `python replay.py`. It reads the two versioned arXiv HTML pages and only the first 14.4 KB of each of 21 public Hubble FITS files. No account or data files are required. It uses the original source tables in memory and prints only derived diagnostics, comparisons, and input version hashes as JSON. Exact replay does not depend on third-party data packaged here.

No third-party HTML, tables, images, FITS headers or datasets are distributed with this note. For the source papers’ reuse terms see the [January preprint license](https://creativecommons.org/licenses/by-nc-nd/4.0/) and [June preprint license](https://creativecommons.org/licenses/by-nc-sa/4.0/). The [public Hubble collection on AWS](https://registry.opendata.aws/mast-hst/) lists its data as public domain.

Based on observations made with the NASA/ESA Hubble Space Telescope, obtained from the Mikulski Archive for Space Telescopes at the Space Telescope Science Institute, which is operated by the Association of Universities for Research in Astronomy, Inc., under NASA contract NAS 5-26555. These observations are associated with program 18152.
