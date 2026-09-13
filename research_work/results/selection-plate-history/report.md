# Pilot-field drilled targets and cohort history

The previous goal turn made progress by recovering the full observed main sample. This follow-up resolves whether its two faint targets were actually part of plate designs. It does not yet establish their inclusion probabilities.

## Direct evidence

Downloaded the eight plateHolesSorted files from the public SDSS platelist repository for field 300+00 and verified each plate/design identifier against the release plate summary. Each contains 250 assigned science fibers: 2,000 assignments in total, representing 1,439 unique science targets. Repeated assignments are not independent stars.

The documented APOGEE2_TARGET1 normal-sample bit identifies 1,404 unique drilled targets. All 1,271 observed EXTRATARG=0 targets occur in this set. Of the remaining 133 drilled normal targets, 110 are absent from this field/telescope's allStarLite rows and 23 have EXTRATARG=1. Absence does not identify a reason: visit completion, reduction acceptance and other metadata require further investigation. We do not relabel the 23 as main targets merely because a design flag says normal.

| Plate | Design | Normal short | Normal medium | Normal long |
|---|---|---:|---:|---:|
| 10176 | 10719 | 100 | 140 | 2 |
| 10177 | 10720 | 103 | 140 | 2 |
| 10178 | 10721 | 108 | 134 | 2 |
| 10179 | 10722 | 107 | 134 | 2 |
| 10180 | 10723 | 107 | 135 | 2 |
| 10181 | 10724 | 109 | 136 | 2 |
| 10182 | 10725 | 109 | 138 | 2 |
| 10415 | 10999 | 110 | 138 | 2 |

The same two long-cohort stars, 2M12262737-6244152 and 2M12272092-6306120, are assigned on every plate. Both retain normal-sample and long-cohort flags. Their 16 assignments have assigned=1 and conflicted=0. They are genuine repeated targeting entries, not unsupported allStar cohort labels.

The separately downloaded apogee2Design file also records [0.5, 0.5, 0] planned cohort fractions for all eight designs, agreeing with allDesign. The drilled allocations demonstrably differ from these nominal fractions, including the long cohort. Why the design algorithm retained these two stars is not derived here. Zero nominal allocation cannot be used to infer zero realized sampling.

## Consequence for the gravity test

Use unique observed main counts and a correctly eligible photometric parent for selection fractions. Do not multiply by eight for repeated long-cohort assignments, use nominal design fractions as completeness, or equate drilled targets with accepted observed targets. Recovering those distinctions prevents survey coverage from being mistaken for a spatial stellar-density signature of companion deposits.

Next reconstruct the eligible photometric parent by footprint, quality, dereddened color and H bounds, with explicit handling of cohort versions. Then model the project's additional sample cuts and the distance-dependent tracer population. These tasks remain necessary before the vertical gravity comparison has an observational likelihood.

## Provenance and verification

- [SDSS platelist description](https://www.sdss.org/dr18/software/packages/svn/) and [path definitions](https://sdss-access.readthedocs.io/en/latest/path_defs.html).
- [Official bit definitions](https://www.sdss4.org/dr17/algorithms/bitmasks/#APOGEE2_TARGET1): short 11, medium 12, long 13, normal sample 14. Signed flag storage is handled with integer bit operations.
- Every downloaded plate URL/hash, plate/design count, exception assignment and missing-main classification is recorded in results.json. Input plate-summary and observed-metadata hashes are included.
- The limited text reader checks every STRUCT1 row against its declared field count, retaining only targeting columns. Identifier-format, uniqueness and membership assertions pass. Run `python research_work/results/selection-plate-history/run.py` to reproduce with cached inputs.
- Two existing FITS header-format warnings (SURVEY and DATE) are retained; they do not affect the columns used.

This is standard survey metadata analysis, not a new physical formula. No held-out kinematics were read, no gravity parameters were changed, and no causal or observational success of the companion hypothesis is claimed. All six goals remain open.
