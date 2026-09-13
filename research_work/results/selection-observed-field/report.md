# Observed membership of the pilot field

This metadata audit advances the Milky Way selection model; it is not a gravitational measurement. The previous turn's priority discussion did not change research state. This run recovers new evidence from the authoritative catalogs.

## Counts and integrity

For field 300+00, telescope lco25m, allStarLite contains 1,343 rows: 1,271 with EXTRATARG=0, 57 with EXTRATARG=1, and 15 with EXTRATARG=5. The 1,271 main-sample identifiers are unique and all match the downloaded parent catalog. J, H, K, AK_TARG and the extinction-method label agree exactly for every matched main target. All 521 previously selected main training stars are contained in this larger set.

The [SDSS selection documentation](https://www.sdss4.org/dr17/irspec/targets/selection-biases/) identifies EXTRATARG=0 as the convenient main-sample flag and requires selection fractions by field, magnitude cohort and color. Main membership is not proof that every object lacks overlapping special-target criteria.

| Recorded H bounds | Recorded color bounds | Main observed count |
|---|---|---:|
| 10–12.2 | 0.5–0.8 | 376 |
| 10–12.2 | Above 0.8 | 362 |
| 12.2–12.8 | 0.5–0.8 | 267 |
| 12.2–12.8 | Above 0.8 | 264 |
| 12.8–13.3 | 0.5–0.8 | 1 |
| 12.8–13.3 | Above 0.8 | 1 |

These groups use stored targeting bounds, not newly reconstructed dereddened colors. Exact floating-point bounds and the open-upper-bound sentinel appear in results.json. All H measurements are inside their recorded inclusive bounds; this is not a verification of targeting endpoint conventions.

## Exception requiring reconciliation

The eight acquired designs record zero long-cohort fraction, yet two observed main targets have APOGEE2_LONG and APOGEE2_NORMAL_SAMPLE flags: 2M12262737-6244152 and 2M12272092-6306120. Their metadata agree with the parent catalog. Do not drop them or assign a zero probability automatically. Actual plate/design histories and the interpretation of planned cohort fractions must be reconciled first.

## What this permits next

The observed counts can now be compared with eligible parent counts after verifying footprint, photometric quality, dereddening, bin boundaries and cohort-version overlap. Neither 521/189661 nor 1271/189661 is the required completeness estimate. Project-specific distance, chemistry, quality and matching cuts then need their own treatment before spatial stellar densities can constrain gravity.

Only identifiers and targeting/photometric metadata were accessed, including metadata outside the training subset. No held-out velocities were read. No physical parameters, sample membership or holdout assignments changed. All six scientific objectives remain open.

Run `python research_work/results/selection-observed-field/run.py` from the repository root. Input hashes, counts and agreement checks are in results.json; the local metadata parquet is ignored by Git. The counting and crossmatching are standard survey bookkeeping, not a new physical formula.
