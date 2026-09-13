# Historical science-input audit

The previous turn made progress by reconstructing photometric eligibility and identifying a failed infrared-error cut. This follow-up traces those exceptions to the original ranked science inputs referenced by the eight pilot plate files.

## Results

The eight input files contain 2,329 rows, representing 1,745 distinct targets; 1,710 distinct targets have the normal-sample bit. Every one of the 2,000 drilled science assignments matches its own design's input list.

All 42 observed stars rejected by the proposed 0.1-mag infrared-error filter occur in the input lists with the normal-sample bit. This rules out an explanation in which normal-sample classification was introduced only in the final allStar catalog. It does not establish which infrared-error values the historical selection algorithm used, nor why the published limit was not reflected in the retained parent metadata. These input files do not carry the mid-infrared error needed to resolve that question.

| Design | Ranked input rows |
|---|---:|
| 10719 | 294 |
| 10720 | 291 |
| 10721 | 292 |
| 10722 | 291 |
| 10723 | 289 |
| 10724 | 289 |
| 10725 | 291 |
| 10999 | 292 |

One design/target pair appears twice in an input list. The first implementation's uniqueness assertion caught this; matching now preserves multiple input entries instead of dropping one. All drilled assignments still match. The duplicate key is recorded in results.json. Repeat entries and assignments must not be counted as additional unique stars.

## Consequences for the stellar-density test

These ranked lists have only about 290 entries per design and are already selected from a larger photometric pool. They cannot replace the eligibility denominator. The evidence now connects the anomalous 42 stars across the photometric parent, historical normal-target inputs, drilled plate metadata, and observed catalog. It does not justify arbitrarily removing the infrared cut from every unseen target or asserting that either previous candidate denominator is exact.

The next useful step is to examine the targeting implementation or quantify a clearly declared selection sensitivity, while retaining camera/footprint and subsequent project-selection effects. No inverse-probability weights are produced from this audit. No new redshift, deposition, motion, or lensing success is claimed; all six goals remain open.

## Reproduction

Run `python research_work/results/selection-science-inputs/run.py`. URLs are taken from each previously verified plate's SCI input reference, resolved under the public SDSS platelist inputs directory. Input hashes, design counts, duplicate keys and all 42 target histories are in results.json. Cached full inputs and extracted metadata remain outside Git. Every row is checked against its declared field count, and matching is performed on design plus normalized identifier. No held-out kinematic columns are read.

This is survey provenance analysis, not a new physics equation or evidence of a companion interaction.
