# First object-level SBF source mapping

The author's digital files provide 280 records in table.good and 19 in table.poor, with 299 unique source names. A fixed 30-arcsecond match to all 164 exposed representative positions finds **110 single positional counterparts and 54 unmatched rows**, with no multiple matches. Of the matches, 103 are in table.good and seven in table.poor. These classifications belong to those original files; they do not determine the quality of CF4's combined distance or justify deleting a current row.

The files are linked in [Tonry et al. 2001](https://arxiv.org/abs/astro-ph/0011223), page 7, and were recovered from the author's [table.good](https://www.ifa.hawaii.edu/~jt/SBF/table.good) and [table.poor](https://www.ifa.hawaii.edu/~jt/SBF/table.poor). We preserve their actual downloaded counts and hashes rather than force them to the paper's advertised 300-galaxy total. Exact release equivalence and the missing-record discrepancy remain unresolved. The EDD endpoint failed with an expired TLS certificate, and the attempted CDS catalog location returned 404; neither failure is evidence that a particular galaxy lacks provenance.

## Source membership and the existing residual pattern

| Positional source match | Number | Existing RMS (km/s) | Existing mean residual (km/s) | Adopted distance range (Mpc) |
|---|---:|---:|---:|---:|
| Tonry candidate | 110 | 429.766 | -43.835 | 10.205-61.461 |
| No match within 30 arcseconds | 54 | 490.633 | +25.839 | 16.398-93.197 |

Residuals are the saved coarse-sky out-of-fold constant-baseline prediction minus observation, expressed as c times the redshift difference. They are not newly fitted velocities. All targets remain in the combined 164-row analysis.

Fifteen of the 16 targets in region 2, the region with the largest negative mean residual, have Tonry counterparts. But every one of the 17 region-3 targets also has a Tonry counterpart, and region 3's mean residual was only +56 km/s. Thus membership alone does not explain the difference. The two source-match categories have different distance distributions and strongly different sky coverage; their RMS difference is not a controlled survey-quality comparison or a measured calibration offset.

This is a useful measurement-provenance lead, not an independent void measurement. No source label was assigned by residual, and no source-specific rate, velocity or distance correction was fitted.

## What the match does and does not prove

A close coordinate match supports likely source membership, but it does not certify a PGC alias or the actual weight given to that source in CF4. CF4 may combine an original Tonry measurement with later data and recalibrate it. A poor-class historical entry may have a better later measurement; an unmatched row may appear under a shifted centroid or in another source. The script records all names, original source positions, separations and classifications so these links can be checked rather than assumed.

Before interpreting a dust revision numerically, resolve those identity links and source contributions and recover the applicable passband/extinction calibration. The 54 unmatched entries need the other optical/infrared catalogs. This mapping does not revise the user's adopted distances or undo existing exclusions in the fresh-catalog audit.

## Reproduction

Run `python research_work/results/redshift-priority/tonry_membership.py /path/to/cache`, with the two downloaded files named tonry-table.good and tonry-table.poor. The parser verifies source hashes and reads only source name, RA and Dec. Masking all remaining source fields leaves the extracted records unchanged. All 164 angular joins were independently checked using Astropy, including the no-match cases. The recorded protocol discloses the initial format-inspection exposure of seven source velocity prefixes; those source entries must not be presented as newly blind observations.

The saved JSON contains all 299 source metadata records, all 164 mapping decisions and regional counts. **Formula provenance:** matching uses established spherical angular separation and the summaries use standard descriptive statistics; no new physical formula is introduced. The full goal remains active, with no fresh predictive improvement claimed.
