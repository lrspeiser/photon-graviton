# ELVES velocity provenance: current metadata and original measurements

Audit date: 9 September 2026. This advances the independent-catalog and observation-model requirements. No candidate distance or velocity values were queried by the metadata script, transformed, predicted or scored. The original publisher table remains the adopted source; a current database record is not proof of its historical provenance.

## What changed

The fixed 30-arcsecond SIMBAD queries succeeded for all 26 pending positions, with no truncated cones. Every position has at least one possible counterpart; five cones contain multiple objects. Name/previously confirmed PGC evidence supports exactly one counterpart for 15 targets. Eleven require further identity work. The saved responses contain 99 measurement-metadata rows across all counterpart records. These are bibliographic/type records, not 99 independent galaxy measurements or 99 certified targets.

Seven of the 15 supported counterparts currently prefer Yu et al. 2022; seven prefer other references, and one has no current basic velocity reference. The latter, dw0930p2143, has a velocity-present flag in the staged publisher features. That discrepancy demonstrates why a present-day preferred record cannot simply stand in for the authors' retrieval. It is not evidence that the published velocity is wrong.

Identity support is narrower than freshness. All 26 remain pending, and zero adopted per-row measurements have been traced conclusively to a particular original source value/conversion. The existing three historical exclusions remain unchanged. We have not used whether a galaxy's redshift would help the model to resolve any identity.

## Source-level findings

| Reference | Supported measurement meaning | What remains to establish for ELVES |
|---|---|---|
| Haynes et al. 2018, 2018ApJ...861...49H | Heliocentric optical-convention velocity, based on the midpoint of the two half-peak edges of the HI profile; the paper specifies a statistical center error of half the tabulated width error. | Whether this particular measurement was used, and how any database conversion was applied. Measurement error does not include galaxy motion or all profile-systematic effects. |
| Yu et al. 2022, 2022ApJS..261...21Y | Catalog central velocity is weighted by HI flux intensity. Its uncertainty is separately tabulated. | Exact velocity-axis convention/frame in this reprocessing and the path to the ELVES adopted value. The source's use of ALFALFA spectra alone is not sufficient evidence of every conversion. |
| SDSS spectral references | SDSS describes heliocentric vacuum-wavelength spectra and spectroscopic redshift extraction. | Verify the applicable release/product and any database transformation for each adopted entry. |
| Other references in the saved summary | Bibliographic identity only in this pass. | Read the source definition and distinguish original measurement from compiled replacement. |

The [Haynes source paper, catalog column 5](https://wiki.physics.wisc.edu/ObsCos/images/2/24/Haynes_ALFALFA_Source_Catalog_ApJ_2018.pdf) explicitly specifies the optical convention despite the radio observing band. Its [CDS catalog documentation](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/ApJ/861/49?format=html&tex=true) also distinguishes high-quality detections from lower-significance detections associated with an already known optical redshift. This selection dependence must be retained, not treated as random redshift availability.

The [Yu catalog documentation](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/ApJS/261/21?format=html&tex=true) identifies the flux-weighted center, carries profile-confusion flags and imports ALFALFA distances. We do not import those distances, dynamical masses or derived environment quantities. A profile midpoint and a flux-weighted center need not coincide for an asymmetric profile; this is a possible measurement-definition contribution, not a quantified explanation of the model residuals.

The [SDSS pipeline documentation](https://classic.sdss.org/dr7/products/general/edr_html/node59.php) provides source-level wavelength/frame context. It does not certify every historical SIMBAD entry or the specific ELVES adopted source.

**Formula provenance: established optical velocity convention, not new research.** For a velocity explicitly tabulated in that convention, z = V_opt/c. A radio observing band does not imply the radio velocity convention. SIMBAD's wavelength code R describes the observing domain, while a generic v type does not by itself establish a conversion formula. Apply the declared multiplicative frame transformation only after the source convention is known; do not correct the already-CMB CF4 data a second time.

## Reproduction and evidence limits

Run `python research_work/results/redshift-priority/summarize_elves_velocity_metadata.py` to replay the saved audit offline. It checks the exact allowed field sets, target completeness, unique target names, absence of retrieval errors/truncation, and association of measurement records to returned object identifiers. It records input hashes and the specific aliases supporting each match. It fails rather than silently omitting an incomplete query. The saved summary contains all 26 per-target decisions and references.

`elves_velocity_metadata.py` is the separate network-acquisition script. Re-running it refreshes a mutable external database and may legitimately change results; it is not required to replay this snapshot. Query strings and raw-response hashes are preserved, but raw HTTP response bytes are not retained, so their hashes are audit identifiers rather than independently reconstructible response archives. The saved parsed payload itself is hashed by the summary.

Free-text aliases/remarks were returned as metadata; the replay check proves that numeric outcome columns were not returned, not that arbitrary remote free text could never mention an outcome. No target distance/redshift values were identified or used in this audit. Source documentation includes aggregate ranges, which are source-level exposure and are not a blind-data success claim.

## Next usable decision

Finish the measurement chain for supported counterparts and resolve the eleven remaining identities, including the two possible old-group aliases. Use primary identifiers and source definitions, never favorable residuals. Record whether the published value's provenance is verified or remains uncertain. If the historical source cannot be recovered, retain that limitation explicitly; do not manufacture a match by choosing the current measurement closest to a desired prediction.

Before opening outcomes for evaluation, freeze one common model, domain, group exclusions, exact measurement/conversion rules and independently justified uncertainties. Current metadata neither supplies independent galaxy motions nor resolves the lack of a predictive environmental input. The smooth distance-dependent extension has not improved exploratory validation; no model is declared successful by this audit.
