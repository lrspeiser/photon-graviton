# ELVES-Field: 29 staged candidates, no new predictions scored

We retrieved the [publisher's confirmed dwarf table](https://content.cld.iop.org/journals/0004-637X/1001/2/244/revision1/apjae4c5ct1_mrt.txt) associated with [Carlsten et al. 2026](https://arxiv.org/abs/2602.16766). The source is useful for investigating nearby distances beyond the lower boundary of the current 10.2–93.2 Mpc sample, not for establishing distant coverage. Its distance column mixes TRGB, SBF and redshift-derived methods. The latter uses a flow-based redshift-distance conversion and is ineligible for independent validation here. Published isolation flags incorporate group/virial and distance assumptions; they are not an independently measured path through voids.

The downloaded machine-readable table has 95 unique identifiers. The extractor reads only names, sky positions, method codes, SBF quality flags and velocity availability, preserving every row and exclusion reason. It does not parse numeric distances or velocities, rank by model agreement, apply isolation cuts or assign independent motions.

| Distance method | All rows | Listed velocity present |
|---|---:|---:|
| TRGB | 22 | 17 |
| SBF | 51 | 19 |
| Redshift-derived | 22 | 22 |

The predeclared criteria retain 29 provisional candidates after restricting to TRGB/SBF with a velocity, excluding ambiguous/failed SBF measurements and removing already exposed appendix examples. Eight example rows appeared automatically in a web lookup before feature staging; their names are explicitly recorded and excluded even where they would fail another rule. No remaining candidate is certified fresh merely because the table is new or its name is unfamiliar.

## Preserved quarantine boundary

- Full raw table remains in the projectless work/catalog-audit directory. Repository output contains allowed metadata only, source/protocol hashes and reasons for exclusion.
- No candidate redshift is predicted, scored or used to choose a model. Distance values remain unmodified; none is substituted into the existing 164-row dataset.
- Reading a source into parser memory and keeping it on disk is disclosed. The claim is that non-example target values were not displayed to this research conversation or scored in this audit, not that the source was never accessed or that all targets are historically unexposed.
- Next: match all 29 candidate names and coordinates against earlier galaxy identities and group/host memberships; inspect velocity convention/frame and distance calibration/selection. Spectroscopic availability and the survey's original distance-confirmation procedure can bias this subset and need explicit treatment.
- Freeze a common formula, fitted parameters, valid range and independently justified uncertainty model before opening labels for any eventual evaluation. Do not tune 29 object-specific gearing ratios or replace missing environmental paths with the survey's isolation flags.

Run elves_field_features.py with the downloaded table path. The script pins the publisher file hash and asserts 95 distinct rows and retention of all disclosed exposed examples. A separate verification replaced every distance value and every nonblank velocity value with nonnumeric placeholders: the exported identity/eligibility records remained identical. This confirms that target magnitudes do not determine the staging decision; it does not establish historical freshness or validate the physical model.
