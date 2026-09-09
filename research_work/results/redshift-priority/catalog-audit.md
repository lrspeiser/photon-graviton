# Candidate catalog and environmental-input audit

ELVES-Field candidate staging is now available: [source audit](elves-field-audit-report.md) and feature-only JSON retain all 95 publisher rows and 29 provisional candidates. Redshift-derived distances and exposed examples are excluded. Candidate aliases/groups, frame conventions, selection and historical freshness remain unresolved; no labels were scored and no environmental proxy adopted.

New source decision: the [2025 TRGB-SBF identity audit](sbf-2025-audit-report.md) resolves 14 calibration targets plus the anchor and Coma object. All 16 share already excluded CF4 groups. The publication recalibrates older distant SBF observations and does not supply a newly certified group-independent test here. Revised calibration remains a sensitivity option; current adopted distances are unchanged.

Update: the [structured identity audit](object-exposure-report.md) now resolves all six maser PGC aliases and finds three exact earlier pilot overlaps. The conservative group exclusion set is 235. The previous name-only scan and 232-group list below are historical minimum snapshots; they no longer represent the latest audit. General aliases, supernova hosts and unstructured historical exposure remain unresolved, so no fresh sample is certified.

9 September 2026. Metadata audit only: no new target-level redshift-distance pairs scored, no new rate fitted, and no fresh-validation success claimed. Previously used full CF4 table is already on disk. It is not automatically an untouched sample merely because some rows were not used in the 164-group fit.

## Candidate decisions

| Source | Useful input | Dependence and decision |
|---|---|---|
| Cosmicflows-4 method-specific distance columns | SBF, Cepheid/TRGB and other published indicators with identifiers | Retain stipulated distances, but audit calibration, group construction and historical exposure. Do not substitute the velocity-to-distance calculator. The existing SBF subset is exposed. |
| Megamaser Cosmology Project geometric distances | Disk-based geometric distances, independent of a target Hubble-law inversion | Candidate external check, pending identifier/group overlap audit and exact measurement conventions. Do not adopt its fitted Hubble constant or flow corrections. Six-object 2020 compilation alone is too small to establish all-distance accuracy. |
| Higher-distance standard-candle compilations | Potential wider-range published distances | Pending individual release audit: distinguish measured distance modulus from redshift-derived coordinate distance, check calibration, selection and shared objects. Not yet eligible for fresh validation. |
| 2M++ reconstructed density/velocity products | Potential environmental information | Do not use as theory-independent void or motion measurements. The audited reconstruction calibration uses cosmological simulations and a growth/velocity interpretation. An observational tracer catalog would need a separate coordinate and selection audit. |

## Primary sources and limits

[Cosmicflows-4](https://arxiv.org/abs/2209.11238) compiles multiple indicators and groups; the local CDS ReadMe and recovered protocol identify exactly which SBF columns and representative selection were used here. A published distance remains stipulated under the user's contract, while calibration assumptions remain visible. Independently predictive redshift testing excludes distances computed from the target redshift itself even if those distances remain allowed fictional facts.

[Pesce et al. 2020](https://arxiv.org/abs/2001.09213) describes six geometric maser distances. Its default cosmological analysis assumes 250 km/s peculiar-motion uncertainty and also examines alternative corrections. This is an assumed nuisance scale, not an independent measurement of each galaxy's motion or a ready-made precision floor for our sample. The abstract's aggregate results were viewed during this audit; a subsequent test of these six galaxies must disclose that source-level exposure even if individual pairs were sealed.

[Hollinger and Hudson 2024](https://academic.oup.com/mnras/article/531/1/788/7651284) calibrates 2M++ velocity-density comparisons with simulation-based mock catalogs. Its quoted velocity errors and cosmological residual-flow conclusions are conditional on that reconstruction framework. Importing them as independent corrections would violate this model's present comparison contract.

## Exposure and quarantine status

The accompanying known-exposure-exclusions.json contains the union of 164 current groups and 68 earlier pilot groups: 232 unique group identifiers. This is a minimum verified exclusion set, not certification that every other object is untouched. Other recovered galaxy experiments exist, and alternate identifiers/group membership require crossmatching. Numeric substring searches are insufficient to certify object identity.

No fresh sample is yet declared quarantined. Before reading/scoring any proposed evaluation pairs: finish the exposure registry; freeze source version/hash, method and quality criteria, identifier matching and group exclusions; separate feature and label outputs; seal predicted values and parameters before label access; then evaluate once. If audit reveals prior outcome exposure, label the sample external/reused instead of retroactively blind.

## Observation accounting

Retain the catalog CMB-frame redshift convention and audit its original optical/radio/relativistic velocity definition before mixing sources. Never assume that every catalog velocity column is c*z. A motion correction derived as observed redshift minus an expansion prediction cannot be used to prove an alternative redshift law. Shared clock/endpoint factors, motion and instrumental errors need a consistent multiplicative redshift-factor model, with no per-object tuning. We have not yet established credible predictive intervals.

## Next action

Build an exact identifier exposure registry from the other recovered observational experiments and audit the maser measurement table and frame definitions under a predeclared external-test protocol. Continue the search for independently calibrated wider-distance samples and independent path tracers. The frozen smooth-rate comparison remains unsuccessful; no additional flexible environmental fit is authorized by missing inputs.
