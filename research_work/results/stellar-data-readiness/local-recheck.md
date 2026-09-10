# Catalog availability recheck — 10 September 2026

All four requested observational products are already on this machine. No repeat network download is needed. `recheck_local.py` reread catalog headers, verified APOGEE and StarHorse against their original downloaded SHA-256 hashes, verified the recovered BRAVA table, and checked prepared-sample identifiers and regional counts. `local-recheck.json` contains exact paths and hashes.

| Product | Verified rows | Use |
|---|---:|---|
| APOGEE DR17 allStarLite synspec_rev1 | 733,901 | Spectroscopic catalog rows, including repeats; radial velocities and chemistry |
| StarHorse APOGEE v2 | 562,424 | Distance estimates and quantiles |
| Gaia DR3 covariance/quality subset | 140,407 | Astrometry and correlated measurement errors for the selected sample |
| BRAVA | 8,585 | Separate bulge line-of-sight velocity check |

BRAVA was also copied, without modification, from the recovered August project into this repository's ignored catalog cache. Large catalogs remain outside Git; scripts and integrity records are tracked. The prepared Gaia-enriched parquet hash remains unchanged. The Gaia subset is not the full Gaia catalog.

The prepared sample contains 140,407 unique integer source IDs and 128,772 candidates after additional Gaia quality cuts. These are selected observations, not a complete census.

| Comparison region | Radius from Galactic axis | Absolute height | Quality candidates |
|---|---|---|---:|
| Plane beneath bulge | 0.5–3.5 kpc | below 0.2 kpc | 717 |
| Above/below bulge | 0.5–3.5 kpc | 0.5–1.5 kpc | 4,085 |
| Disk plane control | 5–9 kpc | below 0.2 kpc | 7,703 |
| Disk off-plane control | 5–9 kpc | 0.5–1.5 kpc | 27,620 |

These counts establish coverage only. They use median-distance positions and combine hemispheres for this inventory; the actual test must keep the two hemispheres distinguishable and propagate distance uncertainty.

## Model and analysis readiness

The local ordinary-matter field implementation already includes a Sormani bar, stellar and gas disks, nuclear components, and a central mass approximation. See `../bar-field-foundation/report.md` and `../rotating-bar-orbits/report.md`. It separates a standard halo benchmark from the companion hypothesis. Published component normalizations retain prior dynamical-model assumptions and need uncertainty treatment; their existence does not make the baseline independently measured or final.

The remaining work is analytical rather than missing bulk catalogs:

1. Resolve or explicitly exclude suspect cross-identifications under a declared rule. An incorrect match combines the distance of one star with another star's motion and can invent an extreme orbit.
2. Use StarHorse input flags and a consistent distance/astrometry likelihood. Gaia parallax must not be counted twice when it already informed the distance estimate.
3. Match chemistry, bar position and survey selection across regions. Otherwise different stellar populations can look like different gravity.
4. Fit orbital populations in the rotating three-dimensional field. A gravitational well predicts acceleration; a population of orbits predicts the distribution of measured velocities. One circular-speed number per bulge star is inadequate.
5. Specify a shared deposit distribution and gravitational response, then compare the same measured quantities in all models. Independently tuning each region would remove the predictive test.
6. Evaluate a declared independent set after fitting. Existing exploratory regions and previously exposed scores cannot be relabeled as fresh blind tests.

No new velocity fit or holdout score was calculated in this availability recheck. The repository's existing lens and disk results do not yet establish a common photon-to-companion mechanism.

## Official provenance

- APOGEE: https://www.sdss4.org/dr17/irspec/spectro_data/
- StarHorse v2: https://data.aip.de/projects/aqueiroz2023.html — the release documents its revised age/mass prior and unchanged distance estimates.
- BRAVA: https://irsa.ipac.caltech.edu/data/BRAVA/index.html
- Gaia: https://www.cosmos.esa.int/web/gaia/dr3

This report introduces no new physics formula. Catalog joins, checksums and coverage counts are standard data-processing operations.
