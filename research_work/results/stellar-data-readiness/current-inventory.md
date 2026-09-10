# Local catalogs and remaining bulge/disk work

All four requested catalog products are present. A fresh execution of `recheck_inventory.py` verified their row counts and SHA-256 hashes against the earlier inventory, and verified the prepared matched table. No duplicate download is needed. This is a data-readiness audit, not a new observational fit.

Local directory: `C:/Users/henry/Documents/Codex/photon-graviton/research_work/data-cache/stellar-catalogs/`.

| Product | Local file | Verified rows | Role |
|---|---|---:|---|
| APOGEE DR17 | allStarLite-dr17-synspec_rev1.fits | 733,901 | Radial velocities, stellar chemistry and quality flags; rows are not all distinct stars |
| StarHorse v2 | APOGEE_DR17_EDR3_StarHorse_v2.fits | 562,424 | Distance posterior summaries and input flags |
| Gaia DR3 targeted subset | gaia-dr3-quality-covariance.fits | 140,407 | Astrometry, quality information and error correlations; not the full Gaia archive |
| BRAVA | brava_catalog.tbl | 8,585 | Additional bulge radial-velocity measurements |
| Prepared joined sample | matched-with-gaia-covariance.parquet | 140,407 | Existing cross-matched analysis input |

The existing quality selection contains 128,772 sources. Preserve its 77,927 training / 24,755 validation / 26,090 test assignments and documented exposure history. This recheck reads file metadata and integrity hashes, not held-out motion outcomes.

The [APOGEE documentation](https://www.sdss4.org/dr17/irspec/dr_synopsis/) describes the spectroscopic measurements. The [StarHorse release](https://data.aip.de/projects/aqueiroz2023.html) confirms the 562,424 APOGEE entries and explains that v2 changed age/mass estimates while leaving distances unchanged. Distances are inferred using stellar models and priors, not direct geometric facts for every star. Inspect the recorded parallax-input flags to avoid counting Gaia parallax twice. No expansion law is needed to turn these local stellar motions into Galactic velocities.

## What already exists on the prediction side

The repository contains a three-dimensional ordinary-matter potential with stellar bar, stellar disks, gas, nuclear components and central mass, plus a conservative extra-potential candidate and orbit integrators. Published mass normalizations and structural choices remain uncertain and can inherit dynamical assumptions; they are not independent measurements of every part of the Galaxy. The empirical extra-potential response is not a first-principles photon-deposition derivation.

## What is still needed, and why

1. **Expand and verify the orbit library.** The existing 72-orbit library covers only 16.87% of the screened training sample jointly in position and velocity at the illustrative 0.5-kpc / 50-km/s matching resolution. Those thresholds are diagnostic choices, not error bars. Missing kinds of orbit cannot be supplied by adjusting their weights. Longer integrations also need numerical and occupation-stability checks.
2. **Fit comparable stellar populations.** Match chemistry and survey fields, handle distance and motion errors, and account for which stars each survey selected. Otherwise differences between sampled populations can look like a gravity effect.
3. **Compare one consistent three-dimensional potential at a time.** Predict the distribution of radial, rotational and vertical motions in the plane beneath the bulge, above and below the bulge, and corresponding disk regions, retaining bar angle and signed height. Extra potential depth creates forces through its spatial slopes; arrival direction need not be the force direction.
4. **Propagate ordinary-matter uncertainty.** Vary justified mass, bar and solar-frame assumptions consistently. Otherwise underestimated ordinary matter can be mistaken for deposited companion gravity. Any standard-halo benchmark is a separate comparison, not an assumed ingredient of the fictional model.
5. **Freeze the procedure before final testing.** Fit training data, use validation for choices, then evaluate reserved tests and a carefully documented BRAVA comparison. Report unsuccessful predictions and uncertainty as well as favorable cases.
6. **Derive the companion link separately.** A successful stellar-motion fit would still require an independently specified photon-to-companion production, transport, capture, support and gravitational response. Placing deposits wherever the fitted extra gravity is needed is a reconstruction, not a prediction.

Additional downloads should be targeted to concrete gaps: survey selection/targeting information and better distance-likelihood or prior information if the chosen statistical implementation needs them. Another copy of these four catalogs will not fix the current prediction gap. The existing likelihood contract records the required selection and distance treatment.

No new claim that the theory matches bulge or disk observations follows from this inventory. The next observational milestone is a sufficiently complete, verified orbit-population prediction with the same assumptions across the compared gravity models.
