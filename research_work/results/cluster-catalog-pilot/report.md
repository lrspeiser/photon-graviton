# Catalog-fed NGC4486-centered supply pilot

## Result

We have applied the photon-to-companion capture calculation to actual cached catalog positions, distances and fitted luminosities. This is a conditional present-luminosity supply estimate, not accumulated gravity or a successful cluster lensing fit.

The receiver is an explicitly illustrative sphere of radius 1 Mpc centered on the catalog position of NGC4486, with alpha=0.0002488993265191759 per Mpc and kappa R=10. The radius and capture strength are stipulated; alpha is the archived empirical redshift coefficient, not a derived interaction constant. No cluster mass, dark-matter profile, expanding geometry or cosmic age is used as a target or fitted input.

## Catalog and exclusions

The cached DustPedia tables contain 875 rows; 814 have usable positive distance and bolometric luminosity in each dust-model version. The 61 omitted names and hashes are preserved in results.json. Columns are decoded from the accompanying byte-by-byte ReadMe, including J2000 coordinates, distance in Mpc and Lbol in solar luminosities. All names are unique among usable rows. Source: Nersesian et al. (2019), catalog J/A+A/624/A80; local files companion_causal_test/data/{themis,dl14}.dat and ReadMe.

Published distances are used as stipulated facts for this pilot, not inferred here from redshift. Luminosities are CIGALE spectral-energy-distribution inferences, NOT raw measured photon counts. Their stellar-population, dust and adopted brightness-law assumptions have not been refitted under our alternative propagation model. The THEMIS/DL14 comparison tests only one modeling choice; it is not a full uncertainty estimate or proof of independence from cosmological assumptions. This limitation prevents treating the calculated power as an established physical supply.

Each galaxy is treated as an isotropic point emitter. Catalog three-dimensional separation alone classifies 17 sources inside the chosen sphere and 797 outside; this is not a dynamical cluster membership measurement. Distance uncertainty and galaxy spatial extent can change the classification near the boundary. The catalog is not a complete inventory of the universe or all emitters around the cluster.

## Calculated deposited power

| Catalog SED model | Internal source power retained (Lsun) | External source power retained (Lsun) | Total (Lsun) |
|---|---:|---:|---:|
| themis | 6.57342e+07 | 1.52196e+08 | 2.1793e+08 |
| dl14 | 6.61015e+07 | 1.53133e+08 | 2.19235e+08 |

External sources contribute 69.84% of the THEMIS pilot's deposited power. The two dust-model totals differ by 0.599%. This supports the practical importance of including outside sources in this calculation, not the physical validity of the assumed capture rate.

## Equations and verification

Known Euclidean coordinate conversion gives each source-to-center separation. The external-source kernel from ../cluster-external-photon-supply/run.py is reused by loading its imports and function definitions without running its standalone driver. Its exact interception angle and coupled transfer solutions are documented there. Internal sources use the actual catalog source radius, not the preceding uniform-emitter approximation: integrate the same coupled photon/companion energy solution over all forward exit lengths L=-r mu+sqrt(R^2-r^2+r^2 mu^2). At the center all exit lengths equal R. Each source contributes Lbol times its retained fraction.

These are established geometry/transport equations applied to the project's optional constant-rate postulates. No equation is claimed unique, and no mechanism for alpha, capture, storage support or gravitational response has been derived. Source-energy fractions were independently verified in the preceding kernel tests. Doubling the angular quadrature for every catalog entry changes retained fractions by less than 1e-9 absolute. The input and kernel hashes, every source contribution and both model runs are retained.

The output is a conditional deposition POWER assuming established illumination. It is not an accumulated energy: applying these present luminosities indefinitely would invent an emission history. To compute a reservoir, specify and causally integrate each source history, include intervening capture competition and any release/work terms, and evolve the deposit. Point-source emission, no intervening removal and fixed opacity are assumptions, not completeness corrections or measured capture efficiencies.

## What remains

Refit or bound photon luminosities with a consistent brightness/time law; establish uncertainty and survey completeness; evaluate source histories and intervening absorbers; calculate where the delivered energy resides and moves; and predict lensing and matter motion with that same reservoir. This pilot makes no comparison with a required cluster mass and cannot establish an energy shortage or surplus. No final observational holdouts were opened; all six objectives remain open.

Run `python research_work/results/cluster-catalog-pilot/run.py`.
