# Stellar uncertainty, spatial holdouts and the proposed height effect

The prepared observations support a real test, but a simple average of stellar speeds is not a reliable gravity measurement. In the training sample, the apparent bulge height effect changes sign after matching broad spatial and chemical populations. We have now allocated spatial holdouts, propagated measurement uncertainty, and evaluated the frozen photon-conversion formula on the stellar distances. A bar-aware orbital prediction remains necessary before this becomes evidence for or against the gravitational mechanism.

The previous goal turn made concrete progress: it acquired and verified the catalogs. This continuation adds executed uncertainty calculations, prospective model-development splits and a training-only population comparison. The full objective remains active; no strong all-observation case is established.

## Spatial holdouts

| Role | Stars | Plane beneath bulge | Off-plane bulge |
|---|---:|---:|---:|
| Training | 77,927 | 375 | 1,687 |
| Validation | 24,755 | 46 | 1,350 |
| Test | 26,090 | 296 | 1,048 |

The split uses 1,720 sky patches rather than randomly interleaving individual stars. Each nested HEALPix level-4 patch belongs to exactly one role. A salted SHA-256 mapping fixes the allocation independently of velocities, model residuals or desired outcomes. `protocol.json`, `split-groups.json` and the hashed per-star allocation record the exact rule. Gaia's source-ID encoding is documented in the [official data model](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html).

**Limits on independence:** previous all-catalog summaries were computed before this split, so these are prospective model-development holdouts, not pristine blind data. Moreover, 1,133 APOGEE field names span multiple sky patches/roles. Geographic separation does not eliminate shared observing-field calibration or selection effects. The small near-plane validation sample is another limitation. Keep this split fixed, account for shared-field effects and use additional external validation before making a strong claim. We have not inspected held-out orbital-fit scores, because no complete orbital prediction has been run yet.

## Uncertainty calculation

We reconstructed the five-dimensional Gaia covariance for all 128,772 candidate stars from its individual uncertainties and ten correlations. Every correlation matrix was positive definite; no invalid matrix was silently repaired. Gaia's position error in right ascension includes cos(declination), which is accounted for when drawing sky coordinates. Positions and proper motions now use the Gaia coordinates rather than mixing them with APOGEE's coordinates. APOGEE's averaged radial velocities remain the spectroscopic velocity input; stellar/binary evolution across observation epochs is not modeled here.

There are four calculations, each with 64 realizations per star: independent distance/astrometry, and three approximate conditional reconstructions with constant parallax biases -0.037, -0.017 and +0.003 mas. These bracket an assumption; they are not a replacement for Gaia's source-dependent zero-point calibration. Common random numbers make the mode comparisons less noisy.

**Known statistical formula, not a theory-specific law:** conditioning a Gaussian astrometric vector y on parallax p gives

\[
E[y\mid p]=y_{\rm Gaia}+\frac{C_{yp}}{C_{pp}}(p-p_{\rm Gaia}),
\qquad
C_{y\mid p}=C_{yy}-\frac{C_{yp}C_{py}}{C_{pp}}.
\]

Here y contains the two tangent-plane position offsets and two proper motions. We mix this conditional distribution over an approximate StarHorse distance posterior, setting p=1/d+b when d is in kpc and p and b are in mas. The distance draw is asymmetric in log-distance and matches the published 16th, 50th and 84th percentiles; its extrapolated tails are assumed, not supplied by the catalog. **We do not multiply by the Gaia parallax likelihood again.** StarHorse already used it.

This is still a sensitivity reconstruction, not the exact joint StarHorse posterior: StarHorse reprocessed parallaxes and errors, and changing the assumed bias here does not rerun its distance inference. Spatially correlated calibration errors, distance-model systematics, abundance errors and survey selection remain outside these realizations. [StarHorse's release description](https://data.aip.de/projects/aqueiroz2023.html) documents the distance and flag products used.

Across all candidate stars, the conditional -0.017 mas calculation gives median one-star uncertainties of roughly **1.0 km/s in radial motion, 1.3 km/s in rotation and 0.9 km/s vertically**. These all-sample medians conceal the more difficult bulge observations:

| Training region | Typical rotation uncertainty per star | Typical vertical uncertainty per star |
|---|---:|---:|
| Plane beneath bulge | 16.0 km/s | 8.5 km/s |
| Off-plane bulge | 8.2 km/s | 2.9 km/s |
| Disk plane control | 0.9 km/s | 0.9 km/s |
| Disk off-plane control | 0.7 km/s | 0.6 km/s |

On average **5.19%** of stars assigned to a target region change region in a distance/astrometry realization. This requires probabilistic region membership in a final likelihood. The 64 draws introduce Monte Carlo noise of order 9% in an individual Gaussian-like standard-deviation estimate; these files are a pilot for the likelihood, not arbitrarily precise error measurements. Changing the constant parallax bias within the tested range has little effect on the all-sample median errors, but that is not a proof that every source or bulge subsample is insensitive.

## What happens to the apparent height effect?

Only training velocities are used below. The broad comparison uses R=0.5–3.5 kpc for the bulge and R=5–9 kpc for disk controls, with |z|<0.2 kpc versus 0.5–1.5 kpc. These are selected populations, not complete samples of every star in each volume.

| Comparison of off-plane minus plane mean rotation | Bulge | Disk control |
|---|---:|---:|
| Unmatched regional averages | +6.0 km/s | -16.1 km/s |
| Common spatial/chemical cells, at least 3 stars per side | -3.2 km/s | -6.7 km/s |
| At least 5 stars per side | -2.1 km/s | -6.6 km/s |
| At least 10 stars per side | -16.8 km/s | -6.6 km/s |

The matching cells use 0.5-kpc radius intervals, 45-degree bar-azimuth intervals, three iron-abundance bins and two alpha-abundance bins. Each retained cell receives a common weight proportional to the smaller of its plane/off-plane counts. This is a **known empirical standardization procedure**, not a new gravity formula. Changing the minimum count changes the common population being described; it is not repeated measurement of one identical quantity.

At the five-star threshold, only 19 bulge cells remain, containing 215 plane and 374 off-plane stars. At the ten-star threshold only six remain. The bulge result is sensitive to population overlap. The disk comparison is more stable over these thresholds, but slower average rotation away from the plane is not automatically weaker gravity: the orbital distribution and random motions also matter. No significance interval or force inference is asserted from these coarse comparisons. They have not been deconvolved for measurement errors, extinction/selection or full abundance uncertainty. See `training-overlap-comparison.json` for all results, including radial and vertical motions.

## What the existing photon formula predicts for these stars

**Known exponential accumulation solution applied to the proposed conversion mechanism:**

\[
z_c=\exp(\alpha D)-1,\qquad
f_{\rm transferred}=\frac{z_c}{1+z_c},\qquad
\alpha=0.0002488993286382367\ {\rm Mpc}^{-1}.
\]

This coefficient is frozen from the earlier galaxy-redshift fit, not retuned using stellar velocities. We evaluated it for every candidate, preserving its holdout role. The median extra spectroscopic shift is equivalent to **0.233 km/s**; the central 68% span approximately **0.106–0.499 km/s**. The parent catalog includes very distant sources outside the four Milky Way target regions, so its extreme tail should not be interpreted as a typical bulge signal.

The per-photon remaining-plus-transferred energy fractions sum to one numerically. This is a bookkeeping identity, **not** a successful total-energy supply calculation or evidence that companions exist. The predicted shift has not been subtracted from the observed velocities in the above tables. Its comparison with stellar motions requires a joint Doppler/endpoint/conversion forward model; stellar radial velocity is not an independent measurement of cosmological redshift. `stellar-frozen-redshift-predictions.parquet` supplies the prediction for that later model.

## What this changes in the next physics test

The currently drafted potential-response kernel gives an acceleration when a deposit distribution is specified. It does not yet give the mixture of orbits at each observed location. Therefore a claimed predicted speed for each bulge star would introduce an unsupported assumption.

The next model must simultaneously specify the bar-aware ordinary-matter potential, the deposit-induced potential and an orbital population/selection likelihood. The requested standard-halo comparison must use the same observational treatment. Fit only training stars; use validation to choose among declared variants; score the reserved test sample only after freezing them. Propagate measurement covariance and region migration, and assess shared-field systematics. A spatially arbitrary correction per cell is not an acceptable substitute for a shared predictive law.

The broader redshift, supernova-timing, spectral, lensing and energy requirements remain in force. These new calculations improve the credibility of the test by preventing a population effect from being mistaken for gravity; they do not complete the theory or remove its earlier unresolved constraints.

## Reproduction

```powershell
python research_work/results/stellar-holdout-audit/run.py
python research_work/results/stellar-holdout-audit/training_comparison.py
python research_work/results/stellar-holdout-audit/verify.py
```

`results.json` records the file hashes and uncertainty summaries. Large per-star allocations, covariance files and frozen photon predictions are in the ignored `research_work/data-cache/stellar-catalogs/` directory. `verification.json` checks pixel-role separation, identifier integrity, covariance-output sanity, an independent Gaussian conditioning experiment and the per-photon energy identity. None of those numerical checks is described as an orbital-model validation.
