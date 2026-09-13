# Frozen retention improves outer motions but not lensing

The fitted retention law was carried into the six existing lens systems with all capture and optical parameters frozen. Only ordinary stellar mass and orbital anisotropy were refitted to inner stellar bins. The result improves outer-motion residuals but worsens lens-angle errors; no joint success is claimed.

## Luminosity mapping and its limits

The available lens products lack directly calibrated rest-frame 3.6-micron luminosities. We therefore executed the declared **conditional mapping**: independently inferred photometric stellar masses, rescaled to the retained optical luminosity-distance law, divided by 0.5 Msun/L_sun to provide a 3.6-micron proxy. Both Chabrier and Salpeter population assumptions were evaluated, without choosing the better result. The proxy is not the dynamically fitted inner stellar mass.

The retained regular optical geometry adds the factor B=1+f/(1+q_opt*f) to the earlier photometric luminosity normalization. All physical radii, apertures, seeing and lens distances use that same geometry. The equivalent disk scale remains Re/1.67834699; that morphology mapping is an assumption. Stellar population, age, dust and band-conversion uncertainties are not fully propagated. The two population choices are a systematic sensitivity bracket, not confidence intervals. Thus this is not a completed observational transfer using measured L_3.6.

## Results with fitted retention q=0.3399907823

| Quantity | Previous intercepted source | Retention, Chabrier proxy | Retention, Salpeter proxy |
|---|---:|---:|---:|
| Inner stellar covariance score | 48.91702 | 42.04076 | 41.43718 |
| Outer conditional residual-square sum | 63.03850 | 44.04208 | 42.42611 |
| Lens fractional RMS | 13.11169% | 13.93772% | 14.10972% |

| Lens | Catalog angle (arcsec) | Chabrier prediction | Salpeter prediction |
|---|---:|---:|---:|
| J0037-0942 | 1.53 | 1.82712 | 1.83317 |
| J1112+0826 | 1.49 | 1.47756 | 1.48095 |
| J1204+0358 | 1.31 | 1.53142 | 1.53365 |
| J1402+6321 | 1.35 | 1.63944 | 1.64170 |
| J1621+3931 | 1.29 | 1.37183 | 1.37337 |
| J1630+4520 | 1.78 | 1.75388 | 1.75890 |

Catalog angles are lens-model summaries; the RMS is descriptive, not a full likelihood with all angle/systematic uncertainty. All six systems were previously exposed. No stellar orbit fit hits its bounds. Mapping retention probabilities range approximately 0.713–0.818 across scenarios. Identical geometry and stellar observations relative to the original regular-optics calculation are verified, as are retention probability bounds.

**Decision:** the galaxy improvement does not yet carry through to a joint lensing/motion improvement. Keep this result as a conditional failure to improve lensing, while recognizing the unmeasured band/morphology mapping. Do not adjust the capture normalization per lens. The candidate remains a galaxy-rotation improvement, not a unifying theory.

## Is 0.340 effectively one-third?

The fitted value differs from 1/3 by 0.00665745, about 2% of 1/3. We tested exact 1/3 independently of the lens run by fixing q and refitting only the other three parameters on the same 89 training galaxies, then freezing them for the 60 comparison galaxies.

| Galaxy quantity | Fitted q=0.3399908 | Fixed q=1/3 |
|---|---:|---:|
| Training log RMS | 0.1384879 | 0.1384891 |
| Validation RMS km/s | 32.52002 | 32.49486 |
| Test RMS km/s | 23.62155 | 23.59079 |
| Validation log RMS | 0.1153834 | 0.1153940 |
| Test log RMS | 0.0909910 | 0.0910583 |

The practical differences are small. The simpler candidate is

    eta(X)=X^(1/3)/(1+X^(1/3)).

This is **known cube-root/logistic mathematics applied as our hypothetical retention law**. The result supports considering the exact fraction for simplicity; it does not establish a fundamental exponent or a three-dimensional geometric cause. No parameter uncertainty or confidence interval was estimated here.

With q fixed, the new training constants are C_half=47,285,892.42 Msun/kpc^3, k0=0.2039029004/kpc and a/R_disk=2.770766589; C_source=2*C_half. The result's boundary flag is true because q is deliberately fixed with identical bounds, not because a free parameter was forced to an edge. All three starts succeed and agree; finer integration changes RMS by at most 0.00113 km/s. The lens table above uses the original fitted-q constants, as declared, and must not be presented as a separate one-third lens test.

Reproduce:

    python research_work/results/isotropic-galaxy-transfer/lensing.py --retention-optics
    python research_work/results/isotropic-galaxy-transfer/run.py --bounded-retention --third-retention

Full lens and outer-star predictions are in `retention-optics-results.json`. Fixed-third galaxy predictions and constants are in `third-radiation-retention-predictions.json` and `third-radiation-retention-results.json`. All six goals remain open.
