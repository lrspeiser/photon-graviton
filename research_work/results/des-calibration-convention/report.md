# Historical DES magnitude-offset convention

The calibration convention is now traced to a pinned SNANA source snapshot and independently checked against the archived DES calibration FITS table. This supplies a usable forward transformation for that historical calibration system; it does not establish compatibility with every later photometry release or apply corrections to observations.

## Source evidence

The pinned source revision is 80866471571978c382000bf387a316c0925f201f, selected as the most recent kcor.c change at or before the archive date. It is not claimed to be the exact revision of the archived executable. The source manifest records immutable URLs and SHA-256 hashes; the archived FITS zero-point table independently agrees with the input offsets.

In [kcor.c](https://github.com/RickKessler/SNANA/blob/80866471571978c382000bf387a316c0925f201f/src/kcor.c#L4456), the AB-system zero-point contribution retains the input reference magnitude. The [synthetic magnitude calculation](https://github.com/RickKessler/SNANA/blob/80866471571978c382000bf387a316c0925f201f/src/kcor.c#L4274) adds that zero point to the spectrum-derived magnitude. This establishes the sign on the synthetic-model side. The retained FITS table contains the same four reference magnitudes and no additional file-based photometry offsets.

## Known calibration algebra

For the input offset delta and physical AB magnitude m_AB,

\[
m_{\rm synthetic}=m_{\rm AB}+\delta,\qquad
F_{\rm native}=10^{0.4(27.5-m_{\rm AB}-\delta)}.
\]

The 27.5 convention is documented for FLUXCAL in the pinned DES release. If native photometry is interpreted in this same synthetic system, its AB-equivalent FLUXCAL is

\[
F_{\rm AB,eq}=F_{\rm native}10^{0.4\delta}.
\]

This is a calibration conversion, not a photon-companion equation. AB-equivalent FLUXCAL is still a conventional flux unit, not directly photons per second. The real-filter/reference-spectrum integral provides the physical count comparison.

| Band | Input delta (mag) | Multiply AB-model FLUXCAL by to predict native FLUXCAL | Multiply native FLUXCAL by for AB-equivalent value |
| --- | ---: | ---: | ---: |
| g | -0.0008 | 1.000737 | 0.999263 |
| r | -0.0116 | 1.010741 | 0.989373 |
| i | -0.0090 | 1.008324 | 0.991745 |
| z | +0.0071 | 0.993482 | 1.006561 |

Prefer predicting the model in native units while retaining the original observations. Applying both a data-side conversion and the same model-side offset would double-correct the comparison. The linear conversion also handles negative noisy fluxes; taking logarithms of those measurements would be inappropriate.

## Shared calibration uncertainty

**Known first-order uncertainty propagation:** if band offsets have covariance V_delta, then at predicted AB-equivalent flux f_i,

\[
J_{ik}=0.4\ln(10)f_i\,1_{\operatorname{band}(i)=k},\qquad
V_{\rm calibration}=J V_\delta J^T.
\]

This is the calibration contribution only; measurement covariance must also be transformed consistently and included. Evaluate this Jacobian at model predictions, not noisy flux measurements. Common band errors correlate different observations and events. They must not be treated as independent errors that vanish by averaging many measurements.

No actual calibration covariance has been supplied in this pass. The test uses a declared artificial two-band covariance to verify the derivative and cross-observation correlations. Nonlinear uncertainty or large offsets would need a fuller treatment.

## Verification and remaining work

All four zero-AB-reference synthetic magnitudes reproduce the archived FITS reference magnitudes to floating-point precision. The first-order covariance agrees with finite-difference derivatives, with maximum absolute difference 7.7e-12 in the artificial test's squared flux units. Negative and zero measurement fluxes remain representable by the linear conversion. No observed light curve is modified.

Still needed: verify the pinned photometry's compatibility with this historical system, acquire the applicable calibration covariance, and constrain source-population evolution for a predictive joint redshift/timing/brightness test. The ongoing artificial timing batch is unchanged. None of the six scientific demonstrations is complete.

```powershell
python research_work/results/des-calibration-convention/fetch_source.py
python research_work/results/des-calibration-convention/verify.py
```

The prior DES calibration acquisition scripts supply the archived FITS input required by verification.
