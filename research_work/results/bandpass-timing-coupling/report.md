# Linking spectral redshift to timing through real passbands

A controlled source model was integrated through the verified DES g/r/i/z response curves. This advances the forward-model requirements for joint redshift, timing and brightness. It does not fit real supernovae or claim that the constructed source spectrum describes them.

## Assumptions and provenance of formulas

**Declared artificial source, not a proposed new physical law:** the emission-frame duration depends on emitted wavelength as W(lambda_e)=30 days (lambda_e/5500 angstrom)^q, with Gaussian temporal profiles whose full width at half maximum is W. The peak spectrum is proportional to lambda_e^-2 over 1000-20000 angstrom and zero outside. All tested filter mappings lie strictly within this finite range. We test q=-0.5, 0, 0.5, 1 to expose the role of source color evolution, not to tune an astronomical fit.

**Our conditional transport postulates:** measured spectral stretch S, event stretch S^b, photon-number conservation and fixed geometric dilution. b=0 is no event stretching; b=1 is common spectral/event stretching.

**Known change of variables and photon-counting radiometry:**

\[
F_{\lambda_o}(t_o)\propto S^{-(b+2)} L_{\lambda_e}(\lambda_o/S,t_o/S^b),\qquad
C_j(t_o)\propto\int T_j(\lambda_o)\lambda_o F_{\lambda_o}(t_o)d\lambda_o.
\]

The omitted h*c and collecting geometry are fixed and cancel in the tested count ratios. Absolute flux calibration is not being claimed.

## Derived fixed-filter result

For this self-similar source and a fixed observer filter, substitution gives

\[
C_{j,S}(t)=S^{-b}C_{j,1}\!\left(t/S^{b-q}\right),
\qquad W_{j,\rm obs}(S)=W_{j,\rm obs}(1)S^{b-q}.
\]

Thus the apparent exponent from a fixed-filter width is b-q, not necessarily b. The derivation does not require a narrow filter: it holds after integrating through the complete response for this specific spectral/temporal family. Real spectra need the full numerical forward model and need not have this simple power-law behavior.

Examples at S=2.2 (z=1.2):

| Imposed event exponent b | Source wavelength exponent q | Fixed-filter width ratio | Event stretch actually imposed |
| ---: | ---: | ---: | ---: |
| 1 | 0 | 2.200 | 2.200 |
| 1 | 0.5 | 1.483 | 2.200 |
| 1 | 1 | 1.000 | 2.200 |
| 0 | -0.5 | 1.483 | 1.000 |

A width measured through the wrong emitted wavelengths can hide real stretching or suggest stretching where none was imposed. This is an observational-model issue, not evidence for or against expansion or companions.

## Matching emitted wavelengths

An ideal comparison uses response T_j(lambda_o/S), keeping the emitted wavelength response fixed. Direct integration then gives W_obs(S)=S^b W_obs(1), independent of q for this source. This is an ideal shifted-filter diagnostic, not a claim that DES possesses a continuously movable filter. Actual comparisons must find overlapping emitted passbands or model the spectra through each actual response. Matching filter centers alone does not establish identical emitted response shapes.

The earlier DES methodology audit already noted emitted-wavelength matching. This experiment does not overturn the published timing result; it explains why that part of the method is necessary and why our own central-wavelength approximation is insufficient for a precision inference.

## Verification and consequence

All 320 integrations cover four filters, two b values, four q values, five S values, and both fixed-observer and ideally matched-emission modes. Half-maximum roots agree with the analytic width scalings to relative error below 1.9e-13. These are implementation checks on an artificial source, not 320 observational successes. Filter file hashes are checked against the acquired archive manifest.

Our 160-case timing calibration remains useful for its declared simplified source shapes and cadence patterns. It does not yet include this full spectral/passband effect. Before real-data interpretation, the source model must predict wavelength-dependent duration and brightness through the verified filters, and the resulting estimator needs injections with those effects. The photon-companion cause and source-evolution ambiguity remain unresolved; all six scientific demonstrations stay open.

```powershell
python research_work/results/bandpass-timing-coupling/run.py
```
