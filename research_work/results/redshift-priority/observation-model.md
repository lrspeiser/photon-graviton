# Observation model and motion sensitivity

Purpose: calibrate the proposed void-time gearing without mistaking galaxy motion or endpoint standards for accumulated time stretch. This pass keeps every adopted distance, observed catalog value and baseline parameter unchanged. It introduces no fitted velocity for any galaxy.

## Factor convention

Established factor multiplication and logarithm identities, applied to this hypothesis; not novel:

S_obs = S_time S_motion S_endpoint S_measurement,

T_obs = T_time + T_motion + T_endpoint + T_measurement,

where S = 1 + z > 0 and T = ln(S). T_time is the modeled integral of the fractional wavelength stretching rate. Measurement is a statistical observation factor, not a physical energy reservoir. Define each contribution using consistent source, intermediate and detector standards; a factor already included by a field/clock calculation must not be applied again as an endpoint correction.

Established special-relativistic longitudinal Doppler formula, an optional local-motion model rather than a proposed gravitational law:

S_motion = sqrt[(1 + beta)/(1 - beta)], T_motion = atanh(beta), beta = v_radial/c.

Positive v_radial denotes recession in the declared local comparison frame. This scalar expression assumes longitudinal relative motion; transverse motion and general frame changes require the full photon-direction and observer-velocity transformation. Constant local c and locally normal clocks alone do not prove this special-relativistic matter law in the fictional universe. Retaining it is a declared comparison assumption.

Derived inversion within that optional observation model; originality not claimed:

v_equivalent = c tanh[ln(1 + z_obs) - T_time].

This is the radial speed that would account for the entire remaining shift if motion alone caused it. It is not an independent motion measurement. Never insert it back as a correction and count the resulting perfect fit as evidence.

## Catalog conventions

The [CF4 Appendix, Tables 2–4](https://arxiv.org/pdf/2209.11238) confirms the individual Vcmb frame and distinguishes it from curvature-adjusted and inferred peculiar velocities. Its method-specific distance moduli share an MCMC-registered scale. The audited appendix and CDS ReadMe do not fully document the original optical/radio velocity convention or precise frame-transform implementation for each source. Current calculations preserve the existing Vcmb/c working convention; final precision validation still needs that provenance. No new CF4 target outcome rows were selected or evaluated in this pass.

Do not import the curvature-adjusted velocity, expansion-derived peculiar velocities, or redshift-based Cartesian coordinates as independent environmental evidence. A CMB-frame label specifies a reference convention; it does not by itself require an expanding-universe interpretation, nor does it establish the physical origin of the CMB dipole in this model.

Established spectral coordinate conventions, from [Greisen et al., Representations of spectral coordinates in FITS](https://www.aanda.org/articles/aa/pdf/2006/05/aa3818-05.pdf), not new formulas:

Optical: v_opt = c z.

Radio: v_radio = c z/(1 + z), hence z = (v_radio/c)/(1 - v_radio/c).

Relativistic radial convention: v_rel/c = [(1 + z)^2 - 1]/[(1 + z)^2 + 1].

The word velocity alone is insufficient to choose a conversion. The implementation requires an explicit convention, rejects invalid domains, and verifies round trips. These conventions describe spectral coordinates; they do not establish that the total measured shift is a physical recession speed.

## Calculation and result

Apply the unchanged earlier baseline alpha=0.0002488993286382367 per Mpc (empirically calibrated here; not derived from time physics) to all 164 exposed rows. Add hypothetical longitudinal motions of plus/minus 100, 300 and 1000 km/s through factor multiplication. These are sensitivity scenarios, not population estimates, correction choices or confidence limits.

For positive 300 km/s, the difference between the exact combined shift and the shortcut z_time + v/c spans 0.914–7.194 km/s when expressed as c times a redshift difference. This includes both the product cross term and the small difference between the exact Doppler factor and v/c. It is much smaller than the roughly 450 km/s exploratory cross-validation scatter. Correct composition matters for precision but does not resolve the principal prediction mismatch.

observation-factor-diagnostics.csv preserves all rows, baseline predictions, signed log-factor residuals, equivalent-motion diagnostics and each hypothetical sensitivity. observation-factor-results.json records source hash and 38 checks, including inverse transformations, factor composition, approaching motion, cancellation and invalid-input rejection. These mathematical checks do not verify a physical time field or observational accuracy.

## Uncertainty still required

Conditional interval propagation, derived from the factor identity; not a new physical law:

If N = T_motion + T_endpoint + T_measurement lies between L and U, then

exp(T_time + L) - 1 <= z_obs <= exp(T_time + U) - 1.

The bounds need independently justified meaning and probability coverage; arbitrary speeds do not supply it. For inverse calibration, T_time lies between T_obs - U and T_obs - L under the same assumption. Correlated motions and calibration errors must remain correlated when predicting a group or sky region. A full probability model must also include common-parameter uncertainty and any physical model discrepancy, without silently refitting the fixed user-adopted distances.

We have not obtained independently justified motion/endpoint bounds or a calibrated joint distribution. Consequently no confidence interval, precision floor, corrected time rate or improved predictive performance is claimed. The next data step remains an independent source/convention/environment audit; additional fit flexibility is not a substitute for those inputs.
