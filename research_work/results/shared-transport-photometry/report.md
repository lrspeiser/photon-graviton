# Shared transport: brightness, spectra, distance calibration and clocks

2026-09-10. Stage 3 conditional predictions using the existing conversion rate; no fitted extra photometric parameter, fresh observation or complete interaction theory. [Protocol](protocol.md), [calculation](run.py), [predictions](predictions.csv), [results](results.json).

The candidate in which light and gravitational waves share the same evolving arrival map survives the earlier relative-delay check because their common propagation delay cancels. This does not prove that the source signals were emitted simultaneously or that their speeds differ. The observed 1.74-second lag includes source timing; later counterpart discovery times are not wavelength-dependent flight delays. That distinction was already included algebraically and is now explicit in the joint audit.

## One stretch factor gives linked photometric predictions

**Proposed specialization:** static Euclidean geometric dilution, unchanged ordinary endpoint standards, conserved photon number, no angular redistribution/lensing, and common achromatic frequency and event stretch S=exp(alpha D). The existing empirical alpha is 0.0002488993286382367 per Mpc. Equal photon/GW transport and these optical assumptions still need a joint causal interaction and matter-clock completion. This is not a demonstrated first-principles model.

**Known change-of-variables and radiometry, conditional on those assumptions:** if an emitted spectral luminosity is L_nu and source time is t_e,

\[
F_{\nu_o}(t_o)=\frac{L_{\nu_e}(S\nu_o,t_e)}{4\pi D^2 S},\qquad
\frac{dt_o}{dt_e}=S.
\]

Each photon has 1/S of its original energy; events last S times as long; an observed frequency bin corresponds to S times its width at emission. Keeping all three factors produces the above spectral flux, whose integral is

\[
F_{\rm bol}(t_o)=\frac{L_{\rm bol}(t_e)}{4\pi D^2S^2}.
\]

Thus, relative to ordinary static geometric dilution, photon arrival rate is reduced by 1/S, instantaneous total-energy flux by 1/S^2, and time-integrated event energy per area by 1/S. The latter preserves the stipulated energy loss per photon; the extra reduction of instantaneous brightness is caused by event stretching, not additional energy destruction. The fate of transferred energy and momentum still requires the receiving sector.

| Static path | Predicted transfer redshift | Instantaneous brightness remaining | Extra bolometric magnitude |
|---|---:|---:|---:|
| About 100 million light-years | 0.0076605 | 98.485% | 0.01657 |
| 40.7 Mpc | 0.0101817 | 97.994% | 0.02200 |
| 100 Mpc | 0.0252023 | 95.144% | 0.05405 |
| 1000 Mpc | 0.2826129 | 60.787% | 0.54048 |

These are model predictions, not measured rows. Extrapolation to 1000 Mpc is not validated by the nearby calibration. Real passbands require integrating the shifted spectrum through the instrument response.

## Spectral shape is a separate test from absolute brightness

**Conditional deduction using the known Planck formula:** for a source thermal spectrum at temperature T,

\[
\frac{B_\nu(S\nu,T)}{S}=S^2 B_\nu(\nu,T/S).
\]

An ideal deterministic rescaling retains a Planck-shaped spectrum at T/S with changed amplitude. It does not broaden a line relative to its center or introduce a frequency-dependent group delay in this specified achromatic map. These statements follow from the assumed map, not from a derived photon-companion scattering amplitude. Quantum receiving states or stochastic interactions can predict additional distortions and need their own calculation.

If source angular size, emissivity and calibration fix the normalization independently, amplitude is an observable and cannot be freely adjusted to declare a pass. A free-amplitude spectral-shape fit alone is a weaker test. This isolated-source flux calculation is not the same as evolving an initially homogeneous radiation density: the earlier FIRAS fixed-volume/number-preserving history used different boundary and transport assumptions. It would be incorrect either to import its failed score into this source-flux calculation or to claim this identity has solved the cosmological thermal-background test.

## Brightness affects the distances used to test the model

**Conditional standard-candle and ideal bolometric SBF bookkeeping:** if a reduction ignores the added S^-2 dimming and its luminosity calibration is unaffected, it infers D_app=D*S. For SBF, the variance/mean ratio of unresolved stellar fluxes scales with the same per-star flux factor when population and passband are fixed. This is an idealization, not a correction to a real published SBF pipeline.

With S=exp(alpha D), the known Lambert W inverse gives

\[
D=\frac{W(\alpha D_{\rm app})}{\alpha}.
\]

For the published 40.7-Mpc SBF host distance used in the GW170817 diagnostic, this idealized sensitivity gives 40.294 Mpc, about a 1% distance effect. We do not replace the previously stipulated distance or claim the actual reduction omitted every relevant correction. The measured distance includes passband, stellar-population and calibration assumptions that require direct auditing. [Host-distance measurement](https://arxiv.org/abs/1801.06080)

If calibration used a geometric anchor with stretch S_a, its inferred luminosity carries that dimming too: the ideal relative relation is D_app=D*S/S_a. Calibration is therefore a shared nuisance, not an independent per-galaxy correction that may be tuned to remove residuals.

## Clock prediction: the surviving missing piece

Let C_e and C_o be ticks per coordinate time of identical atomic standards at emission and observation. Let S_coord be the stretch generated by the coordinate transport map. **Known frequency and clock-ratio algebra**, assuming the source transition follows that clock scale:

\[
1+z_{\rm measured}=S_{\rm coord}\frac{C_o}{C_e},\qquad
\frac{\Delta\tau_o}{\Delta\tau_e}=S_{\rm coord}\frac{C_o}{C_e}.
\]

Equal endpoint clock factors preserve the coordinate stretch; C_o/C_e=1/S_coord cancels it in both observables. The hypothesis must derive this ratio rather than selecting it after seeing the result. Ordinary endpoint clocks do not establish how clocks behave inside a void.

For rods of coordinate length R and clocks ticking at rate C, the measured local propagation speed is v/(R*C). Altering coordinate wave speed while leaving both rods and clocks untouched does not preserve locally measured c. The previous homogeneous universal-metric completion showed that its nonexpanding branch cancels measured redshift, while the branch retaining it makes physical spatial separations expand. That conditional result remains a problem for that completion, not a universal exclusion of every nonexpanding interaction.

This audit therefore does not mark clocks as passed. Separate matter coupling could be proposed, but its local speed, transition frequencies, forces and conservation must then be calculated. Equal photon/GW delays alone do not supply those equations or prove companions retain their energy.

## Verification and test status

Independent quadrature verifies bolometric energy, photon rate and event fluence. The Planck integral matches pi^4/15 within 9e-16, and the transformed spectral shape identity within 7e-16 over the sampled range. Bracketed root finding independently matches the Lambert W distance inversion within 8e-15 Mpc. These are checks of conditional equations, not observations or theory validation.

Next photometric testing needs observer fluxes, passbands, source populations and calibration covariances under the same fixed mechanism. The physical clock/receiver coupling remains unfinished, and map selection remains unresolved for stage 2. No observational holdout may be labeled a full-model validation until these definitions are frozen. All four goal stages remain active and incomplete.
