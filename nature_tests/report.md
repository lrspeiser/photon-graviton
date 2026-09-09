# Tests motivated by natural mechanisms

9 September 2026. This report separates a comparison with published cavity measurements from synthetic field and optical-medium calculations. No raw cavity likelihood, new supernova sample, angular-distance measurement, or completed microscopic interaction is claimed.

## Outcome

The strongest new observational finding is a tension between the unscreened fixed-atom, fixed-length index model and reported long-term silicon-cavity drifts. The field reconstruction conserves energy but is not a robust recovery mechanism for the fitted finite-time history. A lossless resonant medium illustrates how temporal frequency conversion can become nearly achromatic, while also exposing an atomic-binding problem if its dielectric response is applied universally. These results narrow the acceptable mechanisms; they do not validate a complete theory.

## 1. Published cavity drift comparison

Source: Lee et al., Frequency stability of 2.5×10^-17 in a Si cavity with AlGaAs crystalline mirrors, arXiv:2509.13503v1, https://arxiv.org/html/2509.13503v1 (long-term drift section and Figure 3); journal DOI https://doi.org/10.1103/zgrm-cjbb . The quantitative comparison here uses the explicitly identified preprint version. The authors report recent fractional slopes of -2.6e-19, -7.7e-20, -1.0e-19 and +5.2e-20 per second for Si2, Si3, Si5 and Si6. Their optical frequencies were tracked against hydrogen masers or optical clocks; instrument drifts have an unresolved physical origin. Short-term stability plots remove linear drift and are not the quantities used here.

The temporal model predicts -gamma = -7.7315e-11/year = -2.4499645e-18/second for a fixed-length cavity relative to a fixed atomic standard. At 1542 nm this is -476.3 microhertz/second, compared with reported approximate recent absolute drifts of -50, -15, -20 and +10 microhertz/second. Using the published fractional values, the predicted magnitude is 9.42, 31.82, 24.50 and 47.11 times the respective measured magnitude. The positive Si6 drift has the opposite sign.

These are scale comparisons, not sigma values: measurement uncertainties and a raw-data systematic covariance were not supplied to this calculation. The model alone does not reproduce the reported total drift. To retain its universal contribution requires positive additional fractional drifts of 2.19–2.50e-18/second in the four instruments, substantially canceling the predicted term. Such contributions were calculated as requirements, not independently established. With one freely chosen instrumental offset per cavity, a common universal slope is exactly unidentifiable: the four-row design matrix containing a common column and four individual columns has rank four for five parameters. An exclusion therefore needs independent constraints on those offsets, not arbitrary fitting of them.

## 2. Disturbance test without retuning

Use the previous reconstructed potential V/(K gamma^2)=-n^(2p)/2-beta/n, p=0.2639061166 and beta=U0/(K gamma^2). Integrate in tau=gamma t:

n'=v; v'=p n^(2p-1)-beta/n^2+u/n; u'=-vu/n.

Here u is radiation energy density in units K gamma^2. The potential's beta is fixed while radiation or initial velocity is perturbed. Initial n=0.3; the baseline integration ends when its analytic trajectory reaches n=1.5. Cases use beta=0.1 and 1, radiation changes ±10% and ±50%, and velocity changes ±10%. These are illustrative dimensionless strengths, not inferred radiation densities.

At beta=0.1, ±10% radiation changes leave final speeds approximately +2.12% and -2.19% away from the target law v=n^p. At beta=1 the same changes leave +18.73% and -26.21% deviations. Reducing radiation by 50% in the beta=1 case reverses the trajectory and reaches the imposed n=0.05 stopping boundary. Initial velocity changes of ±10% leave approximately ±4.3% speed-law deviations at the end, so some relative sensitivity diminishes but the original trajectory is not recovered over the test interval.

Energy and the photon invariant remain conserved numerically; the largest absolute dimensionless energy drift is 3.56e-9 in the boundary-reaching case. This is not numerical energy loss masquerading as feedback. Distinct conserved invariants prevent a unique full-state attractor that erases all those perturbations. This does not rule out asymptotic agreement of a reduced fractional quantity at much later times, nor feedback in a larger open subsystem. It does show that the current finite-time reconstruction is not a demonstrated robust regulator. Any damping repair must include the degrees of freedom receiving the energy.

## 3. Simple optical-medium test

Motivation: temporal frequency conversion in physical media is real, but its driver must be modeled. See Galiffi et al., Electrodynamics of photonic temporal interfaces, https://arxiv.org/abs/2411.15984 and https://doi.org/10.1038/s41377-025-01947-2 . The calculation below is our own lossless Lorentz-oscillator endpoint diagnostic, not a reproduction of their experiment or a full time-dependent driver simulation.

Assume epsilon(omega)=1+F/(Omega^2-omega^2), mu=1, c=1, conserved wave number k, and adiabatic following of the lower electromagnetic branch. Start at F=0. Set F=3(Omega^2-0.25) at the endpoint so a reference mode with initial frequency 1 ends at frequency 0.5. The lower branch solves omega^4-(Omega^2+F+k^2)omega^2+k^2 Omega^2=0. Test initial frequencies 0.5, 1 and 2 for resonance frequencies Omega=3,10,100.

Final-to-initial frequency ratios are [0.5039597,0.5,0.4846893] for Omega=3; [0.5003520,0.5,0.4985962] for Omega=10; and [0.50000352,0.5,0.49998594] for Omega=100. Moving the resonance far above the band makes the desired factor-of-two shift increasingly uniform. The largest relative deviation from a 0.5 ratio decreases from 3.06% to 0.281% to 0.00281%.

However, the corresponding static permittivity approaches four. If this response also governs local Coulomb binding with fixed charges and masses, the leading optical binding energy scales as epsilon(0)^-2 and approaches 1/16 of its original value. Thus the simple universal dielectric implementation has a large atomic effect even when its traveling-wave shift is nearly achromatic. This atomic inference assumes the same local response law remains valid at atomic spatial scales; a real material's macroscopic dielectric function need not apply there. A proposed vacuum interaction must derive its spatial response rather than assume either extrapolation or protection.

Only endpoint frequencies were tested. No pulse-duration, absorption, temporal-reflection, photon-number, or total driver-energy calculation is claimed. The externally prescribed oscillator-strength change is not yet a self-sustaining cosmic field.

## Next decision

Keep the empirical history as a candidate, but downgrade the simple universal fixed-matter index as an established mechanism. Prioritize a coupling that calculates spatial and temporal electromagnetic response together, with the cavity slopes as quantitative targets. A mechanism that distinguishes traveling radiation, cavities and atomic binding must do so through equations, not independent exemptions. Before adding damping, establish whether the scalar's degrees of freedom and energy transfers can provide recovery without retuning the radiation-dependent potential. Curvature, rotation and CMB construction remain parked.

## Reproduction

Run nature_tests/run.py with Python, NumPy and SciPy. results.json records the complete test grid, source values, scale comparisons and limitations. Cavity inputs are transcribed published summaries, not raw data. Numerical ODE and dispersion calculations are synthetic. No new independent observational validation has been performed.
