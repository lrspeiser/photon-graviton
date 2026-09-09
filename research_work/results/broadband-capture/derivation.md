# Derive the receiver response and test the whole companion spectrum

## Why we need this

The previous calculation found that an absorption cross section proportional to energy^-4 can cancel velocity slowing in an ideal monochromatic bath. That was a specified function, not a material. We now derive a related cross section from an explicit oscillator and test what happens when the incident energies are spread out. Without this, a favorable single-energy example could be mistaken for a receiver that works with the actual companions produced by stars.

The oscillator is an optional scalar-companion comparison, not an ordinary graviton detector or an adopted material. Thermal spectra are used as a diagnostic class; we do not assume that companions are thermal, adopt the Big Bang, or assign them an observed background temperature.

There is related established work on [Einstein–Hopf drag and blackbody friction](https://arxiv.org/abs/1112.5195). Our scalar receiving-channel cross section and broadband identities below are conditional derivations for this project, not results attributed to that electromagnetic paper.

## 1. A cross section from a specified interaction

Reuse the worldline oscillator coupling g q chi from the [earlier oscillator action](../oscillator-response/derivation.md), with canonical scalar field chi. In the receiver rest frame and hbar=c=1,

    q_double_dot + gamma q_dot + Omega^2 q = g chi_in
    gamma = gamma_rad + gamma_rec
    gamma_rad = g^2/(4 pi).

The radiative term follows from the scalar far-field flux. Gamma_rec represents an additional passive receiving continuum, such as internal material modes, whose energy must be tracked. A continuum spectral density approximating constant damping over the band is assumed; it is not derived from an identified material. Finite-size regularization and a cutoff above the useful band remain necessary.

For chi_in=A cos(omega t), the response amplitude is g A/[Omega^2-omega^2-i gamma omega]. The incident scalar energy flux is omega^2 A^2/2, while the two outgoing powers are gamma_rad<q_dot^2> and gamma_rec<q_dot^2>. Therefore

    sigma_rec(omega) = g^2 gamma_rec / [(Omega^2-omega^2)^2+gamma^2 omega^2]
    sigma_sc(omega) = g^2 gamma_rad / [(Omega^2-omega^2)^2+gamma^2 omega^2].

The first transfers energy into the receiving continuum; it does not prove permanent storage. The second reradiates companions. The code verifies cycle-averaged input work equals both powers combined. It also checks the scalar optical theorem using the forward amplitude g^2/[4 pi(Omega^2-omega^2-i gamma omega)]. Treating reradiated energy as deposited would violate this distinction.

Both cross sections have the same spectral shape when the damping parameters are constant. They have an inverse-fourth-power high-frequency tail within the validity of that response approximation. Unlike an unregulated omega^-4 law extending to zero frequency, the response is finite at zero frequency for nonzero Omega. This supplies a candidate origin for the shape; it does not establish sufficient capture strength or usable storage.

For the absorption-only force calculation below, sigma means sigma_rec. Scalar monopole scattering is isotropic in the receiver frame and adds a force contribution with the same spectral weighting in this model. It cannot simply be omitted from an absolute momentum calculation. Its receiving-energy contribution is zero, so the total force per unit stored power also depends on the branching ratio. The normalized coefficients below are not complete force-to-permanent-deposit efficiencies.

## 2. A broadband criterion

Let u(epsilon) be lab spectral energy density per unit companion energy, not per logarithmic interval. For isotropic illumination and slow receiver speed, integrating the previous angular result over energy gives

    P_abs = c integral u(epsilon) sigma(epsilon) depsilon
    dv/dt |_abs = -C P_abs v/(M c^2)
    C = integral u sigma [4+dlnsigma/dlnepsilon] depsilon
        / [3 integral u sigma depsilon].

This assumes full retention of the energy counted in P_abs during this local absorption step. Subsequent radiation, heating, reverse events and recoil support require their own accounts.

Integration by parts supplies an independently useful form:

    C = { [epsilon u sigma]_0^infinity
          + integral sigma(3u-epsilon du/depsilon) depsilon }
        / [3 integral u sigma depsilon].

Endpoint terms matter. If the integrals are finite and the endpoint term vanishes, cancellation depends on the weighted **incident spectral shape**, not only a local cross-section slope. Defining m=dlnu/dlnepsilon, the numerator weights 3-m. An actual companion spectrum must come from the same production and propagation calculation; it cannot be replaced by whichever spectrum makes capture convenient.

## 3. A positive result for thermal spectra: slowing cannot cancel

For a Planck spectral energy density u proportional to epsilon^3/[exp(epsilon/T)-1], with T expressed in energy units, write x=epsilon/T. Direct differentiation gives

    3u-epsilon du/depsilon = u x/(1-exp(-x)).

For nonnegative passive sigma, nonzero absorption, finite integrals and vanishing endpoint term,

    C = weighted_average[x/(1-exp(-x))]/3 > 1/3.

Thus this isotropic thermal-spectrum absorption class has positive velocity slowing whatever regular passive spectral filter is chosen. This is not a universal statement about nonthermal companion baths, anisotropic illumination, active receivers or all field theories. Nor is it a Big-Bang premise: the result concerns the mathematical shape of thermal radiation alone.

The numerical oscillator examples check the direct derivative expression against the manifestly positive expression at two damping ratios and three temperature ratios. Agreement is better than 1e-7 in C. Values range from about 0.353 to 2.210. A constant cross section would instead give C=4/3 for any finite broadband spectrum with suitable endpoints.

These numerical examples extend the ideal response across the integration axis. They do not establish its ultraviolet validity. A physical cutoff and additional channels can change the numerical coefficient; the positive thermal identity still applies to a resulting regular passive cross section under the same assumptions.

The earlier pure sigma proportional to epsilon^-4 counterexample cannot be extended unchanged to the whole thermal spectrum. At low energy, u is proportional to epsilon^2, so u sigma is proportional to epsilon^-2 and its power integral diverges. The code demonstrates that cutoff dependence. A physical regularization must be included in the force calculation; it cannot be ignored while retaining the cancellation.

## 4. Narrow nonthermal bands can behave differently

For the oscillator set e=epsilon/Omega and d=gamma/Omega. Its spectral slope is

    s(e) = -[4e^4+2(d^2-2)e^2] / [(1-e^2)^2+d^2 e^2].

The single-energy low-speed coefficient (4+s)/3 vanishes at

    e^2=2/(2-d^2), for d^2<2.

This is one tuned energy, not exact cancellation at every energy or every velocity. Far above resonance the coefficient approaches zero from a side determined by d; the cross section simultaneously falls steeply, so weak slowing alone does not demonstrate useful capture throughput.

For d=0.1, center a lognormal spectral energy distribution at the single-energy cancellation point e=1.002509. The following results retain the same response and vary only the incident band's standard deviation in log energy:

| Log-energy width | Band-averaged C |
|---|---:|
| 0.001, approximately 0.1% | 0.00185 |
| 0.05, approximately 5% | 0.82469 |
| 0.2 | 1.19676 |
| 0.5 | 1.27666 |

These are specified test spectra, not measured companion distributions. Broadening this particular centered band restores slowing. Other nonthermal bands can still give negative or small coefficients; the saved comparisons include a center at e=10. We therefore do not generalize the thermal bound to every spectrum or declare the entire cancellation idea impossible.

The integration uses a lognormal energy-density distribution over ten standard deviations on each side, with the negligible tails truncated for quadrature, and resolves the oscillator resonance explicitly. The result is a low-speed drift coefficient; random momentum kicks and the receiver's own emission are not a complete thermal/kinetic evolution here.

## What this changes in the research plan

We have moved from an arbitrarily assigned capture function to an explicit passive scalar oscillator response with input, receiving and reradiated energy checked. It can approach the desired spectral shape, but its usefulness depends on the companion spectrum and its capture/storage rates. Its receiving continuum is not yet a long-lived, finite-capacity deposited state.

The next joined calculation should use the spectrum generated by the conversion mechanism, propagate it to the receiver, and apply a physically specified receiving response and all exit channels. It must not combine independent favorable limits while omitting their shared field dynamics or radiation spectrum. Receiver capacity, warming, orbital support and the simultaneous gravity/lensing prediction remain required. The original 20 tasks and 32 observational areas are unfinished.

Reproduce with `python -X utf8 research_work/results/broadband-capture/check_broadband_capture.py`. Saved evidence is in `broadband-capture-results.json`; the full runner includes this diagnostic. Mathematical agreement is not astronomical validation.
