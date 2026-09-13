# A resonance ensemble can flatten the conditional redshift rate

13 September 2026. Constructive effective-response calculation; no completed microscopic or observational claim.

## Result

The previous finite-size interaction left a fractional loss rate approximately proportional to photon energy. An ensemble of independent resonant receiving systems with weight proportional to squared resonance frequency can cancel that dependence over an interior band. In the broad-cutoff examples the normalized fractional rate varies by at most 0.0063%, 0.0630% or 0.1891% across 0.01–100 eV, for damping ratios 0.1, 1 or 3.

This is evidence that a specified energy-dependent response can meet the color-shape target approximately. It is not evidence that companions have these resonances, that the absolute redshift is reproduced, or that the entire electromagnetic spectrum and event timings are explained. The distribution was chosen by inverse scaling analysis, not predicted from the deposition law.

## Formula provenance

**Known damped-oscillator response:** in units hbar=c=1,

\[
\beta(E;\Omega)=\frac{B}{\Omega^2-E^2-i\eta\Omega E},\qquad \eta>0.
\]

The oscillator equation gives this frequency response, including its dispersive and dissipative parts. See [MIT's Lorentz-oscillator notes](https://ocw.mit.edu/courses/6-007-electromagnetic-energy-from-motors-to-lasers-spring-2011/resources/mit6_007s11_lorentz/). Applying this form to the transition polarizability of hypothetical companions is our new assumption. It is not a claim that ordinary atoms or their known polarizabilities are the companion stores.

**Proposed ensemble:** independently scattering receiving systems have a positive weighted resonance distribution dN/dOmega proportional to Omega^p between finite cutoffs. Assume a common numerator B, a common damping ratio eta, and the same small transfer gap and Gaussian spatial overlap as in the preceding calculation. Their incoherent rates add, giving the response factor

\[
R_p(E)\propto\int_{\Omega_{\min}}^{\Omega_{\max}}
\frac{\Omega^p\,d\Omega}{(\Omega^2-E^2)^2+\eta^2\Omega^2E^2}.
\]

This is a sum of squared amplitudes for distinct independent targets, not the squared coherent sum of amplitudes. A coherent medium needs a different calculation, including its forward response.

## Why p=2 works in the scale-separated limit

**Scaling consequence derived here, not a claim of new mathematics:** setting Omega=E*u gives R_p(E) proportional to E^(p-3) times an integral with limits Omega_min/E and Omega_max/E. Far from both cutoffs, for p=2,

\[
R_2(E)\propto\frac1E\int_0^\infty
\frac{u^2\,du}{(u^2-1)^2+\eta^2u^2}
=\frac{\pi}{2\eta E}.
\]

The integral follows by substituting u->1/u in a second copy, adding them, and using t=u-1/u to obtain the integral of 1/(t^2+eta^2). The previous large-size, small-gap result alpha_spatial proportional to E therefore gives

\[
\alpha_{\rm combined}(E)\propto\alpha_{\rm spatial}(E)R_2(E)
\simeq\text{constant}.
\]

The oscillator shape is established mathematics; the weighted population, coupling assumptions and application to inelastic photon-companion transfer are hypotheses. The exponent two was selected to supply the missing inverse energy factor. No novelty claim or first-principles derivation of this population is made.

The infinite limits are used only for the analytic scaling. An actual Omega^2 population extending to arbitrarily high frequencies has divergent total weight and is not adopted. Finite cutoffs, mode abundance and oscillator strengths must be physically supplied. The infinity-limit integral converges for -1<p<3; the p=3 control is explicitly finite-cutoff dependent.

## Numerical comparison

Keep gap Delta=1e-8 eV and Gaussian size a=1e6 inverse eV (approximately 0.197 m). Test p=1,2,3; eta=0.1,1,3; and two cutoff pairs. Each of the 18 cases is normalized at 1 eV and evaluated at 17 logarithmically spaced photon energies between 0.01 and 100 eV. This gives 306 rate-shape evaluations, not 306 observational tests.

| Resonance cutoffs (eV) | Damping ratio | Maximum departure from normalized constant rate, p=2 |
|---|---:|---:|
| 1e-5–1e5 | 0.1 | 0.0063% |
| 1e-5–1e5 | 1 | 0.0630% |
| 1e-5–1e5 | 3 | 0.1891% |
| 1e-3–1e3 | 0.1 | 0.6346% |
| 1e-3–1e3 | 1 | 6.3278% |
| 1e-3–1e3 | 3 | 18.5148% |

The narrower ensemble exposes the cutoff dependence. The other exponents are retained in the output as controls; they do not generically cancel the residual spatial energy scaling. No exponent is fitted separately at each photon energy.

Tightening adaptive integration from 1e-10 to 1e-12 changes every computed response by less than 1e-8 fractionally. The eta=2, p=2 integral agrees with pi/4 to the expected finite-tail accuracy. These checks support the numerical response calculation, not the hypothesized population's existence.

## What this does and does not repair

This escapes the preceding positive-gap-spectrum bound because the squared interaction amplitude now depends on incident energy. It does not contradict that bound or show that a continuum of fixed gaps alone was sufficient. Resonance energies Omega and the inelastic creation gap Delta are distinct assumptions here.

The existing size, gap and absolute-rate costs remain. At a=0.197 m the preceding fixed-count angular diagnostic gave roughly 127 arcseconds of accumulated direction RMS. The resonance factor alone does not change that angular shape in this separable model. This is not an actual image-width prediction, but it prevents calling the present parameter set a complete optical solution. A larger size would require reevaluating the rate and coupling budget.

Damping is not free. Positive eta represents energy leaving oscillator motion. Its destination must be identified as radiation, traveling companions, stored excitation or another accounted sector. Inserting a complex response only into a selected inelastic rate does not supply that complete energy ledger. If radiative damping is involved, the linewidth and coupling cannot be chosen independently without checking scattering and extinction constraints. The same interaction must also define inverse transitions, forward dispersion and occupation evolution. A causal oscillator denominator by itself does not prove that the proposed combined scattering model meets all those constraints.

No absolute normalization, galaxy mode count, source supply, redshift-distance fit, radio/gamma extrapolation, supernova duration stretch, graviton identity, lifetime, retention exponent or lensing prediction follows from this calculation.

## Decision

Retain the resonance ensemble as a constructive candidate for the direct-transfer branch, with its inverse-designed spectrum labeled. The next decisive step is an energy-accounted scattering response: link the damping to actual outgoing channels and compute extinction, inelastic transfer and dispersion using common couplings. If those cannot coexist with the required rate and image preservation, the apparent color success must not be treated as a physical solution. Keep the coherent time-field and gravitational-deposition branches separate until joined by a consistent interaction.

Reproduce with `python research_work/results/companion-extensions/resonant-spectrum.py`. [Source](resonant-spectrum.py), [complete results](resonant-spectrum-results.json).
