# Shared pressure laws: equilibrium does not preserve the fitted speed improvement

13 September 2026. Goal continuation on physical support, with exact-third captured inventory retained. This follows the local ideal-Bose pressure failure.

## Outcome

We tested three shared polytropic pressure laws, calibrated to the required inner pressure of the earlier fitted profiles, then allowed each law to determine its own spherical equilibrium. Four finite mass-matched configurations were found for the 5/3 and 2 exponents. Their predicted stellar speeds are worse than the original capture and fitted-settling controls. The 4/3 exponent's mass-matching candidates failed to reach a zero-pressure surface within the computational domain and were rejected; they are not valid finite equilibrium detections.

The result is a failed test of these pressure calibrations and closures, not an exclusion of all interacting companions. It demonstrates why adding enough pressure to a prescribed profile is different from deriving a profile that both supports itself and matches observations.

## Known equations and proposed interpretation

The trial equation of state is

\[
P=\rho_{\rm ref}v_0^2(\rho/\rho_{\rm ref})^\gamma,
\quad \rho_{\rm ref}=10^7\ M_\odot/\mathrm{kpc}^3,
\quad\gamma\in\{4/3,5/3,2\}.
\]

These are known polytropic forms. Using them as companion support is proposed. The quadratic form is motivated by the contact-interaction condensate equation of state, not invented here; see [Boehmer and Harko](https://arxiv.org/abs/0705.4158). We do not adopt that paper's cosmological origin model. Nor does choosing gamma=4/3 identify our reservoir as radiation: this calculation still treats its density as nonrelativistic gravitating mass.

One v0 per exponent was fitted by mean logarithmic pressure error across both matter baselines at 100 log-spaced radii from 0.1 to 10 kpc. The fit targets are already-derived required pressure, not new observed pressure measurements. The 10-30 kpc interval was excluded from pressure calibration. All these profiles and speed data have previously been inspected.

The total-pressure fits yield v0=187.76, 155.90 and 129.44 km/s respectively. We also fit the positive extra-pressure residual after subtracting the ideal thermal pressure implied by the assigned compact fraction; all results are retained. The equilibrium calculation uses the total effective pressure law. A residual-pressure fit alone cannot define an equilibrium until the thermal and interaction components are coupled consistently.

To determine a stopping radius, we solve the established hydrostatic equations in enthalpy form:

\[
H=\frac{\gamma}{\gamma-1}v_0^2(\rho/\rho_{\rm ref})^{\gamma-1},
\quad \frac{dH}{dr}=-g_b(r)-\frac{GM(r)}{r^2},
\quad \frac{dM}{dr}=4\pi r^2\rho(H).
\]

The surface is H=0. Central density is adjusted only to match the fixed exact-third total mass, not observed stellar speeds. Both the pressure constants and the capture inventory are shared or frozen as declared. The ordinary-matter radial force is spherically averaged for the equilibrium, while the observable circular-speed baseline retains the earlier flattened stellar force and gas monopole. Thus this is not an exact equilibrium in the fully flattened Galaxy.

## Results

RMS errors are km/s over the archived 38 circular-speed bins. I and II are alternative ordinary-matter descriptions of the same Galaxy.

| Model | Equilibrium radius I / II | Inner RMS I / II | Outer RMS I / II | All-bin RMS I / II |
|---|---|---|---|---|
| Original exact-third control | Not a polytropic surface | 3.68 / 11.83 | 8.36 / 9.78 | 6.34 / 10.91 |
| Earlier fitted shared settling | Fitted redistribution | 2.74 / 6.70 | 8.36 / 9.78 | 6.09 / 8.30 |
| Gamma=5/3 equilibrium | 29.14 / 30.48 kpc | 14.81 / 10.05 | 35.88 / 26.24 | 26.93 / 19.47 |
| Gamma=2 equilibrium | 19.64 / 20.05 kpc | 34.59 / 25.26 | 52.04 / 44.99 | 43.73 / 35.98 |

The finite equilibria are too concentrated relative to the useful empirical profile. Retaining the same total mass but putting it inside a smaller radius increases the force over much of the observed range. No speed-based choice among solutions was made.

For gamma=4/3, the central-density scan bracketed mass matches that, upon solving, reached 100,000 kpc without a zero-pressure surface. Those candidates are explicitly stored as rejected. A mass match inside a numerical box is not a demonstrated finite reservoir. The bounded search does not prove that all other boundary conditions, amplitudes or densities fail.

## Why an everywhere increasing pressure-density law cannot reproduce the original shape

Hydrostatic pressure must fall outward wherever density and inward gravity are positive. If pressure depends only on density with dP/d(rho)>0, density must also fall outward. The fitted settled profile contains outward-rising density intervals: approximately 20.3%/0.49% of logarithmically sampled radii between 0.1 and 30 kpc in baselines I/II. These are sampling fractions, not mass fractions; they are not observationally measured shells.

That sign mismatch is independent of the fitted amplitude. It identifies a structural problem for a single positive-compressibility barotropic fluid reproducing the exact prior profile. It does not prohibit entropy gradients, multiple coupled components, anisotropic stress, nonlocal wave stress or time-dependent states. Those alternatives must supply their own equations and conservation accounting.

## Creative next steps with concrete tests

- Orbital support: invert the desired density and potential into a candidate distribution of orbital energies, then check that the distribution is nonnegative. A positive pressure alone is insufficient. Explore anisotropy only with explicit shared rules, not an arbitrary correction at every radius.
- Multiple states: let a broadly distributed population support itself while a smaller bound population carries the central correction. Solve each component's force balance and exchanges rather than assigning all pressure to one fluid.
- Temperature/history-dependent pressure: permit pressure to depend on entropy or stored internal energy as well as density. The required gradients must arise from a transport/cooling equation, not be inserted to reproduce the chosen profile.
- Wave stresses: derive finite-gradient or interaction terms from an energy functional. Check whether the resulting equilibrium and stability preserve the speed fit; the prior Gaussian-wave failures remain relevant controls.

The next distinguishing calculation is an orbital-distribution positivity test. It addresses the support failure without adding an arbitrary pressure term and gives a clear rejection criterion before fitting more data. The propagation, energy-supply and radiation tracks remain active requirements of the full goal.

## Verification and reproducibility

[Protocol](interaction-support-protocol.md), [code](interaction-support.py) and `interaction-support-results.json` retain all amplitudes, pressure errors, 48-point central-density scans for each baseline/exponent, rejected candidates and accepted speed predictions. The integrator reproduces the analytic no-external-field n=1 Lane-Emden radius and mass within 1e-6. Accepted mass roots agree within 1e-6 relative tolerance; tighter integration changes predicted speeds by under 0.00014 km/s. Runs pass with warnings treated as errors. These checks validate the numerical calculation within its assumptions, not dynamical stability.

No local interaction has been derived from photon conversion, and these equilibrium cases do not establish formation, energy supply, cooling rates or lensing. The complete research goal remains open. The prior turn made concrete progress by identifying a pressure deficit; this turn advances that finding to explicit equilibrium solutions and shows that the simplest tested pressure laws do not retain the empirical success.
