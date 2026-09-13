# What could stop companion settling?

13 September 2026. A constrained energy diagnostic, not a time-evolution simulation.

## Finding

The previous fitted contraction, s=0.803193303, is not an energy minimum of the coupled circular-endpoint family. In all six cases, energy falls as contraction deepens at every sampled step from s=1 to 0.4. An energetically downhill direction therefore remains at the fitted radius. This does not prove that an interaction can access that direction or that the reservoir must collapse: orbital isolation or a transport cutoff can prevent further evolution.

In plain language, making an arrangement that balances gravity is not the same as explaining why material stops there. The proposed system still needs either resistance to further compression, a shutoff of the exchange, or a limit on what the outer reservoir can receive.

## Calculation and provenance

The exact-third inventory, donor fractions, both ordinary-matter baselines and receiver bands remain frozen. All 12 contraction samples per case and three finer-resolution samples converge. The donor fraction retains its phenomenological meaning; the earlier ideal thermal-condensate failure still stands. The prior phase mass and contraction were fitted using these observations, so this is not a new prediction test.

**Known energy-extremum criterion, applied here:** along the constrained family, write the computed Newtonian mechanical energy as E(s). An interior stationary point would require E'(s)=0. Local stability against this one coordinate would additionally require E''(s)>0. Neither condition tests arbitrary perturbations or collective stability.

| Ordinary matter | Receiver band (kpc) | E'(s*) (10^51 J) | E''(s*) (10^51 J) |
|---|---:|---:|---:|
| I | 15–30 | 2.299 | -4.853 |
| I | 30–60 | 2.665 | -4.778 |
| I | 60–120 | 2.828 | -4.749 |
| II | 15–30 | 2.266 | -4.919 |
| II | 30–60 | 2.651 | -4.841 |
| II | 60–120 | 2.824 | -4.813 |

s is dimensionless, so both derivatives have energy units. Halving the finite-difference step from 0.02 to 0.01 changes the slope by less than 0.032%; doubling spatial resolution changes it by less than 0.035%. Positive slope at the reference and the sampled monotonic trend are robust at these resolutions. Sampling does not prove global monotonicity at every intermediate radius or beyond the tested interval.

More contraction is not automatically a better speed fit. At s=0.4, all-bin RMS is 9.84–10.42 km/s for I and 9.70–12.84 for II, compared with 6.27–6.28 and 8.10–9.26 at the reference. No contraction is selected anew using these errors.

## Alternative 1: restoring energy

**Hypothetical effective addition using a known power-law form:**

\[
E_{\rm total}(s)=E(s)+A s^{-p},\qquad p\in\{1,2,3\}.
\]

**Algebraic inverse calibration in this work, not a new fundamental law:** to place a stationary point at the prior fitted s*,

\[
A=\frac{E'(s_*)s_*^{p+1}}p,\qquad
E_{\rm total}''(s_*)=E''(s_*)+\frac{(p+1)E'(s_*)}{s_*}.
\]

All 18 calibrated cases give positive A and positive local curvature; this remains true at finer resolution. Choosing A this way builds the desired stop into the calculation. It has not derived the stop from an independently measured interaction.

The energy increase of this restoring component from the uncontracted state must be paid for:

\[
\Delta E_{\rm support}=A(s_*^{-p}-1),\qquad
Q_{\rm remaining}=Q_{\rm mechanical}-\Delta E_{\rm support}.
\]

| p | Added stored energy (10^50 J), across six cases | Remaining release (10^50 J) |
|---|---:|---:|
| 1 | 3.58–4.47 | 0.072–0.328 |
| 2 | 3.23–4.03 | 0.425–0.768 |
| 3 | 2.92–3.65 | 0.731–1.150 |

These are conditional bookkeeping values along the original endpoint family. A real pressure or interaction changes the local force balance and orbital support; the current speed curves cannot simply be retained after adding it. Most of the previous release budget can become stored support energy, so the previous luminosity estimates would change.

**Known scaling, not unique to us:** for a fixed-mass homologously compressed polytropic component with P=K*rho^gamma, internal energy scales as s^[-3(gamma-1)]. Thus p=1,2,3 corresponds to gamma=4/3,5/3,2 under those assumptions. This motivates examples, not their microscopic identification. Repulsive-condensate polytropic models already exist in [Boehmer and Harko](https://arxiv.org/abs/0705.4158). Our earlier full equilibrium tests of simple shared polytropes failed to preserve the successful profiles; this inverse one-coordinate exercise does not overturn that result.

## Alternative 2: shut off transport when states fill

**New hypothesis for this project, using familiar gradient-flow and occupancy mathematics:**

\[
\dot s=-\mu(s)E'(s),\qquad
\mu(s)=\mu_0\max[0,1-\rho_D(s)/\rho_{\rm sat}].
\]

Here mu_0 has units (energy*time)^-1, rho_D is a chosen donor-density measure, and rho_sat is an unspecified capacity. It arrests the allowed coordinate even when lowering energy remains possible. It does not by itself supply radial support; that would still come from orbits or another stress. A discontinuous cutoff, leakage and stability against perturbations would need separate evaluation.

**Known volume scaling:** for a homologous donor parcel, rho_D(s)=rho_D(1)/s^3. Reproducing the old contraction would require rho_sat/rho_D(1)=1.930. This is an inverse target, not a measured density threshold. A single absolute saturation density would not generally contract an initially nonuniform galaxy by a uniform factor. Its cube root comes from volume geometry and is not a derivation of the one-third retention exponent.

## Alternative 3: exhaust the receiver capacity

Bounded angular-momentum acceptance could stop the exchange before an energy minimum is reached. The previous endpoints show finite outward displacement can absorb the required transfer; they do not establish any maximum capacity. A model needs a condition such as escape, loss of coupling, or saturated available states, derived from the receiving population. Arbitrarily imposing the already desired radius would add no explanation.

## Next step and limits

The least redundant next calculation is a density-dependent transport cutoff applied locally, allowing a nonuniform final profile and recalculating the motions. This tests the user's saturation concept without repeating the failed single-pressure-law assumption. Its threshold must be fixed on one declared calibration case and transferred to other cases; different ordinary-matter baselines of one galaxy are sensitivity tests, not independent galaxies.

No new photon supply, redshift, lensing, three-dimensional equilibrium, torque microphysics, age constraint, or universal law has been established here. The evidence is that continued settling is energetically accessible along this assumed family, and that several explicitly costed stopping mechanisms remain candidates.

Reproduce with `python research_work/results/companion-extensions/stopping.py`. [Protocol](stopping-protocol.md), [source](stopping.py), [complete outputs](stopping-results.json).
