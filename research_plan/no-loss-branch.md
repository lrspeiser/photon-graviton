# Optional branch: energy-preserving companions and permanent deposits

User-suggested candidate, recorded 9 September 2026. Energy conservation is mandatory throughout the research program. Zero companion loss and permanent storage are optional physical assumptions, not necessary consequences of conservation and not yet derived properties of gravitons.

## 1. Precise statement of the branch

In a fixed physical volume, omitting spatial flux and pressure work for this illustrative calculation, assume

\[
\dot u_\gamma=j_\star-h_{\rm loss}u_\gamma,\qquad
\dot u_c=h_{\rm loss}u_\gamma-\Gamma_{\rm cap}u_c,\qquad
\dot u_d=\Gamma_{\rm cap}u_c.
\]

Here u are energy densities; j_star is injected energy per volume per time; h_loss and Gamma_cap are nonnegative inverse times. h_loss is a photon fractional-energy-loss rate, **not an expansion rate**. This branch sets both companion secondary loss lambda_c and deposit leakage/decay lambda_d to zero. It also assumes every bit of the specified photon loss enters the companion reservoir.

Adding the equations gives

\[
\frac{d}{dt}(u_\gamma+u_c+u_d)=j_\star.
\]

If the stellar fuel reservoir satisfies u_star-dot=-j_star, their combined energy is constant. With radiative sources outside the volume, injection must instead enter through the boundary flux. In a spatially resolved model add flux divergences; for changing volume/metric include pressure work and covariant exchanges. Include any field doing work on photon frequencies. Its contribution cannot be omitted just because the three-reservoir sum balances.

This proves consistency of the stated bookkeeping. It does not yet prove momentum conservation, a microscopic transfer mechanism, a stable background, or a physically permanent deposit.

## 2. Exact impulse solution

Let h=h_loss>0 and Gamma=Gamma_cap>=0 be constant. Initially place U in photons, with u_c=u_d=0 and j_star=0 thereafter. Then

\[
u_\gamma(t)=Ue^{-ht}.
\]

For Gamma unequal to h,

\[
u_c(t)=\frac{hU}{\Gamma-h}(e^{-ht}-e^{-\Gamma t}),\qquad
u_d(t)=U-u_\gamma(t)-u_c(t).
\]

At Gamma=h the removable singularity has limit

\[
u_c(t)=Uh t e^{-ht},\qquad
u_d(t)=U[1-(1+ht)e^{-ht}].
\]

All three densities are nonnegative. For Gamma>0 all the initial energy eventually reaches deposits in this model, but never more than U. For Gamma=0, u_c=U(1-e^{-ht}) and u_d=0: companions accumulate freely and are never deposited. For h=0 no photon-to-companion conversion occurs.

This differs from the archived shared-secondary-loss calculation: there is no mandatory free-companion peak of U/e and no timing-sector residual in this idealized branch. That is a change of equations, not a violation of energy conservation. It must have a compatible action before it becomes a physical model.

## 3. Continuous injection

For constant j, h>0 and Gamma>0, starting with empty reservoirs, the photon solution is u_gamma=(j/h)(1-e^{-ht}). For Gamma unequal to h,

\[
u_c=j\left[\frac{1-e^{-\Gamma t}}{\Gamma}-\frac{e^{-ht}-e^{-\Gamma t}}{\Gamma-h}\right],\qquad
u_d=jt-u_\gamma-u_c.
\]

For Gamma=h,

\[
u_c=\frac{j}{h}(1-e^{-ht})-jt e^{-ht}.
\]

At late times u_gamma tends to j/h and u_c tends to j/Gamma, while u_d grows approximately jt-j/h-j/Gamma. Thus the free reservoirs can approach steady values, but a permanently growing deposit cannot have a finite stationary density under continuing positive capture.

With finite stellar fuel, time-varying sources or a capture mechanism that switches off, deposition can stop. Those are additional physical histories to calculate. A changing self-gravitating potential could also make a static-volume approximation fail. None is equivalent to quietly adding deposit decay while claiming lambda_d=0.

## 4. What must be recalculated

| Question | Necessary change from the shared-loss checkpoint |
|---|---|
| Source-to-target companion flux | Reintegrate emission histories, geometry, retarded times and survival. Removing one attenuation factor from a final fitted formula is not a general derivation. |
| Gravitational source before capture | Free companions already carry stress-energy; include it rather than making gravity switch on at deposition. |
| Halo assembly | Integrate capture over the full history, with changing sources, potential and finite fuel. |
| Long-term observations | Predict dependence on age and illumination, halo growth, cluster growth and total cosmic stored energy. |
| Tensor/photon propagation | Derive why companion energy is constant in the chosen physical frame while photon energy decreases. If both obey the same time-dependent omega=c k/n for conserved k in that frame, this distinction does not follow. |
| Observed gravitational waves | If companions are ordinary gravitons, derive their relationship to measured GW frequency/amplitude/timing. If they are a different field, state and constrain that explicitly. |
| Thermal and entropy budgets | Specify whether capture heats matter, creates a cold field, traps radiation or transfers to another state. |

No-loss propagation can increase arriving companion energy relative to the same source under secondary losses. It does not exceed the emitted photon energy. Therefore the original 10-billion-year, present-luminosity, local 100%-conversion benchmark still has its approximately 5,428 median energy shortfall for the 26-object sample. A different external source inventory or much older formation history is a distinct budget to calculate; the quoted factor is not a universal bound on those alternatives.

Do not inherit an arbitrary enormous age from the previously fitted history. The recovered positive-n solution with p about 0.264 and gamma about 7.7315e-11/year has a finite past endpoint roughly 17.6 billion years before its reference epoch. A much older universe needs a new or extended physical history that still reproduces the observed low-redshift behavior and passes the clock/background tests. The endpoint is a limitation of that particular ansatz, not an independently measured cosmic age.

## 5. Decision and Codex outcome

Keep this branch alongside shared-secondary-loss and finite-residence alternatives. T20 first derives and tests their ledgers. T03–T05 must supply a compatible interaction and observable transfer law; T16 supplies an admissible history; T06 computes independent supply; T07–T08 calculate capture and growing deposits. Local, lensing and cosmic tests then decide viability.

Acceptance requires all of the following: a conserved full ledger, derived rather than declared propagation/capture differences, adequate independently supplied energy, a stable stored state, and observationally acceptable halo/background evolution. Failure of any one rejects this version or identifies an explicit premise to revise. The equation sum alone cannot establish acceptance.
