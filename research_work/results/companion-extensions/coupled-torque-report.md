# Coupled companion redistribution: circular endpoints

13 September 2026. Conditional calculation, not a derived formation mechanism.

## Result in plain language

An inner population can move inward while passing orbital angular momentum to outer companions. The outer population moves outward; the rest adjusts to the changed gravity. All six declared cases converge and leave positive energy available to escape. Thus letting the receivers react does not destroy this possible route to settling. It also changes the rotation predictions: a receiver population overlapping the measured outer disk can worsen the fit.

This is like moving some material inward in an orbiting system by letting other material carry its orbital motion outward. It does not mean a perfectly spherical gravitational field supplies a torque: an additional interaction or nonspherical disturbance must mediate the exchange. No such interaction, rate, or preferred final radius has been derived here.

## Frozen assumptions and formula provenance

The original exact-one-third retention prescription and accumulated inventory remain fixed. The previous shared fit supplies effective mass 17.7827941 eV and contraction s=0.803193303. Its phase-selection fraction is used only as a phenomenological donor label. The failed ideal-condensate support audit is not reversed by this calculation. Neither parameter is refitted to the stellar speeds.

**Known Newtonian circular-orbit mechanics**, applied to the spherical average of ordinary matter plus a spherical companion reservoir:

\[
j_c^2(r)=r^3g(r)=G r[M_b(<r)+M_c(<r)].
\]

Here j is specific angular momentum. The spherical averaging is an approximation to the ordinary disk and bulge; it is not the actual three-dimensional orbital solution. General orbital-support context is given in [Bovy's Jeans-equation treatment](https://galaxiesbook.org/chapters/I-04.-Equilibria-of-Collisionless-Stellar-Systems_4-The-Jeans-equations.html).

**Our imposed endpoint prescription, not a claimed new fundamental law:** donors D finish at r_f=s r_i. Their scalar angular-momentum decrement is

\[
\Delta J_D=\sum_{i\in D}m_i[j_i-j_c(s r_i)].
\]

Receivers R are the initially extended companions in each declared radial band. Every receiver gains the same specific increment:

\[
j_{f,k}=j_{i,k}+\Delta J_D/M_R\quad(k\in R),\qquad
j_{f,k}=j_{i,k}\quad(k\notin D\cup R).
\]

We solve j_c(r_f)=j_f while recomputing the complete companion gravity. Matched orientation families could balance vector angular momentum while having zero total spin; the computed scalar allocation assumes such matching. It does not establish a realizable interaction between randomly oriented populations.

**Known Newtonian energy accounting:**

\[
E=\tfrac12\sum_i m_i(j_i/r_i)^2+\sum_i m_i\Phi_b(r_i)
-G\sum_i\frac{m_i[M_c(<r_i)+m_i/2]}{r_i},\qquad Q=E_i-E_f.
\]

In the shell sum, enclosed mass excludes the current shell. The half-self factor prevents double counting. Positive Q is the energy that must leave during the transition; it is not a predicted luminosity or a demonstrated cooling channel. Fixed ordinary matter is an external potential. Companion rest mass is held fixed in this Newtonian approximation; a relativistic completion must also include the gravitational effect of emitted binding energy.

## All declared outcomes

Errors below are RMS differences in km/s against the 38 previously used Milky Way circular-speed estimates. These are reused observations, not blind predictions or individual-star velocities. Baselines I and II are alternative ordinary-matter models.

| Baseline | Initial receiver band (kpc) | Receiver mean radius increase | Q (10^50 J) | Inner RMS | Outer RMS | All RMS |
|---|---:|---:|---:|---:|---:|---:|
| I | 15–30 | 5.41% | 3.733 | 3.47 | 8.35 | 6.27 |
| I | 30–60 | 5.59% | 4.471 | 3.47 | 8.36 | 6.28 |
| I | 60–120 | 7.16% | 4.799 | 3.47 | 8.36 | 6.28 |
| II | 15–30 | 6.27% | 3.654 | 6.23 | 11.75 | 9.26 |
| II | 30–60 | 6.44% | 4.432 | 6.23 | 9.78 | 8.10 |
| II | 60–120 | 8.22% | 4.780 | 6.23 | 9.78 | 8.10 |

Matched no-settling all-bin RMS is 6.34/10.91 for I/II. Holding the extended reservoir fixed in the prior imposed contraction gave roughly 6.09/8.30. Allowing a physical orbital response modestly degrades that imposed result for I and improves II when receivers lie beyond the measured disk. No receiver band wins universally, and none is selected by these errors. The nearly identical outer-band predictions arise because their redistribution is outside the measured radii in the spherical approximation.

## Numerical checks and limits

The initial thin-shell crossing representation did not converge: it produced force jumps, with angular-momentum residuals around 0.1%. Those attempts are superseded, not physical failures. The final solver deposits positive packet masses continuously onto neighboring logarithmic grid nodes. Ordinary-matter forces are evaluated directly, avoiding subtraction artifacts below the legacy radial clamp.

The calculation uses 8,192 source bins and 16,384 force nodes, then doubles both. All six cases and displaced initial guesses converge. Main-run transfer balances are within 6e-8 of the transferred angular momentum; unforced specific angular momenta change by less than 9e-9. Zero contraction returns the original radii and energy. Refinement changes any predicted speed by less than 0.084 km/s and released energy by less than 0.039%. Direct circular kinetic energy and the independent shell virial estimate change Q by less than 7e-7 fractionally. These are numerical consistency checks, not observational error bars. The algebraic energy-ledger residual alone would not prove physical conservation.

## Consequence for the theory

Internal angular-momentum exchange plus separate energy release remains feasible at the endpoint level. It avoids requiring the same outgoing rays to carry the entire orbital angular momentum, which was expensive in the prior direct-ray test. It does not resolve capture from light-speed propagation, collective stability, the failed thermal interpretation, the redshift/time mechanism, lensing, photon supply, or a universal stopping law.

The next physical step is to specify an exchange interaction with a calculable rate and show that its evolution approaches a supported profile without imposing s. A dissipative coupling, collective wave torque, and nonspherical gravitational transport should remain separate candidates until their energy and angular-momentum transfers are derived. The receiver bands here are diagnostic choices, not a discovered environmental law.

Reproduce with `python research_work/results/companion-extensions/coupled-torque.py`. See [protocol](coupled-torque-protocol.md), [executable](coupled-torque.py), and [all results](coupled-torque-results.json).
