# Radiation-powered initiation and persistent redshift

## Result

The model now starts without an externally supplied companion seed. It transfers net radiation energy into traveling sectors and deposits. However, this is predominantly **whole-photon removal**, a dimming process; it is not predominantly the energy lost by redshifting surviving photons. The tested capture cases produce a startup redshift, followed by a small blueshift once the field settles.

This establishes a conditional initiation and energy-supply mechanism, not the required persistent distance-redshift law. The same calculation reports the surviving-light timing, wavelength and photon-count predictions separately.

## Explicit extension

**Project postulate, not claimed unique:** allow a whole-photon conversion channel I=epsilon R. R is background photon reference-energy density. Negligible GW probes share the ray law but do not supply this photon-seeding source. T is traveling companion energy, M is traveling temporal-sector energy, and D is deposits. Retain n=1+M/K, b=R/(Kn), d=1-b, v=1/n and Q=T-M in response-time units:

\[
\begin{aligned}
R_t+(vR)_x&=-bM_t-I,\\
T_t+(vT)_x&=-Q-\Gamma T+I,\\
dM_t+(\sigma vM)_x&=Q,\\
D_t&=\Gamma T.
\end{aligned}
\]

Here K=2 and sigma=0.5 in nondimensional units c0=1. M's energy flux is explicitly slower; photons, GWs and T retain local c0. This is neither an identified ordinary-graviton process nor a completed spin/momentum theory. The conversion rate epsilon is per reference time and is **not** the previously fitted alpha per distance.

I removes photons with their energy into T. It does not change the frequency of each survivor. The field's changing n separately changes survivor frequencies through the shared Hamiltonian omega=c0|k|/n. This distinction is essential to both the spectrum and energy budget.

**Conditional conservation algebra:** the total R+T+M+D has flux v(R+T)+sigma vM and no unaccounted source. The radiation ledger includes two distinct transfers:

\[
W=\int bM_t\,dx\,dt,\quad A=\int I\,dx\,dt,
\quad R_{\rm final}+R_{\rm out}+W+A=R_{\rm initial}+R_{\rm in}.
\]

The sum is a reference-energy identity, not a demonstrated physical stress-energy conservation law. Positive W+A means net radiation energy supplies the other sectors.

## Why a settled time field loses accumulated stretching

The stipulated clock rule is d tau_clock=dt/n. The known ray and arrival-Jacobian identities, applied conditionally to that rule, give

\[
\frac{dt}{dx}=n,\quad\frac{d\ln J}{dx}=n_t,
\quad S=1+z=J\frac{n_e}{n_o}.
\]

Along the ray, d ln n/dx=n_x/n+n_t. Subtracting its endpoint integral gives another **conditional derivation using ordinary calculus, not a novelty claim**:

\[
\boxed{\ln S=-\int_{x_e}^{x_o}\frac{n_x(x,t(x))}{n(x,t(x))}\,dx.}
\]

n_x means the partial spatial derivative at fixed reference time. This is not simply an endpoint integral when n evolves along the ray. It makes two limits transparent:

- Spatially uniform n(t): S=1 because all measuring clocks change together.
- Stationary n(x): J=1 and S=n_e/n_o. Interior structure alone gives no accumulated stretch when endpoints have equal n.

A startup can produce a measurable redshift through evolving, nonuniform n. But a stationary maintained energy flow is not the same as a continually evolving time field. Sustaining conversion and capture does not automatically sustain the temporal evolution needed for cumulative redshift. No general exclusion of nonexpanding models follows from these conditional identities.

## Declared spatial test

The domain is x=0..10, run to time 40. Initially R=0.22 and T=M=D=0. Constant incoming photon reference-energy flux is 0.22; no companion or temporal energy enters at the boundary. Two capture regions use Gamma=strength*[exp(-x squared)+exp(-(x-10) squared)], with strength zero or ten. These are imposed capture locations, not an inferred gravity map or a forced reset of clock readings.

This is a one-dimensional one-direction test, not an isotropic cosmic radiation calculation. All parameters and three probe emission times were declared before computation. The length and times have not been assigned Mpc or years, and no coefficient was fitted to observations.

At 320 cells, probes travel from x=0 to10:

| Case | Early measured z, emission 0.1 | Late measured z, emission 20 | Early surviving photon fraction |
|---|---:|---:|---:|
| No initiation | 0 | 0 | 1 |
| epsilon=0.01, capture | +0.00107293 | -0.00085772 | 0.904657 |
| epsilon=0.05, capture | +0.00373003 | -0.00342738 | 0.603870 |

The no-capture initiation control blueshifts at all tested emission times. The capture cases show positive early shifts but do not preserve them into the late regime. Their late coordinate Jacobians are within approximately 6e-10 of one; the residual measured shift comes from different endpoint clock states.

**Conditional number-loss prediction:** photon survival is exp[-epsilon*(arrival-emission)]. This follows from the proposed achromatic removal rate; it is not the spectral factor 1/S. The early lower-rate case loses about 9.53% of its photons while survivors redshift by about 0.107%. The higher-rate case loses about 39.61% while survivors redshift by about 0.373%. These are model outputs, not observed dimming estimates or a fit to an astronomical sample. Surviving waveform intervals follow S; removal adds attenuation rather than the missing interval stretch.

## Energy outcome

| Fine-grid quantity at time 40 | epsilon=0.01, capture | epsilon=0.05, capture |
|---|---:|---:|
| Initial T+M seed | 0 | 0 |
| Net radiation transfer W+A | 0.84674844 | 3.56997825 |
| Whole-photon conversion A | 0.84218223 | 3.55258237 |
| Temporal work W | 0.00456622 | 0.01739588 |
| Final deposits D | 0.71662334 | 3.01846019 |

The deposits are radiation-funded in this reference ledger: there is no supplied seed. Only about 0.54% and 0.49% of net radiation transfer, respectively, is the temporal-work term. That integral includes the entire background's history, not only the early probe; it must not be equated to the energy loss of a selected redshifted signal.

M and T are transported and continually replenished. A nonzero stationary density under continuing inflow does not mean individual energy is permanently trapped in the void; it is an open-flow balance. Permanent D accumulates only where the chosen capture rule acts. No deposit-to-force or lensing response has been established here.

## Evidence and remaining work

[run.py](run.py) regenerates [results.json](results.json) under the [declared protocol](protocol.md):

```text
python research_work/results/radiation-powered-initiation/run.py
```

All twelve 80/160/320-cell runs preserve nonnegative reservoirs and the checked propagation domain. Total and radiation ledger errors are below 7.4e-13. Independent finite-event proper-clock integrations agree with the spectral factors within 8.44e-6. All refinement gates pass; the largest 160-to-320 absolute-z change is 4.83e-5 against a declared 0.003 threshold. These are numerical diagnostics, not observational uncertainties or universal stability proofs.

Next derive a common long-term temporal-field history, rather than a transient chosen for each emitter. It must sustain the desired redshift, predict its change with observing epoch, and account for the accompanying photon loss and energy supply. Increasing epsilon alone increases dimming and does not prevent settling in this tested family. Alternatively, a different per-photon conversion/clock interaction requires a new derivation of both frequency and event timing. Existing alpha, astronomical data and holdouts remain unchanged; all nine goals remain open.
