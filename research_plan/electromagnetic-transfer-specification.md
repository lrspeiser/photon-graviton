# Electromagnetic redshift, signal timing and companion transfer

9 September 2026. This is the current working specification following the request to test both energy conversion and altered time across the electromagnetic spectrum. It supersedes the claim that a single cause has been selected and completed. The [direct-conversion branch](direct-conversion-bound-companions.md) remains the baseline comparison; a time-dependent propagation candidate is reopened. Neither is a completed theory. Published distances remain fixed under the [universe contract](universe-contract.md); expansion, Big-Bang initial conditions and inserted dark-matter halos are not explanatory premises.

The [executed audit](../research_work/results/electromagnetic-audit/report.md) supplies calculations, input hashes, measured-versus-predicted rows and remaining gaps. Synthetic frequencies are not observations. Earlier mathematical failures remain evidence about their stated assumptions, not prohibitions against every modified law.

## 1. What counts as an observed redshift

**Established observational definition.** For identified transition j with today's laboratory frequency nu_lab,j:

\[
1+z_j=\frac{\nu_{\mathrm{lab},j}}{\nu_{\mathrm{observed},j}}.
\]

An emitted line need not have precisely today's laboratory frequency if the proposed field changes atoms. Define r_j=nu_emit,j/nu_lab,j. Under separable propagation, motion and endpoint effects:

\[
1+z_j=\frac{S_E(\nu_j,\text{path})\,S_{\rm motion}\,S_{\rm endpoint}}{r_j}.
\]

**Conditional factorization using established frequency ratios.** Do not count one gravitational effect twice: S_E excludes any endpoint effect already in S_endpoint. A full ray/clock solution replaces the factorization when its assumptions fail. A line imprinted by intervening absorbing gas starts its relevant path at that gas, not at the background quasar. Distinct gas velocities, blends, calibration and extinction need their own likelihood. Redshift-free and blueshifted observations must be retained. Positive conversion can coexist with total blueshift if independently constrained motion or endpoint effects dominate; arbitrary fitted velocities cannot erase every residual.

## 2. The baseline energy rule and a testable color modification

**Proposed partial-transfer law, written as established rate equations:**

\[
\frac{dE_\gamma}{ds}=-\alpha(E_\gamma,s)E_\gamma,
\qquad \frac{dE_c}{ds}=+\alpha(E_\gamma,s)E_\gamma.
\]

The second equation is an isolated packet account before capture, assuming no driver work, recoil energy or companion loss. A mechanism must demonstrate those approximations. It describes surviving photons losing energy, not merely whole photons disappearing into another field.

**Empirical baseline, not first principles:** alpha_0=0.0002488993286382367 per Mpc, from the exposed 164-group fit. Define A=integral alpha ds. If alpha is independent of photon frequency:

\[
S_E=e^A,\quad 1+z_{\rm transfer}=e^A,\quad
E_{\gamma,f}=E_{\gamma,i}e^{-A},\quad
E_c=E_{\gamma,i}(1-e^{-A}).
\]

**Established exponential solution and E=h nu.** Not a new mathematical redshift formula. A common fractional rate is the economical broadband target: equal paths give equal fractional shifts, not equal absolute energy losses. It does not imply all astronomical sources have identical shifts.

**Proposed diagnostic family, not a chosen physical interaction:**

\[
\alpha(E)=\alpha_{\rm ref}(E/E_{\rm ref})^p,
\quad
\frac{E_f}{E_i}=\left[1+pA_{\rm ref}(E_i/E_{\rm ref})^p\right]^{-1/p}.
\]

**Established solution of that postulate**, for a positive bracket; p=0 means its exponential limit. p=-1 gives fixed absolute loss and can exhaust low-energy photons. The radio comparison disfavors the tested large nonzero p values under the common-path assumption. It does not establish p=0 over every unobserved frequency or derive its rate. A tiny conditional fitted radio p is not adopted: it could reflect gas/calibration effects, and fitting two centroids with two effective quantities predicts no independent new object.

## 3. Delay and event stretch must be calculated separately

For a fixed source and receiver in a specified common time coordinate, write a characteristic as dt/ds=a(s,t). If J=partial t(s)/partial t_emit:

\[
\frac{dJ}{ds}=\partial_t a\,J,\qquad J(0)=1.
\]

**Established characteristic-variation identity.** For endpoint clock factors q_emit and q_obs, the measured infinitesimal duration ratio is S_t=(q_obs/q_emit)J. The derivative is at fixed s and fixed physical field, not a different fitted environment for each photon. For stationary a and matching endpoint clocks, J=1 regardless of how much common delay accumulates. A proposed law for d ln S_t/ds alone is a target, not a source for a(s,t).

### Baseline C0: stationary partial energy conversion

**Conditional consequence:** a=1/c gives t_arr=t_emit+D/c and S_t=1. Lower carrier frequency is stipulated by an inelastic interaction, whose phase/coherence properties remain unspecified. This closure does not supply the supernova timing trend with unchanged source populations.

### Candidate T1: an explicit accumulating time field

**Proposed inverse-affine lapse; known mathematical structure, originality unverified:**

\[
q(s,t)=\frac{1}{1+c\kappa(s)(t-t_*)},\qquad
d\tau=q\,dt,\qquad
d\ell^2=ds^2.
\]

Local light speed remains c: d ell/d tau=c. The denominator must be positive. kappa is a physical field profile that must ultimately be defined throughout space, not independently fitted along each observer-source line. s here is the chosen static spatial coordinate along the diagnostic ray. q=1 at emitter and receiver is imposed in the tests by kappa=0 at both endpoints; it is not yet a demonstrated environmental protection mechanism.

**Conditional null-ray consequence:**

\[
\frac{dt}{ds}=\frac{1}{cq}=\frac1c+\kappa(s)(t-t_*),
\qquad
S_t=J=\exp\!\left[\int\kappa(s)\,ds\right].
\]

**Conditional geometric-optics consequence of the established photon Hamiltonian H_gamma=cq|p|:**

\[
\frac{d\ln H_\gamma}{ds}=\frac{\partial_tq}{cq^2}=-\kappa(s),
\qquad S_E=S_t=e^{\int\kappa ds}
\]

at the matching endpoints. Thus the same prescribed background stretches carrier wavelengths and event intervals, independent of carrier frequency in its nondispersive geometric-optics regime. It can leave fixed-path redshift constant between observing epochs even though successive signals have increasing travel time. This resolves the *kinematic* timing distinction without importing a spatial expansion factor. It does not prove that such a background exists.

The audit tests compact smooth profiles kappa=2 alpha_0 sin²(pi s/L), with zero kappa outside the interval. These are alternative illustrative voids, not measured void maps. With L=100 million light-years, their mean kappa equals the earlier fitted alpha, giving z=0.00766048 and a ten-day emitted interval arriving over 10.0766048 days. Calibration of alpha remains empirical. The trace and integrating-factor solution agree; no observation establishes the illustrative q amplitude or epoch.

**Required repairs still open for T1:**

- Derive the field and its forces. In a universal lapse coupling a particle initially at rest has a_local=-c² grad ln q, an established conditional result. Matching endpoint clocks does not remove forces in transitions.
- Supply a consistent gravitational source. The exact fixed-flat-space, zero-shift metric has zero normal energy density under the ordinary Einstein constraint (with zero cosmological constant). The earlier [geometry calculation](../research_work/results/time-only-geometry/report.md) therefore prevents treating this prescribed lapse as already sourced by positive photon/companion energy. A curved nonexpanding geometry or explicit modified gravity must be solved.
- Handle the finite past domain: the tested maximum kappa=2 alpha gives a denominator boundary about 6.55 billion years before t_*. That is a failure of extrapolation, not a measured universe age. A longer history requires a different field evolution with new predictions.
- Derive the energy receiver. Time-dependent Hamiltonian work initially belongs to the driving field. Do not additionally credit companions without a dynamical exchange law.
- Reconcile no-loss companions. Ordinary massless waves traveling on the same metric also redshift by the same factor. Assigning all carriers the same clock law while declaring companions exempt is inconsistent. A distinct companion coupling/dispersion law, or revising the no-loss premise, is necessary; neither change is adopted without a full calculation.

## 4. Spectral flux, photon counts and microwave radiation

For number-preserving deterministic shifts, static Euclidean geometric distance D, isotropic source spectral luminosity L_nu, and a common event stretch S_t:

\[
F_{\nu_o}=\frac{\mathcal P\,L_{\nu_e}(S_E\nu_o)}{4\pi D^2 S_t},
\qquad
F_{\rm bol}=\frac{\mathcal P L_{\rm bol}}{4\pi D^2 S_E S_t}.
\]

**Conditional deductions from photon accounting and the frequency Jacobian.** P is a frequency-independent survival fraction in this formula. Frequency-dependent absorption requires its own kernel. At P=1, C0 dims bolometric flux by S_E^-1; T1 dims it by S_E^-2. A curved geometry must replace 4 pi D² by the derived ray-bundle area. A band flux needs the actual filter throughput, source spectrum, dust and detector response. Neither formula by itself validates a luminosity-distance catalog or the SBF calibration.

**Established Planck law, applied to a proposed fixed-volume homogeneous history:** if every photon in an initially thermal bath loses the same energy factor r while number and usual mode density remain fixed,

\[
u_{\nu,f}(\nu)=r^{-3}u_{\nu,\mathrm{Planck}}(\nu,rT_i).
\]

This has the wrong normalization to be a new ordinary blackbody. It is not the same problem as the directed source flux above. A geometric time field cannot automatically inherit this fixed-volume kinetic formula; its mode density, spatial transport and detector calibration must be derived.

**Optional repair R1, not adopted:** retaining only P=r³ photons restores this particular Planck normalization. A concurrent removal rate 3 gamma, where gamma=c alpha, gives dU/dt=-4 gamma U. Partial shifting receives gamma U and the removed-photon sink receives 3 gamma U. All sink energy must be included. If applied universally along beams as well, it makes T1 flux scale as S_E^-5: at S_E=2, eight times fainter than the P=1 T1 prediction. Matching one microwave spectrum this way does not establish acceptable brightness, heating or a background origin. Thermalization, escape, source replenishment, and changed mode structure are other distinct possibilities requiring their own equations.

## 5. A kernel is needed for line widths, images and timing

**Established transport structure with unspecified proposed interactions:** for photons per frequency and solid angle n(s,t,nu,Omega), use a positive transport/scattering kernel or a derived wave equation. Its moments must yield the selected energy drift, photon removal, time-delay distribution, angular scattering, polarization response and receiving-sector exchange. An independently fitted coefficient for every observable is not one predictive theory.

As one limited noise check, independent jumps E -> (1-epsilon)E with Poisson mean N=A/epsilon give:

\[
\langle E\rangle/E_i=e^{-A},\qquad
\frac{\operatorname{Var}(E)}{\langle E\rangle^2}=e^{A\epsilon}-1.
\]

**Established Poisson moments applied to a proposed jump model.** The methanol width yields an approximate upper allowance epsilon<=6.07e-10 if all its width is assigned to this mechanism. Coherent deterministic conversion is not governed by that bound and needs separate linewidth/image calculations. A stationary linear background usually mixes same-frequency modes; whole-photon conversion is attenuation rather than automatic redshift of survivors.

## 6. Conservation, free travel and capture

**Established continuity equations with proposed exchange terms:**

\[
\begin{aligned}
\partial_tu_\gamma+\nabla\cdot F_\gamma&=j-Q-R,\\
\partial_tu_c+\nabla\cdot F_c&=Q-C,\\
\partial_tu_d+\nabla\cdot F_d&=C,\\
\partial_tu_{\rm other}+\nabla\cdot F_{\rm other}&=R.
\end{aligned}
\]

Here R is any additional removed-photon energy sent to a separately accounted sector. For a time-field model include driving-field work and recoil in the sector equations, with total exchange cancelling; do not assert these four scalar equations are a full gravitational conservation law. In a covariant completion require conservation of total stress-energy, with gravitational energy treated in the chosen consistent framework.

No permanent void storage, no free-companion loss, and capture into long-lived extended well-bound states remain desired postulates. An effective C=Gamma_cap u_c is a placeholder. Its cross-section, reverse processes, momentum, pressure, lifetime and spatial modes remain uncalculated for T1. Under a specifically nonrelativistic low-stress deposit, rho_d=u_d/c² and the Newtonian Poisson equation form a known conditional baseline. Relativistic waves require stresses and a metric predicting rotation and lensing together. An ordinary graviton, a gravitational wave and a new companion field are not interchangeable definitions.

## 7. Next executable goals and acceptance outcomes

**User clarification and current galaxy evidence:** deposits may reshape/deepen the local well; no incoming-direction push or cold-source interpretation is fixed. The [scalar-potential response and bulge test](../research_work/results/joint-galaxy-audit/bulge-local-well.md) gives an explicit kernel formulation and three synthetic geometry comparisons. The [joint galaxy/redshift audit](../research_work/results/joint-galaxy-audit/report.md) supplies calibrated force/redshift benchmarks, not a complete mechanism. Its spherical cold-source calculations are conditional examples. Preserve the shared photon/companion energy and timing requirements when developing any new gravitational kernel.

1. **Complete a shared dynamical action or transport kernel.** Derive alpha/kappa, wavefronts, finite-event timing, receiving work, inverse processes and clock response. Outcome: one mechanism supplying all observables without a photon-by-photon or source-by-source tuned rate.
2. **Test real paired spectra across bands.** Obtain calibrated radio, infrared, optical/UV and X-ray line products with common gas-component associations, response matrices, covariance and independent distance provenance. Outcome: measured centroid/width/flux predictions on held-out systems. High-energy continua alone are not known-rest-frequency lines.
3. **Complete observer-time transient inference.** Repair the existing DES numerical gate, validate recovery and interval coverage on multiple injected shapes, model selection/source evolution and then analyze real fluxes. Outcome: independent timing inference without assuming the desired temporal factor in preprocessing.
4. **Solve the radiation bath and brightness together.** Choose and derive sources, sink or thermalizing physics, mode density and heating. Outcome: microwave residual likelihood plus directed-source photometry from the same law, not a separate adjustment for each.
5. **Protect measured local physics by calculation.** Compute optical/radio clock ratios, gravitational shifts, cavities, Doppler links, pulsar timing and laboratory bounds with actual geometries. Outcome: demonstrate small effects or acceptable signals without confusing propagation loss with altered dimensionless constants.
6. **Derive companion transport and gravity.** Resolve no-loss versus shared-metric propagation, production and retention; solve formation, forces, rotation and lensing. Outcome: supported halos with an energy source. The formation-budget magnitude is still deferred, not satisfied.
7. **Validate independent predictions.** Freeze samples, calibration, nuisance parameters, selection and tolerances before new outcomes. Outcome: comparable predictive performance over the declared observational scope; reusing a fit or failing to find a significant difference does not prove parity.

The [coverage register](../research_work/results/electromagnetic-audit/coverage.md) retains every original R01-R32 requirement and distinguishes numerical calculations, measured summaries, conditional comparisons and missing inputs. No finite pass over available files can honestly establish compatibility with every observation in the universe.
