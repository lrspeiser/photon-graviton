# Prediction checklist: photon energy, companions and gravity

Latest diagnostic: [capture and storage](../research_work/results/capture-storage/derivation.md). A reversible two-step receiver separates capture from long-lived storage in a conditional three-state model. Its energy ledger includes returned companions and shelving radiation; finite capacity and environmental reverse transitions limit accumulation. Receiver identity, absolute rates, spatial support and the full gravity response remain unproved.

Updated 9 September 2026. This checklist develops the user's request to list outcomes and available prediction formulas. It supplements the full roadmap; it does not mark its tasks complete or adopt a microscopic interaction.

## Working starting point

The user has now explicitly authorized a conversion-first research pass without special time/void stretching as the cause. Their original time concept meant that time stretches light and thereby causes photon-energy loss; it did not mean nearby galaxies lack redshift. See [the revised goal](active-goal.md) and [the executed exploratory pass](../research_work/results/conversion-first/report.md). This clarification does not waive timing observations.

Use **partial photon-energy transfer into companion waves as the working explanation of the distance-related redshift component**. The surviving photons must lose frequency, rather than merely disappear. Smooth transfer and many small transfer events remain alternatives. Ordinary source motion and gravitational endpoint shifts must be modeled separately; not every measured redshift is assigned to conversion.

This is a working hypothesis, not a verified physical cause. No interaction-derived rate, companion identity, special void-clock law, no-loss propagation rule, permanent storage rule or enhanced gravity law is locked in. The [universe contract](universe-contract.md) remains authoritative: fixed published galaxy distances, nonexpansion, no assumed dark-matter population or Big-Bang origin, complete energy conservation.

## The first forward formula

Let s be physical path length, E photon energy in an explicitly chosen common energy standard, and alpha the fractional loss per unit length. For deterministic drift:

```
dE/ds = -alpha(s, E, environment, time) E
tau = integral_along_actual_path alpha ds
E_received = E_emitted exp(-tau)
1 + z_conversion = E_emitted/E_received = exp(tau)
```

The integral is along the solution when alpha depends on E. Alpha has units of inverse length and tau is dimensionless. The final equality is an observed spectral prediction only after source/receiver atomic standards and other shifts have been treated consistently. A frequency-independent alpha gives the same fractional shift to different spectral lines in this stipulated model.

For constant alpha and path length R:

```
z_conversion(R) = exp(alpha R) - 1
alpha_inferred = ln(1 + z_conversion)/R
```

The archived empirical coefficient is 7.7315e-5 per million light-years. It is a comparison input, not newly fitted or derived here. Reconstructing an alpha separately for every observed redshift is not a predictive test. A common law must predict new objects with its parameters held fixed.

In a full-transfer branch, before subsequent companion propagation/capture:

```
E_companion_produced = E_emitted - E_received
E_companion_produced/E_emitted = z_conversion/(1 + z_conversion)
```

Examples: a conversion redshift of 0.1 transfers 9.09% of initial energy; 1 transfers 50%; 2 transfers 66.67%. These are conditional energy identities, not measured companion production. If only a fraction eta enters companions, multiply the lost energy by eta and account for the rest in other sectors. Eta here is an energy-partition fraction, not the gravity-response parameter in older reports. Photon number and arrival intervals are separate quantities.

## Available formulas and their limits

Status: **conditional** means calculable after supplying assumptions/parameters; **missing** means no adopted physical law currently predicts the outcome. None of these rows claims new observational validation.

| Outcome | Available formula or calculation | Status / what remains |
|---|---|---|
| Redshift versus fixed distance | z_conversion = exp(integral alpha ds) - 1 | Conditional drift law; derive or globally fit alpha, specify clocks and other shifts, then test unexposed objects. |
| Companion energy supplied | E_c,produced = eta E_emitted z_conversion/(1+z_conversion) | Conditional bookkeeping; eta and the microscopic transfer must be derived. |
| Photon survival | P_gamma = exp(-integral a_remove ds) for independent removal | Conditional; removal alone does not redden surviving photons. Derive both removal and partial transfer. |
| Added spectral width | sigma_E/mean(E) = sqrt(S^epsilon - 1), S=E_initial/mean(E) | Exact for independent Poisson events of constant fractional loss epsilon. Smooth deterministic drift has zero added energy width in its idealized kernel. Neither establishes image fidelity. |
| Event duration | S_t = dt_received/dt_emitted; fixed path and fixed travel speed give S_t=1 | Missing successful joint frequency/timing law. Redshift alone does not imply S_t=1+z. |
| Brightness | F = L P_gamma/(4 pi R^2 S_E S_t), S_E=E_emitted/E_received | Conditional static Euclidean, isotropic, unlensed illustration with common energy/clock definitions; derive geometry and full spectra for real comparisons. |
| Angular size / surface brightness | theta approximately ell/R; I=F/Omega | Conditional small-angle Euclidean geometry and fixed physical size ell; derive ray bundles and rulers in candidate theory. |
| Image blur / polarization | Joint energy-direction-time-polarization transport kernel | Missing interaction-derived kernel; the energy ODE supplies none of these. |
| Clocks and rods | Calculate atomic and cavity frequencies from the same matter/field law | Missing adopted completion; a special void-clock effect is not assumed. |
| Capture probability | P_cap=1-exp(-integral beta ds) | Conditional independent capture along a specified path, beta inverse length. Does not establish permanent trapping or bound-state support. |
| Stored energy | dE_d/dt=P_cap,in - E_d/t_d | Conditional one-zone balance: P_cap,in is captured power, t_d lifetime. Escape/decay energy must enter another ledger sector. Permanent storage is t_d infinite. |
| Total conservation | dE_total/dt=P_external,in-P_external,out | Accounting requirement with stellar fuel, fields, heat, recoil and escaping energy included; not a conversion or gravity law. |
| Radial deposit production | dP/dr=-alpha P; dC/dr=alpha P-beta C; dD/dr=beta C | Conditional outward point-source power model. P photon power, C companion power, D cumulative deposited power; P+C+D constant. Not stored-energy equilibrium. |
| Rotation | v_c^2=r dPhi_motion/dr for weak-field circular orbits | Conditional observable map; need field equation/source distribution. Ordinary benchmark M_d(<r)=integral u_d dV/c^2 is not the adopted response. |
| Lensing | Deflection vector = (1/c^2) integral grad_perp(Phi+Psi) dl | Weak-field static metric ds^2=-(1+2Phi/c^2)c^2dt^2+(1-2Psi/c^2)dx^2 with light following its null rays. Derive potentials and photon coupling before using it. |
| Source history and supply | E_emitted=integral L(t) dt; E_d cannot exceed its accounted source supply | Conditional finite-fuel/source-history calculation, including external sources and any explicitly debited driver. Longer time alone does not supply fuel. |
| Background radiation, structure, abundances and strong gravity | Full coupled evolution equations plus initial/boundary conditions | Missing unified candidate predictions; templates or fitted curves are not a first-principles explanation. |

For Poisson spectral transfer, S is a centroid-energy shift, not automatically an observed line redshift. For lensing, light-modifying interactions can invalidate the assumed metric-ray formula; propagate light with the actual candidate law. In all branches, momentum, reverse reactions and stability must also be checked.

## Concrete next calculations and required outcomes

1. **Declare the redshift forward model.** Specify energy and clock standards, fixed adopted distances, source-motion and endpoint corrections, and either one global alpha or a finite environment-dependent law. Output a parameter/units table and explicit assumptions. Do not assign a free alpha to each object.
2. **Test spectral shift, width and photon counts jointly.** Compare deterministic drift with the existing fractional-event kernel. Output predicted line profiles, transferred energy and count survival using the same parameters. Label any unexplained timing, direction or polarization behavior.
3. **Fit a declared training set and freeze parameters.** Audit distance and redshift provenance, measurement errors and selection effects. Use the archived coefficient only as a baseline. Output residuals and uncertainty intervals; account for uncertainty metadata without changing the adopted fictional distance values.
4. **Predict an unexposed validation set.** Evaluate distance-redshift behavior, frequency dependence, spectral width, brightness and transient duration together. Predeclare metrics and acceptable residuals from observational errors; report both successes and failures. No observational fit has been done by creating this checklist.
5. **Translate validated photon losses into an energy-source map.** Include luminosity history, dust, external radiation, partition and transport. Output companion production with uncertainties; do not equate production with captured energy.
6. **Derive capture and storage.** Predict where energy accumulates, how long it remains and what supports it. Output spatial/time profiles and a complete fuel/field/heat/escape ledger.
7. **Derive motion and lensing together.** Obtain both from the same adopted field interaction and the predicted deposit distribution, calibrated against local gravity. Output rotation and lensing before adding per-object normalization.
8. **Run the entire 32-area coverage register.** Missing predictions stay explicitly missing. New constraints extend the register; favorable results from incompatible theories cannot be combined into one claimed solution.

These are research specifications, not scheduled jobs. A candidate proceeds from a stipulated transport law toward a microscopic cause by deriving its joint transport kernel and gravitational response from one consistent interaction. Agreement with data can constrain or distinguish causes but does not guarantee a unique cause.

## Evidence and full coverage

- [Existing spectral kernels and mathematical checks](../research_work/results/frequency-transfer/conversion-kernels.md)
- [Energy-exchange definitions and conditional gravity bound](energy-exchange.md)
- [Current checkpoint and linked derivations](../research_work/results/RESEARCH-CHECKPOINT.md)
- [Full roadmap with all 32 requirements](theory-development-roadmap.md)

The complete original coverage register is reproduced below for checking that no observation category is lost. Its task IDs refer to the existing backlog. The formula table above supplies only preliminary coverage; rows without a completed candidate prediction remain open.

| ID | Area | Required prediction | Existing tasks |
|---|---|---|---|
| R01 | Operational definitions | Define physical rods/clocks, nonexpansion, photons, gravitons, deposits and the matter metric. Distinguish coordinate changes from observable changes. | T02, T03 |
| R02 | Microscopic conversion | Specify the interaction, allowed channels, energy/momentum/angular-momentum exchange, rate and inverse processes. | T03, T04 |
| R03 | Conservation | Account for photons, free gravitons, deposits, matter, driver, heat and boundary flux. No free energy or double counting. | T04, T06, T08 |
| R04 | Quantum/classical consistency | Gauge constraints, propagating degrees of freedom, positive kinetic energies, causal characteristics, perturbative regime and radiative stability. | T03, T09 |
| R05 | Graviton identity | Determine whether the new radiation is the ordinary massless tensor sector or a new field. Give its mass, spin, dispersion, polarization, interactions and observational consequences. | T03, T11 |
| R06 | Source production | Stellar and AGN luminosity histories, radiation bands, dust reprocessing, obscured sources, completeness and formation of early halos. | T06, T15, T17 |
| R07 | Redshift | Predict observed atomic-line ratios versus independently inferred distance and source/observer conditions. Include peculiar velocities and gravitational endpoint shifts. | T05, T11 |
| R08 | Time dilation | Predict complete transient light curves and spectral aging, not merely carrier-frequency shifts. | T05, T11 |
| R09 | Brightness and photon counts | Predict flux, attenuation, luminosity distance and source calibration from the same transport model. | T05, T11 |
| R10 | Angular distances and surface brightness | Predict source sizes, angular distances and reciprocity/duality relations; identify justified rulers instead of importing an expansion-based ruler. | T05, T11, T17 |
| R11 | Spectral and image fidelity | Line width, chromaticity, scattering blur, polarization/birefringence, coherence and arrival-time dispersion. | T04, T05, T11 |
| R12 | Atomic clocks and cavities | Optical/hyperfine transitions, cavity lengths, fine-structure response, secular and oscillatory drifts, terrestrial versus interstellar conditions. | T05, T10 |
| R13 | Local gravity | Equivalence principle, inverse-square tests, Solar-System ephemerides, gravitational redshift, Shapiro delay, perihelion advance, light deflection, binaries and pulsars. | T09, T10 |
| R14 | Capture and retention | Derive absorption/scattering cross-sections, optical depth, escape, residence time, equation of state and any bound state. | T07, T08 |
| R15 | Halo profiles | Predict cores, radial slopes, outer truncation, total mass, shape and time evolution without fitting a halo separately to each test object. | T08, T12 |
| R16 | Galaxy rotation and scaling | Spiral and dwarf rotation curves, radial-acceleration relation, baryonic Tully–Fisher relation, diversity and scatter. | T12, T18 |
| R17 | Other galaxy dynamics | Pressure-supported dwarfs/ellipticals, satellite orbits, stellar streams, vertical disk gravity, bars and dynamical friction. | T12, T15 |
| R18 | Lensing | Strong and weak lensing, shear, convergence, deflection and time delays from the same fields that predict motion. | T09, T13 |
| R19 | Clusters | Hot gas, hydrostatic bias, pressure support, galaxy velocities, radial profiles and lensing with no cluster-only normalization. | T13, T14 |
| R20 | Merging systems | Predict whether deposits follow galaxies, plasma or neither during collisions; model offsets, stripping, lag and reconnection/reformation. | T07, T14 |
| R21 | Environmental dependence | Isolated versus group galaxies, voids, tidal dwarfs, low-luminosity and unusually low/high inferred-mass systems; dependence on illumination and history. | T06, T12, T14 |
| R22 | Stellar and thermal physics | Stellar structure/lifetimes, solar constraints, cooling, heating, ionization, metallicity and whether capture would overheat gas or evaporate a halo. | T07, T15 |
| R23 | Gravitational waves | Propagation speed, damping, dispersion, polarizations, binary radiation, multimessenger delay and generated stochastic/high-frequency backgrounds. | T04, T10, T11 |
| R24 | Compact objects | Black-hole capture versus long-lived halos, accretion, strong-field boundary conditions, ringdown and consistency with compact-object observations. | T09, T15 |
| R25 | CMB spectrum | Explain the near-blackbody radiation field, its intensity, temperature and permitted spectral distortions without assuming the desired result. | T17 |
| R26 | CMB angular structure | Predict temperature and polarization correlations and lensing from perturbations, with a physical source of the observed angular scales. | T16, T17 |
| R27 | Structure and BAO | Explain galaxy clustering, voids, growth, baryon acoustic features and the matter power spectrum in the chosen history. | T16, T17 |
| R28 | Abundances and cosmic chronology | Light-element abundances, nucleosynthesis alternative/history, stellar ages, reionization, metal production and distant populations. | T15, T17 |
| R29 | Static background and topology | Solve background field equations, distinguish curvature from topology, test stability, redshift drift, age, past/future endpoints and repeated images. | T16 |
| R30 | Global radiation and entropy | Finite fuel, accumulated radiation, dark night sky, infrared/optical/X-ray backgrounds, heat sinks and entropy production. Recycling photons cannot renew their energy. | T06, T16, T17 |
| R31 | Statistical identifiability | Separate gravity response, capture, age, luminosity and reservoir degeneracies. Include selection functions and shared calibration errors. | T01, T06, T18 |
| R32 | Distinguishing predictions | Freeze a small set of tests that separate this theory from baryons alone, fitted halo models and empirical modified-force laws. | T18, T19 |
