# Codex task backlog — current checkpoint

Latest diagnostic: [candidate decisions and priorities](../research_work/results/candidate-review/review.md). The current candidate review separates failed core mechanisms, reusable component calculations and untested extensions. No tested candidate completes the same-action redshift/timing/capture/gravity chain. All 20 tasks and 32 requirements remain represented. Next prioritize an interaction-derived joint kernel and operational light predictions before another galaxy normalization fit; no candidate or new law is adopted.

**Updated active objective:** [conversion-first prediction program](active-goal.md). The user permits leaving special time/void stretching out as the cause, while retaining all timing observations. The [executed pass](../research_work/results/conversion-first/report.md) is preliminary evidence for T02/T04/T05/T07/T08/T18/T20; it does not complete their full scope or supply fresh validation.

**Governing contract:** [published galaxy distances are fixed facts; no assumed dark matter, expansion or Big Bang](universe-contract.md). This overrides earlier open-candidate wording. The distance decision is resolved.

Active research. T01 completed with documented numerical differences; later tasks have varying amounts of preliminary evidence. No full theory is established. Preserve each task’s scope and dependencies. Keep conceptual alternatives open and consult the user before selecting consequential new laws.

## T01 — Reproduce and freeze the recovered baseline

Status: completed_with_documented_numerical_differences

Prerequisites: none.

**Prompt**

Verify the original snapshot and all latest input hashes. Reproduce run.py, followup.py and circulation.py in an isolated copy with recorded dependencies and seeds; compare numerical outputs with archived results using justified tolerances. Preserve the originals. Identify the five missing historical paper-build files without claiming they block current physics. Produce a source/dependency graph and bounded run-time estimates.

**Required deliverables**

- baseline/manifest-audit.json
- baseline/reproduction-report.md
- baseline/environment-lock.txt
- baseline/result-differences.json

**Completion criteria**

- All available original snapshot hashes match or differences are explained.
- Every latest reported baseline quantity is reproduced or a specific discrepancy is isolated.
- No archived scientific output was overwritten.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

## T02 — Define the universe and freeze the evidence protocol

Status: in_progress

Prerequisites: T01.

**Prompt**

Define operational nonexpansion and what observations remain fixed. Separate measured fluxes/angles/spectra/timing from model-derived masses, distances and ages. Build a registry for R01–R32 with units, sources, nuisance parameters and availability. Register reused data as exploratory; define genuinely unused validation data and predeclare metrics, tolerances and stopping rules before new model selection.

**Required deliverables**

- specification/universe-contract.md
- specification/observable-registry.json
- specification/validation-protocol.json

**Completion criteria**

- Every R01–R32 has an observable or derivation requirement and responsible task.
- The material clock/rod convention and allowed physical changes are explicit.
- No standard cosmology-derived quantity is silently treated as a raw measurement.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

Evidence: {"path": "research_work/results/data-audit/data-foundation-report.md", "scope": "32 requirements mapped; 21 selected inputs hashed/present; SPARC distance methods and Pantheon covariance verified. Contract, likelihoods and fresh validation data still open."}

## T03 — Specify gravitons, deposits and candidate actions

Status: preliminary_candidates_recorded_prerequisites_open

Prerequisites: T02.

**Prompt**

Build a bounded candidate registry for ordinary massless tensor gravitons, any explicitly new companion field, and any deposited state. Specify action, symmetries, degrees of freedom, couplings, cutoff and matter frame for each actual candidate. Compare whole-photon mixing, partial energy transfer and field-driven frequency conversion. Identify how each candidate escapes or fails the already derived clock/conformal obstructions. Do not leave arbitrary functions unlimited.

**Required deliverables**

- theory/candidate-registry.json
- theory/action-specification.md
- theory/assumption-map.md

**Completion criteria**

- Every candidate has explicit field content and either a defined action or a verdict that no action has yet been supplied.
- The deposited state and ordinary gravitational-wave sector are identified without double counting.
- Unsupported candidate families are marked incomplete, not accepted.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

Evidence: {"path": "research_work/results/microphysics/action-specification.md", "scope": "Four candidate families; three explicit comparison actions, one action missing. Five check groups verify illustrative kinematics and static mixing, not complete actions. No candidate adopted."}

Evidence: {"path": "research_work/results/scalar-wave/scalar-wave-derivation.md", "scope": "Derived Z(chi)F^2 wave/scalar equations and verified three closed projected-mode examples. Leading adiabatic characteristic remains metric-null. Distinct electric-response operator recorded but not adopted. No full background, stability or cosmological redshift solution."}

## T20 — Solve the alternative energy ledgers analytically

Status: in_progress

Prerequisites: T01, T02.

**Prompt**

Implement and analytically check the optional no-loss branch: photon loss h_loss, companion secondary loss lambda_c=0, capture Gamma_cap, and permanent deposit leakage lambda_d=0. Compare with shared-secondary-loss and finite-residence branches. Include the stellar reservoir, driver where required, heat, volume changes and boundary flux. Derive impulse and continuous-source solutions, Gamma_cap=h_loss limit, zero-capture limit and long-time growth. All branches must conserve energy; none is selected by assumption.

**Required deliverables**

- energy/ledger-derivation.md
- energy/ledger_solver.py
- energy/analytic-checks.json
- energy/branch-comparison.json

**Completion criteria**

- Impulse energy is exactly partitioned and all physical reservoirs remain nonnegative.
- Continuous injection causes the correctly identified unbounded reservoir unless a derived sink/finite source is present.
- The no-loss branch does not inherit shared-loss energy fractions or a cosmic age without rederivation.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

## T04 — Derive conversion rates and spectral consequences

Status: planned

Prerequisites: T03.

**Prompt**

For each specified candidate derive the kinematically allowed process, amplitudes or mixing equations, forward/inverse rates and all recoil/driver terms. Calculate mean photon energy loss, photon survival, angular diffusion, line broadening and polarization response. Distinguish a photon disappearing from a surviving photon redshifting. Obtain the loss coefficient from the coupling, or show why it remains an empirical parameter. Use primary literature for known limits and declare background assumptions.

**Required deliverables**

- theory/conversion-derivation.md
- theory/conversion-kernels.json
- tests/conversion-limits.json

**Completion criteria**

- Energy, momentum, angular momentum and gauge constraints are checked.
- A nonzero required rate is derived in the stated environment or the candidate is rejected.
- Frequency/angle diffusion and inverse processes accompany the mean loss law.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

Evidence: {"path": "research_work/results/frequency-transfer/conversion-kernels.md", "scope": "Exact stipulated energy kernels demonstrate spectral/photon-count non-identifiability of the energy ledger. Fixed travel conditions have no temporal dilation. No microscopic action rate or clock solution derived."}

Evidence: {"path": "research_work/results/scalar-wave/scalar-wave-derivation.md", "scope": "Derived Z(chi)F^2 wave/scalar equations and verified three closed projected-mode examples. Leading adiabatic characteristic remains metric-null. Distinct electric-response operator recorded but not adopted. No full background, stability or cosmological redshift solution."}

## T05 — Derive emission, propagation, clocks and detection together

Status: planned

Prerequisites: T03, T04.

**Prompt**

Derive emitted atomic frequencies, free-wave dispersion, pulse arrival mapping, cavity response, material rods, flux, angular distance and detection from one candidate action. Reproduce the restricted cavity identity and identify any genuinely changed premise. For the no-loss branch calculate companion versus photon dispersion and distinguish companions from ordinary GW if necessary. Do not impose a separate favorable brightness law.

**Required deliverables**

- theory/observable-transfer.md
- theory/clock-response.json
- theory/ray-and-pulse-solver.py

**Completion criteria**

- Redshift, interval stretch, attenuation and angular distance follow from the same model.
- Operational nonexpansion survives in matter units or the candidate is rejected for this project.
- Atom/cavity protection is derived to a stated approximation, not declared.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

Evidence: {"path": "research_work/results/frequency-transfer/conversion-kernels.md", "scope": "Exact stipulated energy kernels demonstrate spectral/photon-count non-identifiability of the energy ledger. Fixed travel conditions have no temporal dilation. No microscopic action rate or clock solution derived."}

Evidence: {"path": "research_work/results/scalar-wave/scalar-wave-derivation.md", "scope": "Derived Z(chi)F^2 wave/scalar equations and verified three closed projected-mode examples. Leading adiabatic characteristic remains metric-null. Distinct electric-response operator recorded but not adopted. No full background, stability or cosmological redshift solution."}

## T09 — Derive gravity and audit consistency

Status: planned

Prerequisites: T03, T04.

**Prompt**

Vary the specified action and derive field equations, stress-energy exchanges, constraints and the weak-field potentials Phi and Psi. Derive motion and lensing, including pressure and any confining sector. Analyze the quadratic perturbations, degrees of freedom, kinetic/gradient signs, characteristic speeds and EFT range. Distinguish allowed gravitational growth from pathological short-wavelength instabilities. Identify the conditions a stored state would have to satisfy.

**Required deliverables**

- theory/field-equations.md
- theory/weak-field-limit.md
- theory/stability-matrix.json

**Completion criteria**

- No extra gravity is inserted without a stress source or explicit modified-gravity term.
- Motion and light deflection have a common derivation.
- Every retained candidate has a stated viable perturbative domain or a scoped failure.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

Evidence: {"path": "research_work/results/gravity-response/motion-and-lensing.md", "scope": "Explicit two-sector conformal scalar action; linear motion and lensing, measured-G normalization and coupling relation derived. Four force cases, eight lens integrals and three identities checked. No full nonlinear stability or observational fit."}

## T16 — Construct a physical nonexpanding history

Status: planned

Prerequisites: T05, T09, T20.

**Prompt**

Solve background equations in the physical matter frame, then test homogeneous and inhomogeneous perturbations in the stated regime. Determine admissible past/future history, source fuel, entropy and topology. If proposing an extremely old no-loss universe, explicitly replace or extend the fitted n history: its finite past positive-n endpoint cannot simply be extrapolated for trillions of years. Compute redshift drift and returning-path limits. Curvature, topology and stability are separate questions.

**Required deliverables**

- background/solutions.md
- background/history-parameters.json
- background/stability-report.md
- background/drift-and-topology-predictions.json

**Completion criteria**

- At least one history solves the equations and respects operational nonexpansion, or the candidate is rejected.
- Age and luminosity histories are compatible inputs to supply calculations.
- A static balance is not mislabeled stable; sinks and driver energy are specified.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

## T06 — Compute independent photon and companion supply

Status: planned

Prerequisites: T04, T05, T16, T20.

**Prompt**

Build time-dependent photon emissivity from stellar/AGN population histories, dust reprocessing, source geometry and catalog completeness, without using target rotation/lensing residuals to normalize supply. Calculate source-to-target transport separately for no-loss, shared-loss and finite-residence branches that remain viable. Use consistent distances, source emission rates and retarded times. Derive uncertainty envelopes and the finite fuel budget; distinguish a required flux from a predicted flux.

**Required deliverables**

- supply/emissivity-history.json
- supply/source-target-flux.csv
- supply/budget-report.md
- supply/systematics.json

**Completion criteria**

- Reproduces the recovered 10-Gyr limits under their assumptions.
- Each new budget is independently normalized with recorded age, source population and capture assumptions.
- Any claimed resolution of the deficit specifies its source, mass, energy and observational costs.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

## T07 — Derive capture and the stored state

Status: planned

Prerequisites: T03, T04, T09, T20.

**Prompt**

Calculate how incident companions can become retained in a weak galaxy potential. Derive cross-sections, optical depths, scattering, focusing, escape/residence and inverse processes. Specify whether capture produces trapped radiation, matter excitations, a condensate or another field; derive its equation of state and heat. For permanent storage, calculate whether detailed balance, saturation, stress and dissipation permit it over the proposed history. A mere rule that deep wells capture more does not complete this task.

**Required deliverables**

- capture/cross-sections.json
- capture/stored-state.md
- capture/retention-and-heating.csv

**Completion criteria**

- A physically specified capture route has a rate and energy destination, or a quantitative obstruction is shown.
- Pressure and confinement are included.
- No incoming radiation is assumed cold or permanently trapped solely because gravity is attractive.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

Evidence: {"path": "research_work/results/radial-capture/capture-shape-findings.md", "scope": "Analytic spatial benchmarks only. Prerequisite interaction, source history and gravity tasks incomplete. No task completion or physical mechanism claimed.", "findings": ["Fast capture of centrally generated companions gives a conditional flat outer deposited contribution.", "Constant absorption of an isotropic external bath gives uniform or edge-weighted deposition in the tested model.", "Fast-capture local supply exactly matches the archived optimistic local budget; no amplitude resolution."]}

## T08 — Solve transport, deposition and self-gravity

Status: planned

Prerequisites: T06, T07, T09.

**Prompt**

Solve the coupled transport/field problem with declared initial and boundary conditions. Derive the diffusion approximation only where justified; otherwise use kinetic/ray transport. Calculate capture backreaction, evolving potential, anisotropy, finite halo boundary and the stress of free companions as well as deposits. Test whether broad interior profiles arise without tuning to target velocity curves. Compare finite-residence and permanent-growth solutions.

**Required deliverables**

- transport/solver
- transport/convergence.json
- transport/halo-profiles.csv
- transport/energy-ledger.json

**Completion criteria**

- Closed energy/momentum ledger and convergent predictions within predeclared tolerances.
- Halo amplitude, core, shape and outer extent follow from inputs rather than test-target fits.
- Permanent reservoirs are evolved in time rather than forced into an unsupported steady state.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

Evidence: {"path": "research_work/results/radial-capture/capture-shape-findings.md", "scope": "Analytic spatial benchmarks only. Prerequisite interaction, source history and gravity tasks incomplete. No task completion or physical mechanism claimed.", "findings": ["Fast capture of centrally generated companions gives a conditional flat outer deposited contribution.", "Constant absorption of an isotropic external bath gives uniform or edge-weighted deposition in the tested model.", "Fast-capture local supply exactly matches the archived optimistic local budget; no amplitude resolution."]}

## T10 — Test laboratory, local gravity and binaries

Status: planned

Prerequisites: T05, T07, T09.

**Prompt**

Derive and compare optical/hyperfine/cavity drifts, equivalence-principle signals, local gravitational redshift, Solar-System dynamics, Shapiro delay and binary/pulsar behavior. Include instrument/systematic uncertainty and source assumptions. Do not turn a frequency-stability number into a drift bound or a scale ratio into a sigma exclusion. Any environment screening must be solved from the same interaction.

**Required deliverables**

- constraints/local-predictions.json
- constraints/local-likelihood-report.md

**Completion criteria**

- A single parameter region is compatible with declared local tests or is excluded with scoped assumptions.
- No independent local exemptions or per-instrument unconstrained offsets count as proof of agreement.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

## T11 — Test distant light and gravitational waves

Status: planned

Prerequisites: T04, T05, T09, T16.

**Prompt**

Forward-model redshift, supernova flux and time dilation, spectral fidelity, angular sizes/surface brightness and gravitational-wave propagation. Include GW emission delays and source calibration. Derive whether no-loss companions are the same measured tensor waves; test the consequences rather than impose identical losses. Compare to the appropriate published likelihoods without importing cosmology-dependent distances unexamined.

**Required deliverables**

- constraints/light-gw-predictions.json
- constraints/light-gw-likelihood-report.md

**Completion criteria**

- Same parameters predict frequency, timing, flux and tensor propagation.
- No separate refit of the redshift history conceals a changed clock or brightness model.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

## T12 — Predict galaxy dynamics and environmental variation

Status: planned

Prerequisites: T08, T10, T11.

**Prompt**

Freeze global parameters on declared training data and predict other galaxies. Include disk geometry, gas, bulges, distance/inclination uncertainty, pressure-supported systems, vertical gravity, satellites, baryonic scaling relations and environmental controls. Compare luminosity/history predictors against constant-bath and shuffled-source controls with the registered metric. Treat recovered SPARC splits as exploratory.

**Required deliverables**

- validation/galaxy-predictions.csv
- validation/galaxy-controls.json
- validation/galaxy-report.md

**Completion criteria**

- No per-test-galaxy halo/gain fits.
- Scatter, outliers, scaling relations and controls are reported, not only mean RMSE.
- A causal illumination advantage is demonstrated or explicitly absent.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

## T13 — Predict independent lensing

Status: planned

Prerequisites: T08, T09, T12.

**Prompt**

Use the frozen dynamical parameters to predict strong/weak lensing, shear, convergence and time delays from Phi and Psi and the model distance relation. Model source geometry, baryons and line-of-sight structure. Select at least one lensing dataset not used to set the halo amplitudes. Declare any remaining new nuisance parameters.

**Required deliverables**

- validation/lensing-predictions.csv
- validation/lensing-report.md

**Completion criteria**

- Lensing is generated from the same action and source distribution as rotation.
- No lensing-only gravity normalization rescues the fit.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

Evidence: {"path": "research_work/results/gravity-response/motion-and-lensing.md", "scope": "Explicit two-sector conformal scalar action; linear motion and lensing, measured-G normalization and coupling relation derived. Four force cases, eight lens integrals and three identities checked. No full nonlinear stability or observational fit."}

## T14 — Test clusters, mergers and capture offsets

Status: planned

Prerequisites: T08, T12, T13.

**Prompt**

Reproduce the conditional eleven-cluster benchmark and expand to radial dynamics and independent lensing with gas/hydrostatic systematics. Model merger transport to predict whether deposited energy follows galaxies, gas or an independent component. Include capture/stripping times, drag, offsets and post-merger survival. Use the galaxy parameter set without a cluster-specific gain.

**Required deliverables**

- validation/cluster-predictions.csv
- validation/merger-evolution.json
- validation/cluster-merger-report.md

**Completion criteria**

- Both equilibrium and collision behavior are explained within one capture/retention prescription or the branch fails.
- Hydrostatic-derived quantities are not mislabeled independent lensing measurements.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

## T15 — Test stellar, thermal and formation consequences

Status: planned

Prerequisites: T06, T07, T08, T09.

**Prompt**

Calculate stellar cooling/lifetimes, gas heating/ionization, disk survival, dynamical friction, compact-object capture and the formation chronology of halos. Address how stars produce the first deposits before a pre-existing large halo is available. Include reemission, metal production and finite fuel. Test whether strong capture locally would overheat or destabilize the systems.

**Required deliverables**

- validation/thermal-stellar-budget.json
- formation/history-report.md
- validation/compact-object-scope.md

**Completion criteria**

- Energy losses/heating are compatible with the stipulated stellar/gas observations under quantified assumptions.
- There is no circular requirement for an already assembled final halo to create its own sole source.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

## T17 — Construct and test full-universe predictions

Status: planned

Prerequisites: T06, T11, T15, T16.

**Prompt**

Derive a radiation/thermal history and perturbation evolution predicting CMB spectrum plus TT/TE/EE and lensing, large-scale structure/BAO, light-element abundances, reionization, source chronology and diffuse backgrounds. Derive any standard ruler in this universe. Add photon/companion/deposit and heat budgets, including the dark night sky and entropy. If only a galactic effective theory survives, report its restricted scope rather than claiming a completed cosmology.

**Required deliverables**

- cosmology/transfer-functions
- cosmology/observable-predictions.json
- cosmology/joint-report.md

**Completion criteria**

- The same background and interaction predict the registered cosmic observables without importing a contradictory expansion history.
- Missing cosmological mechanisms remain explicit and prevent a full-universe completion claim.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

## T18 — Freeze discriminating tests and assess identifiability

Status: planned

Prerequisites: T12, T13, T14, T15, T17.

**Prompt**

Audit parameter identifiability, model-selection history and uncertainty propagation. Choose a bounded set of distinguishing predictions such as illumination at fixed baryons, growth with age, spectral broadening, merger lag, GW losses or redshift drift. Freeze predictions before revealing reserved outcomes. Compare with baryons-only and declared standard/empirical benchmarks using a predeclared joint scoring rule, including complexity and selection effects.

**Required deliverables**

- validation/frozen-predictions.json
- validation/prediction-seal.json
- validation/joint-verdict.md

**Completion criteria**

- At least one nontrivial independent prediction distinguishes the model or observational equivalence is reported.
- All tested candidates and failed predictions remain in the record.
- Unidentifiable coupling/capture/history products are not reported as separately measured constants.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.

## T19 — Assemble the final theory or rejection dossier

Status: planned

Prerequisites: T18.

**Prompt**

Assemble an auditable paper/specification with the explicit action, derivation chain, universal parameters, initial/boundary conditions, every registered requirement, numerical pipeline and observational verdict. Package reproductions and source attributions. If the concept fails, identify the minimal failing assumptions and the evidence without pretending the whole class of alternatives is ruled out.

**Required deliverables**

- release/theory-or-rejection.md
- release/requirements-status.json
- release/reproduction-guide.md
- release/source-manifest.json

**Completion criteria**

- Every R01–R32 is passed, scoped, rejected or blocked with evidence.
- No claim exceeds demonstrated scope.
- A new reader can reproduce the reported conclusions from the released specification and data.

**When a finding conflicts with assumptions**

Record the result under its assumptions, keep the conceptual path open, and ask the user in plain language which physical assumption to revise when a material choice is needed. Do not silently change laws or report an incompatible model as successful.
