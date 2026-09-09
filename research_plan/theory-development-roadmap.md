# Photon-to-graviton deposition: theory-development roadmap

## Governing universe contract

The user has resolved the distance choice: published galaxy distances are fixed fictional facts, even if originally obtained from Hubble flow. Assumed dark matter, cosmic expansion and Big-Bang explanations/initial conditions are outside the active search. Retain observations, but do not import those interpretations or derived rulers, ages and halo profiles. This overrides earlier language that left all cosmological alternatives open. See [universe-contract.md](universe-contract.md).


Prepared 9 September 2026 from the recovered photon-graviton project. This is a proposed research program for a fictional, operationally nonexpanding universe constrained by the same observational records as ours. It is not a completed derivation or a claim that our universe follows these laws.

## 1. The objective and the present position

Develop one explicitly defined interaction in which radiation supplies a gravitating sector, that sector is transported into and retained around galaxies and clusters, and the resulting field predicts the required dynamics and lensing. The same interaction must predict the light received by observers and remain consistent with laboratory physics and the universe's history.

“Same observations” means the measured spectra, fluxes, angles, periods, proper motions, light curves and other selected records. Catalog distances, stellar masses, gas masses, ages and hydrostatic/lensing mass estimates are often derived quantities. Keep their calibration assumptions explicit and recompute them when the new laws change the inference. Initially retain the recovered astronomical inputs as conditional benchmarks; do not silently change them to improve a fit.

The target excludes an independently inserted, arbitrary dark-particle halo. A deposited, nonluminous gravitating field would still be an extra gravitational component; changing its name does not remove its mass/energy accounting. Its distribution and amplitude must follow from sources and dynamics.

There is not yet a successful first-principles theory to reproduce. The equations in section 3 are a derivation of necessary relationships under stated assumptions, followed by the missing steps that would make them a complete theory. A fit, a conservation identity and an interaction derived from an action are different levels of evidence.

### What the recovered work already establishes

| Finding | Meaning for the program |
|---|---|
| 271/276 original snapshot files recovered with matching SHA-256 hashes | The current companion scripts, inputs, outputs and final v9 paper are available. Five historical v9 construction files remain absent; they do not block numerical work. |
| SPARC checkpoint: 149 selected galaxies, 3,150 velocity points, old 89/29/31 split | Useful regression benchmark; these galaxies are already exposed, not a fresh blind holdout. |
| DustPedia overlap: 26 galaxies, using 814 usable catalog sources | Independent luminosity information is available, but completeness and radiation histories remain unsolved. |
| Required/supplied energy median: 5,428 even at 100% conversion of present luminosity for 10 billion years | Ordinary deposited mass-energy cannot meet this benchmark under those assumptions. This is not a theorem about every possible history or external reservoir. |
| Local conversion plus catalogued neighbors: median shortfall 6.85e8 | Modest efficiency adjustments cannot repair this particular supply model. |
| Catalog illumination does not beat the constant-bath control on the declared log-residual metric | Existing fits do not establish the proposed environmental source connection. |
| Thin-edge deposition performs poorly; broad interior profiles perform better | Derive transport through the interior, not just collection at an outer surface. |
| A pointwise empirical depth-dependent force law improves galaxy validation and has promising cluster transfer | A target for derivation, not a solved normalization from photon supply. Cluster tests are not independent lensing tests. |
| Universal index/clock constructions encounter rod, clock and operational-expansion conflicts | A candidate must address these existing obstructions explicitly. |
| A fitted background history can be reconstructed but has sensitivity and stability problems | Reverse-engineering a potential is not proof of a robust physical solution. |

These numerical statements are from `companion_causal_test/report.md`; clock and action obstructions are in `minimal_clock_interaction/derivation.md`, `conformal_action_derivation/derivation.md`, `nature_tests/report.md`, and `interaction_stability_attempt/checks.json`. Paths here and below are relative to the recovered repository, ..

## 2. Everything the theory must answer

This is a coverage register, not a claim that one finite checklist captures every conceivable observation. Each row needs a prediction, input provenance, uncertainty model and pass/fail result. Newly discovered constraints extend the register rather than being ignored.

| ID | Area | Required response / observable | Main tasks |
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

Lensing offsets in merging clusters are a direct reason to include R20; see the [Clowe et al. observations](https://arxiv.org/abs/astro-ph/0608407). The CMB program must cover spectra and correlations rather than importing fitted standard-model parameters as new-universe inputs; see the [Planck observational analysis](https://arxiv.org/abs/1807.06209). These papers identify constraints, not proof that the fictional mechanism is impossible.

## 3. First-principles derivation: what follows, and what still needs a law

### 3.1 Define the degrees of freedom before fitting

Use natural units only for the action below. An ordinary-graviton candidate can begin with the *template*

\[
S=\int d^4x\sqrt{-g}\left[\frac{M_{\rm Pl}^2}{2}R-\frac14F_{\mu\nu}F^{\mu\nu}+\mathcal L_m+\mathcal L_X+\mathcal L_{\rm int}\right].
\]

Here g is the dynamical metric, F is electromagnetism, m denotes ordinary matter, and X denotes any explicitly introduced reservoir/driver needed for conversion or storage. Ordinary gravitons are propagating tensor perturbations of g; do not add an independent copy of their energy on top of the same perturbation contribution. If a separate massive tensor/scalar/vector is proposed, it needs its own action and constraints and must be labeled as that branch.

This is a template with unspecified terms, **not a working interaction**. Choosing L_int, demonstrating its symmetries, and deriving its rates is the central open problem. Variation must give the matter, electromagnetic, reservoir and gravitational equations. Covariance yields on-shell conservation of the total stress-energy, with compatible sector exchanges. Effective gravitational-wave stress-energy requires a stated averaging regime; there is no generally covariant local gravitational-energy tensor that can be assigned arbitrarily.

Minimal electromagnetic coupling to gravity does not by itself supply the required cosmological conversion rate and trapping law. Background-field photon/graviton conversion is an existing comparison mechanism, with assumptions about geometry and magnetic fields that must be carried through; see [Hwang and Noh](https://arxiv.org/abs/2310.04150). Its existence does not establish redshift or halo deposition.

### 3.2 Decide what “photons convert” means

Three distinct mechanisms require distinct calculations:

1. Whole-photon conversion, gamma -> g, possibly coherent mixing in a background. This can remove photons; in a stationary conversion background it does not automatically lower the frequency of the surviving photons.
2. Partial loss, gamma -> gamma-prime + g, with a residual lower-energy photon. Specify what receives recoil and supplies any necessary background energy/momentum. In Lorentz-invariant empty space with standard massless dispersion, p=p-prime+k and p²=p-prime²=k²=0 force p-prime·k=0, hence collinear outgoing momenta. Generic noncollinear decay is unavailable. Whether a nonzero collinear process/rate exists depends on the interaction; this kinematic statement is not a blanket prohibition of background-assisted processes.
3. Frequency exchange with a dynamical medium/field. A time-dependent dispersion can change photon energy, but the work initially belongs to the driver sector unless the action derives transfer into gravitons. It cannot also be counted as independent graviton production without an additional term.

Derive the scattering/mixing kernel, forward and inverse processes, and angular/polarization dependence. As a diagnostic, independent events with fractional loss epsilon and Poisson mean N give mean log-energy loss approximately N epsilon and variance approximately N epsilon² for small epsilon. This ties a redshift mechanism to line broadening; correlations/coherence require their own derivation, not that approximation.

### 3.3 Recover the propagation law conditionally

If the residual photons obey an achromatic fractional-energy loss per geometric length,

\[
\frac{dE_\gamma}{ds}=-\alpha E_\gamma,
\quad E_\gamma=E_e e^{-\alpha s},
\quad S_E=E_e/E_o=e^{\alpha R}.
\]

The recovered local coefficient is alpha = 7.7315e-5 per million light-years. This is fitted input until the interaction predicts it or derives it from independently measured parameters. Observed spectroscopic redshift equals S_E only when the reference transition at emission and observation has the stipulated unchanged response, after other shifts are accounted for.

A stationary path-loss law has t_o=t_e+constant for repeated identical emissions, so dt_o/dt_e=1. Frequency loss alone therefore does not derive signal stretching.

For the separate homogeneous nondispersive history branch, take omega=c k_wave/n(t), conserve k_wave, and hold geometric separation fixed:

\[
R=\int_{t_e}^{t_o}\frac{c\,dt}{n(t)},\qquad
\frac{dt_o}{dt_e}=\frac{n_o}{n_e}.
\]

If atomic reference frequencies are nu_A(t), then

\[
1+z_{\rm spec}=\frac{n_o}{n_e}\frac{\nu_{A,o}}{\nu_{A,e}}.
\]

The interval relation uses the selected coordinate clock; convert to material proper-time intervals in a completed action. Setting n-dot=gamma n^p recovers the empirical history branch, not a microscopic explanation. A candidate must derive the appropriate n, dispersion and material response together.

### 3.4 Derive brightness and angular observables from the same model

For a static Euclidean illustration, let P_gamma be the survival fraction, S_E the photon-energy ratio and S_t the arrival-interval stretch measured in the relevant emitter/receiver clocks. An isotropic source then gives

\[
F=\frac{L P_\gamma}{4\pi R^2S_ES_t},\qquad
D_L=R\sqrt{\frac{S_ES_t}{P_\gamma}}.
\]

If source rulers are fixed and rays follow Euclidean geometry, D_A=R. For P_gamma=1 and S_E=S_t=S, D_L=RS. Whole-photon conversion changes P_gamma, so it changes this calculation. Curvature, lensing, varying rods and nonmetric propagation require a fresh area/ray-bundle derivation. Test surface brightness and distance duality rather than assuming either the usual or this special relation. [More et al.](https://arxiv.org/abs/1612.08784) provide a primary derivation of distance-duality robustness in a specified generalized-electrodynamics class.

### 3.5 Respect the clock/rod obstruction already derived

For a cavity mode nu_C=j c/(2 n L), exact protection of nu_C/nu_A implies n L nu_A=constant. Combined with the preceding redshift relation,

\[
1+z_{\rm spec}=L_e/L_o.
\]

If material lengths are also fixed, this restricted model gives no redshift. If rods shrink, distances grow in rod units; the model may violate the stipulated operational nonexpansion. The recovered universal conformal-metric family reaches a related obstruction. These are statements about specified classes, not every possible interaction.

A new branch must name the assumption it changes and calculate the resulting spectra, ruler behavior, force law and laboratory response. “Cavity photons are exempt” is not an interaction. Exact protection through all history must also be distinguished from the finite sensitivity of present-day observations. [Atom–cavity measurements](https://arxiv.org/abs/2008.08773) motivate a joint clock analysis; an oscillatory limit must not be treated as a secular-drift likelihood.

### 3.6 Keep the entire energy ledger

In a local static-volume illustration, a possible bookkeeping decomposition is

\[
\partial_tu_\gamma+\nabla\cdot F_\gamma=j_\star-Q_{\gamma g}-Q_{\gamma X},
\]
\[
\partial_tu_g+\nabla\cdot F_g=Q_{\gamma g}-Q_{gX}-Q_{\rm cap}+Q_{\rm return},
\]
\[
\partial_tu_d+\nabla\cdot F_d=Q_{\rm cap}-Q_{\rm return}-Q_{\rm heat},
\]
\[
\partial_tu_X+\nabla\cdot F_X=Q_{\gamma X}+Q_{gX},\qquad
\partial_tu_{\rm heat}+\nabla\cdot F_{\rm heat}=Q_{\rm heat}.
\]

The stellar reservoir loses j_star. Sum all sectors and boundary fluxes; exchanges cancel. In a curved/dynamical setting derive the corresponding covariant equations and pressure-work terms. These equations enforce accounting but do not determine the Q terms. The microscopic action must supply them, including signs and inverse reactions. There need not be a globally conserved coordinate energy in an arbitrary time-dependent spacetime; the above scalar ledger is limited to its stated regime.

For one illustrative shared-secondary-loss model, write x=ln S and impose dE_gamma/dx=-E_gamma, dE_g/dx=E_gamma-E_g, with E_g(0)=0. Integration gives

\[
E_\gamma/E_0=e^{-x},\quad E_g/E_0=xe^{-x},\quad
E_X/E_0=1-(1+x)e^{-x}.
\]

The free graviton/companion fraction peaks at 1/e. The last term is an explicit receiving sector, not vanished energy. This reproduces the recovered bookkeeping branch; whether actual gravitons obey these rates remains to be derived. Choosing different secondary losses requires recalculating the whole ledger.

### 3.7 Derive capture instead of assuming “a well stores gravitons”

Ordinary massless radiation in a weak galactic potential is deflected, not automatically converted into a stationary cold halo. A typical galactic escape speed is far below c. Permanent storage requires a demonstrated interaction, trapping configuration, altered dispersion/bound state, or transfer to other degrees of freedom. Strong-field photon orbits near compact objects are not a demonstration of stable galactic retention.

If an absorber density n_a and cross-section sigma_cap are derived, then alpha_cap=n_a sigma_cap and optical depth tau_cap=integral alpha_cap ds. In a simple independent-absorption model, epsilon=1-exp(-tau_cap). Calculate scattering and reemission separately.

For an isotropic bath u_g, group speed v_g, and a spherical collector radius R,

\[
P_{\rm cap}=\epsilon\pi R^2v_g u_g.
\]

The factor follows from inward hemispheric flux v_g u_g/4 multiplied by 4 pi R². Focusing, anisotropic illumination or a different collector changes the effective cross-section and must be solved explicitly. For constant P_cap, zero initial deposit and residence time tau,

\[
\dot E_d=P_{\rm cap}-E_d/\tau,
\quad E_d(T)=P_{\rm cap}\tau(1-e^{-T/\tau})\le P_{\rm cap}T.
\]

Show where escaped or thermalized energy goes. Massless isotropic free radiation has pressure u/3; treating it as pressureless cold matter requires a real transformation or confining stresses. If a massive/condensed deposit is proposed, derive its creation and equation of state. Do not change the particle identity halfway through the calculation.

### 3.8 Infer the gravitational target, then solve the forward problem

Under the restricted cold, spherical-extra-source Newtonian benchmark,

\[
M_X(<r)=\frac{r[v_{\rm obs}^2-v_b^2]}{G},\qquad E_{\rm req}=M_Xc^2.
\]

This determines a diagnostic requirement. It cannot be used as the independent source normalization for the same object's prediction. The disk baryonic force need not be spherical.

A cored template

\[
\rho_X=\frac{\rho_0}{1+(r/r_c)^2}
\]

integrates to

\[
M_X=4\pi\rho_0r_c^2[r-r_c\arctan(r/r_c)],
\quad v_X^2=4\pi G\rho_0r_c^2[1-(r_c/r)\arctan(r/r_c)].
\]

It has a finite center and approximately flat outer velocity contribution, but needs an outer truncation to avoid infinite mass. The task is to derive rho_0, r_c and that truncation from illumination, transport, capture and feedback. The recovered cored fit only establishes a useful shape target.

A diffusion limit may yield partial_t u_d=div(D grad u_d)+q_cap-u_d/tau. D must follow from a mean free path/interaction and this approximation must be valid. With only a boundary source, constant D and a leakage term, the regular stationary interior solution is proportional to sinh(r/ell)/r, ell=sqrt(D tau). It is not automatically the desired cored inverse-square profile. This reproduces why edge-only deposition is a substantive problem.

For full gravity, solve the metric and every sector's stress, including free radiation and confinement. In a weak-field metric with potentials Phi and Psi, slow matter probes Phi and light deflection probes the line-of-sight integral of grad_perp(Phi+Psi)/c². The action fixes the relation between these potentials. Rotation curves alone cannot supply it.

The promising empirical target is

\[
g_X=A a_0(g_b/a_0)^q[\sqrt{2rg_b}/(200\,\mathrm{km/s})]^m,
\]

with recovered approximate A=0.699, q=0.3131, m=0.2651 and a0=1e-10 m/s². None of these is yet derived from photon production. The explicit center/radius is not a general local covariant scalar. Recover this behavior, or a comparably predictive alternative, from the solved fields. If amplified gravitational response is needed, label that branch as modified gravity, retain its driver energy and derive lensing and stability. A fitted amplification is not extra supplied energy.

### 3.9 Close the cosmic background

Specify physical scale factor, spatial curvature, topology, matter frame and source history. Setting a coordinate scale factor constant does not prove physical nonexpansion. For the ordinary Einstein-gravity homogeneous/isotropic subclass, a static physical solution must satisfy both background constraints. With energy density u and pressure p,

\[
Kc^2/a^2=8\pi G u/(3c^2)+\Lambda c^2/3,\qquad
0=-4\pi G(u+3p)/(3c^2)+\Lambda c^2/3.
\]

This is only the Einstein subclass; modified-gravity equations must be derived afresh. Perturb the complete solution, including radiation and matter, to test homogeneous, scalar, vector and tensor modes. Algebraic force balance is not stability.

Closed paths allow repeated encounters, not unlimited energy: a photon retaining fraction s per circuit has transferred at most E0(1-s^N) after N circuits. In the recovered n-dot=gamma n^p branch with p>0, fixed geometric separation and v=c/n, the future geometric path tends to c/(gamma p)=1/(kappa p), approximately 49 billion light-years for the fitted parameters. This is a conditional extrapolation, not a measured horizon or a general theorem about closed universes.

A complete cosmology also needs finite fuel/replenishment, entropy accounting, radiation backgrounds, a chronology and perturbation formation. “The universe is old” is an input to test, not an unlimited energy source.

### 3.10 Optional no-loss, permanent-deposit branch

The user-suggested branch sets companion secondary loss and deposit leakage to zero while transferring all stipulated photon losses into companions. It conserves the three-reservoir energy plus stellar source under its stated fixed-volume assumptions. The companion file `no-loss-branch.md` derives impulse and continuous-injection solutions, the growth consequence of permanent deposits, and the required recalculations. Keep it as an optional candidate; conservation is mandatory in every branch. A shared-loss model with an explicit receiving sector can also conserve energy. T20 benchmarks these alternatives before they are used in a halo fit.

## 4. Research gates and completion standards

| Gate | Required result | What happens on failure |
|---|---|---|
| G0: reliable benchmark | Recover and reproduce the numerical baselines, classify data assumptions and freeze the protocol. | Fix reproduction/data errors before interpreting new fits. |
| G1: a viable mechanism | Define the graviton/deposit; derive an allowed conversion channel, a clock/pulse-consistent transfer law and a plausible capture mechanism. | Reject or explicitly revise the candidate assumptions. Record the scope of the obstruction. |
| G2: supply and dynamics | Close energy/momentum/entropy accounting, solve transport and gravity, and demonstrate local/background viability. | A good rotation fit cannot advance an inconsistent branch. |
| G3: predictive astronomy | Use one frozen parameter set for galaxy dynamics, independent lensing and cluster/environment tests. | Diagnose failure with predeclared alternatives; do not quietly add object-specific gains. |
| G4: complete universe | Explain background radiation, structure, chronology and global stability in the same model. | Retain a limited galactic model only if that is its demonstrated scope. |
| G5: discriminating validation | Freeze and evaluate new predictions with reproducible artifacts and full uncertainty accounting. | Report the falsification or remaining non-identifiability honestly. |

A task can succeed by deriving a viable candidate, by demonstrating that a specified candidate fails, or by proving that additional data are needed. A branch is not accepted merely because all its tasks were executed.

Do not demand that every free constant come from pure mathematics. A first-principles model can have a small number of empirically measured universal constants. What must be derived is how those constants control conversion, spectra, capture and gravity; what must be avoided is a separate fitted knob for each desired observation.

Define numerical tolerances before inspecting new results. Starting engineering proposals are a 1e-8 relative energy-ledger residual in well-conditioned analytic toy tests and less than 1% change in declared observables when resolution is doubled; test suitability and replace them with justified problem-specific tolerances at T01/T09. Observational acceptance uses source likelihoods/covariances and a frozen model-comparison rule, not those engineering tolerances or an arbitrary “looks close” threshold.

## 5. Codex execution contract

Use the accompanying `codex-task-backlog.md` or JSON queue. They contain bounded prompts, prerequisites, deliverables and completion criteria. They are plans; no autonomous runs or recurring jobs have been scheduled.

For every task:

- Work in a new task-specific directory/checkout; preserve original sources and archived outputs. Put regenerated results in a new run directory and record code/data hashes, dependencies, seeds and units.
- Read the prerequisite reports and honor their conditional scope. Label every statement as assumption, derived identity, numerical result, fit, observational comparison or open question.
- Derive before optimizing when the task concerns a mechanism. Include dimensional checks and limiting cases; use independent formulations where a subtle derivation warrants it.
- Never infer input energy from the held-out velocity/lensing target. No per-object capture/gain fit may masquerade as a predicted halo.
- Include negative controls, nuisance parameters, calibration and selection effects. Reused splits are exploratory. Reserve new data before trying candidate variations; use redacted/blinded targets where feasible.
- State an explicit verdict: pass to next gate, reject this branch, or blocked with the exact missing input. Do not repair a failure by silently changing the physical model.
- Do not run downloaded pickle files as code. For archived serialized data, inspect provenance and use an isolated trusted-data conversion workflow before deserializing.
- Publish a short decision record with equations, reproducible command, machine-readable results, assumptions, failures and the next authorized task. A stalled symbolic derivation is not a reason to invent a coupling.

Start with T01–T03 and T20 as their prerequisites allow. Follow with T04/T05/T06/T07 as the specification permits. Expensive halo fitting and cosmological simulation should wait until their mechanism and normalization inputs exist. Independent workstreams can be assigned separately later; this request does not launch them.

## 6. Primary references for the constraint register

The local reports are the source for project-specific numbers. External references checked for this roadmap are examples of primary constraints and methods, not an exhaustive new literature review.

- Photon/graviton conversion and background-field assumptions: [Hwang & Noh](https://arxiv.org/abs/2310.04150).
- Photon/tensor propagation and multimessenger timing: [GW170817 / GRB 170817A](https://arxiv.org/abs/1710.05834). A candidate must model emission delay and distances rather than adopt a speed bound without its assumptions.
- Distance reciprocity: [More et al.](https://arxiv.org/abs/1612.08784).
- Rotation-curve inputs: [SPARC](https://arxiv.org/abs/1606.09251).
- Stellar/dust luminosity inference: [DustPedia analysis](https://arxiv.org/abs/1903.05933). Avoid counting dust-reprocessed starlight twice.
- Merging-cluster lensing offsets: [Clowe et al.](https://arxiv.org/abs/astro-ph/0608407).
- Atom/cavity comparisons: [Kennedy et al.](https://arxiv.org/abs/2008.08773).
- CMB spectra and model-dependent parameter inference: [Planck 2018](https://arxiv.org/abs/1807.06209).

Before a quantitative constraint task begins, retrieve the applicable source tables, methods and current likelihood documentation. Abstract-level source checks are not a substitute for the planned data analysis.

## Active user direction

The user has asked us to begin this program and keep fictional-universe alternatives open. Where the original roadmap says “reject,” treat this as documenting a limitation of the stated assumptions and opening a user decision about revisions, not permanently closing the conceptual path. Ask in plain language before choosing a consequential change of physical laws. Energy conservation and honest reporting remain mandatory.
