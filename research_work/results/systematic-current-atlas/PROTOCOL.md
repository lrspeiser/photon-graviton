# JR-7: systematic companion-current atlas and observable-characteristics census

Date: 21 September 2026, America/Los_Angeles. Baseline inspected: f55c07df0f2da42b6a2d4d49c0b73975ad17d4ec. This is exploratory construction, not exhaustive search of all physics or a discovery of a valid gravity theory. The protocol is declared before the new numerical campaign.

## Theory first

Matter and radiation may generate an extended companion state; structured sources, gas, relative currents and memory can reorganize its gravitational response. Preserve the broader branch portfolio. This campaign extends JR-6's particular positive two-channel transport implementation, not a requirement that all branches use a local Markovian diffusion law or one metric. Keep measured inputs and existing R10 predictions unchanged.

## Physical question

Can uneven gas/source overlap and transport change mean companion occupation as well as its directional pattern, under one reciprocal exchange rule? Which interactions control the sign? Which measured galaxy or cluster characteristics could actually identify those rules, as opposed to unconstrained per-object guesses?

## Coupled spatial construction

In a frame where the declared source and gas patterns have fixed relative orientation, solve

 P_t + Omega_P P_phi = D_P lap(P) - k_plus P + k_minus C - lambda_P P + J,
 C_t + Omega_C C_phi = div(D_C grad C) + k_plus P - k_minus C - lambda_C C.

Cylindrical R,phi,z grid, periodic phi and reflecting radial/vertical boundaries. Angular transport uses a positivity-preserving exponential-fitted finite-volume flux. P and C are coarse-grained energy densities; this is not propagation of real photons at c or a derived entrainment law. Escaping channels and source fuel have explicit energy accounts. This does not yet close recoil, gas dynamics, gravitational binding energy, or relativistic causality.

Freeze total gas normalization and injected luminosity across geometry variations. Four rate families: constant exchange; density-enhanced forward conversion; density-enhanced return; mixed forward/return. All rates are positive, k_plus=0.5 exp(a h), k_minus=0.5 exp(b h^2), h=rho_g/(rho_g+0.15), with (a,b)=(0,0),(2,0),(0,4),(2,4).

Full factorial varies relative source/gas phase 0/45/90 degrees, angular companion drift -0.6/0/+0.6 in the pattern frame, smooth/structured gas, vertical thickness 0.35/0.7, companion loss 0.05/0.3, and vertical/horizontal diffusivity ratio 0.25/4. D_P=0.5, horizontal D_C=0.05, lambda_P=0.5, Omega_P=-0.6; source luminosity is one in toy units. The exact smooth spatial shapes, normalization, grids and all case parameters are written by the executable before each run. 576 primary factorial cases are intended. Select matched contrasts by a deterministic statistic, then refine both members; do not relabel coarse values as converged. Include phase flips, angular averaging, zero-structure controls and independent energy-balance reconstruction. The phase-locked pattern assumption is explicit; unbounded physical phase locking is not proved.

Read the same companion density and both-channel density through one fixed softened 3D potential and its weak-field projected lens integral in representative cases, with no independently fitted lensing multiplier. Distinguish source occupation, mean radial force, quadrupole, tangential force and projected lensing. Additional offset/tilted source controls, including a two-source cluster-like geometry, are synthetic sensitivity cases, not fits to real clusters.

## Observational characteristics

Use immutable JR-1 R10 predictions and the existing 149 SPARC names, raw profiles and catalog. Derive signed galaxy-level residuals, not new force-law parameters. Fit a small fixed collection of regularized/interacting regressions on the 89 training objects, select on the 29 development objects and report the 31 comparison objects only after selection. These galaxies were exposed previously, so no fresh blind validation. Features are catalog/source properties, never observed rotation amplitude or a residual-dependent reconstructed halo. Record correlated predictors, small sample size and nuisance-variable leakage. The exercise asks how much existing static metadata explains the residuals; it is not a new physical law or test of unmeasured currents.

Build a source-backed galaxy/cluster observation registry: available density/light maps, gas phases, line-of-sight velocities, geometry, lensing and uncertainty. Spectral cubes have two sky coordinates plus wavelength/velocity, not directly measured depth. Absence of a measured spin, thickness, clump or phase variable must remain missing, not be supplied from gravitational residuals. Cluster thermodynamic and lensing-derived mass products retain their model assumptions; do not adopt dark matter as an axiom or use inferred dark-halo masses as independent source inputs.

## Interpretation

No new free per-galaxy gravity amplitude or photon-specific factor. Finite numerical possibilities do not establish energy supply, actual matter-light coupling, observed cluster offsets, binary radiation or Solar-system behavior. A negative result restricts the stated implementation and assumptions, not every companion/current/memory mechanism. Save code, parameters, all primary results, finite-grid failures and input hashes. Stop claims at what was actually executed; no promise of background work.
