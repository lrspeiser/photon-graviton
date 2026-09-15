# CC-2, stage 2A: complete-channel formation and supply revision

Declared before execution, 14 September 2026. Baseline: `main` at 299ba1f. The stage-1 review corrections are committed separately, after this protocol.

**Why it is run.** Requested by the project owner's review of 299ba1f. Stage 1 took CF-1's seeded-only density as an input and found that an omitted channel of the same law, collisions between two incoming companions, dominates there. It did not solve again for the density the complete mechanism needs. The owner: "Before fixing stage 2's abundance requirements, recompute supply" with all three collision classes included together.

This stage produces the requirements that a source (stage 2B) must meet. It is not a fit and it proposes no source.

## Question

With incoming–incoming formation, incoming–bound capture and ejection, and bound–bound evaporation acting together under one law, and with mass, energy, orbits, bath depletion and gravity evolving together:
1. How much incident material does each system need to form a given bound population in the benchmark time?
2. Where does that population sit, and what energy state does it reach?
3. Is the approximation still applicable at that supply?
4. Does the seed matter, and does self-gravity?

## Inputs

- **Systems.** CF-1's baryon mass profiles, unchanged: the Milky Way (baryon model I, boundary R_b = 300 kpc), the J1630+4520 lens host (Sérsic stars at the Auger Chabrier mass, 300 kpc), and Coma at both ends of its bracketed gas-and-star inputs (30 Mpc). The baryons are fixed and do not respond to the companions (labeled).
- **The interaction.** CF-1's law, applied to every pair: elastic, equal-mass, isotropic in the centre-of-mass frame, constant σ/m. Trial scales σ/m = 0.1, 1, 10, 100 and 1,000 cm²/g.
- **Incident distributions at R_b** (declared trials, not predictions; the source's own spectrum is stage 2B's job):
  - **S:** a single speed u = 300 km/s, isotropic (CF-1's trial value, kept for continuity);
  - **M:** the isotropic Maxwellian with the same mean kinetic energy, one-dimensional dispersion u/√3 ≈ 173 km/s. This is the distribution a single-speed bath relaxes to under its own elastic collisions.
- **Bath transport.** The unbound population follows the transparent-limit (Liouville) distribution in the current potential, including the centrifugal-barrier limit CF-1 validated. It is depleted to first order by exp(−τ_in(r)), with τ_in the optical depth along the radial inflow path from R_b, counting collisions with both the bath and the confined companions. The distribution is recomputed as the potential changes and is replenished at R_b by the fixed boundary distribution (a steady external supply).
- **Bath gravity.** The bath's excess over the incident density, ρ_bath(r) − ρ_∞, gravitates. The uniform incident density's own gravity is assigned to the cosmological background: a labeled assumption, since in a static, uniform medium that term is not defined locally. The uniform part's enclosed mass is still reported, as an environmental estimate.
- **Initial state.**
  - The bath is present at t = 0 with its steady distribution in the baryonic potential. That is a declared initial distribution; inflow from R_b does not establish it within the benchmark time for Coma (stage 1's transport benchmark uT).
  - The bound population starts empty (primary). The control starts from CF-1's counted seed: 1% of the baryons, distributed like them, with the isotropic Jeans dispersion.
- **Benchmarks** (labeled, not observational requirements):
  - T = 10 Gyr.
  - **B1:** confined companion mass inside R_b equals the baryonic mass inside R_b (CF-1's benchmark).
  - **B2:** companion mass enclosed within the baryonic half-mass radius r_half equals the baryonic mass enclosed there.
  - The cosmic mean matter density for H0 = 70 and Ωm = 0.3 is a comparison unit only.

## Method

**C0. CF-1 regression.** CF-1's growth rule, with its archived kernel tables, must reproduce its archived required exposures A*.

**C1. Fixed-background control** (prescribed shape; labeled a control, not a formation prediction):

dM/dt = S(ρ_∞) + (σ/m)ρ_∞ K M − Λ M²,

where:
- S is the seedless rate from the transparent bath in the baryonic potential;
- K is CF-1's seeded kernel;
- Λ is stage 1's evaporation kernel, for a population shaped like the seed.

C1 is solved for the ρ_∞ that meets B1. It also gives the frozen-coefficient velocity thresholds of the whole process: the sum of all three channels for Maxwellian dispersions of 100, 200, 300, 500 and 1,000 km/s, and for the 300 and 3,000 km/s mixture.

**Coupled evolution.**
- **Tracers.** Monte Carlo companions of equal mass m_p, followed on exact orbits in the evolving spherical potential of baryons, confined companions and bath excess. Each tracer has its own adaptive time step; the integrator is compiled.
- **Collisions** are sampled at a global step Δ, at most 5 Myr for galaxies and 50 Myr for Coma, and reduced so that every collision probability per step stays below 0.2.
  - **Incoming–incoming.** The production rate per radial shell is ½(σ/m)∫ρ_bath²⟨g·P_conf⟩dV, sampled from the local bath distribution. New tracers take the sampled post-collision velocities, and the partner's escaping energy is booked.
  - **Incoming–bound.** Each tracer meets a bath companion drawn from the local distribution with probability (σ/m)ρ_bath·g·Δ. Both fates follow the confinement test, and a confined bath partner becomes a new tracer.
  - **Bound–bound.** Radial-neighbour pairing, valid by spherical symmetry, with probability (σ/m)ρ_c·g·Δ, where ρ_c is the local confined density.
- **Confinement.** "Confined" means unable to reach R_b in the current potential: E < 0, or E ≥ 0 behind the centrifugal barrier. Escaping companions are removed and their energy is booked.
- **Potential.** Updated every global step from the tracers' radii. Each tracer's specific-energy change from the update is booked as work done by the evolving potential.
- **Tracer count.** m_p is set per run so that about 10⁴ tracers represent the benchmark mass. Above 4×10⁴ tracers, half are removed at random and the rest double in mass.

**The required density.** For each system, σ/m and incident distribution:
- zero-seed runs at C1's estimate ρ_1 and at 3ρ_1 and ρ_1/3;
- log–log interpolation of M(T) for B1 and B2, with up to three further runs if a target lies outside the bracket;
- a verification run at the interpolated B1 density ρ_req;
- two controls at ρ_req: the counted seed, and a frozen potential (no self-gravity).

Every run has a wall-clock budget. A run that exhausts it is reported as not completed.

## Outputs (per system, σ/m and incident distribution)

1. **Supply.**
   - ρ_req for B1 and B2, in M_sun/kpc³ and in units of the cosmic mean.
   - A_req = (σ/m)ρ_req.
   - The ratio to CF-1's old requirement A*_old/(σ/m), and C1's estimate.
2. **Regime at ρ_req.**
   - The optical depth τ_in to r_half and to the centre.
   - The bath's self-collision time against T.
   - The bath excess mass within r_half and at test radii (the Milky Way at 25 kpc; Coma at 1, 2 and 3 Mpc), relative to the baryons there.
   - The uniform part's mass inside R_b (environmental estimate).
   - Capture efficiency: confined mass gained per mass of inflow through R_b.
   - The bath's Jeans diagnostic (an instability diagnostic).
3. **Formation history.**
   - M(t), in total and within r_half.
   - Event rates per channel: seedless births, captures, ejections, evaporations.
   - Zero seed against counted seed; coupled against frozen potential.
4. **State at T.**
   - The enclosed-mass profile and its log slope across the tested range.
   - Density, velocity dispersion and anisotropy profiles.
   - The energy distribution, and the speed distribution at r_half, tail included.
   - The companions' acceleration and projected mass within r_half and at the test radii, relative to the baryons.
   - Whether the population is still growing at T (d ln M/d ln t).
5. **Energy.**
   - The confined population's specific-energy ledger: births, removals, collisions and potential work.
   - The energy exported by escaping companions, and its spectrum at R_b.

## Validation (tolerances declared)

- **V1. No interaction.** With σ/m = 0, the counted seed in a static potential for T:
  - no mass is gained or lost;
  - the per-tracer energy error, relative to |Φ(r_half)|, is below 10⁻³ at the 99th percentile;
  - angular momentum is conserved to 10⁻⁹.
- **V2. Bath sampler.** The code's local bath density equals CF-1's orbit-integrated density (`focus.density_by_orbits`) to 10⁻³ at five radii, one of them barrier-limited. Sampled speeds equal w(r), and half the sampled directions point inward within three standard errors.
- **V3. Seedless production.** Where no barrier applies, the sampled production per shell equals stage 1's analytic kernel ½ρ_in²⟨g·P⟩ to 2%.
- **V4. Incoming–bound.** In a frozen potential, with the other channels off, the counted seed's initial net retention rate equals CF-1's (σ/m)ρ_∞KM_seed within 5%.
- **V5. Collisions.** Energy and momentum are conserved in every collision to rounding.
- **V6. Ledgers.** The mass ledger is exact. The specific-energy ledger closes to 10⁻⁹ relative.
- **V7. Potential.** The baryons-only potential equals CF-1's `Potential.phi` to 10⁻⁴ relative. A sampled Plummer sphere's enclosed mass matches the analytic profile within Poisson error.
- **V8. Resolution.** For the Milky Way at σ/m = 1 cm²/g (incident distribution S), doubling the tracer count and halving Δ changes ρ_req(B1) by less than 10%.

## Declared labels

- **"Seedless formation lowers the supply substantially"** at a given σ/m if ρ_req(B1) ≤ 0.1 × A*_old/(σ/m).
- **"Transparent"** if τ_in(r_half) < 0.1 throughout the run. **"Outside the model's regime"** if it exceeds 0.3.
- **"Bath gravity minor at the tested radii"** if the bath excess within r_half is below 10% of the baryons there.
- **"Bath self-relaxed"** if its self-collision time is shorter than T. A single-speed bath is then inconsistent, and case M applies.
- **"Seed-independent"** if the zero-seed and counted-seed M(T) differ by less than 10%.
- **"Self-gravity matters"** if the coupled and frozen-potential M(T) differ by more than 10%.

## Not claimed

- A source, or its spectrum (stage 2B).
- Survival beyond the evolved interval. Numerical verification is not proof of survival.
- A fit to motions or lensing (CR-3).
- Any prescribed final profile. The profile is an output.
- Bath–bath scattering that produces no confined companion enters only through first-order depletion; its thermalization of the bath is represented by case M, not evolved.
