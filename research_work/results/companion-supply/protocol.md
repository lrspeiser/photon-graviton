# CC-2, stage 1: regime and supply audit of CF-1's companion bath

Declared before execution, 14 September 2026. Baseline: `main` at f13c09f.

**Why it is run.** Requested by the project owner's review of f13c09f: "Proceed with CC-2, beginning with a CF-1 regime and supply audit." This stage checks whether the population CF-1 assumed is physically consistent at the supply it requires. It also works out what the omitted channels of the same interaction do.

It produces requirements for stage 2 (a declared source interaction) and stage 3 (a small coupled formation test). It predicts no reservoir.

## Question

At the supply CF-1 requires, can the incoming companions be what CF-1 assumed? That means:
- transparent, single-speed and isotropic;
- not colliding with each other;
- steady out to the exterior boundary R_b.

And what change when two things are added: the collisions CF-1 left out under its own law, and realistic velocity distributions?

## Inputs (unchanged from CF-1)

- **Systems.** The Milky Way and Coma at both ends of its bracketed inputs. The six lens hosts are reported where the cost is small.
- **Archived material.** CF-1's potentials, seed profiles, kernel tables K(u) and required exposures A* = (σ/m)ρ_∞ for retained mass equal to the baryons in 10 Gyr.
- **The interaction.** Elastic, equal-mass, isotropic in the centre-of-mass frame, with constant σ/m. Here it is applied to every pair of companions, as the law implies. Leaving a channel out would need a stated physical reason.
- **Labeled assumptions and comparison units.** The 10 Gyr reference time is an assumption. The cosmic mean matter density is a comparison unit only.

## Audit items

- **A1. Bath scattering.**
  - The characteristic length 1/A*, and the mean free paths for a single-speed bath, 3/(4A*), and a Maxwellian bath, 1/(√2 A*).
  - The bath's optical depth across R_b and across the half-mass radius.
  - The transparent assumption holds only where the optical depth is much less than one. At fixed A*, trading σ/m against ρ_∞ leaves these lengths unchanged.
- **A2. Bath gravity.** For σ/m = 0.1, 1, 10, 100 and 1,000 cm²/g, all trial scales, the implied ρ_∞ = A*/(σ/m) gives:
  - the bath mass inside R_b relative to the system's baryons;
  - the bath's Jeans length u√(π/(Gρ_∞)), compared with R_b.
- **A3. Travel and equilibration.**
  - The distance uT a companion travels in the reference time, compared with R_b.
  - The crossing times of R_b and of the half-mass radius.
  - The collision time 1/(A*u), compared with T.
- **A4. Particle against wave.** The de Broglie wavelength for the archived constituent masses, compared with the half-mass radius and the mean free path, plus the smallest companion mass for which the particle treatment holds at each system and trial speed.
- **A5. Velocity distribution.** CF-1's kernels are linear in the incident distribution, which gives:
  - net retention for mixtures of 300 and 3,000 km/s companions, with the fast fraction that reverses retention;
  - net retention for Maxwellian baths with one-dimensional dispersions of 100, 200, 300, 500 and 1,000 km/s, integrated over the tabulated kernel.

  CF-1's single-speed results are kept as controls.
- **A6. The omitted channels, under the same law:**
  - **Incoming–incoming (seedless).** Two unbound companions arriving at the same point with speed w = √(u² + v_esc²) and independent random directions. The probability that the collision leaves one of them bound. Energy conservation allows at most one. Then the formation-rate kernel ½(σ/m)ρ_in²⟨g·P⟩, integrated over each system with ρ_in = ρ_∞√(1 + v_esc²/u²), compared with CF-1's seed-driven rate at the required exposure.
  - **Bound–bound (evaporation).** The probability that a collision between two seed companions sends one past escape, with each drawn from the seed's Maxwellian truncated at v_esc. Then the evaporation-rate kernel, compared with seed-driven capture.
  - **Reverse processes.** Bath companions scattered by other bath companions (attenuation, from A1), and the ejections already inside K, are both reported.

## Validation (tolerances declared)

- **V1.** For two equal-speed companions, the analytic seedless capture probability
  P(φ) = max(0, 1 − u²/(w² sin φ)),
  with φ the angle between the incoming directions, averaged with the collision-rate weight, is reproduced by Monte Carlo to 1% at five well depths.
- **V2.** A Maxwellian bath narrowed around a speed reproduces CF-1's single-speed kernel there to 10⁻³.
- **V3.** Energy and momentum are conserved in every sampled collision to rounding.

## Outputs

- A regime map for each system: which assumptions of CF-1 hold at the required supply, for which σ/m.
- The velocity-distribution thresholds.
- The seedless formation and evaporation rates against seed-driven capture.
- A list of requirements that a stage-2 source must meet: density, velocity distribution and tail, interaction regime, and time available.

## Not claimed

- A source mechanism.
- Any reservoir, or its survival.
- Any fit.
