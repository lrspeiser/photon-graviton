# CF-1: gravitational focusing and capture of a slow companion population (conditional diagnostic)

Declared before execution, 14 September 2026. Baseline: `main` at e5f6fc3.

**Why it is run.** Requested by the project owner's review of e5f6fc3: "Proceed with CF-1 as a conditional transport-and-capture diagnostic."

## Question

Can one common population of slow companions and one declared interaction law behave differently in shallow galaxy potentials and deep cluster potentials? Both focusing and retention are compared, with no per-object adjustment.

This is a conditional diagnostic:
- the incoming population is a declared trial bath;
- results are reported per unit incident density;
- supply requirements are reported separately.

It does not predict reservoir masses or cluster lensing.

## Four quantities kept distinct

Consider a bath of companions with speed u at the exterior boundary.
- **Entry rate.** The flux reaching radius R. For a single speed, σ_enter = πR²(1 + v_esc²(R)/u²) when no centrifugal barrier outside R blocks entry.
- **Unbound density.** For a transparent potential and an isotropic single-speed bath, ρ_unbound(r)/ρ_∞ = √(1 + v_esc²(r)/u²).
- **Net capture.** Incoming companions that become bound through the declared interaction, minus bound companions it ejects.
- **Retained mass.** Net capture integrated over a declared time, with potential feedback.

## Setup (declared)

The potentials are Newtonian, nonrelativistic and spherical. That is a control, not a law of the fictional universe.

**Potentials.** Each is built from ordinary matter plus an explicitly counted seed, never from the inferred missing mass. Escape speed is measured relative to an exterior boundary R_b, where the bath is specified. R_b is ten times the baryonic extent: 300 kpc for the galaxies and 30 Mpc for Coma. Escape speeds relative to infinity are also reported.
- **The Milky Way.** Baryon model I, spherically averaged, from the archive.
- **Six massive lens hosts.** The SLACS lenses' published Sérsic stars with their Auger Chabrier masses, in the FLRW comparison geometry. These are adopted inputs.
- **Coma.** The repository has no matched gas-and-star model, so the inputs are carried explicitly:
  - **Gas shape.** A β-model with β = 0.75 and r_c = 10.5′ ≈ 296 kpc for h = 0.7 at z = 0.0231 ([Briel, Henry & Böhringer 1992](https://ui.adsabs.harvard.edu/abs/1992A&A...259L..31B/abstract)), with μ_e = 1.17, truncated at 3 Mpc.
  - **Gas normalization, bracketed.** Central electron density 2.5–4.5 ×10⁻³ cm⁻³. The bracket spans the commonly quoted 3.36×10⁻³, not verified here, and the Planck SZ-fitted 2.7×10⁻³ ([Planck Collaboration Int. X 2013](https://ui.adsabs.harvard.edu/abs/2013A&A...554A.140P), Table 2). That paper itself warns its fitted density and temperature profiles are strongly correlated and unlikely to describe the actual cluster.
  - **Stars, a bracketed placeholder.** M* = 0.5–2 ×10¹³ Msun, distributed like the gas. Not a measured input.

  Coma is run at the low end of both brackets and at the high end. Both brackets are to be replaced by a matched gas-and-star model.
- **A separately labeled inverse diagnostic.** An NFW total potential of literature scale (M₂₀₀ = 10¹⁵ Msun, c = 4). It shows focusing and retention in a full well. It is not a formation calculation and not evidence that focusing makes that well.

**Bath.** Isotropic at R_b, with single speed u. The trial values are u = 300, 1,000 and 3,000 km/s. The 3,000 km/s value was chosen after considering the cluster–galaxy contrast, so it is a trial, not a prediction.

**Wave check.** The de Broglie wavelength is compared with each system's half-mass radius. The particle treatment is declared valid only where the wavelength is much smaller. At 3,000 km/s the archived constituent masses give about 3 kpc (1.34×10⁻²⁴ eV) and 40 pc (10⁻²² eV). Otherwise the companion mass is left general, and results are per unit mass.

## Validation (no capture), tolerances declared

- **V1.** Orbits integrated in each fixed potential reproduce σ_enter at the half-mass radius to 10⁻³ relative, for each trial speed. The critical angular momentum is found by bisection on the integrated pericentre.
- **V2.** The time-weighted density from the orbits reproduces √(1 + v_esc²/u²) to 10⁻³ at 0.1, 0.3, 1, 3 and 10 half-mass radii.
- **V3.** With capture disabled, sampled incoming orbits in the stationary potential all leave, so the retained fraction is zero.
- **Interaction check.** For stationary targets, the per-event probabilities are capture = v_esc²/w², ejection = u²/w² and net = (v_esc² − u²)/w², with w² = u² + v_esc². The Monte Carlo reproduces these to 1%, with energy and momentum conserved to rounding.

## One declared interaction

**The law.** Elastic, equal-mass, companion–companion scattering, isotropic in the centre-of-mass frame, with a constant cross-section per unit mass σ/m. It acts between incoming companions and an explicitly counted bound seed.

**The seed.** In every system it is 1% of the baryonic mass, distributed like the baryons, with the isotropic Jeans dispersion in the host potential. That one rule holds for all systems.

**Per event,** energy and momentum are conserved, and the targets' velocities are sampled. Tallied weighted by relative speed:
- capture of the incomer;
- ejection of the target;
- both, and the net;
- the energy the bound population gains and the energy that escapes.

**Kernels.** Per unit (σ/m)ρ_∞:
- the growth kernel K = ∫ (ρ_unbound/ρ_∞)⟨|v_rel|·net⟩ ρ_seed dV / M_seed;
- the matching heating kernel for the bound population's energy.

**Potential feedback.** The retained mass keeps the seed's shape, so the potential scales by (1 + q), with q the retained mass over the baryonic mass. All speeds then scale by √(1+q), and exactly K(q; u) = √(1+q) K(0; u/√(1+q)). The retained mass is integrated over a reference time of 10 Gyr, labeled an assumption, for a range of exposures (σ/m)ρ_∞T. Growth, erosion or runaway is reported.

**Energy, and what receives it.** Elastic scattering conserves the companions' total energy, and each incomer arrives with +u²/2 per unit mass relative to R_b. The energy the bound population gains per net retained companion is reported against the seed's specific binding energy. Without a declared energy sink, net capture heats the reservoir; the report states whether a reservoir could stay bound. CF-1 adds no sink.

## Declared outputs

- For each system and trial speed:
  - escape speeds;
  - entry and density factors against radius;
  - per-event capture, ejection and net probabilities;
  - the growth and heating kernels;
  - the retained-mass history against exposure.
- The contrast in K between Coma, at both bracket ends, and the Milky Way and lens hosts.
- **Supply, reported separately.**
  - The exposure needed for the retained mass to equal the baryonic mass in 10 Gyr, which is a reference scale only.
  - That exposure converted to an incident density at σ/m = 1 cm²/g, a trial scale, in units of the cosmic mean matter density for H0 = 70 and Ω_m = 0.3. The unit is for comparison only.
- The inverse diagnostic, labeled.

## Not claimed

- That such a population exists or is supplied (CC-2).
- That focusing makes the cluster excess.
- Any lensing prediction.
- Any fit. Nothing is adjusted to data.
