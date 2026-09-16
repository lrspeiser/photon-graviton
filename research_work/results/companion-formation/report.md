# CC-2 stage 2A report: seedless formation lowers the supply, but it builds an outer envelope, not an inner reservoir

All numbers are at 10 Gyr, for the labeled benchmarks B1 (confined companion mass inside R_b equal to the baryons there) and B2 (companion mass inside the baryonic half-mass radius r_half equal to the baryons there). "Within the regime" means the bath's optical depth to r_half stayed below 0.3; that holds in 33 of the 47 combinations that reached B1, all at σ/m of 0.1–10 cm²/g.

- **Supply.** With all three collision classes of CF-1's law acting together, the incident density that forms B1 falls well below CF-1's seeded-only requirement.
  - With the bath's gravity omitted, within the regime, it is 29–580 times lower in the Milky Way, 2.8–44 times in J1630 and 77–1,600 times in Coma.
  - The declared label, at least ten times lower, holds in 35 of the 47 combinations.
  - The supply is still large: 3,900–42,000 times the cosmic mean density in the Milky Way, 2,500–29,000 in J1630 and 19–220 in Coma. The cosmic mean is a comparison unit only.
- **Where the mass sits.** What forms is an extended envelope.
  - At B1 the companions inside r_half are 0.4–3% of the baryons there with the bath's gravity omitted, and 0.02–0.9% with it on.
  - Across the tested radii the enclosed companion mass rises as r^1.5 to r^2.7.
  - At the outermost test radius (25 kpc in the Milky Way, 21 kpc in J1630, 3 Mpc in Coma) companions add 0.5–6% of the baryons' enclosed mass, which is also their share of the acceleration, and 3–13% of the projected mass.
  - The population is still growing at 10 Gyr, with d ln M/d ln t of 1.0–1.4.
- **The inner benchmark.**
  - Making the interior companion-dominated (B2) takes 2.2–3.9 times the B1 density.
  - At that density the total confined mass is 16–47 times the baryons in Coma, 87–130 times in J1630 and 490–670 times in the Milky Way.
- **The bath's own gravity.**
  - **The static limit.** With its focused excess gravitating, the bath has a stable static state only below a limit density.
  - **Galaxies.** B1 lies below that limit in nine galaxy combinations at 0.1–10 cm²/g, and there it needs 0.16–0.67 of the gravity-omitted density.
  - **Coma.** No combination within the regime forms B1 with the bath's gravity on, because every such run loses its static bath: the boundary model's own limit.
  - **What that means for the Coma numbers.** With the bath's gravity omitted, its focused excess inside Coma's r_half is 0.26–5.5 times the baryons, so the Coma results rest on an omission that is not small.
- **Controls.**
  - Every combination within the regime is seed-independent: the counted seed changes M(T) by −1% to +5%.
  - Self-gravity matters in 29 of the 33: a frozen potential forms 0.59–0.98 of the coupled mass.
- **What leaves.** Each seedless birth's partner escapes. Those partners carry 1.0–1.7 times the confined mass away at 300–500 km/s. Ejections by incoming companions remove a further 0.01–0.75 of it in the galaxies and under 0.4% in Coma.

Protocol: [protocol.md](protocol.md), declared in 1ccd762 before execution. Code: [mc.py](mc.py) (engine), [formation.py](formation.py) (driver) and [checks.py](checks.py) (suite job). Results: [formation-results.json](formation-results.json).

**Why it was run.** The owner's review of 299ba1f. Stage 1 found that collisions between two incoming companions dominate at CF-1's seeded-only density, but it did not solve again for the density the complete mechanism needs. This stage does. All three collision classes of CF-1's law act together, and mass, energy, orbits, bath depletion and gravity evolve together.

**What is held fixed.**
- CF-1's baryon profiles and boundaries: the Milky Way and J1630 at 300 kpc, Coma at 30 Mpc.
- CF-1's interaction law.
- The incident distributions S (a single speed of 300 km/s) and M (the Maxwellian with the same mean energy).
- The benchmarks B1 and B2.

These are labeled trials and benchmarks, not observational requirements.

## How the calculation works

- **Tracers.** Companions are followed on exact orbits in the spherical potential of the baryons, the confined companions and the bath's focused excess, updated every step.
- **Bath.** The incoming population follows CF-1's transparent (Liouville) distribution in the current potential, including its centrifugal barrier. It is depleted to first order along the inflow and replenished at the boundary.
- **Collisions.** All three classes use the same law:
  - incoming–incoming collisions create confined companions, and the escaping partner's energy is recorded;
  - incoming–confined collisions capture or eject;
  - confined–confined collisions exchange energy and can evaporate a companion.

  Each collision conserves energy and momentum.
- **Controls.** CF-1's own growth rule (C0), a fixed-background rate model (C1), a counted-seed start, and a frozen potential.
- **Finding the requirement.** For each system, incident distribution and σ/m, a ladder of runs brackets the density that meets each benchmark. Verification runs at the interpolated density are repeated, up to three times, until they land within 25% of B1. The seed and frozen-potential controls run at the first verification density.

## Validation (declared tolerances)

| Check | Declared tolerance | Result | Passed |
|---|---|---|---|
| V1, no interaction (10 Gyr) | energy error < 10⁻³ of \|Φ(r_half)\| at the 99th percentile; no mass change; angular momentum 10⁻⁹ | energy 3.9×10⁻⁴; escapes 0; angular momentum 6.6×10⁻¹² | yes |
| V2, bath sampler | density 10⁻³ at five radii, one barrier-limited | 10⁻⁵ (barrier-limited radius 164 kpc) | yes |
| V3, seedless production | 2% of the analytic kernel | 8.1×10⁻³ | yes |
| V4, incoming–bound kernel | 5% of CF-1 | 0.014 (CF-1 recomputed); 0.03 (archived) | yes |
| V5, collisions | rounding | kinetic 1.3×10⁻¹⁵, momentum 3.3×10⁻¹⁶ | yes |
| V6, ledgers | mass exact; energy 10⁻⁹ | energy 3.3×10⁻¹⁴, mass 4.2×10⁻¹³ | yes |
| V7, potential | Φ 10⁻⁴; Plummer within Poisson error | Φ 4.2×10⁻⁵ (to 0.9 R_b); Plummer 4.5×10⁻³ against 0.012 | yes |
| V8, resolution | ρ_req(B1) within 10% | 0.075 | yes |
| C0, CF-1 regression | archived exposures | 1.6×10⁻⁶ | yes |

The resolution check V8 compares the Milky Way's B1 density (bath S, 1 cm²/g, bath gravity omitted) with a ladder run at twice the tracers and half the step. The two differ by 7.5%, inside the declared 10%.

240 kept runs: 38 thinned (largest mass change from thinning 3.9×10⁻¹⁶ of the final mass); largest resampling mass change 0.27; collisions capped for at least 1% of the mass, or pair probabilities above 1, in 10 runs: J1630 | M | 100 cm2/g | bath gravity omitted, J1630 | M | 1000 cm2/g | bath gravity omitted, J1630 | M | 1000 cm2/g | bath gravity on, J1630 | S | 100 cm2/g | bath gravity omitted, MW | M | 100 cm2/g | bath gravity omitted.

## The bath's own gravity: how dense a stable static bath can be

Before any tracer run, each combination was probed for the largest incident density at which the bath, with its focused excess gravitating, still has a stable static state (bisection to 2%). The table gives the limits in M_sun/kpc³ for S / M. The letter after each limit is what fails above it:
- **B, barrier onset.** The excess pushes the circular speed at R_b above the incident speed, and a centrifugal barrier empties the interior.
- **D, deficit.** The net enclosed mass turns non-positive somewhere: the bath's deficit, from depletion at high σ/m or from barrier shielding, outweighs everything inside it.
- **R, runaway.** The largest positive eigenvalue of the focusing response reaches one.
- **F, no fixed point.** The iteration found no static state. It is treated as having none, whether the cause is a fold or a numerical failure.

| σ/m (cm²/g) | Milky Way | J1630 | Coma, low bracket | Coma, high bracket |
|---|---|---|---|---|
| 0.1 | 1.9×10⁵ D / 2.2×10⁶ R | 4×10⁵ B / 2.2×10⁶ R | 657 B / 227 R | 657 B / 227 D |
| 1 | 1.9×10⁵ D / 10⁶ D | 2.1×10⁵ B / 1.7×10⁶ R | 674 B / 227 R | 641 B / 184 R |
| 10 | 7.5×10⁴ D / 1.3×10⁵ F | 1.8×10⁵ D / 4.5×10⁵ F | 1,105 B / 219 R | 539 B / 219 R |
| 100 | 2.2×10⁴ D / 2.3×10⁴ D | 4.7×10⁴ D / 5.1×10⁴ D | 302 B / 173 R | 302 B / 175 R |
| 1,000 | 7,242 D / 7,332 D | 1.6×10⁴ D / 1.7×10⁴ D | 55 B / 60 F | 102 R / 67 R |

In Coma the bath's focused excess alone uses up the boundary's allowance, u²R_b/G. For S at 300 km/s that allowance is 6.3×10¹⁴ M☉, 4.4 times the baryons (item 13 below). The ladder runs with bath gravity on made 58 stops: 50 lost their static bath as companions formed, and 8 had none to begin with. None stopped for runaway capture, for its wall-clock budget, or on an error.

## Supply

Each system's table gives C1's estimate and the densities that meet B1 (M_sun/kpc³). The ratio to CF-1's old requirement is in brackets. Where B1 lies above the static limit, the entry gives what a run at 0.9 of the limit formed. The **check** columns give the last verification's M/M_b for B1 and M(<r_half)/M_b(<r_half) for B2; a B2 density whose check lies outside 0.8–1.25 is uncertain within its ladder. Labels: **T**, transparent; **O**, outside the regime; **U**, collisions not fully resolved.

**Milky Way**

| Bath | σ/m | C1 | ρ_B1, bath gravity on | ρ_B1, omitted | B1 check (on / omitted) | ρ_B2, omitted (check) | τ max (on / omitted) | labels (on / omitted) | seed/zero (omitted / on) | frozen/coupled (omitted / on) |
|---|---|---|---|---|---|---|---|---|---|---|
| S | 0.1 | 2.2×10⁶ | above the limit; 0.9× limit forms 0.17 M_b | 1.7×10⁶ (2.4×10⁻³) | – / 0.84 | 3.9×10⁶ (1.1) | – / 0.015 | – / T | 1.04 / – | 0.714 / – |
| S | 1 | 8.5×10⁵ | 1.6×10⁵ (2.3×10⁻³) | 6×10⁵ (8.3×10⁻³) | 0.92 / 0.95 | 1.3×10⁶ (1.1) | 0.015 / 0.051 | T / T | 1.04 / 1.03 | 0.697 / 0.807 |
| S | 10 | 7.3×10⁵ | above the limit; 0.9× limit forms 0.048 M_b | 2.4×10⁵ (0.034) | – / 0.97 | 5.3×10⁵ (1.4) | – / 0.19 | – / – | 1.03 / – | 0.587 / – |
| S | 100 | 7.3×10⁵ | above the limit; 0.9× limit forms 0.021 M_b | 2.7×10⁵ (0.37) | – / 0.91 | 8.1×10⁵ (0.024) | – / 1.2 | – / O | 1.05 / – | 0.504 / – |
| S | 1,000 | 7.3×10⁵ | above the limit; 0.9× limit forms 7.6×10⁻³ M_b | not reached | – / – | – (–) | – / – | – / – | – / – | – / – |
| M | 0.1 | 1.6×10⁶ | 2×10⁵ (2.8×10⁻⁴) | 1.2×10⁶ (1.7×10⁻³) | 1.1 / 0.85 | 3.3×10⁶ (1.3) | 2.4×10⁻³ / 0.011 | T / T | 1.02 / 1.01 | 0.686 / 0.976 |
| M | 1 | 5.9×10⁵ | 1.1×10⁵ (1.5×10⁻³) | 4.3×10⁵ (5.9×10⁻³) | 0.95 / 0.96 | 1.1×10⁶ (1.1) | 0.011 / 0.038 | T / T | 1.05 / 1.04 | 0.691 / 0.914 |
| M | 10 | 5×10⁵ | 7.5×10⁴ (0.01) | 1.6×10⁵ (0.022) | 0.82 / 1 | 3.8×10⁵ (1.2) | 0.07 / 0.14 | T / – | 0.992 / 1 | 0.597 / 0.747 |
| M | 100 | 5×10⁵ | above the limit; 0.9× limit forms 0.052 M_b | 9×10⁴ (0.12) | – / 0.87 | 2.7×10⁵ (11) | – / 0.63 | – / O | 0.97 / – | 0.525 / – |
| M | 1,000 | 5×10⁵ | above the limit; 0.9× limit forms 0.016 M_b | not reached | – / – | – (–) | – / – | – / – | – / – | – / – |

**J1630**

| Bath | σ/m | C1 | ρ_B1, bath gravity on | ρ_B1, omitted | B1 check (on / omitted) | ρ_B2, omitted (check) | τ max (on / omitted) | labels (on / omitted) | seed/zero (omitted / on) | frozen/coupled (omitted / on) |
|---|---|---|---|---|---|---|---|---|---|---|
| S | 0.1 | 1.5×10⁶ | above the limit; no static bath | 1.2×10⁶ (0.03) | – / 0.84 | 2.8×10⁶ (1.2) | – / 0.011 | – / T | 1.02 / – | 0.685 / – |
| S | 1 | 5.3×10⁵ | 1.5×10⁵ (0.037) | 3.9×10⁵ (0.096) | 0.89 / 0.82 | 8.6×10⁵ (0.86) | 0.016 / 0.036 | T / T | 1.04 / 1.03 | 0.683 / 0.831 |
| S | 10 | 4.1×10⁵ | 9.8×10⁴ (0.24) | 1.5×10⁵ (0.36) | 0.81 / 0.86 | 2.9×10⁵ (0.42) | 0.096 / 0.14 | T / – | 1.05 / 1.05 | 0.637 / 0.725 |
| S | 100 | 4.1×10⁵ | above the limit; 0.9× limit forms 0.15 M_b | 8.6×10⁴ (2.1) | – / 0.89 | 1.8×10⁵ (1.5) | – / 0.7 | – / O | 0.986 / – | 0.531 / – |
| S | 1,000 | 4.1×10⁵ | above the limit; 0.9× limit forms 0.031 M_b | not reached | – / – | – (–) | – / – | – / – | – / – | – / – |
| M | 0.1 | 1.1×10⁶ | 3×10⁵ (7.5×10⁻³) | 9.1×10⁵ (0.023) | 1.1 / 0.88 | 2.4×10⁶ (0.92) | 3.9×10⁻³ / 9.5×10⁻³ | T / T | 1.02 / 1.02 | 0.702 / 0.973 |
| M | 1 | 3.8×10⁵ | 1.4×10⁵ (0.034) | 2.9×10⁵ (0.073) | 1 / 0.88 | 7.7×10⁵ (1.1) | 0.016 / 0.031 | T / T | 1.03 / 1.03 | 0.681 / 0.911 |
| M | 10 | 2.8×10⁵ | 6.8×10⁴ (0.17) | 10⁵ (0.26) | 0.98 / 0.89 | 2.4×10⁵ (1.1) | 0.079 / 0.11 | T / – | 1.03 / 1.03 | 0.643 / 0.788 |
| M | 100 | 2.8×10⁵ | above the limit; 0.9× limit forms 0.6 M_b | 4.5×10⁴ (1.1) | – / 0.8 | 10⁵ (0.5) | – / 0.49 | – / O | 1.03 / – | 0.688 / – |
| M | 1,000 | 2.8×10⁵ | above the limit; 0.9× limit forms 0.069 M_b | 7.7×10⁵ (191) | – / 0.88 | 8.2×10⁵ (0.024) | – / 4.1 | – / OU | 5.68 / – | 0.543 / – |

**Coma, low bracket**

| Bath | σ/m | C1 | ρ_B1, bath gravity on | ρ_B1, omitted | B1 check (on / omitted) | ρ_B2, omitted (check) | τ max (on / omitted) | labels (on / omitted) | seed/zero (omitted / on) | frozen/coupled (omitted / on) |
|---|---|---|---|---|---|---|---|---|---|---|
| S | 0.1 | 1.1×10⁴ | above the limit; static bath lost | 8,835 (7.1×10⁻⁴) | – / 0.89 | 2.5×10⁴ (0.97) | – / 0.01 | – / T | 1.05 / – | 0.717 / – |
| S | 1 | 3,453 | above the limit; static bath lost | 2,821 (2.3×10⁻³) | – / 0.9 | 8,264 (0.9) | – / 0.032 | – / T | 1 / – | 0.705 / – |
| S | 10 | 1,086 | above the limit; static bath lost | 939 (7.6×10⁻³) | – / 0.94 | 3,053 (1) | – / 0.1 | – / – | 1.05 / – | 0.703 / – |
| S | 100 | 338 | above the limit; static bath lost | 339 (0.027) | – / 1 | 1,620 (1) | – / 0.33 | – / O | 0.996 / – | 0.677 / – |
| S | 1,000 | 102 | 49 (0.039) | 145 (0.12) | 0.99 / 0.93 | 6,648 (0.78) | 0.52 / 0.98 | O / O | 1.03 / 1.01 | 0.702 / 1.06 |
| M | 0.1 | 8,914 | above the limit; static bath lost | 7,640 (6.2×10⁻⁴) | – / 0.97 | 2.3×10⁴ (1) | – / 9.5×10⁻³ | – / T | 1.01 / – | 0.741 / – |
| M | 1 | 2,813 | above the limit; no static bath | 2,440 (2×10⁻³) | – / 0.97 | 7,709 (0.97) | – / 0.03 | – / T | 1.01 / – | 0.732 / – |
| M | 10 | 884 | above the limit; static bath lost | 795 (6.4×10⁻³) | – / 0.97 | 2,917 (1.2) | – / 0.095 | – / T | 1.01 / – | 0.761 / – |
| M | 100 | 274 | above the limit; static bath lost | 280 (0.023) | – / 1 | 1,602 (1.2) | – / 0.3 | – / O | 1.04 / – | 0.744 / – |
| M | 1,000 | 82 | above the limit; 0.9× limit forms 0.56 M_b | 120 (0.097) | – / 0.99 | 5,740 (0.45) | – / 0.92 | – / O | 1.02 / – | 0.737 / – |

**Coma, high bracket**

| Bath | σ/m | C1 | ρ_B1, bath gravity on | ρ_B1, omitted | B1 check (on / omitted) | ρ_B2, omitted (check) | τ max (on / omitted) | labels (on / omitted) | seed/zero (omitted / on) | frozen/coupled (omitted / on) |
|---|---|---|---|---|---|---|---|---|---|---|
| S | 0.1 | 9,415 | above the limit; static bath lost | 7,519 (1.2×10⁻³) | – / 0.9 | 2.4×10⁴ (1) | – / 0.01 | – / T | 1.01 / – | 0.713 / – |
| S | 1 | 2,968 | above the limit; static bath lost | 2,415 (4×10⁻³) | – / 0.9 | 7,755 (1.1) | – / 0.032 | – / T | 1.04 / – | 0.711 / – |
| S | 10 | 929 | above the limit; static bath lost | 795 (0.013) | – / 0.92 | 2,796 (1) | – / 0.1 | – / – | 1.04 / – | 0.705 / – |
| S | 100 | 285 | above the limit; static bath lost | 286 (0.047) | – / 1 | 1,403 (0.88) | – / 0.33 | – / O | 1.04 / – | 0.704 / – |
| S | 1,000 | 84 | above the limit; static bath lost | 122 (0.2) | – / 0.97 | 6,048 (0.59) | – / 0.98 | – / O | 1 / – | 0.684 / – |
| M | 0.1 | 8,435 | above the limit; static bath lost | 7,441 (1.2×10⁻³) | – / 0.99 | 2.3×10⁴ (0.9) | – / 0.01 | – / T | 1 / – | 0.788 / – |
| M | 1 | 2,657 | above the limit; static bath lost | 2,373 (3.9×10⁻³) | – / 0.97 | 7,714 (1) | – / 0.032 | – / T | 1.01 / – | 0.795 / – |
| M | 10 | 830 | above the limit; static bath lost | 779 (0.013) | – / 0.99 | 2,848 (1.1) | – / 0.1 | – / – | 1.04 / – | 0.79 / – |
| M | 100 | 253 | above the limit; static bath lost | 273 (0.045) | – / 1 | 1,573 (1.2) | – / 0.32 | – / O | 1.01 / – | 0.776 / – |
| M | 1,000 | 73 | above the limit; 0.9× limit forms 0.45 M_b | 118 (0.19) | – / 1 | 5,666 (0.39) | – / 0.98 | – / O | 1.02 / – | 0.789 / – |

- **C1 as a guide.** The fixed-background control C1 uses a prescribed seed-shaped population. At 0.1–1 cm²/g it overestimates the coupled B1 density by a factor of 1.2–1.4. In the galaxies at 10 cm²/g its evaporation term holds it on a plateau, and it overestimates by about 3.
- **The bath's gravity lowers the requirement.** Where a static bath exists, letting its excess gravitate deepens the potential, and B1 then needs 0.16–0.67 of the gravity-omitted density.
- **Outside the regime.** At 100–1,000 cm²/g the bath becomes opaque inside r_half (τ up to 4.1), so those rows lie outside the model's regime. The one combination labeled U, J1630 M at 1,000 cm²/g, also reached the sub-step cap: with 18% of its verification run's mass and 65% of its seed control's. Its seed/zero ratio of 5.7 reflects that, not the physics.

## Where the formed companions sit

Within the regime, at B1: the companion mass inside r_half over the baryons there; the enclosed companion mass's logarithmic slope across the tested radii; the growth rate d ln M/d ln t at 10 Gyr; and companions over baryons at each test radius. Enclosed ratios come first; projected (lensing) ratios are in brackets. In spherical symmetry the enclosed ratio is also the companions' share of the acceleration relative to the baryons. The last two columns are the bath's focused excess inside r_half over the baryons there, and the fraction of the r_half speed distribution above 0.9 of the local escape speed.

| Combination | inside r_half | slope | d ln M/d ln t | companions/baryons at the test radii, enclosed (projected) | bath excess inside r_half | tail |
|---|---|---|---|---|---|---|
| Milky Way S 0.1, omitted | 4.7×10⁻³ | 1.9 | 1.4 | 5 kpc 3.6×10⁻³ (0.01); 10 kpc 9×10⁻³ (0.023); 25 kpc 0.034 (0.075) | 0.02 | 0.3 |
| Milky Way S 1, omitted | 5×10⁻³ | 1.9 | 1.3 | 5 kpc 3.7×10⁻³ (0.011); 10 kpc 9.4×10⁻³ (0.025); 25 kpc 0.036 (0.083) | 6.2×10⁻³ | 0.22 |
| Milky Way S 1, on | 5.1×10⁻⁴ | 2.3 | 1.3 | 5 kpc 3.8×10⁻⁴ (2.5×10⁻³); 10 kpc 10⁻³ (6.7×10⁻³); 25 kpc 6.9×10⁻³ (0.03) | 2.3×10⁻³ | 0.12 |
| Milky Way S 10, omitted | 6×10⁻³ | 1.7 | 1 | 5 kpc 5.5×10⁻³ (0.012); 10 kpc 9.6×10⁻³ (0.024); 25 kpc 0.036 (0.079) | 1.5×10⁻³ | 0.11 |
| Milky Way M 0.1, omitted | 4×10⁻³ | 2 | 1.3 | 5 kpc 2.9×10⁻³ (9.7×10⁻³); 10 kpc 8.7×10⁻³ (0.023); 25 kpc 0.036 (0.078) | 0.025 | 0.27 |
| Milky Way M 0.1, on | 1.7×10⁻⁴ | 2.7 | 1 | 5 kpc 1.4×10⁻⁴ (2.6×10⁻³); 10 kpc 9.2×10⁻⁴ (7.3×10⁻³); 25 kpc 5×10⁻³ (0.034) | 6.9×10⁻³ | 0.23 |
| Milky Way M 1, omitted | 4.2×10⁻³ | 2.1 | 1.4 | 5 kpc 2.8×10⁻³ (0.01); 10 kpc 8.3×10⁻³ (0.024); 25 kpc 0.039 (0.086) | 8.2×10⁻³ | 0.23 |
| Milky Way M 1, on | 8.1×10⁻⁴ | 2 | 1.1 | 5 kpc 5.6×10⁻⁴ (2.9×10⁻³); 10 kpc 1.5×10⁻³ (7.4×10⁻³); 25 kpc 7.6×10⁻³ (0.031) | 2.7×10⁻³ | 0.49 |
| Milky Way M 10, omitted | 4.5×10⁻³ | 2 | 1.1 | 5 kpc 3.1×10⁻³ (0.01); 10 kpc 8.7×10⁻³ (0.024); 25 kpc 0.038 (0.082) | 2.4×10⁻³ | 0.25 |
| Milky Way M 10, on | 1.3×10⁻³ | 2.2 | 1.2 | 5 kpc 7.4×10⁻⁴ (4.3×10⁻³); 10 kpc 3.1×10⁻³ (0.011); 25 kpc 0.014 (0.042) | 1.4×10⁻³ | 0.28 |
| J1630 S 0.1, omitted | 0.011 | 2.2 | 1.4 | 3 kpc 1.3×10⁻³ (6.7×10⁻³); 10 kpc 0.01 (0.027); 21 kpc 0.028 (0.066) | 0.031 | 0.32 |
| J1630 S 1, omitted | 0.01 | 1.7 | 1.3 | 3 kpc 3.6×10⁻³ (7.9×10⁻³); 10 kpc 9.8×10⁻³ (0.026); 21 kpc 0.028 (0.065) | 9.6×10⁻³ | 0.23 |
| J1630 S 1, on | 1.8×10⁻³ | 2.3 | 1.2 | 3 kpc 2.9×10⁻⁴ (2.1×10⁻³); 10 kpc 1.6×10⁻³ (9.3×10⁻³); 21 kpc 6.9×10⁻³ (0.028) | 4.5×10⁻³ | 0.081 |
| J1630 S 10, omitted | 0.013 | 1.5 | 1.2 | 3 kpc 5.4×10⁻³ (0.01); 10 kpc 0.013 (0.029); 21 kpc 0.032 (0.068) | 2.9×10⁻³ | 0.18 |
| J1630 S 10, on | 6.8×10⁻³ | 2 | 1.3 | 3 kpc 1.4×10⁻³ (4.6×10⁻³); 10 kpc 6.5×10⁻³ (0.018); 21 kpc 0.017 (0.045) | 2.2×10⁻³ | 0.18 |
| J1630 M 0.1, omitted | 0.011 | 2.2 | 1.4 | 3 kpc 1.3×10⁻³ (7.3×10⁻³); 10 kpc 0.011 (0.029); 21 kpc 0.031 (0.072) | 0.038 | 0.24 |
| J1630 M 0.1, on | 2.1×10⁻³ | 2.4 | 1 | 3 kpc 3.6×10⁻⁴ (2.6×10⁻³); 10 kpc 1.7×10⁻³ (0.013); 21 kpc 8.3×10⁻³ (0.038) | 0.019 | 0.07 |
| J1630 M 1, omitted | 9.3×10⁻³ | 2.4 | 1.3 | 3 kpc 1.2×10⁻³ (6.4×10⁻³); 10 kpc 8.6×10⁻³ (0.028); 21 kpc 0.03 (0.072) | 0.012 | 0.23 |
| J1630 M 1, on | 3×10⁻³ | 2.6 | 1.1 | 3 kpc 2.9×10⁻⁴ (2.7×10⁻³); 10 kpc 2.7×10⁻³ (0.013); 21 kpc 9.4×10⁻³ (0.037) | 6.6×10⁻³ | 0.28 |
| J1630 M 10, omitted | 0.015 | 1.9 | 1.2 | 3 kpc 3.2×10⁻³ (10⁻²); 10 kpc 0.014 (0.033); 21 kpc 0.037 (0.077) | 3.6×10⁻³ | 0.22 |
| J1630 M 10, on | 8.6×10⁻³ | 2 | 1.2 | 3 kpc 1.6×10⁻³ (5.5×10⁻³); 10 kpc 7.7×10⁻³ (0.021); 21 kpc 0.019 (0.053) | 2.7×10⁻³ | 0.078 |
| Coma, low bracket S 0.1, omitted | 0.019 | 2.4 | 1.4 | 300 kpc 4.3×10⁻³ (0.016); 1000 kpc 0.011 (0.031); 2000 kpc 0.026 (0.057); 3000 kpc 0.04 (0.088) | 4.2 | 0.22 |
| Coma, low bracket S 1, omitted | 0.018 | 2.5 | 1.4 | 300 kpc 2.3×10⁻³ (0.015); 1000 kpc 0.011 (0.03); 2000 kpc 0.024 (0.055); 3000 kpc 0.037 (0.084) | 1.3 | 0.2 |
| Coma, low bracket S 10, omitted | 0.017 | 2.4 | 1.4 | 300 kpc 4.6×10⁻³ (0.015); 1000 kpc 0.01 (0.029); 2000 kpc 0.023 (0.054); 3000 kpc 0.037 (0.084) | 0.38 | 0.19 |
| Coma, low bracket M 0.1, omitted | 0.026 | 2.6 | 1.3 | 300 kpc 3.8×10⁻³ (0.02); 1000 kpc 0.015 (0.041); 2000 kpc 0.034 (0.074); 3000 kpc 0.054 (0.11) | 5.5 | 0.18 |
| Coma, low bracket M 1, omitted | 0.026 | 2.5 | 1.3 | 300 kpc 3.8×10⁻³ (0.02); 1000 kpc 0.016 (0.041); 2000 kpc 0.034 (0.072); 3000 kpc 0.051 (0.11) | 1.7 | 0.19 |
| Coma, low bracket M 10, omitted | 0.023 | 2.4 | 1.3 | 300 kpc 5.2×10⁻³ (0.019); 1000 kpc 0.015 (0.038); 2000 kpc 0.031 (0.068); 3000 kpc 0.049 (0.11) | 0.51 | 0.22 |
| Coma, high bracket S 0.1, omitted | 0.018 | 2.6 | 1.4 | 300 kpc 2.7×10⁻³ (0.014); 1000 kpc 0.011 (0.03); 2000 kpc 0.025 (0.055); 3000 kpc 0.038 (0.085) | 2.8 | 0.19 |
| Coma, high bracket S 1, omitted | 0.019 | 2.5 | 1.4 | 300 kpc 4.1×10⁻³ (0.015); 1000 kpc 0.011 (0.031); 2000 kpc 0.025 (0.057); 3000 kpc 0.039 (0.088) | 0.88 | 0.21 |
| Coma, high bracket S 10, omitted | 0.017 | 2.5 | 1.3 | 300 kpc 3.1×10⁻³ (0.015); 1000 kpc 9.9×10⁻³ (0.03); 2000 kpc 0.024 (0.055); 3000 kpc 0.038 (0.085) | 0.26 | 0.22 |
| Coma, high bracket M 0.1, omitted | 0.033 | 2.5 | 1.2 | 300 kpc 4.9×10⁻³ (0.025); 1000 kpc 0.021 (0.05); 2000 kpc 0.044 (0.088); 3000 kpc 0.064 (0.13) | 4.1 | 0.17 |
| Coma, high bracket M 1, omitted | 0.031 | 2.6 | 1.2 | 300 kpc 5.4×10⁻³ (0.023); 1000 kpc 0.02 (0.048); 2000 kpc 0.041 (0.084); 3000 kpc 0.063 (0.13) | 1.3 | 0.18 |
| Coma, high bracket M 10, omitted | 0.03 | 2.4 | 1.3 | 300 kpc 6.9×10⁻³ (0.023); 1000 kpc 0.019 (0.047); 2000 kpc 0.04 (0.083); 3000 kpc 0.061 (0.13) | 0.38 | 0.2 |

- **An outer envelope.** The companions form an envelope that is dilute where the baryons are. At B1 its mass is mostly beyond the tested radii, so it adds little to motions or lensing there.
- **The Milky Way's shape.** Across 5–25 kpc the Milky Way's companion mass grows as r^1.7–2.7. CR-2 found that the measured curve needs roughly r^1.3 there.
- **The inner benchmark.** Moving the interior to B2 takes 2.2–3.9 times the density. By then the total confined mass is tens to hundreds of times the baryons (table above).

## Energy exported

Within the regime, with the bath's gravity omitted: mass exported per unit of confined mass at 10 Gyr, by channel, with the speed range that carries most of it. The energy ledger closes to 3×10⁻¹⁴ relative in every run (V6).

| Combination | seedless partners | ejections | evaporation | orbit escapes |
|---|---|---|---|---|
| Milky Way S 0.1 | 1 (300–500 km/s) | 0.05 (200–300 km/s) | 0 | 0 |
| Milky Way S 1 | 1.2 (300–500 km/s) | 0.18 (200–300 km/s) | 3.8×10⁻⁵ (50–100 km/s) | 0 |
| Milky Way S 10 | 1.7 (300–500 km/s) | 0.75 (200–300 km/s) | 1.2×10⁻³ (100–200 km/s) | 0 |
| Milky Way M 0.1 | 1 (300–500 km/s) | 0.035 (100–200 km/s) | 3.2×10⁻⁵ (100–200 km/s) | 0 |
| Milky Way M 1 | 1.1 (300–500 km/s) | 0.11 (100–200 km/s) | 9.1×10⁻⁵ (100–200 km/s) | 0 |
| Milky Way M 10 | 1.4 (300–500 km/s) | 0.43 (100–200 km/s) | 1.1×10⁻³ (50–100 km/s) | 0 |
| J1630 S 0.1 | 1 (300–500 km/s) | 0.039 (200–300 km/s) | 5×10⁻⁵ (50–100 km/s) | 0 |
| J1630 S 1 | 1.1 (300–500 km/s) | 0.12 (200–300 km/s) | 7.7×10⁻⁴ (200–300 km/s) | 0 |
| J1630 S 10 | 1.3 (300–500 km/s) | 0.46 (200–300 km/s) | 6.3×10⁻³ (100–200 km/s) | 0 |
| J1630 M 0.1 | 1 (300–500 km/s) | 0.029 (100–200 km/s) | 2.4×10⁻⁵ (200–300 km/s) | 0 |
| J1630 M 1 | 1.1 (300–500 km/s) | 0.09 (200–300 km/s) | 5.9×10⁻⁴ (100–200 km/s) | 0 |
| J1630 M 10 | 1.2 (300–500 km/s) | 0.31 (100–200 km/s) | 7×10⁻³ (100–200 km/s) | 0 |
| Coma, low bracket S 0.1 | 1 (300–500 km/s) | 1.5×10⁻⁴ (200–300 km/s) | 0 | 0 |
| Coma, low bracket S 1 | 1 (300–500 km/s) | 7.6×10⁻⁴ (200–300 km/s) | 1.5×10⁻⁶ (500–1000 km/s) | 0 |
| Coma, low bracket S 10 | 1 (300–500 km/s) | 3.4×10⁻³ (200–300 km/s) | 1.1×10⁻⁵ (500–1000 km/s) | 0 |
| Coma, low bracket M 0.1 | 1 (300–500 km/s) | 2.3×10⁻⁴ (100–200 km/s) | 0 | 0 |
| Coma, low bracket M 1 | 1 (300–500 km/s) | 8.1×10⁻⁴ (300–500 km/s) | 0 | 0 |
| Coma, low bracket M 10 | 1 (300–500 km/s) | 3×10⁻³ (100–200 km/s) | 0 | 0 |
| Coma, high bracket S 0.1 | 1 (300–500 km/s) | 4.2×10⁻⁴ (300–500 km/s) | 0 | 0 |
| Coma, high bracket S 1 | 1 (300–500 km/s) | 9.9×10⁻⁴ (300–500 km/s) | 0 | 0 |
| Coma, high bracket S 10 | 1 (300–500 km/s) | 2.6×10⁻³ (200–300 km/s) | 0 | 0 |
| Coma, high bracket M 0.1 | 0.99 (300–500 km/s) | 5.5×10⁻⁴ (500–1000 km/s) | 0 | 0 |
| Coma, high bracket M 1 | 0.99 (300–500 km/s) | 1.3×10⁻³ (300–500 km/s) | 0 | 0 |
| Coma, high bracket M 10 | 0.99 (300–500 km/s) | 4×10⁻³ (300–500 km/s) | 4.6×10⁻⁵ (200–300 km/s) | 0 |

Every seedless birth sends its partner out, so the process drives a companion wind of about the confined mass at the incident speed. Ejections matter in galaxies and are negligible in Coma. Evaporation is negligible at these σ/m.

## The fixed-background control and the incident velocities

C1's B1 density at 1 cm²/g (M_sun/kpc³) for Maxwellian baths of increasing mean speed, and for S with a fraction of its mass moved to 3,000 km/s. C1 is a control with a prescribed shape, not a formation result. These are the frozen-coefficient velocity thresholds of the whole process that the protocol asked for.

| Incident distribution | Milky Way | J1630 | Coma, low bracket | Coma, high bracket |
|---|---|---|---|---|
| Maxwellian, mean speed 100 | 3.2×10⁵ | 2.2×10⁵ | 2,015 | 2,176 |
| Maxwellian, mean speed 200 | 6.8×10⁵ | 4.5×10⁵ | 3,181 | 2,905 |
| Maxwellian, mean speed 300 | 1.1×10⁶ | 7.1×10⁵ | 4,573 | 4,014 |
| Maxwellian, mean speed 500 | 4×10⁶ | 1.6×10⁶ | 7,355 | 6,429 |
| Maxwellian, mean speed 1000 | 2×10⁷ | 1.2×10⁷ | 1.4×10⁴ | 1.2×10⁴ |
| fraction at 3,000 km/s 0 | 8.1×10⁵ | 5.2×10⁵ | 3,494 | 2,969 |
| fraction at 3,000 km/s 0.0034 | 8.2×10⁵ | 5.3×10⁵ | 3,506 | 2,979 |
| fraction at 3,000 km/s 0.01 | 8.3×10⁵ | 5.3×10⁵ | 3,528 | 2,999 |
| fraction at 3,000 km/s 0.1 | 10⁶ | 6.2×10⁵ | 3,822 | 3,251 |
| fraction at 3,000 km/s 0.34 | 2.5×10⁶ | 1.2×10⁶ | 4,998 | 4,249 |
| fraction at 3,000 km/s 0.5 | 5.1×10⁶ | 2.3×10⁶ | 6,361 | 5,299 |

Slower incident distributions need less. From 100 to 1,000 km/s the Milky Way's requirement rises about 60-fold, and a third of the mass at 3,000 km/s triples it.

## Deviations from the protocol, and implementation choices

The protocol's method was written before the code. Where the implementation departs from it, the departure and its reason are listed here. Items 1–4 replace declared rules; the others fill in details the protocol left open.

1. **Tracer masses** (protocol: equal masses m_p). Tracers born in a production shell carry a mass proportional to the square root of that shell's initial production rate, within a factor of 10 of the heaviest, so that more tracers resolve the inner region. A shell whose first pool holds no confined event takes the mass of the nearest shell that has one. Seed tracers are sized within the same range.
2. **Bound–bound pairing** (protocol: radial-neighbour pairing with probability (σ/m)ρ_c·g·Δ). Kept, with two additions that unequal masses and sub-cycling require:
   - **Collisional zones.** Where a tracer expects at least one bound–bound collision within the run, its radial zone is collisional. There the tracers are resampled to the zone's mean mass before colliding (unbiased; the mass and energy changes are booked). Equal masses then collide, as the protocol's equal-mass tracers would, and nothing is split. Where collisions are rare, the heavier tracer of a colliding pair is split so that equal masses collide, which conserves energy and momentum exactly.
   - **Partners.** Within a sub-cycled step, partners are redrawn among the eight nearest radial neighbours in every pass.
3. **Collision step control** (protocol: the global step shrinks until every collision probability per step is below 0.2). The orbit step stays at 5 Myr (galaxies) or 50 Myr (Coma), and collisions are sub-cycled with the same 0.2 limit per sub-step, up to 2¹⁰ sub-steps per step: per tracer for incoming–bound, per radial zone for bound–bound. Dense seed cores at 1,000 cm²/g would have needed global steps of about a thousand years.
   - Positions are frozen during the sub-steps, so heat conduction is not resolved where the mean free path is short.
   - The mass fraction at the cap and the number of pair attempts with probability above one are reported. An added label, "collisions resolved", fails if at least 1% of the mass reached the cap or any pair attempt had a probability above one.
4. **Tracer-count control** (protocol: above 4×10⁴ tracers, half are removed at random and the rest double in mass). Radial neighbours are merged in pairs instead. A merged tracer carries both masses and takes the position and velocity of one of the two, chosen with probability proportional to its mass: the mass is exact, and energy is conserved in expectation and booked. With unequal masses, random removal could delete a large share of the total at once. Separately, if more than 1/20 of the tracer limit would be born in one step, the mass of future births is doubled.
5. **Production pools.** Seedless production is sampled from one pool of events per logarithmic shell (40 shells). A pool is regenerated every 10 steps, or whenever the potential at r_half changes by more than 1%. Each pool aims at 100 confined events, with 1,000–20,000 draws set from the previous pool's confined fraction, so that every pool is an unbiased estimate.
6. **The static-bath test.** The protocol said the bath's excess gravitates but did not say what to do if no static state exists. The engine continues from a weightless bath in 12 density steps. At each step:
   - the fixed point must converge (Anderson mixing);
   - the net enclosed mass must stay positive;
   - the largest positive eigenvalue of the focusing response must stay below one.

   Otherwise the combination has no stable static bath at that density.
7. **A bath-gravity-omitted ladder** runs for every combination, as a labeled sensitivity. The gravity-on ladder runs only below each combination's static limit, plus one point at 0.9 of it.
8. **Local density.** The radial neighbour window used for collision rates is at least 2% of the radius. Resampled clones and repeated pool births share a radius exactly, and they had produced spurious collision probabilities.
9. **Stops.** A run stops as "not completed: runaway" if one step's captures exceed twice the tracer count. The ladder is not pushed to higher density once the bath is opaque inside r_half (τ > 3).
10. **Validation details.**
    - V2's barrier-limited case uses the slowest component of case M, in the Milky Way at 164 kpc.
    - V7's pointwise comparison stops at 0.9 R_b, where Φ → 0 makes relative errors meaningless; the error relative to the central depth is also given.
    - V8 compares bath-gravity-omitted ladders, which exist at every density.
11. **Seed truncation.** CF-1's counted seed is cut at E < 0. Here the cut is E < −10⁻³|Φ(r_half)|, the integrator's energy tolerance. In a pre-run 10 Gyr test, one of 10,000 J1630 seed companions, sampled as merely confined, sat on an E ≥ 0 orbit trapped by the boundary barrier and left through R_b. The excluded seed mass is negligible.
12. **A static bath can be lost during a run.** A run with the bath's gravity on stops as "not completed: static bath lost" if the net enclosed mass turns non-positive anywhere. By then the bath's deficit, from the forming companions' depletion or from barrier shielding, outweighs everything inside it, and the declared static-bath model no longer holds.
13. **The boundary model's own limit.** The bath enters R_b at its fixed trial speed. Once the net mass inside R_b exceeds u²R_b/G, the circular speed there exceeds u, and a centrifugal barrier empties the interior. For Coma that allowance is 6.3×10¹⁴ M☉ at 300 km/s, 4.4 times its baryons, and at its static limit the bath's focused excess alone uses it up. This limits the boundary condition, not nature: in a real environment the enclosed mass would accelerate incoming companions before they reached R_b.
14. **First-order depletion and gravity.** In the opaque regime the depletion factor also removes scattered companions, which are still unbound and still have mass, from the gravitating excess. That underestimates the bath's mass and is one more reason runs with τ(r_half) > 0.3 lie outside the model. In the transparent regime the difference is of order τ.
15. **Fully depleted shells and kernel guards.** Where depletion empties a production shell, its sampling weight had underflowed and inverse-CDF sampling returned an infinite radius. Now:
    - such shells produce nothing;
    - sampled radii are clipped to their shell;
    - the compiled grid-position and confinement kernels check for infinite or NaN radii before any integer conversion, which would otherwise index outside the tables.

### How the reported run was reached

There were four canonical launches, and only the fourth is reported.
- **The first** stopped after three minutes on a NaN birth rate (item 15).
- **The second** ran to completion in 6,790 s, and every validation passed. Its counted-seed controls were nevertheless wrong: in the Milky Way and J1630 they formed between 0.003 and 6.9 times the zero-seed mass at the same density, where a seed of 1% of the baryons should change little. Two defects were behind this, and the second affected every run:
  - **Splitting and thinning.** Splitting at every collision of unequal masses multiplied tracers in the dense seed cores, with 2–9×10⁵ splits per run. The repeated random-removal thinning that followed moved up to 135 times the final mass, and each thinning doubled the birth masses, so births nearly stopped. Items 2 and 4 replace both.
  - **Switched-off shells.** A production shell whose first pool happened to hold no confined event was switched off for the whole run (item 1). The Milky Way's outermost shell carries a quarter of the production. Twelve runs at one density fell into clusters at 4.7, 6.6 and 8.1×10¹⁰ M☉, exactly according to which of the two outermost shells had been switched off.
- **After both fixes**, four zero-seed runs at that density (Milky Way, S, 1 cm²/g, bath gravity omitted) agree within ±1%, and four seed controls lie 3% above them. The suite's checks passed, and the regression anchor was regenerated.
- **The third** canonical run started from scratch and took 4,878 s on 16 workers. Its verification rounds 2–4 had drawn their seeds in the order runs finished, so they depended on timing.
- **The fourth** reran from scratch after that was fixed (a53f45c), in 8,462 s on 7 workers while other runs shared the machine. Its rounds 0 and 1 reproduced the third run's line for line. Rounds 2–4 changed with their seeds, and the ranges above moved slightly: for example the Milky Way's shape across 5–25 kpc went from r^1.9–2.7 to r^1.7–2.7, and the gravity-on requirement from 0.16–0.71 to 0.16–0.67 of the gravity-omitted density.

## What this does and does not show

**Shown.**
- **The supply falls.** Under CF-1's law, with every collision class and the mass, energy, orbits, depletion and gravity evolving together, the incident density that forms a baryon-equal population in 10 Gyr falls far below CF-1's seeded-only requirement at 0.1–10 cm²/g.
- **The seed and self-gravity.** The result does not depend on a seed. It does depend on the companions' own gravity.
- **Where the mass goes.** The population is an extended envelope, still growing, with most of its mass beyond the tested radii. An inner reservoir needs a total confined mass tens to hundreds of times the baryons.

**Not shown.**
- **A source.** That is stage 2B, after the source-side radiation check.
- **Survival beyond 10 Gyr.** The numerical checks verify the calculation, not the reservoir's durability.
- **Any fit to motions or lensing.**
- **Coma with a gravitating bath.** Under this boundary model no in-regime Coma run keeps a static bath with its gravity on, and omitting it is not a small approximation there.
- **A bath with no stable static state.** Its gravitational collapse is outside this model.

## Requirements handed on

To the source-side radiation budget and spectrum check, and then to stage 2B:
- **Density.** Within the regime (0.1–10 cm²/g), a source must maintain at each system's boundary, for 10 Gyr:
  - 7×10⁴–1.7×10⁶ M_sun/kpc³ around the galaxies (1,700–42,000 times the cosmic mean);
  - 780–8,800 M_sun/kpc³ around Coma (19–220 times).

  The incident distribution is S or M, at 300 km/s mean energy.
- **Throughput.** The mechanism sends out about as much companion mass as it confines, as seedless partners at 300–500 km/s, plus ejections in the galaxies.
- **The inner benchmark.** An interior dominated by companions (B2) needs 2.2–3.9 times those densities.
- **The envelope's shape.** The Milky Way's measured curve needs companion mass that grows more slowly with radius than this envelope. That is a requirement on CR-3 and on any source spectrum, not a property a source can assume.

## Reproduce

```sh
python research_work/results/companion-formation/formation.py     # about 80 minutes on 16 otherwise idle workers (CC2A_WORKERS)
python research_work/results/companion-formation/checks.py        # suite job: V2, V3, V5-V7, C0 and the regression anchor
```

The driver compares its output with the archived results and overwrites them only with `--canonical`. Wall-clock timings are excluded from the comparison. Everything else, including the suite job's short-run anchor, reproduces exactly on the same machine and thread settings. Since a53f45c the verification rounds take their seeds in the order the tasks were issued, so the output no longer depends on the number of workers or on the order in which runs finish.
