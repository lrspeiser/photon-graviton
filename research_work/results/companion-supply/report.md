# CC-2 stage 1 report: at CF-1's seeded-only normalization, the bath would be opaque and heavy

**Result.** At the supply CF-1's seeded-only capture requires, its incoming companions cannot be the transparent, single-speed, gravitationally negligible background CF-1 assumed. For every trial cross-section from 0.1 to 1,000 cm²/g:
- **Opaque.** The bath is optically thick across the declared boundary region of the Milky Way and Coma, with optical depth 4–8.
- **Heavy, as an environmental estimate.** A uniform bath at that density filling the declared boundary region would hold 10¹–10⁷ times the system's baryons. That is the boundary-volume inventory, not the companion mass at the radii where motions and lensing are measured.
- **Unstable by the Jeans criterion.** In Coma, at every trial scattering strength, the Jeans criterion permits gravitationally growing perturbations on scales smaller than the region. That is an instability diagnostic, not a simulated collapse.

At that normalization the bath's own gravity and collisions could not be neglected.

**The old normalization is not a requirement.** Under the same interaction law, collisions between two incoming companions can leave one bound. The channel is kinematically allowed and its probability is verified. At the old normalization its initial rate exceeds CF-1's seed-driven capture by 10²–10⁷. Because that rate grows as the square of the incident density, it can lower the required supply substantially, so CF-1's density is not a lower bound for the complete mechanism. Stage 2A recomputes it.

**Velocity thresholds of the seeded coefficient.** Averaged over a Maxwellian bath, CF-1's seeded net-retention coefficient changes sign at a dispersion of about 210 km/s in the Milky Way, 330 km/s in a massive elliptical and 510–700 km/s in Coma. These thresholds belong to that coefficient. The complete process, with seedless formation, is evaluated in stage 2A.

This stage predicts no reservoir. The owner's review of 299ba1f narrowed several statements; see the [review qualifications](#review-qualifications-added-after-the-owners-review-of-299ba1f).

Protocol: [protocol.md](protocol.md), declared in 5381c68 before execution. Code: [cc2a.py](cc2a.py). Results: [cc2a-results.json](cc2a-results.json).

Inputs are CF-1's potentials, seeds, kernels and required exposures A* = (σ/m)ρ_∞, for retained mass equal to the baryons in 10 Gyr at a 300 km/s bath. The 10 Gyr span is a labeled assumption; the cosmic mean density is a comparison unit only.

## Validation (declared tolerances)

| Check | Result | Tolerance |
|---|---|---|
| V1: seedless capture probability, analytic against Monte Carlo, at five well depths | 1.1×10⁻³ relative; no collision ever left both bound | 10⁻² |
| V2: a Maxwellian narrowed to one speed reproduces the single-speed kernel | 2.0×10⁻⁵ | 10⁻³ |
| V3: energy and momentum per sampled collision | 2.4×10⁻¹⁵, 4.0×10⁻¹⁶ | rounding |

For two incoming companions with speed w = √(u² + v_esc²), the probability that a collision leaves one bound is P(φ) = max(0, 1 − u²/(w² sin φ)), where φ is the angle between their paths. Averaged over collisions it rises steeply with depth: 0.09, 0.35, 0.72, 0.91 and 0.98 for v_esc/u = 0.5, 1, 2, 4 and 8.

## A1. The bath scatters on itself

| At the required supply | A* (kpc⁻¹) | Scattering length 1/A* | Optical depth across R_b | Across the half-mass radius | Bath–bath collision time |
|---|---|---|---|---|---|
| Milky Way (R_b 300 kpc) | 1.5×10⁻² | 66 kpc | 4.5 | 0.10 | 0.22 Gyr |
| J1630 lens host (300 kpc) | 8.4×10⁻⁴ | 1.2 Mpc | 0.25 | 0.009 | 3.9 Gyr |
| Coma, low end (30 Mpc) | 2.6×10⁻⁴ | 3.9 Mpc | 7.8 | 0.40 | 12.6 Gyr |
| Coma, high end (30 Mpc) | 1.3×10⁻⁴ | 7.9 Mpc | 3.8 | 0.20 | 25.6 Gyr |

- **Mean free paths.** For a single-speed bath they are 3/4 of the scattering length, and for a Maxwellian 1/√2 of it.
- **Trading cross-section for density changes nothing.** At fixed A*, raising σ/m while lowering ρ_∞ leaves every entry in this table unchanged.
- **Where the bath is still transparent.** Only across the half-mass radius, and across the boundary for the lens host.

## A2. A uniform bath at that density would outweigh the baryons, and the Jeans criterion permits it to fragment

| σ/m (cm²/g) | Milky Way: bath mass in R_b / baryons; Jeans length / R_b | J1630 | Coma low | Coma high |
|---|---|---|---|---|
| 0.1 | 8×10⁵; 0.03 | 1×10⁴; 0.14 | 1×10⁷; 0.002 | 3×10⁶; 0.003 |
| 1 | 8×10⁴; 0.10 | 1×10³; 0.43 | 1×10⁶; 0.008 | 3×10⁵; 0.011 |
| 10 | 8×10³; 0.32 | 112; 1.4 | 1×10⁵; 0.024 | 3×10⁴; 0.035 |
| 100 | 812; 1.0 | 11; 4.3 | 1×10⁴; 0.08 | 2.6×10³; 0.11 |
| 1,000 | 81; 3.2 | 1.1; 13.5 | 988; 0.24 | 258; 0.35 |

- **A far higher density.** The incident density equals A*/(σ/m). In units of the cosmic mean density, that runs from 1.8×10⁷ down to 1.8×10³ for the Milky Way, and from 3×10⁵ down to 15 for Coma.
- **CF-1's premise does not hold at this supply.** CF-1 put a fixed baryonic potential around a bath of negligible gravity. At the supply CF-1 itself requires, a uniform bath over the boundary region would dominate the mass inventory, and in Coma the Jeans criterion permits it to fragment.
- **What these two columns measure.** The mass ratio fills the whole declared boundary region (for Coma, 30 Mpc) with the uniform incident density: it is an environmental estimate, not the companion mass or acceleration at the tested radii. The Jeans length u√(π/(Gρ_∞)) is an instability diagnostic for that uniform background under a simple support prescription; it does not evolve collapse, fragmentation or heating.

## A3. Travel and time

| | Distance travelled in 10 Gyr / R_b (at 300 / 1,000 / 3,000 km/s) | Crossing time of R_b at 300 km/s | Crossing time of the half-mass radius |
|---|---|---|---|
| Milky Way | 10 / 34 / 102 | 1.0 Gyr | 0.02 Gyr |
| J1630 | 10 / 34 / 102 | 1.0 Gyr | 0.03 Gyr |
| Coma | 0.10 / 0.34 / 1.0 | 98 Gyr | 5.0 Gyr |

- **A transport benchmark, not a source horizon.** A companion moving at a constant 300 km/s covers about 3 Mpc in 10 Gyr (the distance uT). A newly supplied 300 km/s population therefore cannot be assumed to establish CF-1's steady 30 Mpc bath within the 10 Gyr benchmark.
- That does not exclude a population already present nearby, one produced locally, or one with a longer declared history. Gravity and scattering also change transport relative to the constant-speed estimate. The initial distribution and production history have to be specified.

## A4. Particle against wave

For the particle treatment, the de Broglie wavelength must be at most a tenth of the half-mass radius. At 300 km/s that needs a companion mass of at least:

| Milky Way | J1630 | Coma |
|---|---|---|
| 6.3×10⁻²³ eV | 3.8×10⁻²³ eV | 2.6×10⁻²⁵ eV |

The archived 1.34×10⁻²⁴ eV constituent fails in galaxies but passes in Coma. The 10⁻²² eV constituent passes everywhere. For the archived masses, the wavelength is at most 0.6 of the mean free path (the light constituent in the Milky Way).

## A5. The velocity distribution decides the sign of the seeded coefficient

These use CF-1's seeded kernels, which are linear in the incident distribution.

| | Fast fraction at 3,000 km/s that reverses 300 km/s retention | K with a 1% fast tail (km/s) | K for a Maxwellian bath of 1-D dispersion 100 / 200 / 300 / 500 / 1,000 km/s |
|---|---|---|---|
| Milky Way | 0.34% | −20 | +552 / +47 / −235 / −654 / −1,531 |
| J1630 | 12.7% | +391 | +1,500 / +537 / +94 / −455 / −1,431 |
| Coma, low | 34% | +1,390 | +3,839 / +1,744 / +901 / +28 / −1,191 |
| Coma, high | 53% | +2,897 | +7,349 / +3,560 / +2,119 / +759 / −825 |
| Full well (inverse, labeled) | 95% | +19,480 | positive throughout |

- **Seeded retention windows.** The seeded coefficient reverses at a bath dispersion of about 210 km/s in the Milky Way, 330 km/s in the lens host, and 510 and 700 km/s for Coma's two ends. A Maxwellian's mean speed is 1.6 times its one-dimensional dispersion.
- **Depth selectivity of the seeded coefficient survives.** A bath with a dispersion of 250–450 km/s gives the Milky Way a negative seeded coefficient and Coma a positive one.
- **Scope of these thresholds.** They belong to CF-1's incoming–bound coefficient with a seed-shaped population. They are not the thresholds of the whole formation process, which also has the seedless channel and bound–bound losses. The 0.34% fast fraction belongs to the specific 300 and 3,000 km/s mixture; it is not a general limit on a source's tail.
- The table covers the kernel from 30 to 10,000 km/s. The probability outside that range is below 0.7% for every dispersion, and below 0.03% at 300 km/s and above.

## A6. The channels CF-1 left out

At the required exposure, the table gives two rates relative to CF-1's seed-driven capture, (σ/m)ρ_∞ K M_seed:
- seedless formation, ∫½(σ/m)ρ_in²⟨g·P⟩dV;
- bound–bound evaporation, ∫½(σ/m)ρ_seed²⟨g·P_esc⟩dV.

| σ/m (cm²/g) | Milky Way: seedless / seed; evaporation / seed | J1630 | Coma low | Coma high |
|---|---|---|---|---|
| 0.1 | 3×10⁶; 7×10⁻⁴ | 1×10⁴; 3×10⁻⁴ | 2×10⁷; 1×10⁻⁸ | 8×10⁶; 3×10⁻⁸ |
| 1 | 3×10⁵; 7×10⁻³ | 1×10³; 3×10⁻³ | 2×10⁶; 1×10⁻⁷ | 8×10⁵; 3×10⁻⁷ |
| 10 | 3×10⁴; 0.07 | 123; 0.03 | 2×10⁵; 1×10⁻⁶ | 8×10⁴; 3×10⁻⁶ |
| 100 | 3×10³; 0.73 | 12; 0.29 | 2×10⁴; 1×10⁻⁵ | 8×10³; 3×10⁻⁵ |
| 1,000 | 316; 7.3 | 1.2; 2.9 | 2×10³; 1×10⁻⁴ | 771; 3×10⁻⁴ |

- **At the old normalization, seedless formation's initial rate exceeds seeded capture everywhere,** by 10²–10⁷, except in the lens host at the largest cross-section, where they are comparable.
  - The owner's kinematic example holds in general: two incoming companions can leave one bound while the other carries the energy away.
  - The channel is itself depth-keyed, through both the capture probability and the focused density.
  - "Works" here means the channel is allowed and its probability is verified. It does not yet mean a durable reservoir of the right size forms: that needs the new companions' orbits, residence times and later collisions (stage 2A).
- **Evaporation is negligible in Coma,** but significant for galaxies at σ/m ≳ 100 cm²/g.
- **Caveat.** These rates use the transparent bath of CF-1 and are computed at its initial distribution. A1 and A2 show that at this supply the bath is itself opaque and self-gravitating, so they are indicators, not predictions.

## Requirements at the old normalization (superseded as source requirements)

After the owner's review of 299ba1f these are no longer requirements for a source: they were derived at CF-1's seeded-only normalization. Stage 2A recomputes the supply with all channels together, and the requirements handed to a source come from it. They are kept as the diagnosis at the old normalization:
1. **A slow distribution,** with a one-dimensional dispersion below about 200 km/s for positive seeded retention in galaxies and below 500–700 km/s in Coma. For the tested two-speed mixture, the tail at 3,000 km/s had to stay below about 0.3% for galaxies and below 30–50% for Coma.
2. **Exposure** A* of 10⁻⁴–10⁻² kpc⁻¹. At that supply the population is opaque on the boundary scale, and a uniform bath over the boundary region outweighs the baryons.
3. **A declared time, region and initial distribution.** The transport benchmark uT is about 3 Mpc for a 300 km/s population in 10 Gyr.
4. **A companion mass** of at least about 6×10⁻²³ eV for the particle treatment in galaxies, or a wave treatment otherwise.

## What this does and does not show

**Shown.**
- CF-1's bath is internally inconsistent at the supply CF-1's seeded-only capture requires: it would be opaque, and as a uniform environment it would be heavy and Jeans-unstable.
- The same interaction's seedless channel is allowed, and at that normalization its initial rate outruns the seeded one.
- The seeded coefficient's depth selectivity survives Maxwellian baths, and its sign depends on the fast tail.

**Not shown.**
- The supply the complete mechanism requires (stage 2A).
- A source (stage 2B).
- Any reservoir or its survival.
- That the medium this implies would produce the observed halos or cluster excess.

## Review qualifications (added after the owner's review of 299ba1f)

- **No new premise decision is needed.** The first version of this report asked the owner whether a dominant, self-interacting companion medium is compatible with the premises. The [universe contract](../../../research_plan/universe-contract.md) already answers it:
  - an independently inserted halo or bath, or dark matter under another label, is prohibited;
  - a gravitating companion population derived from the permitted energy supply and interactions, with its abundance, energy, motion and gravity calculated, is a permitted hypothesis to test.

  A derived population would still be an additional gravitating component; its explained origin and linked predictions are what would distinguish it, not its name.
- **The old density is not a lower bound.** The stage takes CF-1's density as an input and finds that an omitted channel dominates there. It does not solve again for the density the complete mechanism needs. Seeded capture grows as ρ_∞M_seed and seedless formation as ρ_∞², so the old requirement cannot be carried forward as a minimum. The defensible statement is: *the seeded-only requirement lies outside CF-1's assumed regime; the supply the seedless-inclusive mechanism needs is undetermined* until stage 2A.
- **Labels.** Jeans results are instability diagnostics, not simulated collapse. Boundary-volume masses are environmental estimates for a uniform bath, not companion mass at the tested radii. uT is a transport benchmark, not a source horizon.
- **Thresholds keep their scope.** The Maxwellian sign changes and the 0.34% fast fraction belong to CF-1's seeded coefficient and the specific distributions that produced them. A population can gain bound companions through incoming–incoming collisions even where incoming–bound collisions have negative net retention.
- **CC-1 stays a separately labeled branch.** Its clock field belongs to the co-scaling branch, in which separations counted in atomic rulers grow. A source compatible with operational nonexpansion remains the active objective; a CC-1-based source can only be a labeled comparison.

## Reproduce

```sh
python research_work/results/companion-supply/cc2a.py   # about 3.5 minutes; suite job
```
