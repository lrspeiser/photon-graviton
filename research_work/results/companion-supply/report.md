# CC-2 stage 1 report: the bath CF-1 needs would be the dominant, self-gravitating medium

**Result.** At the supply CF-1 requires, its incoming companions cannot be the transparent, single-speed, gravitationally negligible background CF-1 assumed. For every trial cross-section from 0.1 to 1,000 cm²/g:
- **Opaque.** The bath is optically thick across the declared boundary region of the Milky Way and Coma, with optical depth 4–8.
- **Heavy.** Its own mass inside that region exceeds the system's baryons by 10¹–10⁷.
- **Unstable.** In Coma, at every trial scattering strength, it is Jeans-unstable on scales smaller than the region.

It would be the dominant self-gravitating, self-interacting component, closer to a self-interacting dark-matter medium than to a trace bath.

**The seed stops mattering.** Under the same interaction law, collisions between two incoming companions leave one bound often enough to beat CF-1's seed-driven capture by factors of 10²–10⁷, and they are also depth-keyed.

**Depth selectivity survives realistic distributions.** Retention needs a Maxwellian bath dispersion below about 210 km/s in the Milky Way, 330 km/s in a massive elliptical, and 510–700 km/s in Coma.

This stage predicts no reservoir. It turns CF-1 into requirements for a source (stage 2) and for a coupled test (stage 3), and that test must include the bath's own gravity and collisions.

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

## A2. The bath would outweigh the baryons, and could collapse

| σ/m (cm²/g) | Milky Way: bath mass in R_b / baryons; Jeans length / R_b | J1630 | Coma low | Coma high |
|---|---|---|---|---|
| 0.1 | 8×10⁵; 0.03 | 1×10⁴; 0.14 | 1×10⁷; 0.002 | 3×10⁶; 0.003 |
| 1 | 8×10⁴; 0.10 | 1×10³; 0.43 | 1×10⁶; 0.008 | 3×10⁵; 0.011 |
| 10 | 8×10³; 0.32 | 112; 1.4 | 1×10⁵; 0.024 | 3×10⁴; 0.035 |
| 100 | 812; 1.0 | 11; 4.3 | 1×10⁴; 0.08 | 2.6×10³; 0.11 |
| 1,000 | 81; 3.2 | 1.1; 13.5 | 988; 0.24 | 258; 0.35 |

- **A far higher density.** The incident density equals A*/(σ/m). In units of the cosmic mean density, that runs from 1.8×10⁷ down to 1.8×10³ for the Milky Way, and from 3×10⁵ down to 15 for Coma.
- **CF-1's premise does not hold at this supply.** CF-1 put a fixed baryonic potential around a bath of negligible gravity. At the supply CF-1 itself requires, the bath dominates the mass and, in Coma, fragments on its own.

## A3. Travel and time

| | Distance travelled in 10 Gyr / R_b (at 300 / 1,000 / 3,000 km/s) | Crossing time of R_b at 300 km/s | Crossing time of the half-mass radius |
|---|---|---|---|
| Milky Way | 10 / 34 / 102 | 1.0 Gyr | 0.02 Gyr |
| J1630 | 10 / 34 / 102 | 1.0 Gyr | 0.03 Gyr |
| Coma | 0.10 / 0.34 / 1.0 | 98 Gyr | 5.0 Gyr |

- A 300 km/s companion travels only about 3 Mpc in 10 Gyr. Coma's declared 30 Mpc boundary therefore cannot hold a steady bath in that time, and the region that could supply Coma within 10 Gyr is about a tenth as large.
- A longer history is allowed in the fictional universe, but it would change the exposure and the source's own evolution. It has to be declared.

## A4. Particle against wave

For the particle treatment, the de Broglie wavelength must be at most a tenth of the half-mass radius. At 300 km/s that needs a companion mass of at least:

| Milky Way | J1630 | Coma |
|---|---|---|
| 6.3×10⁻²³ eV | 3.8×10⁻²³ eV | 2.6×10⁻²⁵ eV |

The archived 1.34×10⁻²⁴ eV constituent fails in galaxies but passes in Coma. The 10⁻²² eV constituent passes everywhere. For the archived masses, the wavelength is at most 0.6 of the mean free path (the light constituent in the Milky Way).

## A5. The velocity distribution decides the sign

These use CF-1's kernels, which are linear in the incident distribution.

| | Fast fraction at 3,000 km/s that reverses 300 km/s retention | K with a 1% fast tail (km/s) | K for a Maxwellian bath of 1-D dispersion 100 / 200 / 300 / 500 / 1,000 km/s |
|---|---|---|---|
| Milky Way | 0.34% | −20 | +552 / +47 / −235 / −654 / −1,531 |
| J1630 | 12.7% | +391 | +1,500 / +537 / +94 / −455 / −1,431 |
| Coma, low | 34% | +1,390 | +3,839 / +1,744 / +901 / +28 / −1,191 |
| Coma, high | 53% | +2,897 | +7,349 / +3,560 / +2,119 / +759 / −825 |
| Full well (inverse, labeled) | 95% | +19,480 | positive throughout |

- **Retention windows.** Retention reverses at a bath dispersion of about 210 km/s in the Milky Way, 330 km/s in the lens host, and 510 and 700 km/s for Coma's two ends. A Maxwellian's mean speed is 1.6 times its one-dimensional dispersion.
- **Depth selectivity survives.** A bath with a dispersion of 250–450 km/s erodes the Milky Way but fills Coma.
- **The source must deliver the whole distribution,** its fast tail included, not a mean speed. The table covers the kernel from 30 to 10,000 km/s. The probability outside that range is below 0.7% for every dispersion, and below 0.03% at 300 km/s and above.

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

- **Seedless capture dominates everywhere,** by 10²–10⁷, except in the lens host at the largest cross-section, where it is comparable. The seed is not needed.
  - The owner's kinematic example holds in general: two incoming companions can leave one bound while the other carries the energy away.
  - The channel is itself depth-keyed, through both the capture probability and the focused density.
- **Evaporation is negligible in Coma,** but significant for galaxies at σ/m ≳ 100 cm²/g.
- **Caveat.** These rates use the transparent bath of CF-1 and are computed at its initial distribution. A1 and A2 show that at this supply the bath is itself opaque and self-gravitating, so they are indicators, not predictions.

## Requirements for a stage-2 source

The source must deliver, with its energy taken from the clock field and the field's back-reaction included:
1. **A slow distribution,** with a one-dimensional dispersion below about 200 km/s for any galaxy retention and below 500–700 km/s for Coma. Its tail above 3,000 km/s must stay below about 0.3% for galaxies and below 30–50% for Coma.
2. **Exposure** A* of 10⁻⁴–10⁻² kpc⁻¹. At that supply the population is opaque on the boundary scale, outweighs the baryons and can collapse on its own, so stage 3 must evolve its gravity and collisions together with any retained mass.
3. **A declared time and region.** A 300 km/s population reaches Coma from only about 3 Mpc in 10 Gyr.
4. **A companion mass** of at least about 6×10⁻²³ eV for the particle treatment in galaxies, or a wave treatment otherwise.

## What this does and does not show

**Shown.**
- CF-1's bath is internally inconsistent at the supply CF-1 requires: it would be opaque, heavy and self-gravitating.
- The same interaction's seedless channel outruns the seeded one.
- Depth selectivity survives realistic velocity distributions, but the fast tail decides galaxies.

**Not shown.**
- A source (stage 2).
- Any reservoir or its survival (stage 3).
- That the medium this implies would produce the observed halos or cluster excess.

This audit also raises a question for the programme: the bath this mechanism needs behaves like a self-interacting dark-matter medium in its own right. Whether that is compatible with the project's premises is the owner's decision.

## Reproduce

```sh
python research_work/results/companion-supply/cc2a.py   # about 3.5 minutes; suite job
```
