# Driven collective reservoir (CR-1)

**Result.** In this first implementation, no declared channel supports a collective reservoir.

- **Supply.** In 10 Gyr the derived photon source adds 1.6×10³ Msun to a 10⁹ Msun seed of the light constituent (1.34×10⁻²⁴ eV), and 0.46 Msun to one of the heavy constituent (10⁻²² eV). That is 1.3×10⁸ and 4.5×10¹¹ times short of the 2.09×10¹¹ Msun reference inventory.
- **Stationary control.** Without a bound source the seed is stationary, with excess energy below 6×10⁻⁷.
- **Required-rate diagnostic.** This diagnostic assumes the supply problem is solved, which lets support be judged separately from supply.
  - **Growth into the occupied mode is adiabatic** (excess energy at most 1.5×10⁻⁴).
    - The light condensate grows 123-fold but reaches only 59% of the inventory, because it contracts and collects from a shrinking volume.
    - The heavy one contracts below the grid within 1.04 Gyr. Its analytic endpoint is a stable 1.6 pc soliton at 25% of the Kaup limit, whose rotation contribution the Milky Way excludes.
  - **Growth that follows the photons is worse.**
    - It heats the heavy condensate to 7.7% of the ground-state energy before it too contracts below the grid.
    - For the light condensate it runs away: the collected microwave-background energy grows with the condensate's own radius. The condensate reaches the edge of the box at 0.35 Gyr with positive total energy.

Supply, support and the collection region fail separately. None of this tests a condensate with self-interaction, a population of modes, or another history.

**Protocol and files.**
- Protocol: [protocol.md](protocol.md), declared at `main` 660a4f3 and committed as e2d712e before any seed or source run.
- Code:
  - [gpp.py](gpp.py), the spherical Gross–Pitaevskii–Poisson solver;
  - [cr1.py](cr1.py), the driver;
  - [checks.py](checks.py), the suite job.
- Results: [cr1-results.json](cr1-results.json). The first execution's record is kept in [first-run-results.json](first-run-results.json).

## Assumptions of this first implementation

The review asked for these to be labeled as assumptions, not universal restrictions on collective reservoirs.

- **One coherent mode.** A single spherical, ℓ = 0 wavefunction: no population of modes and no angular momentum.
- **No self-interaction (g = 0).** Self-gravity and the baryonic potential are included. A repulsive interaction, one proposed support mechanism, is a different model: CR-2 in the [model contract](../../../research_plan/model-contract.md).
- **A counted seed.** 10⁹ Msun, which is 1% of the Milky Way variant I baryons. How the seed formed is not addressed.
- **Prescribed illumination.** Present-day starlight, using the L_bol = L[3.6] upper bound, plus the microwave background. It is collected only inside r₉₉, the radius holding 99% of the current condensate.
- **A 10 Gyr reference duration.** This is a benchmark, not a limit on the universe's age.
- **A nonrelativistic solver.** Runs stop when the condensate becomes unresolved. That is a statement about the solver, not a simulation of collapse to a black hole.

## Execution history, with every defect kept visible

**The first execution** of the declared driver completed, but only after 9.8 hours.
- **Overflow.** The imaginary-time relaxation drove a factor exp(−V·dt/2κ) to overflow. The resulting NaNs filled the overheating diagnostic of the light constituent's photon-shaped required-rate run, and that run kept going long after its condensate had left the box.
- **Silent failure.** The driver recorded the failure quietly, and nothing surfaced until the end. The project owner rightly objected to waiting hours to discover a NaN.

The record stays in first-run-results.json. Its supply numbers, and every run apart from that diagnostic, agree with the final run.

**Corrections.** All are numerical; the declared physics is unchanged.
- **No overflow.** The potential is shifted by its minimum inside the imaginary-time exponent.
- **Second-order relaxation.** The potential is frozen within each imaginary-time step, which leaves a second-order rather than a first-order bias in the relaxed state.
- **Resolved-domain stop.** The protocol's clause "until the solution leaves the resolved domain" is implemented as a stop when r₉₉ reaches the absorbing layer.
- **Growth cap.** Mass growth is capped at 0.2% per step.
- **Fail-fast driver.** Any overflow or non-finite value raises at once. Progress is logged with flushing, and each run has a 45-minute wall-clock budget.
- **A time step calibrated on energy.** At the original step of 0.1 rad of potential phase, the fail-fast monitoring exposed an artifact in the light constituent's runs. Strang splitting conserves a modified energy that, for that extended seed, sits about 1.4% from the true one. It showed up as a 1.4% "excess energy" and a slow leak into the absorber. Each mass now uses the largest phase fraction whose source-free evolution conserves energy to 10⁻⁴ over 0.3 Gyr.
- **A ground-state bound.** The overheating diagnostic compares the state's energy with the ground-state energy at the same mass. That reference is now the lower of two variational upper bounds: the current density relaxed in imaginary time, and the self-gravitating soliton of that mass placed in the baryonic potential. Relaxing only the current density had failed for the light photon-shaped run: its 457 kpc condensate never relaxed to its 0.4 kpc ground state, which gave an impossible −0.995. The soliton bound gives +1.09.

## Validation (declared tolerances)

| Check | Result | Tolerance |
|---|---|---|
| V1: Schrödinger–Poisson ground state against the archived values | E 1.1×10⁻⁵, μ 1.9×10⁻⁵, r_half 1.5×10⁻⁵; virial 1.2×10⁻⁵ | 10⁻⁴ |
| V2: real-time stationarity over 100 dynamical times | energy 1.1×10⁻¹⁰; density 4.6×10⁻⁶ | 10⁻⁸; 10⁻³ |
| V3: absorbing boundary | 99.9995% absorbed; 3.1×10⁻⁶ reflected | ≥ 99.9%; < 10⁻³ |
| V4: source into the mode | mass 2.7×10⁻¹³; excess energy −1.4×10⁻⁵; ledger 8×10⁻¹² | 10⁻⁸; 10⁻³; 10⁻⁶ |

| Constituent | Seed stationary residual | Calibrated step (rad of potential phase) | Energy drift over 0.3 Gyr |
|---|---:|---:|---:|
| 1.34×10⁻²⁴ eV | 1.3×10⁻³ | 0.00625 | 2.4×10⁻⁶ |
| 10⁻²² eV | 6.5×10⁻⁸ | 0.1 | 2.0×10⁻⁹ |

The light seed's residual costs only about 10⁻⁶ in energy, since the error enters quadratically.

## Supply

| Constituent | Seed half-mass radius | Collection rate | Starlight share | Collected in 10 Gyr | Short of 2.09×10¹¹ Msun by |
|---|---:|---:|---:|---:|---:|
| 1.34×10⁻²⁴ eV | 11.4 kpc | 161 Msun/Gyr | 53% | 1.6×10³ Msun | 1.3×10⁸ |
| 10⁻²² eV | 0.13 kpc | 0.046 Msun/Gyr | 99.9% | 0.46 Msun | 4.5×10¹¹ |

The compact heavy seed collects only the intense bulge starlight. The extended light seed reaches into the halo, where the microwave background supplies about half.

## Runs

The Milky Way RMSE column adds the reservoir in quadrature to the archived baryon speeds at the 38 Eilers radii. It is reported, not scored; baryons alone give 52 km/s.

| Constituent | Channel | Outcome | Final mass (Msun) | Largest excess energy | Escaped fraction | Milky Way RMSE |
|---|---|---|---:|---:|---:|---:|
| 1.34×10⁻²⁴ eV | A, no bound source | Stationary for 10 Gyr | 1.0×10⁹ | 5.7×10⁻⁷ | 4×10⁻⁹ | 52.0 |
| 1.34×10⁻²⁴ eV | B1, into the mode | Growth 1.6×10⁻⁶ | 1.0×10⁹ | 5.7×10⁻⁷ | 4×10⁻⁹ | 52.0 |
| 1.34×10⁻²⁴ eV | B2, at the local phase | Growth 1.6×10⁻⁶ | 1.0×10⁹ | 5.7×10⁻⁷ | 4×10⁻⁹ | 52.0 |
| 1.34×10⁻²⁴ eV | C, into the mode | Grows 123-fold while contracting from 11.4 to 7.4 kpc | 1.23×10¹¹ (59%) | 1.5×10⁻⁴ | 1×10⁻⁶ | 28.9 |
| 1.34×10⁻²⁴ eV | C, following the photons | Runs away. r₉₉ reaches the absorbing layer at 0.35 Gyr with positive total energy | 4.5×10¹² | 1.09 | 4×10⁻⁶ | 46.3 |
| 10⁻²² eV | A, no bound source | Stationary for 10 Gyr | 1.0×10⁹ | 2.4×10⁻⁷ | 1×10⁻¹⁰ | 51.6 |
| 10⁻²² eV | B1, into the mode | Growth 5×10⁻¹⁰ | 1.0×10⁹ | 2.4×10⁻⁷ | 1×10⁻¹⁰ | 51.6 |
| 10⁻²² eV | B2, at the local phase | Growth 5×10⁻¹⁰ | 1.0×10⁹ | 2.4×10⁻⁷ | 1×10⁻¹⁰ | 51.6 |
| 10⁻²² eV | C, into the mode | Contracts below the grid (r_half < 0.049 kpc) at 1.04 Gyr | 6.7×10⁹ | 1.3×10⁻⁵ | 1×10⁻⁹ | 46.4 |
| 10⁻²² eV | C, following the photons | Contracts below the grid at 0.36 Gyr | 5.5×10⁹ | 0.077 | 2.5×10⁻⁵ | 47.4 |

Mass ledgers close to 4×10⁻¹⁰ or better in every run. The first execution did not stop at the box edge in the light photon-shaped run; it continued to 10 Gyr, by which time 81% of the injected mass had escaped.

## Analytic endpoints at the reference inventory (2.09×10¹¹ Msun)

| Constituent mass | Self-gravitating half-mass radius | Compactness GM/(r c²) | Fraction of the Kaup limit | Reservoir speed at the Eilers radii | Milky Way RMSE with baryons (baryons alone: 52 km/s) |
|---|---:|---:|---:|---:|---:|
| 1.34×10⁻²⁴ eV | 8.9 kpc | 10⁻⁷ | 0.3% | 170–236 km/s | 57 km/s |
| 10⁻²² eV | 1.6 pc | 0.006 | 25% | 190–413 km/s | 110 km/s |

Neither endpoint is a black hole: the heavy constituent's core is 80 Schwarzschild radii across. Both, however, are single solitons whose shape the rotation curve rules out.

## Assessment against the declared rule

The rule required all four of the following:
- growth to the required inventory within 10 Gyr;
- no overheating (excess energy below 0.1);
- no dispersal (escaped mass below 10%);
- staying in the nonrelativistic domain.

Channel C tests only dynamics and cannot support formation on its own.

- **Channel A** has no bound growth, by construction.
- **Channels B1 and B2** grow the seed by 1.6×10⁻⁶ (light) and 5×10⁻¹⁰ (heavy), so they fail on growth.
- **In channel C**, results differ by channel and constituent:
  - the light condensate's growth into the mode is well behaved but falls short, at 59%;
  - its photon-shaped growth overheats (1.09) and leaves the box;
  - the heavy condensate leaves the resolved domain in both forms, the photon-shaped one close to the overheating limit (0.077).

**Verdict:** no channel supports a collective reservoir in this implementation.

## What this does and does not rule out

- **Ruled out, for this implementation.** A single coherent, non-interacting spherical condensate fed by the conversion law cannot be the galaxy reservoir.
  - Supply fails by a factor of 10⁸ to 5×10¹¹.
  - Even with supply solved, a coherent mode contracts. The heavy constituent heads toward a parsec soliton. The light one grows along the ground-state sequence into a soliton more compact than the reference reservoir: 7.4 kpc at 59% of the inventory, and analytically 8.9 kpc at the full inventory, which over-rotates the Milky Way (RMSE 57 km/s).
- **Failures kept separate, as the review asked.**
  - **Supply:** the derived source.
  - **Support:** with no self-interaction there is only quantum pressure, and the endpoints are compact.
  - **Collection region:** the r₉₉ rule drives the light photon-shaped runaway through the uniform microwave background.
  - **Duration:** the 10 Gyr span is a reference. A longer span cannot close a 10⁸-fold gap without a stated history.
- **Not tested.**
  - a repulsive self-interaction, which gives Thomas–Fermi support (CR-2 in the model contract);
  - a population of modes, meaning a warm or incoherent halo with angular momentum;
  - other illumination histories;
  - any channel in which something other than the products absorbs the photon's momentum. Channel A records the kinematics for when the products do.
- **Unresolved is not collapsed.** A nonrelativistic solver stopping on an unresolved core is not a simulation of black-hole formation. The heavy endpoint is a stable soliton at 25% of the Kaup limit.

## Reproduce

```sh
python research_work/results/collective-reservoir/cr1.py      # about 50 minutes; fails fast and logs progress
python research_work/results/collective-reservoir/checks.py   # suite job: V1-V4 and both seeds against the archive, about 3 minutes
```
