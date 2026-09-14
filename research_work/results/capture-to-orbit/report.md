# Capture to orbit: receiver-assisted threshold production (RB-1) is not promoted

13 September 2026. Protocol declared before execution: [protocol.md](protocol.md). Baseline `main` at 3884b4f. All data were previously exposed; no parameter was fitted or selected. This is a conditional feasibility test in a fixed potential, not a formation simulation or a stability calculation.

**Revised the same day after review** ([protocol-addendum.md](protocol-addendum.md), results in `revision-results.json`). The first evaluation mixed two incident spectra, assumed an isotropic drag coefficient and left the receivers' energy change out of the energy balance. All three are corrected below. The conclusion is narrowed to the tested prescription.

## Outcome in plain language

We specified one local interaction: an incoming companion strikes a star or gas atom and creates one new particle of fixed mass. Conservation laws and the standard near-threshold rule then fix every outcome. We followed each product on its actual orbit instead of placing it on a circular one.

The interaction builds a finite, positive, extended reservoir, with a convergent tail (density falling as r^-4) and increasingly radial orbits outward. At the required reference inventory it fails three requirements:

1. **Momentum.** The receivers pay for the motion of every particle created in their frame. Even if every captured quantum were retained, ordinary matter would give up a median **34 times its own angular momentum**. This fails in 148 of 149 galaxies.
2. **Supply.** Only about **one part in 10^10** of the captured energy stays bound. The median supply multiplier is 4.1×10^9 when the spectrum extends below threshold and 8.2×10^9 when it is cut at threshold. Nearly everything leaves as fast new particles.
3. **Shape.** The population is too concentrated toward the centre. With the same inventory it overpredicts rotation speeds at all radii and scores worse than the original model on every SPARC partition.

**The tested RB-1 prescription does not produce a viable reservoir at the required reference inventory.** Its interaction and illumination assumptions imply severe supply and receiver back-reaction problems, and its normalized orbital population worsens the main rotation benchmarks. This does not exclude other capture interactions.

## Declared interaction and its consequences

The reaction is c + R → R + X, with X of fixed mass m and threshold E'_th = m c^2 (1 + m/2M_R) in the receiver frame. Near threshold the s-wave term dominates (Wigner threshold law): the cross section scales as p_f/p_i and emission is isotropic in the centre-of-momentum frame. Receivers are the archived SPARC stellar disks (SBdisk × 0.5), bulges (Vbul, × 0.7) and exponential gas disks (catalog 1.33 M_HI, scale fitted to Vgas only), or the archived Milky Way baselines. The incident field is I(E, n) = S(E) A(r, n): S is a flat number spectrum, and A is the reference model's direction-dependent attenuation (paper equation 6).

The kinematics, the threshold law, radiation drag and orbit mechanics are established physics. Applying them to companions is our hypothesis. Two spectrum controls share the same density per unit energy and are always reported separately:

- **S1, threshold-cut:** E in [E_th, 2 E_th].
- **S2, extends below threshold:** E in [0.5 E_th, 2 E_th].

1. **Bound production.** For S2, bound products are uniform inside the galaxy-frame escape ball |v| < v_esc(r), whatever the receiver's velocity. For S1, the Doppler shift of a moving receiver removes half of the slow products, because only companions overtaking the receiver can reach threshold. The rate falls by 0.500 at 200/500 km/s and by 0.511 at 0.01c/0.03c (exact-kinematics Monte Carlo gives 0.515). Deterministic quadrature over the actual field matches exact-kinematics Monte Carlo within 1.1% for both spectra.
2. **Efficiency.** Bound rest energy per absorbed energy is about (v_esc/c)^3/(3 I_E) for S2 and half that for S1, with I_E = 1.0736. The cubic scaling holds for any smooth spectrum.
3. **Momentum drawn from receivers.** A receiver moving at β through the field loses momentum κ (E/c^2) v per absorbed energy E. For a threshold cross section and an anisotropic field with pressure tensor P, κ = 1 + (1+χ) p⊥ for tangential motion, where p⊥ = P⊥/u and χ = arccosh(2)/I_E = 1.227 is the power-weighted d ln σ/d ln E. The isotropic, gray limit is the familiar 4/3; the isotropic threshold value is 1.742. The attenuated field has p⊥ = 0.26–0.36 across the receiver sites, so κ ranges over about 1.58–1.80.
4. **Line spectrum.** A bath-frame line with threshold inside the receivers' Doppler band gives a bound fraction v_esc^3/[2(δ+β)]^(3/2), confirmed within 0.6%. Tuning the threshold to the band edge makes products co-move with their receiver, but only receivers within (v_esc − v)^2/2c of one speed can react, about 0.15 km/s. The result is a fine-tuned ring, not a reservoir.

## Population, orbits and density

Populations are steady, phase-mixed orbit averages in a spherical monopole potential.

- **Primary potential:** receiver baryons plus the reference exact-third reservoir, held fixed.
- **Labeled sensitivities:** a baryons-only potential, and a stationary self-consistent fixed point in which the population's own gravity replaces the reference reservoir. The fixed point converged in all 167 systems.
- **Normalization:** retained inventory fixed at the reference total.

The revision changes the per-site bound-production factor uniformly: 0.500 for S1 (maximum deviation 0.1%) and 1.000 for S2 (0.005%). Population shapes and rotation predictions are therefore unchanged.

SPARC medians:

| Quantity | Median (10th to 90th percentile) |
|---|---|
| Mean orbital circularity L/L_c(E) | 0.41 (0.36 to 0.45) |
| Anisotropy β at 0.35, 1.4, 2.8, 11 and 45 R_d | 0.18, 0.34, 0.47, 0.77, 0.92 |
| Bound mass beyond 10 R_d / beyond 100 R_d | 47% / 6.9% |
| Outer density slope d ln ρ/d ln r | −4.02 (−4.18 to −3.68) |
| Half-mass radius / receiver half-mass radius | 4.0 (reference 3.7) |
| Plunging, pericenter < 0.05 R_d | 0.19% (max 14.7%) |
| Bound mass on orbits with radial period > 10 Gyr | 9.3% (3% to 19%) |

Production weighted by v_esc^3 favours the deep centre, so the enclosed mass inside the disk is too large. Predictions are sampling-invariant: requesting different radii changes enclosed fractions by at most 2.2×10^-15.

## Frozen transfer results

SPARC: 149 galaxies, 3,150 radii, original 89/29/31 split.

| Model | RMSE train / validation / test (km/s) | log RMS train / validation / test |
|---|---:|---:|
| Ordinary matter | 52.565 / 58.219 / 47.772 | 0.2904 / 0.2704 / 0.2561 |
| Simple MOND, archived fitted a0 | 19.890 / 26.876 / 16.398 | 0.1092 / 0.0954 / 0.0783 |
| Original exact-third reference | 29.025 / 32.495 / 23.591 | 0.1385 / 0.1154 / 0.0911 |
| MOND-guided mixture (f frozen) | 20.468 / 27.136 / 16.241 | 0.1192 / 0.1042 / 0.0832 |
| **RB-1, primary potential** | **49.876 / 55.556 / 45.397** | **0.1524 / 0.1414 / 0.1281** |
| RB-1, baryons-only potential | 209.750 / 133.465 / 247.631 | 0.3063 / 0.2743 / 0.3541 |
| RB-1, stationary self-consistent | 130.577 / 83.255 / 104.475 | 0.2128 / 0.1837 / 0.2261 |

Mean equal-galaxy errors by radial bin (r/R_d < 1, 1 to 3, ≥ 3):
- original reference: −6.73, −7.49, +10.41
- MOND-guided mixture: +1.91, −4.99, −4.83
- RB-1 primary: +23.94, +21.85, +15.75

RB-1 beats the reference in 61 of 149 galaxies and in 6 of the 21 inventory-capped galaxies.

Milky Way, fiducial R_d = 2.6 kpc and luminosity factor 1, RMSE in km/s:
- baseline I: reference 6.761, MOND-guided 9.259, RB-1 24.981
- baseline II: reference 10.410, MOND-guided 11.870, RB-1 7.760

RB-1 improves on the reference in 7 of the 18 archived scenarios.

## Closed energy and momentum ledger

All values are per unit retained rest energy, M_ret c^2: SPARC medians, with 10th to 90th percentiles in brackets.

| Entry | S1, threshold-cut | S2, extends below |
|---|---:|---:|
| Absorbed incident energy (supply multiplier) | 8.21×10^9 [7.8×10^8, 3.4×10^11] | 4.10×10^9 [3.9×10^8, 1.7×10^11] |
| Bound rest plus kinetic energy | 1 + 2.8×10^-7 | 1 + 2.8×10^-7 |
| Receiver energy change | −1984 [−7461, −790] | −992 [−3731, −395] |
| Escaping, rest plus kinetic | absorbed − bound − receiver change | same |
| Receiver angular momentum lost / held, at this supply | 5.1×10^11 | 2.5×10^11 |
| Push history for net radial push < 10% of gravity | ≥ 4.2×10^13 Gyr | ≥ 2.1×10^13 Gyr |

Receivers lose kinetic energy through the drag, at κ β^2 per absorbed energy. That is negligible against the throughput but about 10^3 times the retained rest energy, and it now enters the escaping energy.

The ideal case assumes every quantum is retained. With the field's own anisotropy and the threshold cross section, receivers would lose a median **33.7** times their angular momentum (6.6 to 129; above 1 in 148 of 149). With gray anisotropic drag the median is 25.6, against the archived isotropic 25.3, so anisotropy moves this ratio by only about 1%. The Milky Way fiducials give 5.1 and 5.7, or 3.9 and 4.3 gray. Receiver recoil is at most q/2 of the absorbed energy, with q = E'/(M_R c^2).

An exact-kinematics aggregate check closes incident = bound + escaping + receiver to 5×10^-11 for both spectra. In that control the receivers' energy loss is 40 times the bound energy.

## Promotion assessment

| Requirement | Result |
|---|---|
| Numerical verification | All 11 unit checks pass. They cover the kinematics, spectrum controls against Monte Carlo, drag against the pressure-tensor formula on four control fields, aggregate ledger closure, orbits, sampling and the Plummer relation. The 0.5 km/s quadrature-doubling gate still fails in 18 of 167 systems (max 2.38 km/s, UGC03546). Passing the verification jobs does not mean the full numerical model has converged in those systems |
| Fitted agreement | None. No parameter was fitted or selected |
| Transfer | Fails. Worse than the original reference on all SPARC splits. Milky Way baseline I worsens and II improves; 7 of 18 scenarios improve |
| Physical support | Collisionless orbits support the population in the fixed potential by construction. The stationary self-consistent version is more concentrated and worse. No collective-stability test was run |
| Demonstrated source supply | Fails for the tested prescription. The supply multiplier is 4–8×10^9, and receivers' momentum debt exceeds 1 in 148 of 149 galaxies even at ideal efficiency |

**Verdict: not promoted.** The original reference and MOND-guided branches are unchanged. RB-1 is archived as a failed candidate.

## What this establishes, and what it does not

For this interaction, the momentum debt is fixed by conservation. A receiver that stops a companion and emits the product isotropically in its own frame must supply the product's share of its own motion. That debt scales with the absorbed energy, not the retained energy. For a threshold reaction, a smooth spectrum makes bound efficiency scale as (v_esc/c)^3, and a sharp line needs fine tuning.

These statements apply to receiver-mediated two-body threshold production with heavy receivers and isotropic emission in the centre-of-momentum frame. They do not exclude collective receivers, multi-step absorb-then-emit processes, non-threshold interactions, or fields whose rest frame co-moves with the receivers. A uniformly boosted bath that co-moves with the receiver produces no drag; this is verified exactly.

## Next distinguishing calculations

1. **Self-illumination pilot.** Does a radiation field derived from a declared rotating source geometry reduce the total receiver torque enough to build the required reservoir? The account must include the energy and angular momentum lost by the emitting stars and carried away by escaping particles. Build the field from the emitters, transport it to the receivers and apply the same reaction. Assign no co-rotation factor and no new inventory normalization. Reduced drag would not solve the threshold-efficiency problem, which has to be evaluated alongside.
2. **Receiver back-reaction.** Evolve disk angular momentum under the drag, to bound the reservoir mass any receiver-mediated channel can create before measurable contraction.
3. **Two-step absorb-and-emit** with a fixed internal Q-value, carrying the same closed ledger.

Net rotation of the reservoir is not a success criterion; what matters is a supported orbital distribution.

## Reproduction

```sh
python research_work/results/capture-to-orbit/checks.py      # 11 unit checks, about 20 s
python research_work/results/capture-to-orbit/runner.py      # populations and predictions, about 7 min on 8 workers
python research_work/results/capture-to-orbit/revision.py    # consistency revision of the ledger, about 20 s
```

Scripts write to fresh directories under `research_work/generated/`; `--canonical` rewrites the archived `results.json`, `capture-to-orbit.png` or `revision-results.json` here. Modules:
- `kinematics.py`: exact two-body kinematics and threshold laws
- `incident.py`: unified incident field, moments, bound-production quadrature, exact reaction-weighted channels
- `supported_profile.py`: potentials, orbit averaging, population profiles, Plummer closure
- `ledger.py`: archived first ledger
- `inputs.py`: frozen receivers and comparison rows
- `runner.py`, `revision.py`: drivers

The v1.5 PDF predates this result.
