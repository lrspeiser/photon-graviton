# Capture to orbit: receiver-assisted threshold production (RB-1) is not promoted

13 September 2026. Protocol declared before execution: [protocol.md](protocol.md). Baseline `main` at 3884b4f. All data were previously exposed; no parameter was fitted or selected. This is a conditional feasibility test in a fixed potential, not a formation simulation or a stability calculation.

## Outcome in plain language

We specified one local interaction: an incoming companion strikes a star or gas atom and creates one new particle of fixed mass. Conservation laws and the standard near-threshold rule then fix every outcome. We followed each product on its actual orbit instead of placing it on a circular one.

The interaction does build a finite, positive, extended reservoir. Its outer tail is convergent (density falling as r^-4) and its orbits are increasingly radial outward. It fails three requirements:

1. **Momentum.** The receivers must supply the momentum of every particle created in their frame. Even if every captured quantum were retained, ordinary matter would have to give up a median **25 times its own angular momentum**. This fails in 147 of 149 galaxies.
2. **Supply.** Only about **one part in 10^10** of the captured energy stays bound (median supply multiplier 4.1×10^9). Nearly all of it leaves as fast new particles.
3. **Shape.** The population is too concentrated toward the centre. With the same inventory it overpredicts rotation speeds at all radii: +24, +22 and +16 km/s mean error in the inner, middle and outer bins. It scores worse than the original model on every SPARC partition. It does better in 61 of 149 individual galaxies and in 7 of the 18 Milky Way scenarios.

These failures follow from conservation and phase space, not from a tuned parameter. They point the next mechanism away from ordinary-matter receivers and threshold windows.

## Declared interaction and its consequences

The reaction is c + R → R + X, with X of fixed mass m and threshold E'_th = m c^2 (1 + m/2M_R) in the receiver frame. Near threshold the s-wave term dominates (the Wigner threshold law): the cross section scales as p_f/p_i and emission is isotropic in the centre-of-momentum frame. Receivers are the archived SPARC stellar disks (SBdisk × 0.5), bulges (Vbul, × 0.7) and exponential gas disks (catalog 1.33 M_HI, scale fitted to Vgas only), or the archived Milky Way baselines. The local reaction rate is proportional to ρ_b J, where J is the reference model's attenuated mean intensity (paper equation 6).

The kinematics, the threshold law, radiation drag and orbit mechanics are established physics. Applying them to companions is our hypothesis. The declared reference spectrum is a bookkeeping choice.

1. **Injection law (derived, then verified).** Take a spectrum smooth across threshold over widths larger than the receivers' Doppler widths. The s-wave rate is then constant per unit velocity-space volume near threshold. Bound products are therefore uniform inside the galaxy-frame escape ball |v| < v_esc(r), whatever the receiver's velocity. The local bound production rate is proportional to ρ_b J v_esc^3. Exact-kinematics Monte Carlo at v_esc = 0.03c, with the receiver moving at 0.01c, reproduces the bound fraction v_esc^3/(3N) within 0.4% (N = √3 − π/3). The speed distribution matches a uniform ball (KS 0.0035), and the mean product velocity is 0.006 v_esc against a receiver velocity of 0.33 v_esc.
2. **Efficiency.** Bound rest energy per absorbed energy is (v_esc/c)^3/(3 I_E), with I_E = 1.0736 for the declared spectrum (flat from threshold to twice threshold). At v_esc = 500 km/s this is 1.4×10^-9. The cubic scaling holds for any smooth spectrum.
3. **Momentum drawn from receivers.** An absorber moving at v through an isotropic bath loses momentum κ (E/c^2) v per absorbed energy E. This holds whether it keeps the energy or emits the product isotropically in its own frame. κ = 4/3 for gray absorption, the known Compton-drag result. When the cross section depends on energy, κ = 4/3 + χ/3, with χ the power-weighted d ln σ/d ln E. The declared spectrum gives χ = arccosh(2)/I_E = 1.227, so κ = 1.742. Exact-kinematics quadrature gives 1.33332 in the gray limit (M_R/m = 10^7) and 1.7401 at β = 3×10^-4, approaching 1.7422 as √β. The receivers therefore pay for the created mass's motion. A receiver population can fund at most about three-quarters of its own mass in co-moving product before losing all its momentum.
4. **Line spectrum.** Consider a bath-frame line whose threshold lies inside the receivers' Doppler band. The bound fraction is v_esc^3/[2(δ+β)]^(3/2), and exact kinematics agree within 0.6%. Tuning the threshold to the band edge makes products co-move with their receiver, but then only receivers within (v_esc − v)^2/2c of one speed can react. That is about 0.15 km/s for v_esc ≈ 500 km/s and v ≈ 200 km/s. The result is a fine-tuned ring, not a reservoir.

## Population, orbits and density

Populations are steady, phase-mixed orbit averages in a spherical monopole potential.

- **Primary potential:** receiver baryons plus the reference exact-third reservoir, held fixed.
- **Labeled sensitivities:** a baryons-only potential, and a stationary self-consistent fixed point in which the population's own gravity replaces the reference reservoir. The fixed point converged in all 167 systems, within 34 iterations.
- **Normalization:** retained inventory fixed at the reference total, as in the MOND-guided branch.

SPARC medians for the primary population:

| Quantity | Median (10th to 90th percentile) |
|---|---|
| Mean orbital circularity L/L_c(E) | 0.41 (0.36 to 0.45) |
| Anisotropy β at 0.35, 1.4, 2.8, 11 and 45 R_d | 0.18, 0.34, 0.47, 0.77, 0.92 |
| Bound mass beyond 10 R_d / beyond 100 R_d | 47% / 6.9% |
| Outer density slope d ln ρ/d ln r | −4.02 (−4.18 to −3.68) |
| Half-mass radius / receiver half-mass radius | 4.0 (reference 3.7) |
| Plunging, pericenter < 0.05 R_d | 0.19% (max 14.7%) |
| Bound mass on orbits with radial period > 10 Gyr | 9.3% (3% to 19%) |

The tail is convergent, as expected for marginally bound radial orbits in a Keplerian exterior. Overall the population is about as extended as the reference. However, production weighted by v_esc^3 favours the deep centre. The enclosed mass inside the disk is therefore larger, and the rotation curve rises too steeply. The last 9.3% of the mass is not phase-mixed over a 10 Gyr history. It is reported, not removed.

The profile lives on a model-defined grid from 10^-3 to 10^5 kpc. Requesting different radii changes predicted enclosed fractions by at most 2.2×10^-15, so the prediction is sampling-invariant.

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

Mean equal-galaxy errors by radial bin (r/R_d < 1, 1 to 3, ≥ 3). Negative means too slow.

| Model | Inner | Middle | Outer |
|---|---:|---:|---:|
| Original reference | −6.73 | −7.49 | +10.41 |
| MOND-guided mixture | +1.91 | −4.99 | −4.83 |
| RB-1 primary | +23.94 | +21.85 | +15.75 |
| RB-1 self-consistent | +57.13 | +59.01 | +35.42 |

Galaxy by galaxy, RB-1 beats the reference in 61 of 149 galaxies and the MOND-guided mixture in 38. Among the 21 inventory-capped galaxies it beats the reference in 6. Their RMSE is 31.44/27.35/14.91 km/s, against 29.03/24.45/13.99 for the reference and 25.42/19.93/14.09 for the MOND-guided mixture.

Milky Way, 38 Eilers bins, fiducial R_d = 2.6 kpc and luminosity factor 1. RMSE in km/s, with mean bias in parentheses.

| Model | Baseline I | Baseline II |
|---|---:|---:|
| Original reference | 6.761 (+1.01) | 10.410 (−7.18) |
| MOND-guided mixture | 9.259 (+2.33) | 11.870 (−7.01) |
| RB-1 primary | 24.981 (+21.64) | 7.760 (+2.85) |
| RB-1 self-consistent | 51.101 (+46.23) | 16.814 (+13.80) |

Baseline I overshoots by 30 km/s in the inner bins; the lighter baseline II benefits. Across all 18 archived scenarios, RB-1 beats the reference in 7 (primary) and 3 (each sensitivity), with primary RMSE ranging from 7.4 to 46.9 km/s.

## Energy and momentum ledger

All values are per unit retained rest energy, M_ret c^2, and are SPARC medians unless stated.

| Entry | Value |
|---|---|
| Absorbed energy, the supply multiplier | 4.10×10^9 (10th to 90th percentile 3.9×10^8 to 1.7×10^11; Milky Way fiducials 1.03×10^9 and 1.31×10^9) |
| Escaping as fast X particles, rest plus kinetic | essentially all of the absorbed energy |
| Bound kinetic energy at injection | 2.8×10^-7 |
| Receiver recoil | ≤ q/2 of absorbed energy, q = E'/(M_R c^2); negligible for eV-scale quanta on nucleons |
| Receivers' angular momentum lost / held, ideal (every quantum retained, κ = 4/3) | **25.3** (5.0 to 96.6); below 1 only for NGC3741 (0.39) and UGC05721 (0.86); Milky Way 3.9 and 4.3 |
| Same, declared spectrum and κ | 2.5×10^11 |
| Receivers' kinetic energy drained / held, declared | 4.6×10^11 |
| Illumination history keeping the bath's net radial push below 10% of gravity | ≥ 2.1×10^13 Gyr (minimum 9.1×10^9 Gyr); universe age left free |
| Reference inventory / receiver baryon mass | 12.5 (3.2 to 44) |

No dissipation is needed for the bound population, which is collisionless and born on its orbits. The captured energy leaves as a flux of fast massive particles. That is an invisible outgoing channel, not a heat signature, and recoil heating is negligible. The drag's physical consequence is that receivers lose orbital angular momentum and spiral inward. That consequence is not included in the fixed potential; it is the reason for the verdict.

## Promotion assessment

| Requirement | Result |
|---|---|
| Numerical verification | All 8 unit checks pass. The archived reference is reproduced to 0.027 km/s, predictions are sampling-invariant, and every fixed point converges. The 0.5 km/s quadrature-doubling gate fails in 18 of 167 systems, all bulge-heavy massive galaxies. The worst is UGC03546 at 2.38 km/s, 0.4% of a ~580 km/s prediction that already overshoots the observed 262 km/s peak. This is retained as a failure and changes no conclusion |
| Fitted agreement | None. No parameter was fitted or selected |
| Transfer | Fails. Worse than the original reference on all SPARC splits (RMSE and log). Milky Way baseline I worsens and II improves; 7 of 18 scenarios improve |
| Physical support | Collisionless orbits support the population in the fixed potential by construction. The stationary self-consistent version is more concentrated and worse. No collective-stability test was run |
| Demonstrated source supply | Fails. Supply multiplier about 4×10^9, and receivers' momentum debt exceeds 1 in 147 of 149 galaxies even at ideal efficiency |

**Verdict: not promoted.** The original reference and MOND-guided branches are unchanged.

## What this teaches

The calculation exposes a general trilemma for converting isotropic light-speed quanta into slow bound mass through a local two-body process:

- **Smooth spectra** make the bound efficiency scale as (v_esc/c)^3.
- **Sharp lines** let products co-move with their receiver, but Doppler detuning confines reactions to a receiver-speed window of order 0.1 km/s.
- **Any receiver that stops the quanta** must supply the new mass's momentum, κ (E/c^2) v with κ between 4/3 and 1.74. Receivers at rest instead give products no net angular momentum.

The median SPARC inventory is 12.5 times the baryons. An extended, MOND-like reservoir therefore cannot be assembled by ordinary-matter receivers from an external isotropic bath. Three escape routes remain open and would need their own calculations:

- a companion bath that co-moves with the receivers, for example companions generated from the galaxy's own rotating starlight;
- a collective receiver whose momentum budget is not that of the visible matter;
- a two-step absorb-then-emit process with an internal Q-value. This removes the (v/c)^3 penalty but not the drag.

## Next distinguishing calculations

1. **Momentum budget of self-illumination.** Compute the mean momentum per energy of companions produced from a rotating disk's own starlight at each receiver. A bath that co-rotates on average reduces the drag on disk receivers; an external bath cannot. This is cheap and reuses the disk models.
2. **Receiver back-reaction.** Evolve disk angular momentum under the drag. This bounds the reservoir mass any ordinary-matter receiver channel can create before the disk contracts measurably.
3. **Two-step absorb-and-emit with a fixed Q-value,** carrying the momentum accounting from item 1.

The universal-K Plummer check (a ∝ M^2 in the isolated limit) passes and remains the analytic anchor for the queued baryon-coupled Coma support calculation.

## Reproduction

```sh
python research_work/results/capture-to-orbit/checks.py      # unit checks, about 15 s
python research_work/results/capture-to-orbit/runner.py      # full comparison, about 7 min on 8 workers
```

The runner writes to a fresh directory under `research_work/generated/`; `--canonical` rewrites `results.json` and `capture-to-orbit.png` here. `results.json` records the protocol, source and input hashes, quadrature, every prediction, per-galaxy orbit diagnostics and ledgers, all 18 Milky Way scenarios and the numerical checks. Modules: `kinematics.py` (exact two-body kinematics, drag, injection laws), `supported_profile.py` (potentials, orbit averaging, population profiles, Plummer closure), `ledger.py` (energy and momentum bookkeeping), `inputs.py` (frozen receivers and comparison rows). The v1.5 PDF predates this result.
