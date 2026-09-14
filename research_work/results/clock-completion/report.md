# CC-1 report: the co-scaling completion is internally consistent

**Branch label.** CC-1 is a separately labeled co-scaling branch, run as a consistency calculation at the project owner's direction.
- It does not meet the original nonexpanding premise: separations counted in atomic rulers grow.
- It does not establish a nonsingular past.
- It does not exclude other completions. The fixed-ruler, nonexpanding branch stays open.

**Result.** Matter and light are coupled universally to one metric, g_m = −c²dt² + n(t)²dx², with every mass, charge and dimensionless constant fixed in local units. All six declared checks pass (below). From the coupling alone, the measured redshift, clock ratios, pulse durations, cavity and ruler behavior, and the timing and amplitude of light and gravitational waves follow without double counting. The complete energy exchange closes.
- **Distances.** They depend on the history n(t), which comes from the field's own dynamics, not from the rulers.
- **The V = 0 dynamics.** With the potential set to zero, the history is coasting and turns around at a finite n.
- **What a hot past requires.** If the field's energy gravitates, a turnaround as hot as recombination needs a field energy density of at least 5.5% of 3H²c²/8πG. Nucleosynthesis temperatures need 1.5×10⁴ times that.

Protocol: [protocol.md](protocol.md), declared in 6ce4c7a before execution. Code: [atoms.py](atoms.py) (N1–N3), [ledger.py](ledger.py) (N5, O8) and [cc1.py](cc1.py) (driver, N4, N6, O1–O7). Results: [cc1-results.json](cc1-results.json).

## The coupling

- **Matter.** Each particle has the action −mc²∫dτ in g_m.
- **The Coulomb field.** Maxwell's equations are conformally invariant, so the field of a charge at rest is the flat Coulomb field in conformal time η = ∫dt/n; measured locally it is 1/(proper distance)².
- **An atom's Lagrangian, in each frame:**
  - field frame (t, comoving x): L = (m/2)n²|ẋ|² + k/(n|x|);
  - conformal frame (η, x): the metric is static and the mass is proportional to n, L̃ = (mn/2)|dx/dη|² + k/|x|.

  These are one action (L dt = L̃ dη).
- **Photons.** Maxwell's equations in g_m are PF-1's matched medium with ε = μ = n, so PF-1's wave law and its propagation results carry over unchanged.
- **Gravitational waves.** They are tensor perturbations of g_m with a constant Planck mass in local units: d/dt(n³ḣ) + nk²h = 0. PF-1's postulate P2, which gave them the photons' law, is compared.
- **The field.** L_n = (M/2)ṅ² − V(n), as in PF-1's postulate P3. The toy sets V = 0; no other potential is declared.

## Declared checks

| Check | Tolerance | Result |
|---|---|---|
| N1. Classical Coulomb atom, field frame, while n rises from 1 to 2 over 3,223 orbits | 10⁻⁶ | Orbital frequency in t constant to 4.8×10⁻⁸; proper size to 2.1×10⁻⁸. The coordinate size halves (0.5000000) |
| N2. The same atom in the conformal frame | 10⁻⁶ | ω_η/n̄ constant to 4.8×10⁻⁸; ω_η doubles. Its crossing times, converted to t, match N1's to 9×10⁻⁹ |
| N3. Quantum clock: 1-D soft-Coulomb atom, Schrödinger equation on the comoving grid | 10⁻⁶ | ω₁₀ and ω₂₀ constant in t to 6.8×10⁻⁹ and 2.0×10⁻⁸; their ratio 1.3125838 to 1.3×10⁻⁸ |
| N4. Measured redshift in each frame | 10⁻⁶ | 1+z = n_o/n_e: field frame 7.0×10⁻⁹ (classical) and 6×10⁻¹¹ (quantum); conformal frame 7.0×10⁻⁹ |
| N5. Energy ledger, homogeneous toy, V = 0 | 10⁻⁸ (energy, ledger); 10⁻⁶ (pair, turnaround) | Energy conserved to 1.6×10⁻¹⁴, ledger closed to 8.8×10⁻¹⁵, bound pair constant to 6.6×10⁻⁹, turnaround at the predicted n to 5.4×10⁻¹⁴ |
| N6. Cavity of fixed proper length while n doubles, ṅℓ/nc ≤ 10⁻⁷ | 10⁻⁶ | Round trip equals 2ℓ/c to 3.3×10⁻¹⁵ at all 21 epochs |

The residual drifts in N1–N3 are the tidal term (n̈/n)r of a changing scale factor, suppressed by the slow ramp. In the real universe it is of order (H/ω_atom)² ~ 10⁻⁶⁶.

## What an observer measures (O1–O8)

- **O1. Redshift.** 1+z_measured = (n_o/n_e)·ν_atom(t_o)/ν_atom(t_e). The coupling fixes the second factor: atoms tick at a constant rate in t (N1, N3). So the field frame puts the whole shift into propagation (a factor of 2.000000 across the ramp) and the conformal frame puts it into the clocks (1.99999999). Their products agree.
- **O2. Clock ratios.** Constant, because every dimensionless constant is fixed. The two quantum transitions keep their ratio to 1.3×10⁻⁸.
- **O3. Pulse durations.** Two emissions arrive stretched by n(t_o)/n(t_e) = 2.0000000 in the observer's atomic time: event stretch equals 1+z.
- **O4. Cavities and rulers.**
  - Light speed in local units is c at every epoch.
  - A bound cavity's round trip stays 2ℓ/c.
  - A laboratory whose ends are fixed in the field's coordinates would see its round trip double as n doubles. This is the sense in which ruler counts grow.
- **O5. Gravitational waves.**
  - Light and gravitational waves share the null cones of g_m under either law, so they arrive together.
  - Tensor-law amplitude falls as 1/n (0.5000000 across the ramp), so the siren distance equals D_L.
  - Under P2 the amplitude does not fall (1.0000000): the siren distance would be D_L/(1+z).
- **O6. Beam geometry.**
  - D_L = (1+z)D_M and D_A = D_M/(1+z), so D_L/((1+z)²D_A) = 1 exactly.
  - For a linear n (coasting), D_M = ln(1+z)/α, which is PF-1's luminosity distance.
  - A fixed-ruler reading, in which an object of fixed size at fixed distance D subtends size/D, would instead give D_L/D_A = 1+z.
- **O7. Energy exchange,** per comoving volume, all into the field:
  - **Radiation** loses (ṅ/n)E_rad exactly. The quantity nE_rad is invariant to 4×10⁻¹⁴ in N5, not just adiabatically.
  - **Free matter** loses (ṅ/n)(p²c²/n²)/E. This is its pressure work.
  - **Bound systems** exchange nothing on average; the instantaneous exchange ∝ (2K + V) averages to zero by the virial theorem.
  - **Rest mass** exchanges nothing.
  - **Gravitational waves** lose (ṅ/n)E_GW: adiabatically under the tensor law, exactly under P2.
  - **In the toy,** from n = 1 to 2 the field gains 2.35954, equal to the radiation's loss (1.00000), the relativistic particle's (1.35950), the slow particle's (3.7×10⁻⁵) and the pair's (3×10⁻¹³).
  - **Today,** the radiation-to-field rate is (ṅ/n)u_γ ≈ 9.5×10⁻³² W/m³, using the supernova value 70.48 km/s/Mpc for ṅ/n.
- **O8. History.** The declared V = 0 dynamics give M n̈ = E_rad/n + (matter terms) > 0, so n has a minimum.
  - The field's energy today equals the radiation energy it absorbed since then: z_max = E_field/E_rad and q₀ = −E_rad/(2E_field).
  - **Recombination.** A turnaround at least as hot as 3,000 K needs E_field ≥ 1,100 E_rad. Then q₀ ≥ −4.5×10⁻⁴, a coasting history, whose D_L PF-1 scored at Δχ² = +35 against flat FLRW. If the field's energy gravitates, that is at least 4.6×10⁻¹¹ J/m³, 5.5% of 3H²c²/8πG.
  - **Nucleosynthesis (8×10⁸ K).** It needs 2.9×10⁸ E_rad: 1.2×10⁻⁵ J/m³, 1.5×10⁴ times that unit.
  - **The turnaround belongs to the toy.** A V(n) falling faster than −W/n as n → 0 would remove it. The turnaround's temperature and density, in atomic units, are those of the corresponding epoch of a hot early universe; the conformal ("static") frame changes no dimensionless ratio.

## What CC-1 settles and what it does not

**Settled.** Within co-scaling there is an explicit coupling that produces the measured redshift, fixed clock ratios, local c, event stretch, simultaneous arrival of light and gravitational waves, and distance duality, with a closed energy ledger and no double counting.

**Not settled:**
- **The history.** Only the V = 0 toy was run. Any other history needs a declared V(n) or other dynamics, and only then a brightness comparison.
- **Gravity.** There is no complete gravitational action, including whether and how the field's energy gravitates.
- **Reservoir supply.** The field's energy is not shown to become companions or reservoirs.
- **The original premise.** Co-scaling does not meet the nonexpanding premise.

## Numerical notes

- **N1 and N2.** The ramp is n = 1.5 − 0.5 cos(πt/T), with T = 2×10⁴ orbital units for the classical atoms, padded with 20 orbits at n = 1 and at n = 2.
- **N3.**
  - **Setup.** T = 4×10⁴; 2,048 points over 160 Bohr radii; Strang steps of 0.05. The beat frequencies are fitted per window with the known beat ripples removed.
  - **Strang offset.** Strang splitting lowers ω₁₀ and ω₂₀ by 1.2×10⁻⁵ and 1.5×10⁻⁵ relative to the continuum eigenvalues. The offset is the same at every n, because the comoving step is a dilation of the proper-frame step.
- **N5.** DOP853 with rtol 10⁻¹².
- **N6.** The cavity times are integrals over each short interval.
- **The first execution** stopped at N6 on a root-finder tolerance below its allowed minimum. N1–N5 were identical in both executions.

## Reproduce

```sh
python research_work/results/clock-completion/cc1.py    # about 3 minutes; compares with cc1-results.json
```
