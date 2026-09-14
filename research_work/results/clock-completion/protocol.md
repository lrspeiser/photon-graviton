# Clock completion CC-1 (the co-scaling branch): what an observer measures

Declared before execution, 14 September 2026. Baseline: `main` at 3e2cf0e.

**What this branch is.** It is separately labeled and opened at the project owner's direction after review. It investigates the co-scaling completion of PF-1 as a **consistency calculation**.

**What it does not claim:**
- that the original nonexpanding requirement is satisfied (counted in atomic rulers, separations grow in this branch);
- that a beginning is eliminated;
- that other completions are excluded.

The fixed-ruler branch stays open and conceptually separate.

## Model (declared)

- **Coordinates.** (t, x), called the field frame. The clock field is n(t) > 0, with L_n = (M/2)ṅ² − V(n), as in PF-1's postulate P3. V is not specified.
- **Matter and electromagnetism couple universally** to g_m = −c²dt² + n(t)²dx². This is the r = 2 member of the repository's [clock–ruler family](../clock-ruler-completion/report.md).
  - Masses, charges and every dimensionless constant are fixed in local (proper) units.
  - Maxwell's equations in g_m are identical to PF-1's matched medium (ε = μ = n), so PF-1's photon propagation is kept exactly.
- **The coupling sets the local reference; nothing assumes it.** In the formula
  1+z_measured = (n_o/n_e)·ν_atom(t_o)/ν_atom(t_e),
  the matter coupling determines the second factor.
- **Gravitational waves** are tensor perturbations of g_m with a standard kinetic term. PF-1's postulate P2, which made gravitational waves obey the photons' scalar law, is replaced here; its consequences are compared.
- **Background gravity is open.** The background dynamics of g_m's conformal factor come from L_n, not from Einstein's equations. A complete gravitational action is an open requirement and is stated as such.

## Outputs (derived, not fitted)

- **O1.** The received frequency against the same transition at the observer. It is computed in the field frame and in the conformal frame (static metric, particle masses proportional to n). The products of the propagation factor and the clock factor must agree, so no frequency change is counted twice.
- **O2.** Ratios between different atomic clocks.
- **O3.** Pulse durations, meaning event stretch.
- **O4.** Cavities and rulers: the local light speed in local units at emitter and observer, and ranging within bound systems.
- **O5.** Arrival times of electromagnetic and gravitational waves, and the siren distance against D_L, under the tensor law and under P2.
- **O6.** Beam geometry: D_A, D_L and the distance-duality ratio.
- **O7.** The energy ledger: field, radiation, free matter, bound systems and gravitational waves, with every exchange rate.
- **O8.** The history. Dimensionless ratios as n becomes small (temperature and number density in atomic units), whether the declared field dynamics turn around, and what V(n) would have to do. No history is fitted or imported, and no supernova brightness fit is made.

## Numerical checks (tolerances declared now)

Units are c = 1 and, for the atoms, ħ = m = k = 1. Where n rises from 1 to 2 it follows n(t) = 1.5 − 0.5 cos(πt/T), with T long compared with the internal periods (stated in the script).

- **N1. A classical Coulomb atom in the field frame.**
  - The Lagrangian follows from g_m and Maxwell's Coulomb field in g_m: L = (m/2)n²|ẋ|² + k/(n|x|).
  - Over whole orbits, the orbital frequency in t and the proper orbit size n|x| must stay constant within 10⁻⁶, while the coordinate size falls as 1/n.
- **N2. The same atom in the conformal frame.** The metric is static and the mass is proportional to n: L̃ = (mn/2)|dx/dη|² + k/|x|. The orbit-averaged orbital frequency in conformal time η, divided by n, must stay constant within 10⁻⁶.
- **N3. A quantum clock in the field frame.** A one-dimensional soft-Coulomb atom obeys iφ̇ = −φ''/(2n²) + U(nx)φ in comoving x.
  - It starts in a superposition of the three lowest states.
  - The beat frequencies ω₁₀ and ω₂₀ are measured in t from projections onto the instantaneous eigenstates.
  - Each must stay constant within 10⁻⁶ while n doubles, and so must their ratio, which is a clock-comparison ratio.
- **N4. The measured redshift.** The frequency of a free mode is computed in each frame's own time and combined with that frame's clock factor from N1–N3. Both frames must give 1+z = n_o/n_e within 10⁻⁶.
- **N5. An energy ledger for a homogeneous toy.**
  - The toy has the field with V = 0 (declared for the toy only), three radiation modes, two free particles (one relativistic, one not) and one bound pair on a circular orbit.
  - Total energy must be conserved within 10⁻⁸ of the total.
  - The ledger must close: the field's gain equals radiation's and free motion's loss, minus the bound pair's change, within 10⁻⁸ of the total.
  - The bound pair's orbit-averaged proper energy must stay within 10⁻⁶ of its binding energy.
  - Integrated backward, the dynamics must turn around at the minimum n predicted analytically from energy conservation, within 10⁻⁶.
- **N6. A cavity.** A light round trip between the two ends of a system of fixed proper length must equal 2ℓ/c within 10⁻⁶ at every sampled epoch while n doubles, with ṅℓ/(nc) ≤ 10⁻⁷.

## Not claimed

- That co-scaling satisfies the original nonexpanding premise.
- A nonsingular past or any particular history.
- Supernova distances, reservoir supply, galaxy gravity or lensing.
- A complete gravitational action.

## Assessment rule

CC-1 passes as an internally consistent completion if N1–N5 pass and O1–O8 are derived without contradiction. That validates the completion's internal structure only; it says nothing about its history or its reservoir supply.
