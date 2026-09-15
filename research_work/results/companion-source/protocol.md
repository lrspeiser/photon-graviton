# CC-2 stage 2B-F1: a coherent field decaying into companion pairs, as a source and phase-space test

Declared before execution, 15 September 2026. Baseline: `main` at 6a3b14b.

**Why it is run.** It follows the owner's review of 68eb17c. Supply is no longer the principal obstacle; where the mass goes is.
- **Stage 2A's envelope.** Stage 2A's elastic mechanism builds an outer envelope, M(<r) ∝ r^1.7–2.7, while the Milky Way needs about r^1.3. A larger source with the same phase space would only build a heavier envelope.
- **What 2B must now prove.** A phase-space distribution as well as an energy supply.
- **The owner's first toy** is a coherent field decaying into companion pairs. Pair production also answers RC-1's momentum objection, that one photon cannot become one slow massive companion.

## The source (a declared toy)

- **The field.** A homogeneous coherent field χ, at rest in each system's frame (labeled; a bulk flow of the system through the field is a later variant). Its rest-energy density is ρ_χc² and its lifetime 1/Γ is much longer than the 10 Gyr span. It decays as χ → C + C.
  - Pairs are produced at a constant rate q = Γρ_χ (companion rest mass per unit volume and time).
  - ΓT ≪ 1 is therefore a labeled assumption, and the field's depletion is reported.
- **Kinematics.** A decay at rest gives two companions back to back, each with energy m_χc²/2.
  - The mass defect ε = m_χ/(2m_C) − 1 fixes their common speed, v_d = c√(1 − (1+ε)⁻²) ≈ c√(2ε).
  - The decay is isotropic in the field frame.
- **Per unit of field energy.**
  - Companion rest mass takes a fraction 1/(1+ε), kinetic energy ε/(1+ε).
  - The pair carries zero momentum.
  - Companion mass is produced at a rate Γ/c² per unit of field energy.
  - The field loses Γρ_χc² per unit volume and time, the equal-and-opposite energy.
- **Trial decay speeds.** v_d = 3, 10, 30, 100 and 300 km/s (ε from 5×10⁻¹¹ to 5×10⁻⁷). Each is universal across systems. The field density, through q, is the one normalization.
- **The companion interaction.** CF-1's elastic law, unchanged, at σ/m = 0 (collisionless), 0.1 and 1 cm²/g: stage 2A's regime plus the collisionless limit.

## How the population is built

Stage 2A's machinery, with one added source channel and a growing bath:
- **Production inside R_b.** Pairs are born uniformly in volume at rate q, with speed v_d in isotropic directions.
  - A born companion that is confined (E < 0, or behind the centrifugal barrier) becomes a tracer.
  - An unconfined one leaves with the far-field bath; its energy is booked as exported.
- **The incident bath.** Companions produced outside R_b arrive there with speed v_d, at a density that grows as ρ_∞(t) = qt. Inside R_b they follow stage 2A's transparent Liouville distribution with first-order depletion.
  - Unbound companions produced locally add a density of order q times the crossing time. Against qt that is negligible, so it is omitted (labeled).
- **Everything else as in stage 2A.**
  - All three collision classes, exact orbits and self-gravity.
  - The bath's focused excess gravitates, with an omitted-gravity sensitivity. The uniform part is assigned to the background (labeled).
  - Stage 2A's collisional zones, pools and numerical controls are unchanged.
- **The zone of influence.** R_b is 300 kpc for the Milky Way and J1630 and 30 Mpc for Coma, CF-1's boundaries. The Milky Way at 150 kpc is a sensitivity. Infall from beyond R_b is not modeled (labeled).

## Scoring

1. **The Milky Way's profile gate.**
   - **The target** is the extra mass the observed speeds require at the 38 Eilers radii, M_req(<r) = r(v_obs² − v_b²)/G, with the archive's baseline-I baryon speeds, as in CR-2 and its interior diagnostic.
   - **The model's speeds** are √(v_b² + G·M_extra(<r)/r). M_extra is all gravitating non-baryonic mass at T: the confined companions plus, where it gravitates, the bath's focused excess.
   - **The normalization.** For each decay speed and σ/m, a ladder in q brackets the value that minimizes the RMSE, and a verification run sits at the interpolated q.

   Three gates, all declared here:
   - **G1, the fit** (CR-2's declared rule, unchanged). RMSE of at most 20 km/s over the 38 bins, and the inner 20 bins no worse than baryons alone.
   - **G2, the shape.** The companions' enclosed-mass log slope over 8–20 kpc lies within 0.3 of the required 1.31, the statistic of CR-2's interior diagnostic.
     - G1 alone does not test shape. In a design check made before this protocol, any power law up to r^2.4 met 20 km/s with its best normalization, and even r^2 scored 15 km/s.
   - **G3, the cost.** At that normalization, all gravitating non-baryonic mass inside R_b is at most 10 times the baryons. That mass is the confined companions plus, where it gravitates, the bath's focused excess. The threshold is a labeled benchmark: the order of the Milky Way's total mass out to about 200 kpc that the motions of its halo stars and satellites indicate.
     - A lenient 20 times is also reported.
     - For comparison, stage 2A's inner benchmark needed 410–670.

   **The profile gate passes** when G1, G2 and G3 all pass.
2. **Energy supply.**
   - **The field energy the gate requires:** qTc² per unit volume of companion rest energy, plus the kinetic fraction, so ρ_χ ≥ qT. It is reported in units of the cosmic mean density (a comparison unit only) and of the microwave background's energy density.
   - **The field's energy-loss rate.**
   - **The Jeans diagnostic of the uniform companion medium** that density implies. A medium of density qT at speed v_d is Jeans-unstable on scales above λ_J; its growth time is reported against T. This is an instability diagnostic. The medium's collapse is outside the model, and belongs to RC-2 and the global background.
3. **Universality.** The field is homogeneous, so q is the same everywhere. With the Milky Way's best-scoring decay speed, σ/m and q, J1630 and Coma run unchanged. Their companion masses at their test radii are reported against their baryons, as predictions, not gates.

## Expectation, stated before running

- **Cold decays** (v_d ≪ v_esc). Companions born nearly at rest fall on almost radial orbits.
  - **Profile.** Their time-averaged density is ρ ∝ 1/(r²·v_esc(r)), so M(<r) ∝ r/v_esc(r), about r^1.4 in the Milky Way's baryonic potential. That would meet G2.
  - **Cost.** The confined mass is then set by the zone: roughly (π/3)(R_b/r)·M_req(<r)/M_b. That is about 16 baryon masses at R_b = 300 kpc and about 8 at 150 kpc, marginal against G3.
  - **The bath.** A cold incident bath's focused density has the same 1/(r²·w) form.
- **Faster decays.**
  - They put pericenters, r_p ≈ r₀·v_d/v, outside 5 kpc, which cores the inner profile.
  - Above the local escape speed nothing is born bound.
- **Why run it.** Self-gravity, the bath's gravity and collisions can move all of this.

## Validation (tolerances declared)

- **F1, decay kinematics.** Energy and momentum are conserved in every sampled decay to rounding. Directions are isotropic within three standard errors.
- **F2, born-bound mass.** With collisions off and the baryons-only potential frozen, the confined mass after T equals qT∫P_conf(r)·4πr²dr within three standard errors, at v_d = 100 and 300 km/s. The integral is a quadrature of the confinement probability at speed v_d over isotropic directions.
- **F3, the cold-birth density.** With collisions off, a frozen potential and v_d = 1 km/s, the density at T in each radial bin across 5–25 kpc matches the radial-orbit quadrature within three standard errors or 5%, whichever is larger. The quadrature counts every passage of each orbit born at rest from r₀ within the finite span T.
- **F4, the growing bath.** With a frozen potential, only incoming–incoming collisions, and σ/m small enough that the bath is transparent, seedless production after T equals ½(σ/m)q²(T³/3)∫W g_n dV from stage 2A's pools, within three standard errors.
- **Stage 2A's V1–V8 carry over.** V6's energy ledger includes the births from the field.

## Declared labels

- **"Profile gate passed"** (G1, G2 and G3), with G1, G2 and G3 also reported separately, together with the lenient cost.
- **"Collisions resolved,"** **"transparent"** and **"outside the model's regime,"** as in stage 2A.

## Not claimed

- **Any fit to lensing or clusters.**
- **A global energy or flux budget.** That is RC-2.
- **The collapse of the uniform companion medium, or infall from beyond R_b.**
- **More than one field toy.** A bulk flow is a variant for later.
- **A microscopic model of χ.** The small mass defect is stated, not explained.

## Pre-declared follow-ups if the gate fails

Each needs its own protocol; the report says which apply:
- a bulk flow of the systems through the field frame;
- decay stimulated by already occupied companion states (Bose enhancement);
- a field whose density is sourced by the baryons.

## Files

`field.py` (decay kinematics and the new source channel), `f1.py` (the driver), `checks.py` (the suite job), `f1-results.json` and `report.md`, with stage 2A's `mc.py` extended.
