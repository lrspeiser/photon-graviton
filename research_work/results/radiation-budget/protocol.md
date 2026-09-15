# RC-1: the radiation budget and spectrum on the source side

Declared before execution, 15 September 2026. Baseline: `main` at 7620668 (CC-2 stage 2A).

**Why it is run.** It follows the owner's direction after 2c9a110 ([persistent-background plan](../../../research_plan/persistent-microwave-background.md), "Shared constraints", item 1). Before stage 2B selects an interaction, every energy transfer into companions must predict what it does to the photons: their energy, their number and their spectrum. The target is the FIRAS spectrum and its residuals. It is the queue item "a source-side radiation budget and spectrum check" in the [model contract](../../../research_plan/model-contract.md).

**What is already established, and not repeated here.**
- **Frequency shift at fixed photon number** ([thermal conversion](../thermal-conversion/derivation.md)). In a fixed volume, such a shift turns a Planck bath into q⁻³·B_ν(qT). The FIRAS diagnostic χ² is 45.02 at q = 1, 86.00 at 0.9999 and 4,339.94 at 0.999. Keeping a Planck shape needs photon removal, dN/dt = −3hN with dU/dt = −4hU, with the energy split 1:3.
- **What a source or absorber would have to do.** Background replenishment and the thermalizer establish this. Neither is a mechanism.
- **Removing photons everywhere** (the electromagnetic audit) keeps FIRAS, but dims every beam by the same survival factor.
- **Supply comparisons under other interpretations.** The optical-growth reading falls 6×10⁵ short. Starlight over 10 Gyr for 26 galaxies falls a median 5,428 short.

None of these compared the radiation with the supply the complete formation mechanism needs, which [stage 2A](../companion-formation/report.md) now gives, and none asked which transfer law and history could deliver it.

## Question

Suppose photons supply the companions, as the project's postulate says.
1. **Energy.** At stage 2A's required incident densities, how much energy do the companions carry per unit volume, compared with the radiation present today?
2. **Spectrum.** For each declared transfer law, what does taking that energy from the microwave background do to its spectrum, and does the FIRAS comparison allow it?
3. **History and cost.** For the law that keeps a Planck shape, what past temperature and duration would the transfer need, and what does the same law do to every other beam?

## Inputs

- **Stage 2A's requirement.** The densities ρ_req for B1 and B2 in every combination labeled within its regime (τ_in(r_half) ≤ 0.3), with and without bath gravity, read from `companion-formation/formation-results.json`. Also the incident distributions: S, and M with the same mean energy, so ⟨v²⟩ = (300 km/s)² in both.
- **The companions' energy.** Their gravitating mass is what the benchmarks count. The energy to be supplied per unit volume is therefore ρc²(1 + ⟨v²⟩/2c²), rest energy included: the universe contract requires the reservoir to be derived from the permitted energy supply. The conversion uses CF-1's constants (M_sun = 1.98847×10³⁰ kg, kpc = 3.0856775814913673×10¹⁹ m).
- **Radiation today.**
  - The microwave background at 2.72548 K (Fixsen 2009), u = aT⁴.
  - The UV-to-mm background at 45, 100 and 170 nW m⁻² sr⁻¹. This is Hauser and Dwek's historical range, as in the optical-growth test; u = 4πI/c.
- **FIRAS.** `time_revision/data/firas.txt`, the LAMBDA monopole product: 43 channels of residuals, diagonal uncertainties and the Galaxy template.
  - Each fit profiles the colour temperature and the Galaxy coefficient, as in the thermal-conversion pass.
  - There is no covariance or calibration model, so scores are diagnostic comparisons, not official limits.
- **The archived loss coefficient** α₀ = 2.488993286382367×10⁻⁴ Mpc⁻¹. It sets the rate h = α₀c where a rate is needed.

## Transfer laws

Each law is declared with its full energy and photon ledger:
- **(a) Frequency shift at fixed photon number** (the archived loss law). Each photon keeps a fraction q of its energy.
- **(b) Removal of whole photons at fixed shape.** A fraction 1 − s of photons is removed at every frequency.
- **(c) The Planck-preserving combination.** dN/dt = −3hN and dU/dt = −4hU; the spectrum stays Planck at a falling temperature.

In every law, the energy the photons lose goes into the companion sector of the plan's bookkeeping and nowhere else. None of the three is a mechanism. Each is a declared limit on what a mechanism could do.

## Method

1. **Energy.** Compute E_req = ρ_req c²(1 + ⟨v²⟩/2c²) for B1 and B2 in each in-regime combination. Report it as a ratio to the microwave background alone, and to the microwave plus UV-to-mm background at each of the three levels.
2. **Laws (a) and (b) against FIRAS.**
   - For each law, find the largest fraction of the microwave background's energy it can transfer before the diagnostic χ² rises by 4 above no transfer, refitting the colour temperature and Galaxy coefficient.
   - Report the companion density that fraction could supply, as a ratio to ρ_req.
3. **Law (c).**
   - **With the microwave background as the donor**, report three quantities:
     - the past temperature T_i = T₀(1 + E_req/u₀)^(1/4);
     - the fraction of photons removed, 1 − (T₀/T_i)³;
     - the duration ln(T_i/T₀)/h at h = α₀c.

     No age is imposed; the duration is reported, not judged.
   - **Applied at the same rate to every beam**, the law removes photons from starlight too. Report the extra dimming, a photon-number factor e^(−3α₀D), in magnitudes at 100 Mpc, 1 Gpc and 3 Gpc.
   - **The microwave-only alternative.** Report the ratio of optical to microwave removal rates that keeps that dimming below 1% at 1 Gpc. That is the frequency selectivity a mechanism would need.
4. **The plan's two first tests.**
   - **Preservation.** A thermal input is a preservation control; only law (c) is expected to keep a Planck spectrum.
   - **Thermalization.** A nonthermal input is a thermalization test: a diluted 5,800 K blackbody as a starlight proxy, under each law. Each law maps a spectrum's shape onto a shifted or scaled copy, so none is expected to make it thermal. The check computes this rather than assuming it.

## Validation (tolerances declared)

- **W1.** The numerically integrated number and energy densities of a Planck spectrum equal (2ζ(3)/π²)(kT/ħc)³ and aT⁴ to 10⁻¹⁰.
- **W2.** Law (a) conserves photon number and scales energy by q to 10⁻¹⁰.
- **W3.** The FIRAS diagnostic χ² at q = 1, 0.9999 and 0.999 reproduces the thermal-conversion pass's 45.02, 86.00 and 4,339.94 to 0.01.
- **W4.** Law (c) keeps a Planck spectrum at T·e^(−ht) to 10⁻⁹ over ten e-folds of temperature, transporting the spectrum along its characteristics. Its time-integrated energy splits 1:3 between frequency shift and removed photons to 10⁻⁹.

## Declared labels

For each in-regime combination, at B1:
- **"Present radiation can hold the requirement"** if u_CMB + u_UVmm(170) ≥ E_req.
- **"FIRAS permits law (a)"** or **"(b)"** if the energy fraction that law can transfer within Δχ² ≤ 4 covers E_req.

Once, for the laws themselves:
- **"Law (c) needs frequency selectivity"** if, applied to all beams, its extra dimming at 1 Gpc exceeds 1%.

The report then states which radiative routes remain for stage 2B: the law and history each would need, or none.

## Not claimed

- A thermalization mechanism, or the origin of the microwave background.
- Official FIRAS limits. There is no covariance, calibration or foreground likelihood.
- A source interaction or its kinematics. Those belong to stage 2B; so does momentum conservation, which forbids one photon from becoming one slow massive companion.
- Any change to the existing redshift or brightness fits. Dimming is reported as a cost, not refitted.

## Files

`rc1.py` (the calculation, which is also a fast suite job), `rc1-results.json`, and `report.md`.
