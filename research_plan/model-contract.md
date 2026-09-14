# Model contract: one clock field for redshift, gravity and the reservoir

**Status.** Proposal of 14 September 2026, written after PF-1, RPG-1 and CR-1. It is a plan, not a result.
- Nothing here is fitted or promoted.
- Decision D1 belongs to the project owner.
- Each experiment in section 5 needs its own declared protocol.

## 1. What the first implementations tested

Each branch's general hypothesis is broader than its first implementation. Each implementation failure below is preserved as declared, and none rules out the broader hypothesis.

| Branch | General hypothesis | Restrictions of the first implementation | Declared outcome |
|---|---|---|---|
| PF-1 ([report](../research_work/results/propagation-field/report.md)) | An evolving propagation field changes frequency, timing, photon number and brightness together, and exchanges energy with the light | Homogeneous n(t); matter left uncoupled by postulate; linear history; one law for light and gravitational waves | The wave law is consistent: one simulation gives all three stretches, and the energy exchange closes. A measured redshift still needs a matter coupling. The linear history's brightness scores Δχ² = +35 against flat FLRW |
| RPG-1 ([report](../research_work/results/radiation-polarized-gravity/report.md)) | Gravity's response depends on acceleration and on a radiation state s | s omitted; μ = x/(1+x) with the archived a*; lensing Φ = Ψ; PF-1 distances | Rotation matches algebraic simple MOND. The Milky Way vertical force is too strong at 1.1 kpc. SLACS Einstein radii are 0.35–0.74 of observed. With FLRW distances and published masses they are 0.48–0.89 (post hoc) |
| CR-1 ([protocol](../research_work/results/collective-reservoir/protocol.md)) | A collective reservoir, fed by photons, supports the extra gravity | One coherent spherical mode; no self-interaction; a 10⁹ Msun seed; present-day illumination collected inside r₉₉; a 10 Gyr benchmark | The derived supply is 10⁸–10¹¹ short within the seed's collection volume. Growing at the required rate, the heavy constituent contracts to a compact core, while photon-shaped growth runs away and disperses. See the CR-1 report |

## 2. The obstruction, stated once

Three results in the repository combine into a single constraint:
- the clock–ruler family ([clock-ruler-completion](../research_work/results/clock-ruler-completion/report.md));
- the tested matter couplings ([matched-wave](../research_work/results/matched-wave/derivation.md) §4, [shared interaction test](../shared_interaction_test/shared_interaction_report.md), [matter-clock-closure](../research_work/results/matter-clock-closure/report.md));
- PF-1's equivalence with coasting FLRW.

**The constraint.** In a homogeneous universe, a redshift that survives comparison with local atoms, with every clock ratio fixed, requires the number of atomic rulers between galaxies to grow. Any coupling that avoids this either changes clock ratios or local light speed, which local precision measurements exclude by orders of magnitude, or cancels the redshift.

**The only remaining alternative is inhomogeneity.** The field would evolve only where there is no matter (screening). But an index that rises away from galaxies bends light away from them. That is the opposite of lensing, and the contrast needed for the redshift is of order unity. The [void-screening](../research_work/results/void-screening/report.md) report already found a sign problem in the clock gradient.

## 3. The proposal: one clock field does all three jobs

A single field χ, the propagation field of PF-1, with n = e^χ, plays every role. It could be the phase of the companion condensate.

- **C1. Redshift (PF-1, completed).**
  - **The completion.** Light propagates with index n, exactly as in PF-1. Matter's rulers and clocks co-scale with the same field: the r = 2 member of the clock–ruler family, or equivalently static space with universally evolving particle masses (Wetterich, arXiv:1303.6878).
  - **Local physics.** Every dimensionless constant, every clock ratio and every locally measured light speed is fixed, so the redshift survives measurement: 1+z = n_o/n_e.
  - **Energy.** The light's lost energy flows into the field at (ṅ/n)u_γ, as PF-1 derived. The field's quanta are the companions.
  - **Geometry.** Distances become FLRW-like (D_L = (1+z)²D_A). That alone removes about 30% of RPG-1's lens deficit. The supernova brightness then fixes the history n(t), which must be FLRW-like (between coasting and de Sitter). That history must come from a declared field dynamics, not a fitted curve.
- **C2. Galaxy gravity (RPG-1, derived rather than postulated).** Around matter, the same field carries a static perturbation.
  - **Kinetic function.** Take a Lagrangian that is linear at large gradients and cubic near the cosmological rolling, of the k-essence or superfluid type (compare Berezhiani & Khoury, arXiv:1507.01019).
  - **The static equation.** The static field equation is then RPG-1's AQUAL equation. The transition sits where the static gradient matches the rolling rate, so a* = γ c |d ln n/dt|. RPG-1's ξ = 0.118 becomes γ, a constant of the kinetic function.
  - **Prediction.** a*(z) ∝ H(z), with H(z) the same history that the supernovae fix. That turns the a0 ~ cH0 coincidence into a test.
- **C3. Light bending.**
  - **The requirement.** Weak lensing shows the same radial-acceleration relation as rotation curves, so at low acceleration lensing must equal dynamics. A conformally coupled scalar bends no light ([motion-and-lensing](../research_work/results/gravity-response/motion-and-lensing.md)).
  - **The construction.** The field must therefore enter the photon light cone through its own rolling frame, the timelike direction of its gradient. This is the structure TeVeS and AeST (arXiv:2109.13287) use with a separate vector. It is a requirement to be derived and checked; nothing assumes it holds.
- **C4. The reservoir (CR-1, reassigned).** The field's condensate supplies mass only where a response cannot.
  - **Galaxy cores.** At SLACS Einstein radii g_N ≈ 10 a*, so any acceleration-keyed response adds only 15–18%. With FLRW distances, the lenses need 2.1–2.7× Chabrier (1.2–1.5× Salpeter), comparable to the archived masses fitted to the same lenses' inner stellar motions.
  - **Clusters.** They need about 2× more than MOND gives.
  - **The job.** A compact core supported by self-interaction or pressure, not a galaxy-wide halo. CR-1 shows that a single coherent mode cannot be the halo.
- **C5. The radiation state s.** The field's local state responds to both of its sources: radiation through u_γ/n, and matter through the co-scaling coupling.
  - **Baryons dominate locally.** Baryon rest-energy density exceeds starlight energy density by about 10⁶–10⁷, so the galaxy-scale field is keyed to baryons, as the radial-acceleration relation requires.
  - **A pure starlight key fails.** RPG-1 showed it would shift the transition 4.6-fold between gas-rich and gas-poor galaxies.
  - **Radiation's role is cosmological.** It sets the clock, not a local switch. Any local s-dependence must be a small, declared correction.

**Provenance.** Every ingredient exists in published work:
- static-frame cosmology;
- k-essence and superfluid MOND;
- TeVeS/AeST lensing;
- a rolling (ghost) condensate.

The proposal is the combination: the redshift field, the galaxy-gravity field and the companion condensate are one field, with a* tied to the redshift rate. Its originality is unverified.

## 4. Decision for the owner

**D1. What does "nonexpanding" mean operationally?**
- **(a) Co-scaling.** Galaxies stay at fixed positions in the propagation field's coordinates while matter's rulers shrink. Counted in rulers, separations grow. Observationally this is identical to an expanding universe with scale factor n(t), and there is no Big Bang singularity in this frame.
- **(b) No ruler growth.** Ruler counts must not grow. Then C1 is excluded, and the first experiment becomes the screening lensing-sign test, which is likely fatal.

Recommendation: (a), because it is the only option that the repository's constraints leave consistent.

## 5. Experiment queue (each needs its own declared protocol)

1. **CC-1, the clock completion.**
   - **Setup.** The minimal action: the field χ with PF-1's photon coupling and co-scaling matter.
   - **Outputs.**
     - the received frequency against the local atomic reference;
     - the ratios between atomic clocks;
     - the drift of local light speed and cavities;
     - event durations;
     - arrival times of light against gravitational waves;
     - D_L, D_A and the distance-duality ratio.
   - **Brightness last.** Only then compare supernova brightness at equal freedom, with n(t) from a declared field dynamics.
2. **CG-1, clock-field gravity.** The static galaxy solution for a declared kinetic function.
   - **Derived quantities.** The interpolating function, a* = γ c H, and the lensing-to-dynamics ratio at low acceleration.
   - **Tests.**
     - SPARC rotation;
     - the Milky Way's rotation and vertical force, including the lighter-disk baseline;
     - SLACS with the C1 geometry;
     - the data needed for the a*(z) test (high-redshift rotation curves).
3. **CR-2, the condensate core.** CR-1 with a repulsive self-interaction and a declared thermal population.
   - **Its job.** The SLACS inner mass and the cluster excess.
   - **Separations.** Supply is reported separately from support. The formation history and the collection region are labeled assumptions.
4. **Population comparisons,** only after 1–3.

## 6. What would stop the proposal early

- **CC-1.** Co-scaling cannot fix every measured ratio while a declared history fits the brightness at equal freedom, or D1 is (b).
- **CG-1.** No kinetic function gives both the radial-acceleration relation and lensing equal to dynamics at low acceleration with a* = γ c H. Or the Milky Way vertical force overshoots for every declared baryon model.
- **CR-2.** A core that supplies the SLACS mass spoils the inner rotation curves.
