# Model contract: the co-scaling branch and the open couplings

**Status.** Revised 14 September 2026 after the project owner's review. It replaces the proposal of the same day (15bd1f6, corrected in 3e2cf0e).
- **The owner's decision.** Co-scaling is investigated as a separately labeled branch (CC-1). That is not a declaration that the original nonexpanding requirement is satisfied, nor that other alternatives are excluded. A genuinely fixed-ruler, nonexpanding branch stays open and conceptually separate.
- **Nothing here is fitted or promoted.** Results are marked as results; everything else is a plan. Each experiment needs its own declared protocol.

## 1. What the implementations tested

Each branch's general hypothesis is broader than its implementation. Each failure below is preserved with its stated scope and does not rule out the broader hypothesis.

| Branch | General hypothesis | Restrictions of the implementation | Declared outcome |
|---|---|---|---|
| PF-1 ([report](../research_work/results/propagation-field/report.md)) | An evolving propagation field changes frequency, timing, photon number and brightness together, and exchanges energy with the light | Homogeneous n(t); matter left uncoupled by postulate; linear history; one law for light and gravitational waves | The wave law is consistent: one simulation gives all three stretches, and the energy exchange closes. A measured redshift needs a matter coupling. The linear history's brightness scores Δχ² = +35 against flat FLRW |
| RPG-1 ([report](../research_work/results/radiation-polarized-gravity/report.md)) | Gravity's response depends on acceleration and on a radiation state s | s omitted; μ = x/(1+x) with the archived a*; lensing Φ = Ψ; PF-1 distances | Rotation matches algebraic simple MOND. The Milky Way vertical force is too strong at 1.1 kpc. SLACS Einstein radii are 0.35–0.74 of observed |
| CR-1 ([report](../research_work/results/collective-reservoir/report.md)) | A collective reservoir, fed by photons, supports the extra gravity | One coherent spherical mode. Labeled assumptions: zero self-interaction, a 10⁹ Msun seed, illumination collected inside r₉₉, a 10 Gyr span | Derived supply 10⁸ to 5×10¹¹ times short in the seed's collection volume. At the required rate, the heavy constituent contracts to a compact core; photon-shaped growth runs away or disperses. The ground-state reference used for its excess energies is checked in the [report](../research_work/results/collective-reservoir/report.md) |
| CC-1 ([report](../research_work/results/clock-completion/report.md)) | Co-scaling: matter and light share one metric whose scale factor is the propagation index | Homogeneous; universal coupling to g_m = −c²dt² + n²dx²; the index's dynamics from L_n = (M/2)ṅ² − V(n) | **Consistent.** 1+z = n_o/n_e in both frames, all clock ratios fixed, D_L = (1+z)²D_A, the ledger closes. With V = 0 the history is coasting, and ruler counts grow. See section 3 |
| CR-2 ([report](../research_work/results/supported-reservoir/report.md)) | A repulsively supported reservoir supplies the lenses' extra mass | Thomas–Fermi n = 1 polytrope with one shared constant; Newtonian gravity; six SLACS lenses fitted jointly (motions plus exact lensing) with identical geometry and mass conventions; Milky Way check | **Fails.** The best shared core is about 80 kpc in both geometries. χ² 113.6 against free NFW's 85.3 (FLRW); Milky Way RMSE 22.5 against a limit of 20 km/s. Supply is 10⁶–10⁹ short. The lenses need 1.5–3 times Chabrier stellar mass |
| CF-1 ([report](../research_work/results/gravitational-focusing/report.md)) | One slow companion population is focused and retained differently in shallow and deep potentials | Conditional; potentials from baryons plus a 1% seed; one declared interaction (elastic companion scattering); trial bath speeds; results per unit incident density | **Depth-keyed, conditionally.** At 300 km/s Coma retains 140–290 times faster than the Milky Way and stays bound. The Milky Way overheats. Faster baths erode every baryonic well. Supply is extreme |
| CG-0 ([report](../research_work/results/clock-gradient/report.md)) | The clock field's gradient can serve as the rolling frame for galaxy gravity and lensing | One scalar; three formulations of its static response | **Fails as a single field.** The gradient turns spacelike at every SLACS Einstein radius, inside 2 kpc in the Milky Way and in the Solar System. Only a two-field (bi-potential) form keeps the frame defined |

## 2. The constraint within the tested family

Three results combine:
- the clock–ruler family ([clock-ruler-completion](../research_work/results/clock-ruler-completion/report.md));
- the tested matter couplings ([matched-wave](../research_work/results/matched-wave/derivation.md) §4, [shared interaction test](../shared_interaction_test/shared_interaction_report.md), [matter-clock-closure](../research_work/results/matter-clock-closure/report.md));
- PF-1's equivalence with coasting FLRW.

**Within the homogeneous universal-metric family tested**, a redshift that survives comparison with local atoms, with every clock ratio fixed, requires the number of atomic rulers between galaxies to grow. Couplings in that family that avoid this either change clock ratios or local light speed, which precision measurements exclude, or cancel the redshift.

**Co-scaling is one completion compatible with the kinematics of the tested family. Other completions have not been exhaustively classified.** Families not classified include:
- inhomogeneous or screened fields (the [void-screening](../research_work/results/void-screening/report.md) report found a lensing-sign problem in one version);
- non-universal couplings;
- couplings outside a single metric, including frequency- or photon-number-changing processes, which are the original companion mechanism.

## 3. The co-scaling branch (CC-1 results)

CC-1 is a consistency calculation, not a brightness fit. Protocol 6ce4c7a; [report](../research_work/results/clock-completion/report.md).

- **The coupling decides the clock factor.** In 1+z = (n_o/n_e)·ν_atom(t_o)/ν_atom(t_e), the clock factor is derived from the atoms' own dynamics in g_m. Classical and quantum atoms tick at a constant rate in t, to 5×10⁻⁸ and 2×10⁻⁸. The conformal frame puts the whole shift into the clocks instead. Both frames give 1+z = n_o/n_e to 7×10⁻⁹, so nothing is counted twice.
- **Local measurements.**
  - Clock ratios are constant.
  - Light speed in local units is c.
  - A cavity's round trip is 2ℓ/c to 3×10⁻¹⁵ while n doubles.
  - Pulse durations stretch by 1+z.
- **Gravitational waves.**
  - They arrive with the light under either law, because both laws share the metric's light cones.
  - As standard tensor waves of g_m, their amplitude falls as 1/n, giving a siren distance equal to D_L.
  - Under PF-1's postulate P2, the amplitude does not fall, giving D_L/(1+z).
- **Energy ledger.** The field gains exactly what radiation and free motion lose. For the homogeneous toy, energy is conserved to 2×10⁻¹⁴ and the ledger closes to 9×10⁻¹⁵. Bound systems exchange nothing on average. Rest mass exchanges nothing.
- **History, derived rather than imported.**
  - The distance history follows from n(t), which the field's dynamics must supply. Shrinking rulers select no particular history, and the flat-FLRW distances with H0 = 70 and Ωm = 0.3 used in RPG-1's post-hoc diagnostic are an adopted comparison, not a consequence.
  - With V = 0, the declared dynamics turn around at a finite n_min with z_max = E_field/E_radiation. This is a property of that toy, not an established nonsingular past. A potential falling faster than −W/n as n → 0 would remove the turnaround. In atomic units the turnaround is as hot and dense as the corresponding epoch of a hot early universe.
  - A turnaround at least as hot as recombination needs E_field ≥ 1100 E_radiation. That makes today's history coasting (|q₀| ≤ 5×10⁻⁴), whose D_L(z) is the one PF-1 scored at Δχ² = +35.
  - If the field's energy gravitates, the implied density is at least 5.5% of 3H²c²/8πG for recombination, and 1.5×10⁴ times it for nucleosynthesis.
- **What CC-1 does not settle.** It gives no complete gravitational action, no history other than the V = 0 toy, no reservoir supply, and no claim about the original nonexpanding premise, which co-scaling does not meet: separations counted in atomic rulers grow.

## 4. The unification goal (a goal, not a result)

**Goal.** One field sets the redshift, the galaxy-scale gravity and the reservoir, with shared parameters. There is no complete action yet, and the parts tested so far do not connect:
- **Galaxy gravity.** [CG-0](../research_work/results/clock-gradient/report.md) shows one scalar cannot serve as both the rolling frame and the galaxy force. The frame's normalization X = (χ̇/c)² − |∇χ|² crosses zero at 8.5 a*, just above the transition. A completion needs one of:
  - an independent timelike vector, with the AeST structure (Skordis & Złośnik, arXiv:2007.00082);
  - a bi-potential construction sourced through the Newtonian field.

  Its interpolating function must also meet the Solar System bound; the simple ν does not.
- **Lensing.** It needs a derived coupling to the light cone. RPG-1's post-hoc diagnostic compared FLRW-geometry lens predictions with archived dynamical masses fitted in a different, regular-optics geometry. That comparison is a reason for a joint fit, not evidence that lensing and dynamics agree.
- **The reservoir.** CR-1's failures stand with their stated scope.
  - CR-2's supported condensate, with one shared scale of about 80 kpc, fits the lenses worse than free halos and misses the Milky Way threshold.
  - Its supply from local starlight fails by 10⁶–10⁹.
  - If reservoirs exist, they need a cosmological supply (CC-2). In the Milky Way they also need an interior more concentrated than a Thomas–Fermi core with the lenses' scale: the core is too uniform across 5–25 kpc, and its edge plays no role (CR-3, after CF-1).
- **The radiation state s.** Whether matter or radiation controls the local field state must come from the coupling equations. The ratio of baryon to starlight energy density alone does not establish it.

**Provenance.** Ingredients exist in published work:
- static-frame cosmology (Wetterich, arXiv:1303.6878), whose representation is a field redefinition that changes no dimensionless ratio, so its statements about singularities do not transfer without the dynamics;
- k-essence and superfluid dark matter (Berezhiani & Khoury, arXiv:1507.01019);
- TeVeS/AeST lensing;
- rolling condensates.

## 5. Experiment queue (each needs its own declared protocol)

1. **CC-2, cosmological supply.** Given a declared field state and interaction, what traveling population is produced, and how does extracting that energy change the field's evolution?
2. **CF-2, galaxies' energy sink.** CF-1 found that galaxies' captured companions overheat without one. Declare a sink or a second interaction, and derive the equilibrium profile a retained cluster population reaches, for comparison with lensing under identical conventions.
3. **CR-3, a revised supported reservoir,** designed from what CF-1 and CC-2 require. Support, supply and formation history are reported separately.
4. **CG-1, clock-field gravity,** only with a named vector or second field. It tests the Milky Way radial and vertical forces first, states a Solar System–safe interpolating function, and then SPARC and SLACS.
5. **The co-scaling history.** A declared V(n) or other dynamics comes before any brightness comparison. No distance history is imported.
6. **The fixed-ruler branch.** Its first experiment remains whether a screened field can bend light toward galaxies.

## 6. What would stop each line early

- **Co-scaling.** A declared field dynamics that reaches the recombination and nucleosynthesis temperatures cannot fit the brightness at equal freedom, or requires field energy that gravitating tests exclude.
- **CG-1.** No vector or bi-potential completion gives the radial-acceleration relation, lensing equal to dynamics at low acceleration, the Milky Way vertical force and the Solar System bound together.
- **CR-3.** No shared-dispersion two-phase law fits the lenses jointly within the declared allowance of free halos while meeting the Milky Way threshold. (CR-2's single-phase law already failed this.)
- **CF-2.** No declared energy sink lets a galaxy-scale well keep the companions it captures without an energy budget that its baryons or radiation could not carry.
