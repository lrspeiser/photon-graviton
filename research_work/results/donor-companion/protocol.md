# RC-2, part 1: a time-dependent donor-and-companion calculation in an open region

Declared before execution, 15 September 2026. Baseline: `main` at 70ca39d.

**Why it is run.** The owner's reviews of ed96b00 and 9232e07 made this the highest-priority physical experiment.
- **What 2B-F1 found.** Stage 2B-F1's cold decays pass the Milky Way's profile gate. The robust part is the inner shape. At 3 km/s, without scattering and without bath gravity, the rotation error is 8.0 km/s and the 8–20 kpc slope 1.32. Low-angular-momentum material falls inward from an extended production region. The owner's review of 9232e07 makes cold collisionless decay the main physical reference case.
- **What rests on assumptions.** The outer mass distribution rests on a zone of influence R_b, with an incident static bath at the boundary standing in for the companions born outside it. The result depends strongly on R_b: with a 150 kpc zone, no sampled rate passed.
- **The donor, counted once.** With bath gravity on, 2B-F1's gravitating extra density is the bound companions plus the incident bath minus ρ∞ = qt. For a homogeneous, stationary donor that already reads as daughter mass minus the donor's depletion. So this calculation does not add a second −qt. It removes the bath together with its subtraction, and puts in their place the daughters followed in time, bound and unbound, minus the donor's depletion.

This calculation drops the zone and the bath. It evolves the donor's depletion and every daughter together, in an open region. It asks: with the source law and its parameters held fixed, do the predicted masses and velocities at fixed observational radii converge as the computational region grows, or does the boundary choose them?

It is the first part of RC-2 as the status queue (item 5) defines it: the finite-time open region, with the computational boundary kept apart from the observational apertures. The compact topologies are the second part.

## The system, the source and the donor

- **The system.** The Milky Way's baryons: the archive's baseline I, as in 2B-F1.
- **The source.** 2B-F1's field: homogeneous, at rest in the host's frame, decaying into companion pairs at rate q per unit volume. Each companion leaves with speed v_d in an isotropic direction.
  - The source has no physical edge. R_comp is a numerical boundary, never a source extent. A source with a physical extent belongs to the next item.
- **The donor** is stationary.
  - **Its initial density** is uniform and belongs to the static background under the universe contract's global law, which this calculation does not model.
  - **What gravitates** is the change from that initial state: the baryons, every daughter, and the donor's depletion. The depletion is −q t per unit volume, because conversion removes donor mass wherever a companion is born.
- **The parent's stock** is not modeled. The report gives the parent density that depletion fractions ΓT of 1%, 10% and 100% would imply.

## The region, the apertures and the physics

- **An open region of radius R_comp.**
  - Companions are born at rate q everywhere inside it.
  - **Every birth is followed, bound or unbound.** 2B-F1 removed births that were unbound at once, but at 3 km/s such a companion needs about 16 Gyr to move 50 kpc, so it stays in the region.
  - A companion that crosses R_comp outward leaves, and its mass and energy are booked as escaped. Nothing enters from outside.
  - There is no incident bath and no zone of influence.
- **Region sizes:** R_comp = 150, 300, 600, 1,200 and 2,400 kpc.
- **Fixed observational apertures.**
  - The 38 Eilers radii (5–25 kpc), scored with 2B-F1's G1 and G2.
  - Outer apertures at 50, 100, 150, 200 and 300 kpc, wherever they lie inside R_comp. Their circular speeds follow from the net enclosed mass.
- **Span:** 10 Gyr, with snapshots at 2, 5 and 8 Gyr.
- **Collisionless** (σ/m = 0).
- **Decay speeds.** 3 km/s is the reference case; 10 km/s is the comparison.
- **The engine.** Stage 2A's orbits in spherical symmetry, with the self-gravity of the whole perturbation (baryons + daughters − depletion).
- **Resolution, held fixed as the region grows.**
  - **The grid and the birth shells** keep the log spacing the engine uses for a 300 kpc region: 1,024 nodes and 40 shells there.
  - **Tracer masses** follow one function of radius in every region. It is 2B-F1's rule: a shell's tracer mass is proportional to the square root of its production, within a factor of 10 of the heaviest. The rule is evaluated for a 300 kpc region in which births total 16,000 over the span, and beyond 300 kpc the square-root rule continues. So a larger region adds births outside without thinning those inside, and the engine never has to thin the population.
  - **The report** gives the tracer count inside each aperture.
- **Budget.** Each run has a wall-clock budget of 4 hours. A run that exceeds it is recorded as not completed, and a region with a missing seed cannot enter the convergence test.
  - A cost probe before this declaration, at a rate that is not declared here, ran 10 Gyr in 1, 4.5 and 6 minutes at 300, 1,200 and 2,400 kpc. It removed unbound births at once, as 2B-F1 did. The machine was loaded, and no probe thinned its population.

## Normalization

- **The fixed rate.** The rate is held fixed across regions and is never refitted for a region. It is 2B-F1's best sampled rate for the matching combination: the same v_d, collisionless, with bath gravity on, since that is the bookkeeping (daughters minus depletion) that this calculation makes explicit.
  - 1,437.8 M☉ kpc⁻³ Gyr⁻¹ at 3 km/s, where 2B-F1's RMSE was 6.80 km/s.
  - 1,466.1 at 10 km/s, where it was 8.50 km/s.
  - These are also the centres of the 2B-F1 revision's brackets for those combinations.
  - The reference case's rate without bath gravity (1,361.3, RMSE 8.0 km/s) lies inside the search's bracket.
- **At 1,200 kpc, a direct search as in the revision.**
  - Rates q₀·3^(e/8) for e = −4, −2, 0, 2 and 4, three seeds each. The fixed-rate runs serve as e = 0.
  - Where the best lies at an edge, up to two outward extensions of two steps each.
  - Then the best rate's two odd neighbours.
  - The selected rate has the lowest seed-mean RMSE, and the gates are applied to its seed means.

## Scoring (declared)

1. **Boundary convergence.** Three seeds per region at the fixed rate. For each doubling (R_comp to 2R_comp), two seed means agree when they differ by less than two standard errors of the difference, or by less than a floor, whichever is larger:
   - the 38-bin RMSE: 1 km/s;
   - the 8–20 kpc slope: 0.05;
   - the companion mass inside each outer aperture that both regions contain: 5%. The donor's depletion inside an aperture is the same in every region, so this also tests the net contrast and the circular speed.

   R_conv is the smallest region from which every later doubling agrees. R_conv,profile is the same for the RMSE and the slope alone. The boundary test passes if R_conv exists within the tested range.
2. **The profile gate at 1,200 kpc.** G1 and G2 as in 2B-F1, on the search's seed means.
   - There is no cost gate: the three quantities below are reported instead.
   - The outer net mass is compared with the roughly 10¹² M☉ inside about 200 kpc that the Milky Way's halo tracers indicate. That is a comparison, not a gate.
3. **Three quantities, reported separately at every aperture.**
   - **Positive inventory:** the companions inside.
   - **Gravitational mass contrast:** the companions minus the depletion, with the baryons listed beside it.
   - **Source energy:** the field's power density, and the parent stock at ΓT = 1%, 10% and 100%.

   Also reported:
   - the donor's depletion inside each aperture;
   - the net inflow across each aperture at each snapshot, which is the companions' flux integrated over time;
   - the mass and energy that escape through R_comp;
   - the fraction of births bound at birth, and the fraction still inside the region at T.
4. **Where the contrast comes from.** At each outer aperture the net contrast splits exactly into two parts:
   - the births inside the aperture minus the donor's depletion inside it, which is zero apart from sampling noise, now that every birth is followed (control C2);
   - the net inflow of companions across it: the companions inside minus the births inside.

   So the net contrast is the mass carried in from outside the aperture: the identity in the owner's review. The births inside each aperture are counted as they are drawn.

## Controls (declared)

- **C1, homogeneous production with no central perturbation.** The same source in the 2,400 kpc region with no baryons: three seeds at 3 km/s and the fixed rate.
  - **What it tests:** whether the method (its boundary, and its tracer noise) manufactures a central mass excess.
  - **It passes** if two things hold at T, against the Milky Way runs at the same region and rate:
    - the magnitude of its net contrast inside each outer aperture stays below 10% of the companion mass inside that aperture;
    - its companion mass inside 25 kpc stays below 10% of theirs.
  - **Its growth over time** is reported either way.
- **C2, conversion without spatial redistribution.** In every run, at every aperture and snapshot, the births drawn inside the aperture must match the donor's depletion inside it within four Poisson standard errors. If the daughters stayed where they were born, the net contrast would vanish: converting the donor does not double-count local gravitating mass.
- **C3, the source switched off after a declared production episode.** At 3 km/s in the 1,200 kpc region, production runs at the fixed rate for 5 Gyr and then stops, so the depletion stops growing too. The run is evolved to 10 Gyr, with three seeds.
  - **Reported:** the inner RMSE and slope and the aperture contrasts at 5 and 10 Gyr, beside the continuous runs' at the same times.
  - **It asks** whether the profile persists, relaxes, or depends on continuing infall. It is not a gate.

## Expectation, stated before running

- **The depletion cancels the uniform part of the daughters' gravity.** A daughter medium born at rest and uniform, with its donor depleted uniformly, exerts no net force. What pulls is the baryons, plus the daughter mass already carried inward.
- **The daughters' own instability.** Over a stationary donor, the daughter medium's density fluctuations still gravitate.
  - At 3 km/s its Jeans length is about 12 kpc.
  - At the fixed rate, fluctuations on larger scales grow by about 6 e-folds over 10 Gyr.
  - So tracer noise may seed a radial collapse, and C1 measures it. If C1 fails, part of the spherical result is noise-seeded collapse, and the three-dimensional test becomes decisive.
- **The infall horizon.** Matter at rest reaches 25 kpc by T from 334 kpc in the baryons' potential alone. It would come from about 720 kpc if 10¹² M☉ were inside from the start. The inflow builds up over the span, so the horizon lies between.
- **Convergence.**
  - The 5–25 kpc profile should converge by 600 kpc, or by 1,200 kpc at the latest.
  - An outer aperture needs a region beyond its own horizon: up to about 800 kpc for the 300 kpc aperture. The 2,400 kpc region tests whether the outer apertures converge.
  - 2B-F1's 300 kpc zone lay inside the horizon and its 150 kpc zone far inside it, which would explain its zone dependence.
- **The inner profile** at the fixed rate should resemble 2B-F1's collisionless runs with bath gravity: an RMSE of 6.8 km/s at 3 km/s and 8.5 km/s at 10 km/s.
- **The outer contrast.** In 2B-F1's bath-gravity accounting, the net non-baryonic contrast reached about 9×10¹¹ M☉ inside 150 kpc, then fell to about 1×10¹¹ M☉ at the 300 kpc zone's edge: the mass carried inward left a deficit behind it. In an open region, infall from beyond 300 kpc should partly fill that deficit as the region grows.

## Validation

- **D1, ledgers,** in every completed run.
  - The companions inside R_comp plus those that escaped equal the births, to 10⁻⁹.
  - The donor's depletion inside R_comp equals the mass the field converted, qT times the region's volume, to 10⁻⁹.
  - The births drawn agree with their expectation within four Poisson standard errors.
  - The energy ledger closes to 10⁻⁹, as in 2B-F1's V6.
- **D2, the depletion's gravity.** The grid potential of the uniform depletion matches the analytic uniform sphere to 10⁻⁴ of the central depth, which is the grid's own trapezoid error. A check made before this declaration found 4×10⁻⁵ at both 300 and 1,200 kpc.
- **D3, cold infall.** The revision's F3 covers the field channel's births and orbits in a frozen potential. It is cited here, not repeated, as the statistical check the revision's amendment describes.
- **D4, the code path.** Run with the depletion off, with unbound births removed at once, and with 2B-F1's mass rule, the donor model and 2B-F1's field model without collisions or bath gravity must give bit-identical results from the same seed and rate.

## Reported predictions

J1630 and Coma run unchanged, at the selected rate of the better-fitting decay speed.
- **Their regions** are max(2, R_conv/r_h,MW) times each system's own infall horizon r_h: matter at rest falling onto the system's baryons, taken as a point mass, within T. The horizons are 332 kpc for the Milky Way, 529 kpc for J1630, and 3.7 and 4.6 Mpc for the two Coma models.
- **Their tracer masses** follow the same rule, referenced to the smaller of the region and their 2B-F1 zone. Coma's zone is 30 Mpc, so a region of a few megaparsecs referenced to the zone would hold only an eighth of its births.
- **Their companion mass and net contrast** at stage 2A's test radii are reported against their baryons, as predictions, not gates. For Coma this replaces the old question, whether a static bath exists, with a new one: what distribution develops from the same source.

## Not claimed

- **The donor's own dynamics.** What keeps the donor distributed as assumed, how its energy and stresses gravitate, and what supports the nonexpanding background all remain open, as the owner's review of 9232e07 notes.
- **Stability or fragmentation in three dimensions.** Spherical symmetry keeps only the radial mode of the daughter medium's instability.
- **A common field frame.** The source is at rest in the host's frame. Host–field relative velocities are a later test.
- **A localized or baryon-triggered source.** That is the next item, after a source-response feasibility calculation.
- **Compact topologies** (S³, T³) and the global background's law. Recirculation, and the crossing time of material that returns, belong to that second part of RC-2. In an open region, what leaves never comes back.
- **A microscopic donor,** or its parent's equation of state.

## Files

In a new directory, `donor-companion/`:
- `protocol.md`;
- `open_region.py`, the model: stage 2A's engine with the depletion term and an open boundary;
- `rc2a.py`, the driver;
- `checks.py`, the suite job;
- `rc2a-results.json`;
- `report.md`.

## Amendment 1 (15 September 2026), after the owner's review of 1585c9f

Declared before any RC-2a result was read. The first canonical attempt started at 12:48 and was stopped during its first round to apply this amendment. It wrote nothing, and its partial log is not used.

### B1. Scope and order

- **One decay speed:** 3 km/s, the reference case. The 10 km/s comparison is dropped.
- **Three stages, each declared now.**
  1. **Controls, and the reference at the fixed rate.**
     - C1a–C1c (B2), C3 and C3b (B5), and D4.
     - The region-size comparison at 1,437.8 M☉ kpc⁻³ Gyr⁻¹ over 150–2,400 kpc, with three seeds per region.
     - C2 is evaluated in every run.
  2. **The direct search at 1,200 kpc,** as declared, but only if stage 1's controls are interpretable: C1a passes, and the reference demonstrates convergence across resolutions (C1c). Otherwise RC-2a stops after stage 1, and the representation is fixed first.
  3. **At the selected rate, frozen,** the region-size comparison over 600, 1,200 and 2,400 kpc, with three seeds each and the rules of B3. The rate is not refitted. Convergence at the old rate does not carry over, because a different source strength changes the mass profile and the distances from which material arrives.
- **J1630 and Coma move to the item after RC-2a,** at the frozen rate, and only if stage 3 converges.

### B2. The representation, and three tests of sampling noise

Over a stationary donor, the daughter medium is Jeans-unstable above about 12 kpc, so sampling noise could seed structure. The tracers should represent the proposed source. Their large individual masses should not silently set its fluctuation spectrum.

- **The quiet representation, used for the reference.**
  - Each birth shell receives its converted mass exactly, every step. That mass is a whole number of tracers, with the remainder carried forward, rather than a Poisson draw.
  - Radii are stratified in volume within the shell.
  - Births come in pairs with opposite radial velocities and the same tangential speed. Their direction cosines are stratified within each step.
- **Which runs use it.** The reference, C1a, C1c, C3, C3b and stages 2 and 3 use the quiet representation. C1b and D4 use 2B-F1's Poisson draws.
- **C1a, the numerical null.** The quiet representation, homogeneous production and no baryons, in the 2,400 kpc region, with three seeds.
  - **It passes** if two things hold at T, against the reference runs in the same region at the same rate:
    - its net contrast inside each outer aperture, taken as a root-mean-square over the seeds (a noise excursion can have either sign, so a seed mean could hide it), stays below 10% of the companion mass inside that aperture;
    - its companion mass inside 25 kpc stays below 10% of theirs.
- **C1b, the numerical-noise test.** The same setup with ordinary random sampling (Poisson counts, random radii and directions), at 16,000 and 64,000 reference births, with three seeds each.
  - **Reported,** as root-mean-squares over the seeds: the net contrast inside each aperture, and the excess inside 25 kpc over the uniform expectation.
  - **How they change with resolution.** Growth seeded by noise should fall as the square root of the number of births, so four times the births should give a ratio near 0.5.
- **C1c, a physical perturbation.** The Milky Way's baryons, declared and fixed, are the perturbation. The reference configuration runs in the 1,200 kpc region at 16,000 and 64,000 reference births, with three seeds each. The two resolutions must demonstrate convergence under B3's rules.
- **A delta-f representation** is not implemented here. It would treat the uniform component analytically and sample only the deviations.

### B3. Convergence of the gravitational signal, with a precision requirement

These rules apply to each doubling in stages 1 and 3, and to C1c's two resolutions.

- **Reported separately:**
  - **Positive inventory:** the companion mass inside each aperture. Its tolerance is 5%, as declared; it is reported but does not decide convergence. At 300 kpc a 5% inventory error would be about 70% of the net contrast.
  - **The net mass contrast inside each aperture,** with an absolute tolerance: the mass that would change the aperture's circular speed by 2 km/s, δM = 2 r v_c δv / G. This stays meaningful when the contrast is near zero.
  - **The circular speed at each aperture,** v_c = √(G (M_b + net) / r). Tolerance: 2 km/s.
  - **The rotation speed predicted at each of the 38 Eilers radii.** Tolerance: 2 km/s at every radius.
  - **The 38-bin RMSE and the 8–20 kpc slope,** with tolerances of 1 km/s and 0.05, as declared.
- **Two verdicts for each quantity:**
  - **Agree.** The seed means differ by less than two standard errors of the difference, or by less than the tolerance, whichever is larger. This is the declared rule.
  - **Resolved.** Two standard errors of the difference are below the tolerance, so the comparison could have detected a difference of that size.
- **Demonstrated convergence** for a doubling requires every signal quantity to agree and be resolved: the circular speeds, the net contrasts, the rotation speeds, the RMSE and the slope. A comparison that agrees without being resolved is reported as "no statistically resolved difference", never as convergence.
- **R_conv (demonstrated)** is the smallest region from which every later doubling demonstrates convergence. It needs an actual comparison with a larger completed region, so the largest region never qualifies by default. R_conv (statistical), the declared definition, is reported beside it.

### B4. Where the mass came from

Every tracer carries four values: its birth radius, birth time and birth specific angular momentum, and the time it first came within 25 kpc.
- **For the companions inside 10, 25, 100 and 200 kpc at each snapshot**, the report gives:
  - mass-weighted quantiles (10%, 50%, 90%) and histograms of birth radius (logarithmic bins) and birth time (1 Gyr bins);
  - the distribution of birth angular momentum.
- **The empirical infall horizon** is the 99th percentile of birth radius among the companions inside 25 kpc at T. It is compared with the pre-stated 334–720 kpc, which is a conditional estimate, not a bound on the self-gravitating run.
- **Where the tags live.** Stage 2A's engine carries them as an optional per-tracer array that follows every birth, escape and merger. The array is absent in stage 2A and 2B-F1 runs, which must stay bit-identical; their smoke runs and the suite check that.

### B5. A matched-budget comparison (C3b)

C3 stays as declared. C3b converts the same total mass as the continuous run, but early: 2q₀ for 5 Gyr, then nothing. It runs at 3 km/s in the 1,200 kpc region, with three seeds, and is evolved to 10 Gyr.

- **Reported beside C3 and the continuous run** at 2, 5, 8 and 10 Gyr: G1 and G2, the RMSE and slope, and the aperture contrasts.
- **Descriptive readings, declared now:**
  - **A persistent reservoir** if C3b passes G1 and G2 at 10 Gyr.
  - **Continuing infall** if only the continuous run passes at 10 Gyr.
  - **A transient match** if a history passes at an earlier snapshot but not at 10 Gyr.

### B6. Scope

Unchanged. The donor stays prescribed and stationary, and its initial density belongs to an unmodeled background. Its dynamics, the parent's equation of state and the nonexpanding background's law remain open. Passing RC-2a would establish a transport result, not a complete source or cosmological solution.
