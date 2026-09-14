# Capture-to-orbit feasibility: receiver-assisted threshold production

Declared before execution, 13 September 2026. Reviewed baseline: `main` at 3884b4f. This is a conditional feasibility test in a fixed potential, not a formation simulation, a self-consistent equilibrium or a collective-stability calculation.

## Question

Given the reference incident companion field and one declared local interaction, what bound population, orbital distribution and density develop? No observed rotation speed, fitted halo, MOND target or observed profile is used to construct the population. Every unhelpful outcome is retained: escape, plunging orbits, momentum drawn from ordinary matter, energy that leaves the galaxy, and failure to form an extended reservoir.

## Declared interaction (RB-1)

- **Receivers.** Ordinary baryons: stellar disk surface density from SPARC `SBdisk` times the archived disk mass-to-light ratio 0.5; spherical bulge mass from `Vbul` times 0.7; gas as a thin exponential disk of catalog mass 1.33 M_HI whose scale length is fitted to the SPARC `Vgas` contribution alone. That gas reduction uses ordinary-matter inputs only. The Milky Way uses its two archived baryon baselines. The reaction probability per unit baryonic mass is uniform. The local rate is proportional to rho_b(r) J(r), where J is the reference model's attenuated mean intensity (paper equation 6). Receiver rest energy greatly exceeds the companion energy.
- **Reaction.** c + R -> R + X, where X is a new species of one fixed mass m. Exact special-relativistic two-body kinematics, with threshold E'_th = m c^2 (1 + m/2M_R). The mass is never retuned per packet.
- **Near-threshold rule.** An s-wave matrix element (Wigner threshold law) with constant |M|^2: the cross section is proportional to p_f/p_i, and emission is isotropic in the centre-of-momentum frame.
- **Spectrum.** (S) The primary spectrum is smooth across threshold over fractional widths much larger than the receivers' Doppler widths. For the efficiency ledger only, the declared reference spectrum is flat in energy from threshold to twice threshold. (L) A monochromatic bath-frame line is treated analytically and checked with exact kinematics; it is not a separate population run.
- **Consequence to verify, not assume.** Under (S), bound production should be uniform within the galaxy-frame velocity ball |v| < v_esc(r) at each production site, independent of receiver velocity. Exact-kinematics Monte Carlo must confirm this before it is used.

## Potential and population

- **Primary potential.** The spherical monopole of the receiver baryons plus the reference exact-third companion reservoir, held fixed.
- **Labeled sensitivities.** (i) A baryons-only monopole. (ii) A stationary self-consistent fixed point in which the population's own gravity replaces the reference reservoir in the potential.
- **Steady state.** A phase-mixed orbit average of injected particles over (E, L) in the spherical potential.
- **Classification.** Bound means E < 0; escaping means E >= 0; plunging means pericenter < 0.05 R_d. Report radial periods and the bound mass on orbits with radial period above 10 Gyr, which cannot be phase-mixed over such a history. Universe age stays free.
- **Normalization.** The fixed retained inventory is the total bound mass, equal to the reference total inventory M_total, the same convention as the MOND-guided branch. The implied supply multiplier (absorbed over retained energy, under the declared reference spectrum) is reported, never used to rescale.
- **Prediction.** v^2 = v_b^2 + G M_pop(<r)/r at the measured radii, using the spherical monopole like the other companion models. No parameter is fitted and nothing is selected. Densities live on a model-defined radial grid with an explicit outer tail. They must not depend on which radii are requested (sampling-invariance check).
- **Optional labeled sensitivity.** Disk-launched orbits have a known orbital-plane distribution. If time permits, the in-plane speed of the resulting oblate density may be computed with the existing axisymmetric multipole solver.

## Frozen comparison runner

Compare side by side: ordinary matter, simple MOND (archived fitted a0), the original exact-third reference, the MOND-guided mixture (f = 0.9235570945357928, frozen) and capture-to-orbit (primary plus the two sensitivities).

- **SPARC.** All 149 galaxies with the original 89/29/31 splits; equal-galaxy RMSE (km/s) and log RMS; radial bins r/R_d < 1, 1 to 3, and >= 3.
- **Inventory-limited set.** The 21 inventory-capped galaxies as an explicit diagnostic set.
- **Milky Way.** Both fiducial baselines (R_d = 2.6 kpc, luminosity factor 1) and all 18 archived sensitivity scenarios, over the 38 Eilers bins.

## Energy and momentum ledger

Report each of the following:

- absorbed companion energy
- retained rest energy
- bound kinetic energy at injection
- orbital binding energy
- escaping rest plus kinetic energy
- receiver recoil energy
- momentum drawn from moving receivers, from the radiation-drag identity F = -(4/3)(P/c^2) v, which must be verified with exact kinematics
- receiver angular-momentum loss as a fraction of the disk angular momentum, at ideal efficiency and under the declared reference spectrum
- the minimum illumination history for which the bath's net radial push stays below 10% of local gravity

## Checks

1. Exact kinematics: per-event energy and momentum balance; the 4/3 drag coefficient at beta = 0.01 to 0.1; the uniform-ball bound injection and bound-fraction formula at artificially large v_esc/c.
2. Orbit averaging: the analytic Kepler time fraction; a stationary Plummer distribution function reproducing the Plummer mass profile.
3. Doubling every quadrature changes predicted speeds by less than 0.5 km/s.
4. Sampling invariance.
5. The reference reservoir reproduces the archived companion speeds.
6. The Plummer universal-K identity K = (G/6)(4 pi/3)^(1/5) M^(4/5) a^(-2/5), so fixed K implies a proportional to M^2. This checks the equations of the latest Coma support candidate.

## Promotion rule

The candidate is promoted only if all three conditions hold:

- the numerical checks pass;
- frozen predictions improve on the original reference in validation and test without degrading both Milky Way fiducials;
- the ledger closes with an available receiver momentum budget and a stated, bounded supply.

Otherwise it is reported as not promoted, naming the failing requirement. A numerical pass is not a physics pass.
