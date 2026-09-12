# Can captured deposits remain in place?

The depleted-capture profiles can be supplied with a conditional circular-orbit support construction. A strongly shielded profile cannot instead be supported by an isotropic population of bound massive particles. A positive pressure integral alone would miss this failure. These results concern a frozen spherical snapshot; they do not explain how capture assembles it or establish collective stability.

## Model scope and provenance

Retain p6,C1,R100 and opacity normalizations0.1,10,1000 from the [depleted-supply calculation](../depleted-capture-feedback/report.md). Reconstruct each deposited density using the same source equations. At exposure1, switch off incident supply for a snapshot-support diagnostic. This does not model continuous capture together with settling. Traveling companions are not assigned massive-particle behavior here; only the retained state is conditionally postulated to behave as nonrelativistic bound particles.

All units remain G=M_b=a=Wstar=1. Density, squared speeds and mechanical energies have units M_b/a^3, G M_b/a and G M_b^2/a respectively. The ordinary spherical Plummer potential remains fixed, and deposited gravity is included. No observed-data fit, independently tested galaxy, or dark-halo inference is introduced.

The statistical and orbital formulas below are known mechanics. The [galpy Eddington documentation](https://docs.galpy.org/en/latest/reference/dfeddington.html) describes the standard isotropic distribution-function framework and relative binding energy. Applying these constraints to our generated profiles is a conditional project calculation, not a claim to have invented the equilibrium relations. We do not perform a full Eddington inversion here; a hard computational cutoff and density interpolation require care before such an inversion could validate a complete distribution function.

## Isotropic particles: a necessary condition

Let Psi=-Phi be relative potential depth, with zero potential at infinity, and binding energy epsilon=Psi-v^2/2. For a nonnegative isotropic bound-particle distribution f(epsilon), known velocity-space integration gives:

`rho(Psi)=4 pi sqrt(2) integral_0^Psi f(epsilon) sqrt(Psi-epsilon) d epsilon`,

`d rho/d Psi=2 pi sqrt(2) integral_0^Psi f(epsilon)/sqrt(Psi-epsilon) d epsilon >=0`.

Positive enclosed mass makes Psi decrease outward. Thus rho cannot increase outward in such an equilibrium. The high-opacity profile violates this by orders of magnitude: shielding leaves little central deposited density while capture is larger farther out. This rules out the isotropic bound-particle interpretation of that profile, not all possible field states or anisotropic particle distributions. The low-opacity cases pass this necessary sampled-density check; that alone does not prove that their full f is nonnegative.

## Why a pressure integral is insufficient

Formal isotropic Jeans support, with pressure set to zero at the numerical boundary, is:

`P_J(r)=integral_r^R rho_D(u) g(u) du`, `sigma_r^2=P_J/rho_D`.

The integral is positive for every positive density, including the hollow profile. But isotropic bound particles must also satisfy `3 sigma_r^2 <= 2 Psi`, since their mean squared speed cannot exceed the squared escape speed. The hollow profile requires a much larger central velocity dispersion than binding permits. This is a second independent necessary-condition failure, not merely a preference for a different fit.

The boundary condition is part of this diagnostic, not a demonstrated physical edge. Full phase-space positivity, boundary behavior and stability remain open even when these inequalities pass.

## Circular orbits: an explicit alternative snapshot

In the frozen spherical potential, place each shell's particles on circular orbits, with orbit planes and phases uniformly distributed. This is anisotropic in velocity: radial velocity is zero, and each of the two tangential components has mean squared velocity v_c^2/2. Spherical symmetry and zero net angular momentum do not imply isotropic velocities.

Known circular-orbit relations give:

`v_c^2=r g`, `j_c=r v_c`, `P_r=0`, `P_theta=P_phi=rho_D v_c^2/2`.

The anisotropic radial force balance is exactly `dP_r/dr + (2P_r-P_theta-P_phi)/r = -rho_D g`. Thus the deposited radial profile can be maintained as an ideal circular-orbit ensemble in the frozen mean potential. Its circular binding energy `Psi-v_c^2/2` is positive in all sampled shells. The radial test-orbit frequency satisfies:

`omega_rad^2=[M_total(r)+r dM_total/dr]/r^3 >0`.

Positive frequency establishes linear radial stability of individual circular test orbits in this fixed potential. It does not establish collective stability of a self-gravitating ensemble, resilience to a Galactic bar, or survival during ongoing capture. Singular circular-orbit velocity distributions are an idealization, not a formation mechanism.

## Energy accounting required by support

The required orbital kinetic energy is `K_c=1/2 integral rho_D r g dV`. The formal isotropic Jeans energy is `K_J=3/2 integral P_J dV`. Independent shell quadrature verifies that these agree for the chosen zero-pressure boundary. Their agreement is the scalar virial identity, not proof that the isotropic distribution exists.

We also calculate deposited self-potential energy `U_self=-4 pi integral rho_D M_D(r) r dr` and interaction energy in the fixed ordinary potential `U_external=integral rho_D Phi_b dV`. U_external is counted once; U_self already contains the correct self-interaction factor. A formation budget must compare K+U_self+U_external with the initial state and any outgoing energy. Binding energy can contribute to the budget, so K is not automatically an unexplained extra energy deficit.

The previous absorption ledger counted incoming energy as deposited mass-equivalent density. If that density is reinterpreted as particle rest mass, its rest energy and these mechanical contributions must be reconciled. No physical value of c in these dimensionless velocity units or assembly history is chosen here. We do not claim that the kinetic energy, angular-momentum redistribution or dissipation needed for circularization has been supplied.

## Reproduction and numerical gates

The [protocol](protocol.md) predates calculation. Initial256/512-shell comparisons fail the1-percent gate for all three cases because of sampled outer pressure/density accuracy. Those failures remain in initial-results.json. Each affected case is refined to1024 shells, retaining its earlier comparison in results.json. The source evolution uses the unchanged fine ODE tolerances and8 impact nodes per annulus; inherited angular-systematic limitations are documented in the previous study. The prior weak-source angular refinement supports these parameters but does not replace a complete accuracy proof.

`run.py` reconstructs profiles and computes support diagnostics. Cached profiles are retained as compressed NPZ files; use `--recompute` to regenerate them from the capture equations. No star catalog is opened. The pressure and energy integrals use12-node quadrature within each uniform-density shell, including the exact shell contribution to gravitational acceleration.

## Implication for the active goal

There is a concrete candidate for retained particles: a bound, tangentially supported population whose gravity is computed from its density. It remains a snapshot construction. A coupled formation model must evolve positions and momenta as energy arrives, rather than deposit mass permanently at the capture point. Alternatively, a retained field state must supply an explicit stress and energy law. The hollow profile has ruled out one simple isotropic closure, not the user's broader companion hypothesis. Redshift generation, joint observables, the physical source budget and withheld validation remain incomplete.

## Final checkpoint

| kappa0 | Deposited mass | Required orbital K | U_self | U_external | Mean squared speed / escape squared at r0.1, formal isotropic closure |
|---|---:|---:|---:|---:|---:|
| 0.1 | 0.2318055 | 0.0397110 | -0.0109619 | -0.1111265 | 0.34655 |
| 10 | 0.2132688 | 0.0359508 | -0.0090283 | -0.1004521 | 0.36286 |
| 1000 | 0.0496456 | 0.0061318 | -0.0002638 | -0.0137818 | 6511.6 |

The high-opacity core ratio is approximately6512, whereas any isotropic population consisting entirely of bound particles requires a ratio no larger than1. Its density increases by approximately6803-fold between r0.1 and r3. These discrepancies are much larger than the numerical changes. The lower-opacity profiles pass the sampled necessary inequalities; their complete distribution functions have not been validated.

All three final512/1024-shell comparisons pass the unchanged1-percent gate; the largest difference is 0.2902 percent. The independent kinetic-energy integrals agree to floating-point precision. This tests the support quadrature and virial identity, not physical formation or collective stability.

Specific circular angular momenta and speeds at r0.1,0.3,1,3,10 are retained in results.json. [comparison.csv](comparison.csv) provides energy and necessary-condition diagnostics; NPZ files retain the generated profiles. Mechanical energies are in G M_b^2/a units and cannot be read as a photon-supply budget without the missing physical scale and assembly history.
