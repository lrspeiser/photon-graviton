# Ordinary bar torques can deflect initially radial material

The existing ordinary bar component generates angular momentum for initially radial probes and prevents central entry on their tested first approach. Thus zero angular momentum throughout the motion was a consequence of the previous spherical force model, not a universal consequence of cold capture. This does not establish a stable deposited population or resolve all central behavior.

## Field and source provenance

This controlled experiment uses only the cached three-dimensional ordinary bar density. It does not include the rest of the Milky Way: no disk, nuclei, black hole, deposited self-gravity or dark halo contributes to these trajectories. The spherical control is the l=0 part of the same bar expansion, preserving its radial monopole while removing nonspherical forces. The experiment is therefore about the torque from this component, not a new complete Galactic model.

The density implementation comes from the pinned [Sormani–Vasiliev AGAMA bar example](https://raw.githubusercontent.com/GalacticDynamics-Oxford/Agama/f302756b8af2b763db58e278e30478517dc8eea3/py/example_mw_bar_potential.py). Its bar approximates a previously fitted Galactic model; the upstream full potential also contains a separately defined halo. Here only the bar-density cache is loaded. Excluding the halo force does not make the bar's fitted normalization an independent measurement, nor remove its historical model dependence. This is a conditional reference mass distribution, not proof of the companion hypothesis from observations.

The rotating case adopts the project's earlier37.5km/s/kpc pattern speed. The stationary bar and spherical control use zero pattern speed. These are specified physical alternatives; rotation is not assumed to leave the particle's inertial energy constant.

## Probe design

Eight probes launch at radii1 and3kpc, azimuths pi/6 and pi/3, and z/r=.2 and.6, each with inward radial speed50km/s and zero initial angular momentum. These deliberately selected directions are not observed stars or an isotropic source sample. The speed is a probe choice, not a recalculated three-dimensional photon/companion capture prediction. Initial Cartesian states are identical between fields.

Integration stops at the first radial turning point, entry inside0.1kpc, or time0.05kpc/(km/s). The inner-radius event is a reporting boundary, not a sticky center or a physical absorption prescription. A spherical-control trajectory stopped on entry has not had its true pericenter measured. Likewise, a bar trajectory that turns outside that radius on its first approach may enter later; no claim of permanent avoidance follows.

## Known equations, conditional application

The equations are ordinary Newtonian test-particle mechanics in a prescribed bar potential:

`dx/dt=v`, `dv/dt=-grad Phi(x,t)`, `dL/dt=x cross (-grad Phi)`.

Nonspherical gravity need not be parallel to position, so it can change angular momentum even when the injected velocity was purely radial. No additional companion force or initial tangential kick is inserted in this test.

For a stationary bar, energy is conserved while angular momentum changes. For a rigidly rotating potential `Phi(R,z,phi-Omega t)`, known Hamiltonian mechanics gives:

`dE/dt=Omega dLz/dt`, hence `E-Omega Lz=constant`.

The computed rotating trajectories obey this Jacobi invariant. The bar is an angular-momentum and energy recipient/donor, not an unexplained source of free energy. Its backreaction is neglected in this test-particle approximation; a finite captured population would require accounting for that exchange in the ordinary matter. These equations are established, not unique project formulas. Their use with the chosen launch probes is a conditional mechanism test.

## Interpretation and limits

Ordinary nonspherical structure provides an existing physical route to the angular momentum missing from purely spherical radial-flow models. The first-pass result therefore motivates a three-dimensional capture-and-motion calculation before inventing a physical core solely to repair a spherical numerical singularity.

It does not establish circularization, an isotropic equilibrium, a permanently core-avoiding orbit family, or the required galaxy/cluster density. Only one component of a conditional ordinary-matter model was used. Full sky directions, capture weighting in the nonspherical radiation field, long trajectories, other Galactic components and changing deposited gravity remain necessary. In particular, the earlier high-rate causal/relativistic breakdown is not repaired by this test.

## Numerical verification

The [protocol](protocol.md) precedes the integrations. DOP853 compares relative tolerances2e-9 and2e-11 with a common maximum step. Accepted event radii agree within1e-4kpc, endpoint velocities within.01km/s, and invariant drift divided by220^2 remains below1e-5. Full-bar cases additionally compare harmonic orders40 and64; the event-radius gate is.01kpc with unchanged event classification. All comparisons are retained, including the spherical zero-torque control.

Run `run.py spherical`, `run.py static`, and `run.py rotating`, then `export.py`. Source-cache hashes are recorded. No survey or withheld sample is read. The numerical checks establish first-pass integration accuracy within this representation, not precision knowledge of the actual Galactic bar.

## Next requirement

Extend the source and moving-deposit calculation to nonspherical ordinary gravity, retaining momentum and energy exchange and testing long-term behavior. The full nine-stage goal remains open: these trajectories do not derive the photon/time interaction, calibrate redshift, fund deposition, or pass joint astronomical observations.

## Final checkpoint

| Field | Probes entering0.1kpc before first turn | First-turn radius range (kpc) |
|---|---:|---:|
| spherical | 8/8 | Not measured: stopped on entry |
| static | 0/8 | 0.113462–0.423135 |
| rotating | 0/8 | 0.123885–0.444172 |

All integration and harmonic-order event checks pass. The largest L40/L64 first-turn radius change is 0.000137322kpc; the largest retained invariant drift divided by220^2 is 3.116e-08. These gates concern the numerical bar-component experiment, not its astrophysical normalization.

Representative probe: launch radius3kpc, azimuth pi/6, z/r=.2:

| Field | Event | Radius (kpc) | Final Lz (kpc km/s) | Energy change ((km/s)^2) | Omega delta Lz ((km/s)^2) |
|---|---|---:|---:|---:|---:|
| spherical | core_entry | 0.100000 | 0.000000 | 0.000277 | 0.000000 |
| static | first_turn | 0.423135 | -111.975648 | 0.000073 | -0.000000 |
| rotating | first_turn | 0.322947 | -54.149296 | -2030.598721 | -2030.598590 |

The eight directions are deliberately chosen probes, not an isotropic sample or a measured capture fraction. The result establishes first-pass deflection and angular-momentum exchange in this conditional bar field. Later passages, symmetry-axis trajectories and a continuously supplied population remain untested.

[comparison.csv](comparison.csv) retains all24 preferred-order results; the JSON files preserve both field orders, tolerances, events and source-cache hashes.
