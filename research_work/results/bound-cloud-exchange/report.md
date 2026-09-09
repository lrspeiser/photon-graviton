# Bound-cloud energy exchange: excitation without population growth

The density-coupled candidate can transfer photon energy into a preexisting bound cloud, but its exact phase symmetry conserves companion occupation. It therefore does not yet explain production of the store. The [derivation](derivation.md) specifies the interaction, units, symmetry, seed energy and scope; [protocol.json](protocol.json) freezes the cases and numerical gates.

Five cases, each at two integration resolutions, evolve a Gaussian cloud and symmetric photon pairs together in one Hamiltonian. No external energy is supplied. All photons have crossed beyond eight final cloud radii by the end. The numbers below are dimensionless, in the declared binding-energy units; the initial cloud binding energy is -0.5 and its separately accounted seed rest energy is 6666.67.

| Case | Initial photon energy | Final cloud energy gain | Minimum / maximum radius |
| --- | ---: | ---: | ---: |
| No light | 0 | 0 | 1 / 1 |
| One pair | 0.04 | 2.51326e-7 | 0.999292 / 1.000709 |
| 16 pairs | 0.64 | 6.15055e-5 | 0.989031 / 1.011215 |
| 64 faint pairs | 0.256 | 4.02898e-6 | 0.997168 / 1.002847 |
| 64 pairs | 2.56 | 4.11897e-4 | 0.972043 / 1.029550 |

The cloud contracts and expands while remaining bound within this approximation and finite history. Its net excitation is paid for by photon energy loss. This energy can still gravitate; unchanged occupation does not mean unchanged total gravitational energy. However, no new companions appear, no incoming companions are captured, and the seed rest energy has not been explained by this calculation.

Across all ten runs, maximum absolute total-energy drift is 5.29e-9 (gate 1e-8), independently integrated work disagrees with cloud energy change by at most 3.61e-14 (gate 1e-8), and final radius, momentum and transferred energy change by at most 1.06e-9 under refinement (gate 1e-7). Occupation conservation follows analytically from symmetry, not from a numerical test that allowed occupation to vary. Machine-readable histories and source hashes are in [the saved results](bound-cloud-exchange-results.json).

This is a coupled energy-exchange test of a single Gaussian mode, not a stable halo solution, ordinary-graviton theory, formation calculation, late-time attractor or observational validation. The next step is a specified production interaction with occupation/charge and momentum accounting, followed by its actual emitted spectrum and access to the bound state. Source/detector clocks, void coupling, sustained cumulative redshift, energy supply and joint rotation/lensing remain open within the full 20-task/32-area goal.
