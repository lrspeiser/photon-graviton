# Orbit checks for the full-bar empirical gravity response

**The new conservative potential has now been used to integrate the same six synthetic stellar launch states under ordinary matter, an axisymmetric extra potential, and the full-bar extra potential.** This is a forward-model check, not a fitted distribution of observed stars. No stellar velocities, validation outcomes or test outcomes were loaded.

All integration, rotating-frame invariant, frame-comparison and sampled-path field-refinement gates, including any recorded individual-path repair, pass: **True**. Each check has the limited domain described below. Passing them does not establish a successful gravity theory; a failed check must be resolved before using this orbit operator for inference.

## What is held fixed

The launch states are copied exactly from the earlier rotating-bar diagnostic: radii 1, 3 and 8 kpc, initial heights 0.1 and 0.8 kpc, and azimuth 30 degrees. Initial radial and vertical velocities are 20 km/s; the original tangential launch velocities were 0.8 times an ordinary-matter circular-speed diagnostic. These identical physical states are used in all three new cases. They are not equilibrium distributions or six measured stars.

Each integration lasts 0.25 kpc/(km/s), about 244 million years, sampled at 501 times. The prescribed bar rotates at 37.5 km/s/kpc, with the earlier prograde sign convention. No pattern speed, initial condition, companion coefficient or population weight is fitted. The separate inertial-frame check covers 0.05 kpc/(km/s), about 49 million years, with 101 times.

The ordinary field matches the field used to generate the [full-bar completion](../full-bar-completion/report.md): disk harmonic expansion, full bar through l=64, axisymmetric nuclear components and central mass. The independent Legendre-versus-normalized-harmonic check agrees to 2.45e-14 in relative force at 60 synthetic points. This check prevents changing the ordinary baseline inadvertently when attaching the extra potential. It does not establish its observational mass accuracy.

The axisymmetric-extra case uses the circular reference correction. The full case adds the nonaxisymmetric correction, including the nonlinear change in its azimuthal mean. Both still include the full ordinary bar. These are not the earlier synthetic deposit rings, and no dark halo is added.

## Known mechanics and limits of the energy check

Known rotating-frame Hamiltonian equations, not a new companion formula:

`H_J=|p|^2/2+Phi(x)-Omega*L_z`,

`xdot=p-Omega cross x`, `pdot=-grad(Phi)-Omega cross p`.

Here p is the inertial velocity expressed in rotating axes. The constant Jacobi quantity H_J monitors test-particle motion in a rigidly rotating prescribed potential. A star's inertial energy and angular momentum can change separately through the bar torque. This is not the energy ledger of the photon/companion/universe system; source backreaction and bar evolution are not integrated.

## Integration and invariant results

DOP853 tolerances follow the earlier declared sequence 2e-9, 2e-11, then 2e-13 when needed. Each adjacent comparison must differ by less than 0.0001 kpc in position and 0.01 km/s in velocity. The maximum Jacobi drift divided by (220 km/s)^2 must remain below 1e-5. Failed initial comparisons remain in the recorded attempt list.

| Model | Retained relative tolerance | Position refinement difference (kpc) | Velocity refinement difference (km/s) | Maximum scaled Jacobi drift |
|---|---:|---:|---:|---:|
| ordinary | 2e-13 | 4.7766e-06 | 0.0027933 | 4.1198e-09 |
| axisymmetric_extra | 2e-13 | 0.00012022 | 0.052378 | 6.1226e-09 |
| full | 2e-11 | 1.2968e-05 | 0.0063089 | 1.0256e-07 |

### Recorded individual-path refinement

The table above preserves the original batch integrations, including failures. The circular-reference failure was localized by recomputing the preceding tolerance and checking each trajectory separately. The following extension integrates only failing launch states independently at relative tolerances 2e-13 and 2.3e-14. Each particle obeys the same fixed potential without interacting with the other probes, so independent integration does not change the physical model. The position and velocity thresholds remain unchanged.

- axisymmetric_extra, probe 0: position difference 3.2291e-06 kpc; velocity difference 0.001454 km/s; passes: **True**.

The separately saved refined trajectory replaces a numerical integration only when these checks pass; it does not replace an observed star, fit a velocity or alter a gravity coefficient. The original paths and failing comparisons are retained. The illustrative occupation table below uses the original ordinary and full cases, which did not require this individual repair.


The domain checks apply at every field evaluation, including intermediate integration stages, and raise errors rather than extrapolating. Full-case minimum/maximum evaluated spherical radii are 0.0484915 / 8.10295 kpc. This is coverage of these integrations, not proof of coverage for an eventual stellar orbit library.

The independent inertial-versus-rotating integration differs by at most 6.7913e-09 kpc and 5.3893e-06 km/s over its shorter interval. This checks the rotation signs and momentum convention as well as the numerical evolution.

## Field resolution along the actual computed paths

At all 3,006 saved positions on the full-field trajectories, changing the three-dimensional correction from the fine to finer grid changes the extra force by at most **0.36423%**, with median **0.004181%**. The predetermined target is 1% of the extra force. This is not normalized to a larger ordinary field that could hide an unresolved addition.

This check is distinct from time-integration refinement and from the earlier 60 fixed spatial probes. It does not sample every integrator stage or every possible orbit. The reference-order and outer-boundary checks were performed in the preceding field study, not repeated at all these positions. Disk and inner-boundary convergence, mass-model uncertainty and population selection remain separate requirements.

## Illustrative orbital occupations

The following fractions are equal-time samples on individual trajectories. They are not survey-weighted stellar counts or predicted equilibrium populations. A changed fraction says how this particular launched path changes, not that a corresponding fraction of real bulge stars should move there.

| Initial radius (kpc) | Initial height (kpc) | Full-field radial range (kpc) | Full-field maximum height (kpc) | Ordinary time fraction within 0.2 kpc of plane | Full-field time fraction within 0.2 kpc |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.1 | 0.035–1.208 | 0.496 | 0.525 | 0.673 |
| 1 | 0.8 | 0.029–1.594 | 0.874 | 0.397 | 0.301 |
| 3 | 0.1 | 0.567–3.133 | 0.243 | 0.810 | 0.884 |
| 3 | 0.8 | 0.656–3.080 | 1.062 | 0.162 | 0.120 |
| 8 | 0.1 | 2.022–8.053 | 0.269 | 0.341 | 0.711 |
| 8 | 0.8 | 2.021–8.056 | 0.938 | 0.108 | 0.154 |

## Consequence for the observational test

This calculation evaluates the specific orbit operator needed for the bulge/plane comparison. It cannot replace a shared stellar population model. The next inference work must use a sufficiently broad orbit library or distribution function, propagate position/distance and velocity errors, account for survey selection and chemical populations, and keep questionable cross-identifications from becoming false high-speed evidence.

No per-region gravity multiplier or orbit weight was adjusted here. The empirical additional field remains a known mathematical gravity construction with fitted galaxy-response coefficients; its photon origin, creation/capture relation, energy retention, observed timing and lensing completion remain unproved. The full scientific goal remains active and the total source-energy budget remains deferred.

## Reproduction

Run `run.py`, `support_check.py`, `refine_individual.py axisymmetric_extra` for the recorded failing-path extension, and `report.py`. `results.json` stores the original attempts; the separately named individual-refinement JSON records the repair. Large trajectories remain in the ignored data cache. The historical initial states, old stellar splits and previously frozen observational predictions are unchanged.
