# Three-dimensional field foundation: provisional

We implemented a conservative three-dimensional potential for the ordinary-matter reference and separate additions for the comparison halo and the three previously drafted companion-deposit geometries. This is a numerical foundation, **not an orbital population fit or evidence that the theory succeeds**. No new stellar validation/test scores were evaluated.

## Components and provenance

The ordinary-matter reference combines the Sormani bar, two holed stellar disks, two gas disks, a nuclear stellar disk, a nuclear star cluster and a central black-hole approximation. The component layout and numerical parameters come from the [official Hunter24 AGAMA example](https://github.com/GalacticDynamics-Oxford/Agama/blob/f302756b8af2b763db58e278e30478517dc8eea3/py/example_mw_potential_hunter24.py). The central disk holes provide space for the bar. We do not also import the older Sormani example's disk or its central mass concentration.

The bar density routines are adapted from [Sormani and Vasiliev's AGAMA implementation](https://github.com/GalacticDynamics-Oxford/Agama/blob/f302756b8af2b763db58e278e30478517dc8eea3/py/example_mw_bar_potential.py), with attribution and the upstream license included. Only the density routines are retained; there is no implicit halo import. The nuclear spheroid density profile was checked against the upstream spheroid implementation. These are published models, not independent mass measurements: their previous dynamical fitting and assumed halo affect the interpretation of their parameters. No assertion that this fixed normalization is the unique ordinary-matter baseline is made.

The published spherical Einasto halo is a separate comparison component. The companion cases instead add the previously specified equatorial ring, paired upper/lower rings or shell, using the same softened inverse-distance response and amplitude 1,000 (km/s)^2 kpc. Those placements and amplitudes are **synthetic diagnostics**, not a derived deposition history or a fit to the measured stars. No spiral-arm perturbation or rotating-frame orbit model has been added in this foundation.

## Calculation and validation

**Known gravitational mathematics:** the bar/nuclei potential is calculated by spherical-harmonic integration of the density. Forces are derivatives of that same potential, so they are not an independently adjusted acceleration field. The disks use galpy's DiskSCF method. Units are kpc, km/s and solar masses. The scalar potential has units (km/s)^2 and acceleration has units (km/s)^2/kpc.

There are 72 deliberately chosen spatial positions from 0.5 to 20 kpc in cylindrical radius, heights 0.1/0.5/1.1 kpc, and three azimuths. Five component combinations produce **360 field predictions**. They are spatial predictions, not 360 observed stars. `field-predictions.json` retains each vector acceleration and potential value.

An independent analytic Plummer model is reproduced to maximum relative force error **0.0081%** over the 60-point reference check. Differencing the computed bar potential agrees with its analytic gradient to about 0.0000023 (km/s)^2/kpc at the three selected gradient-check points. These verify the numerical method in their stated scope; they do not validate astrophysical parameters.

Increasing the bar angular order from 16 to 24 changes its force by up to 4.43% on the probe grid; going from 24 to 40 changes it by up to 2.29%, concentrated near the thin outer bar. The latter is up to **1.70% of the total reference force**. The median bar-force change is much smaller. Increasing disk expansion order from 24 to 40 changes the disk-only force by up to 2.42% on its separate grid. The 360 stored pilot predictions use bar order 40 and disk order 24. They are **not certified for precision orbit inference** merely because an analytic spherical check passed.

A higher-order bar experiment encountered overflow in intermediate radial powers. It terminated without producing a valid high-order bar cache. The integration was rewritten to accumulate scaled radius ratios, and an order-64 analytic Plummer test now remains finite and matches the known force to 0.059%. The full order-64 bar convergence test has not been rerun. This distinction is recorded in `stable-integration-check.json`. galpy emits warnings while probing its holed-disk helper at zero radius; the stored off-axis field evaluations are finite. Those warnings are not treated as evidence of convergence.

The integrated bar mass approaches about 1.8254e10 solar masses; nuclear stellar components total about 7.0088e8 solar masses in this reference, with the black-hole component added separately. These are numerical integrals of the chosen model, not new observational mass determinations.

## Why further refinement is not the immediate priority

The user's broader-priority instruction prompted a [theory-level review](../theory-priority-review/report.md). That review confirms an exact degeneracy between the gravity-fit amplitude and the proposed photon-loss scale. Refining this field solver cannot establish that photon conversion is the physical source of the inferred extra gravity. The current foundation is saved for later use; a derived transfer/capture/gravity connection and a timing/energy-consistent propagation rule now take priority over additional numerical resolution.

Before precision stellar scoring, finish force convergence, impose a declared bar pattern speed/orientation, verify rotating-frame orbit invariants, and specify the orbit population and selection likelihood. The previously frozen spatial roles remain intact. A field value at a star's position is not a prediction of that star's measured speed.

Run `python research_work/results/bar-field-foundation/build.py` to reproduce the pilot using the cached or rebuilt density expansions. Large coefficient caches stay in `research_work/data-cache/bar-field/`. The source implementation, component definitions, outputs and checks are tracked.
