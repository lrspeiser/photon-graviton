# SE-B evolving bundle checkpoint

7 of seven declared cases have complete archives.

| Case | Central bend (model radians) | Arrival offset (model time) | Area gain | Bundle derivative error | Numerical gates |
|---|---:|---:|---:|---:|---|
| emission-only | -0.011270736 | 0.0107100034 | 1.03328373 | 0.000397 | True |
| X | -0.011272043 | 0.0107104647 | 1.03328549 | 0.000397 | True |
| mixed | -0.0112774617 | 0.0107160817 | 1.03329593 | 0.000397 | True |
| Y | -0.0112771959 | 0.010718004 | 1.03329702 | 0.000397 | True |
| time | -0.0112773522 | 0.0107181099 | 1.03329646 | 0.000397 | True |
| space | -0.0103254399 | 0.0105957948 | 1.02441056 | 0.000141 | True |
| probe-radius | -0.0082780675 | 0.00993502136 | 1.02026907 | 0.000104 | True |

![Candidate trajectories](bundles.png)

The detector is x=1.5; rays launch at x=-1.5 near y=1 with parallel +x
momentum. The detector map includes the evolving background, not a frozen
lens snapshot. Arrival offsets are coordinate times relative to flat flight.
Area gain refers to the parallel-ray transport map, not yet an astronomical
magnification. No galaxy or cluster measurements appear in this chart.

The massive-body chart subtracts free flight from the same initial position
and velocity, derived from the zero-field candidate Hamiltonian. These short
test-body paths are not stable orbits, circular speeds or a rotation curve.

Independent audit: audit_bundles.py reconstructs detector crossings and maps
from every-step traces and background energy from raw endpoint states.
Probe cone extrema are recorded by the runner; the archive does not contain
the full field at every step for independent reconstruction of those extrema.

Comparative emitter effects must wait for matching completed cases. Time,
space and probe-radius comparisons are still required before interpreting
a point-ray result. Positive source enhancement alone is not lens enhancement.
Hamiltonian optics and numerical bundle derivatives are established methods;
the shared response law is a candidate assumption. All twelve goals remain active.

## Completed campaign verdict

All seven executions finished. Overall declared campaign pass: False.
The time comparison passes, but the spatial comparison fails. Thus the
small differences between emitter choices are not resolved accurately.

| Comparison | Bend relative difference | Arrival difference | Matrix difference | Passed |
|---|---:|---:|---:|---|
| time | 1.38570673e-05 | 1.05896955e-07 | 2.78627023e-07 | True |
| space | 0.0921758251 | 0.000122209229 | 0.00721826252 | False |

Changing the probe kernel radius changes the central bend by 26.5946%.
This radius dependence is a modeling/regularization sensitivity, not an
observed photon-size effect. A point-ray interpretation is not established.
