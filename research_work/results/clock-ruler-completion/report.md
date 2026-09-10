# Do clocks and rulers preserve the wave model's observable redshift?

**In the homogeneous universal metric family tested here, a nonexpanding spatial geometry gives no measured redshift. Retaining the full reference-frequency redshift instead makes physical separations grow.** The latter does not satisfy the requested nonexpanding universe merely because it can be written with changing coordinate rulers.

This is a conditional classification of matter completions for the latest wave candidate. It is not a theorem against all nonexpanding photon-conversion theories. No observation, galaxy parameter or holdout is fitted.

## Why this test differs from the earlier clock checks

The existing atomic-line analysis considered an optical medium with epsilon=mu=n but fixed matter limiting speed and charge/mass parameters. It found blueshift and unequal responses among atomic transitions. The universal-clock analysis instead treated a common time multiplier with fixed spatial rulers and found cancellation of a homogeneous redshift.

Here clocks and rulers are treated together, with a common matter/photon limiting speed. This tests whether allowing atomic lengths to change can restore observed redshift and local measured c without introducing spatial expansion. It is a different optional matter completion, not a correction that invalidates either earlier result.

The one-dimensional transverse-wave Hamiltonian does not itself specify electrostatic forces, particle masses, atomic standards or a universal metric. The assumptions below are additional. Their consistency at the level of clock/ruler kinematics does not supply a complete interacting action.

## Known metric and electromagnetic identities

Use reference units c0=hbar=m0=1 and a spatially homogeneous isotropic metric

`ds^2 = -N(t)^2 dt^2 + A(t)^2 dx^2`.

N is the clock lapse, and A converts coordinate separation into physical separation on a constant-time slice. These are standard relativistic concepts, not project inventions; [Gourgoulhon's 3+1 lecture notes](https://arxiv.org/abs/gr-qc/0703035) provide the established geometric framework.

Known Maxwell-action algebra for this metric gives electromagnetic coefficients `epsilon=mu=A/N`. Thus the wave candidate's homogeneous optical index n fixes

`A/N=n`, while `dx/dt=N/A=1/n`.

It does not determine A and N individually. A standard relativistic test-particle Hamiltonian in this metric is

`H_m = N*sqrt[1+p^2/A^2]`.

Its low-momentum kinetic term is `p^2/(2*m_coordinate)`, where `m_coordinate=A^2/N`. This is a coordinate kinetic coefficient, not a claim that the locally measured rest mass changes.

Parameterize the remaining freedom by `m_coordinate=n^r`. Solving the two relations gives

`N=n^(r-2)`, `A=n^(r-1)`.

These are conditional algebraic deductions using known metric mechanics. They are not newly discovered redshift equations or a physical derivation of exponent r.

## Atomic standards and locally measured speed

Identical local atomic clocks tick at coordinate frequency proportional to N, while identical physical rods have coordinate length proportional to 1/A. Thus

`clock_frequency_reference proportional to n^(r-2)`,

`rod_length_coordinate proportional to n^(1-r)`.

The same leading scalings follow from a quasistatic Coulomb comparison with strength proportional to 1/n, common limiting speed 1/n, and coordinate inertial coefficient n^r: transition frequencies scale as inertial coefficient divided by n squared, while atomic lengths scale as n divided by that coefficient. Under this particular common-speed assumption the dimensionless Coulomb coupling remains constant. This does not substitute for a full atomic/QED completion of the directional reservoir model.

All r in this family give the same locally measured light speed:

`c_measured = (1/n) / [(1/A)*N] = 1`.

Preserving this local result alone therefore does not determine whether an observer measures redshift, or whether physical separations change.

## What the observer measures

For the wave candidate n(t)=1+a*t, travel over fixed coordinate distance d gives `S_ref=n_o/n_e=exp(a*d)`. A line emitted at the source's atomic transition frequency arrives at reference frequency reduced by S_ref. Comparing it with the identical transition at the detector yields

`1+z_measured = S_ref*(N_o/N_e) = (n_o/n_e)^(r-1) = A_o/A_e`.

The received-to-emitted proper event-duration ratio is the same factor. This is the familiar metric redshift relation obtained within the assumed family, not a distinct explanation derived from photon capture.

The code independently integrates ray travel times and source/receiver proper-time intervals for 54 combinations of r, distance and emission time. The finite event intervals obey the same factor exactly for this power-law family; an infinitesimal approximation is not used to assert the finite-interval result. Hamiltonian forms and local speed are checked separately.

For an illustrative optical-index ratio n_o/n_e=1.4:

| Completion | Detector/source clock rate | Coordinate rod-length ratio | Physical separation ratio | Measured 1+z |
|---|---:|---:|---:|---:|
| r=0 | 0.510204 | 1.4 | 0.714286 | 0.714286 |
| r=1 | 0.714286 | 1 | 1 | 1 |
| r=2 | 1 | 0.714286 | 1.4 | 1.4 |

Every row gives constant locally measured c. The first is a contracting metric with blueshift. The second changes the lapse but keeps physical distances fixed, canceling the homogeneous redshift. The third preserves the full frequency shift with fixed clock rates, but fixed-coordinate laboratories become physically farther apart.

## Consequence for the nonexpanding requirement

In this universal homogeneous family, choosing shrinking coordinate rods is not a way to retain the full redshift while keeping measured separations fixed. The physical separation is A times the coordinate separation. More identical rods fit between fixed-coordinate laboratories as A grows; that is the spatial change responsible for the measured redshift in this interpretation.

This does not assume a Big Bang, solve the Einstein equations or import dark matter. It identifies the kinematics of this particular completion. It is not adopted as the explanation for the fictional nonexpanding universe.

The energy transfer and local current in the prior directional-wave calculation remain valid within that postulated classical Hamiltonian. However, interpreting its coordinate frequency shift as the required measured redshift now requires a matter completion outside the nonexpanding homogeneous universal-metric branch. The reservoir sector's stress-energy and matter work would also have to be included; the metric analogy does not automatically provide a conserved full interacting model.

## What remains open

Possibilities outside this calculation include a nonuniversal light interaction with unchanged local matter standards, or a dynamically generated nonuniform optical environment with physically matched source and detector regions. Neither is adopted merely because it could avoid this restricted cancellation. Each needs a common wave, matter, energy and momentum calculation, including local measured light speed and observable frequency/event timing.

The earlier matched-endpoint ray examples were prescribed kinematics, not a derived boundary mechanism. Reusing them without coupling the new reservoirs and matter would not close the gap. Likewise, adding a separate event-duration factor to an energy-loss law would not explain the common observations.

The next model should therefore be tested against these observable standards before more gravity fitting is interpreted as photon-origin evidence. No holdout data was opened and no claim of a strong case is made.

## Reproduction

Run `run.py`. `results.json` contains every tested exponent/distance/emission combination, the three representative completions, Hamiltonian and interval checks, scope limitations and hashes of the relevant previous work. No existing wave or matter calculation is modified.
