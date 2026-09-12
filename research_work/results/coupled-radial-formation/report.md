# Moving deposits and capture feedback from zero

This experiment combines depleted incoming supply, momentum-derived radial injection, moving deposits and their changing gravity. A lower-rate example remains within the slow-motion regime and passes the declared numerical mass and energy checks over the tested interval. A higher-rate example leaves the physical approximation and fails resolution checks. Even in the lower-rate case, the central potential remains sensitive to numerical regularization.

## Model and formula provenance

Start with zero deposits. The only prescribed gravitational background is ordinary Plummer matter, held fixed: G=M_b=a=1 and Phi_b=-1/sqrt(1+r^2). There is no old deposited profile held in place. This improves on the earlier test-population response, but the response of ordinary matter is still omitted.

The project capture postulate is `kappa=.1 [W/(1+W)]^6`, with W the current ordinary-plus-deposit depth. Two density-rate amplitudes Cdot=.1 and1 imply incident B=Cdot/.1. These are new explicit rate choices per dynamical time, not the old exposure-integrated C silently reinterpreted as a measured luminosity. Isotropic rays are attenuated along full chords within R30. Particles leaving that domain remain in the simulation and continue to gravitate; no escape mass is discarded.

The candidate local-merger postulate gives `beta=|F|/J`, `v_new=-c beta` and newly retained rest mass `dm=dE/c^2 sqrt(1-beta^2)`, with c/v0=1000. Known relativistic composition motivates this local energy/momentum relation. The subsequent trajectories use nonrelativistic mechanics, an approximation acceptable only while speeds stay small. The method is not a relativistic graviton or stored-field theory.

Moving spherical cohorts interact through a regularized shell Hamiltonian:

`U_D=-(1/2) sum_ij m_i m_j / sqrt(max(r_i,r_j)^2+epsilon^2)`.

Its radial force on an ordered shell includes interior mass plus half its own mass, multiplied by `r/(r^2+epsilon^2)^(3/2)`. This is a numerical regularization of spherical Newtonian shell gravity, not a new physical core law. A separate finite-difference gradient check validates the force implementation. All shells follow signed radial trajectories through the center; there is no sticking, imposed pressure or circularization.

Gravity changes capture, and newly captured rest mass changes gravity. Rays are recomputed at each injection update. Radiation transport and gravity propagation remain quasistatic; the method does not solve time-dependent wave propagation or a relativistic metric. No formula is claimed unique: mechanics, shell geometry and attenuation are known; capture and merger are explicit project hypotheses.

## Separate accounting checks

The ray calculation balances incoming energy against absorbed and transmitted energy. The local merger separately balances absorbed energy against rest energy and relativistic bulk kinetic energy. Retained rest mass equals the accumulated injected rest mass, including shells outside R30.

For moving shells the discrete mechanical Hamiltonian is `H=K+U_D+sum_i m_i Phi_b(r_i)`. Injection changes it by added nonrelativistic kinetic energy plus the exact potential-energy difference caused by inserting the new cohort. Between injection events the integrator should conserve H; the numerical residual subtracts all recorded insertion changes.

This does not prove a globally funded physical energy budget. Gravitational insertion work is tracked explicitly, not supplied by an identified radiation/recoil/field channel. In the low-rate fine epsilon0.05 case at T3 it is approximately-0.0269113 in G M_b^2/a units. The local incoming-energy convention, binding release, the fixed ordinary source's response and energy in transit still need reconciliation. Relativistic bulk kinetic energy used in the local conversion identity also differs from the subsequent nonrelativistic K at high speed. A balanced local ledger cannot validate that high-speed evolution.

## Physical and numerical outcomes

The lower rate generates a finite moving deposit population from zero over three dynamical time units. Maximum speed stays near0.0011c and the opacity-change indicator remains small. The deposited population is already centrally concentrated; no stationary halo or long-term stability is established. Total mass and selected enclosed fractions are comparatively insensitive to the tested resolution and core choices.

The potential at the exact center does not share that convergence. Reducing epsilon from0.05 to0.00625 raises the central deposited depth from approximately0.10057 to0.14447 while total mass barely changes. Thus outer enclosed-mass agreement is insufficient to claim a reliable central clock or gravity prediction. Physical angular momentum, a finite retained-state structure or another explicit central treatment is still required; choosing a convenient epsilon would conceal this missing physics. This finite sequence is evidence of sensitivity, not a proof of the exact zero-epsilon asymptote.

At the higher rate, feedback drives the solver into speeds approaching0.5c and very rapid opacity changes. Late masses and profiles cannot be used as predictions of this nonrelativistic, quasistatic model. Some mechanical numerical gates fail; one passes despite the physical breakdown. That distinction is essential: solving the assumed equations accurately does not make their approximation applicable.

After seeing that breakdown, descriptive flags were added at speed0.03c and logarithmic opacity change0.1 per diameter light crossing. These are warning markers, not sharp laws or predeclared observational thresholds. In the strong fine cases, the opacity flag first occurs near T1.05 and the injection-speed flag near T2.5. The indicator is `max |delta ln kappa|/delta t * (2R/c)`, not a full error estimate for causal transport. Low-rate runs do not cross these warning levels.

## Reproduction and scope

The [protocol](protocol.md) records the original grid and the later warning/core extension separately. Run `run.py 0.1`, `run.py 1`, and `run.py 0.1 --core` in this directory with OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. The first two commands retain coarse/fine runs at epsilon0.05 and the fine epsilon0.025 control; the third adds epsilon0.0125 and0.00625 only for the lower rate. `check_hamiltonian.py` tests the force gradient; `export.py` compares checkpoints without hiding failed cases.

The main numerical gates were mass agreement within2 percent, enclosed-fraction changes below0.02 absolute, and mechanical residual below0.5 percent of K+|U|. Changing epsilon is a change of regularized equations, distinguished from mesh/time refinement. The capture boundary remains untested in this fully coupled model. No astronomical data, photon-redshift mechanism, lensing fit, gravitational-wave timing fit or withheld evaluation is supplied here.

## Consequence for the nine goals

There is now a coupled candidate calculation for one particle interpretation, rather than permanently placing energy at arrival positions. It identifies a slow, lower-rate regime worth further study and a stronger regime requiring different dynamics. The low-rate center still lacks a physical closure; the high-rate branch needs causal and relativistic transport before extrapolation. Global energy/momentum and the original photon/time source remain open. This is progress on formation and gravity, not completion of the nine-stage objective.

## Final checkpoint

Lower-rate results at fine resolution; these are dimensionless synthetic calculations, not observed galaxy masses:

| Time | Deposited rest mass, epsilon0.05 | Fraction inside r0.1 | Fraction inside r1 | Fraction inside r3 |
|---|---:|---:|---:|---:|
| 1 | 0.0164317 | 0.4910% | 27.6512% | 77.1877% |
| 2 | 0.0339771 | 4.5453% | 35.7147% | 77.9171% |
| 3 | 0.0528323 | 4.3424% | 42.7951% | 78.7854% |

Core sensitivity at T3, all at the same fine resolution:

| Epsilon | Deposited mass | Central deposited potential depth |
|---|---:|---:|
| 0.05 | 0.05283225 | 0.100568 |
| 0.025 | 0.05283513 | 0.115750 |
| 0.0125 | 0.05283587 | 0.130406 |
| 0.00625 | 0.05283606 | 0.144465 |

Strong-source status, with late values explicitly outside the approximation:

| Source shells | Epsilon | Final maximum speed/c | Maximum opacity log change per crossing | Mechanical numerical gate |
|---|---:|---:|---:|---|
| 64 | 0.05 | 0.4904 | 7.988 | fail |
| 128 | 0.05 | 0.4829 | 11.969 | pass |
| 128 | 0.025 | 0.4901 | 12.002 | fail |

All low-rate mass/fraction resolution comparisons and all low-rate mechanical gates pass. The stronger rate fails the mass-resolution gate at T2 and T3; it also leaves the nonrelativistic and slowly changing radiation assumptions. Its enormous late mass values remain in the raw results for audit, not as physical predictions. All local radiation/rest-mass identities remain numerically balanced, which does not rescue those physical failures.

The [comparison table](comparison.csv) and comparisons.json retain every checkpoint and failed gate. A force-gradient finite-difference check of the shell Hamiltonian agrees within2.8e-9 absolute; it tests the implemented regularized force, not a microscopic gravity theory.
