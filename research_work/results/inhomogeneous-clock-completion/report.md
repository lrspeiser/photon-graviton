# An explicit clock rule can preserve redshift without spatial expansion

An evolving, spatially nonuniform clock field provides an operational counterexample to the earlier homogeneous cancellation. The ray/clock submodel is now specified and tested. It does not yet derive that field from radiation or establish agreement with astronomical observations.

## The postulates and their mathematical consequences

Postulate fixed spatial rods and clocks at rest obeying d tau=dt/n(x,t). Both messenger types follow omega=c0|k|/n and dx/dt=c0/n. Then a local observer using that clock measures dx/d tau=c0. This is a local identity, not proof of a finite laboratory's behavior or an independent experimental test.

These choices are the standard lapse-metric kinematics ds^2=-c0^2 dt^2/n^2+dx^2+dy^2+dz^2. Their use here is a proposed measurement rule, not a new metric formula. Spatial distances between fixed coordinate positions do not expand. Clocks keep their ordinary local operation but accumulate different elapsed times when compared across regions. A fixed-position observer may need force to stay there; a freely falling observer is a different trajectory.

Known phase and arrival-map derivatives imply

    J = dt_o/dt_e = exp[integral partial_t n dx/c0]
    S_measured = J*n_e/n_o
    Omega_o/Omega_e = 1/S_measured.

Thus wave frequency and event timing have the same measured stretch, with endpoint clock changes included rather than ignored. Spatially homogeneous n gives cancellation; a nonuniform evolving field need not. Spatially varying but time-independent n has only endpoint frequency effects in this specialization, not arbitrary accumulated redshift between equal clock environments.

Frequency shifts from evolving gravitational environments, including effects on gravitational waves, have known precedents: [Laguna et al., Integrated Sachs-Wolfe Effect for Gravitational Radiation](https://arxiv.org/abs/0905.1908). We do not adopt that paper's cosmological background or claim its equations validate our model. This calculation cannot be advertised as discovery of a previously unknown mathematical redshift mechanism.

## Prescribed example and results

Take n=1+beta*t*sin^2(pi*x) with unit spatial period and c0=1, over t>=0. Its profile is a diagnostic postulate, not a radiation-generated void solution, protection screen or measured gravitational map. Nodes have n=1 for all time; other positions have explicitly changing clocks. beta is not fitted to astronomical data.

For node-to-node paths of integer length D, n_e=n_o=1 and S=exp(beta*D/2). This recovers an exponential relation conditionally; alpha_eff=beta/2 for those paths only. Generic endpoints require the clock factor and do not follow that simplified formula.

| beta | Path | Measured redshift |
|---|---|---:|
| 0.02 | 0 to 1 | 0.01005017 |
| 0.02 | 0 to 2 | 0.02020134 |
| 0.02 | 0 to 5 | 0.05127110 |
| 0.02 | 0.25 to 1.25 | 0.00001566 |
| 0.02 | 0.25 to 1 | 0.01922433 |
| 0.02 | Homogeneous control, length 2 | 0 |

All emissions are at coordinate time 1. Both messengers have identical arrival times under the imposed common ray rule. This does not predict a nonzero source lag. Ordinary observer clocks at the node endpoints see a redshift that does not disappear on reception. Generic endpoints show why selecting only special endpoints would overstate universality: the same unit path length can have a markedly different measured shift.

## The gravity consequence is part of the model

The prescribed metric requires holding acceleration a_hold=c0^2*grad ln(1/n) for fixed-position observers. An initially stationary freely falling body's local acceleration is opposite that. The time field therefore has force consequences even before a separate deposited-energy potential is added. It cannot be claimed to have no gravitational effect in the void merely because both waves see the same local c.

For beta=0.02, maximum holding-acceleration magnitudes along the integer paths of lengths 1, 2 and 5 are about 0.110, 0.171 and 0.352 in units c0^2 per chosen length unit. These are demonstration scales, not fits to real forces. There is no justification for discarding them or adding another gravity term to cancel them after inspecting data.

The lapse profile is prescribed. Ordinary Einstein field equations with a strictly flat spatial slicing and zero extrinsic curvature would impose additional source constraints; we have not shown compatibility with positive temporal-field energy using those equations. The fictional model can postulate another gravitational response, but must supply it explicitly. Neither coordinate kinematics nor the earlier open-cell energy sum supplies it automatically.

## What is completed and what remains

Completed here: one fully explicit ray/clock measurement rule, its observable frequency/event relation, homogeneous control, generic-endpoint comparison and the associated holding-force consequence. This advances roadmap goals 1–3 at the kinematic level. It establishes that explaining the ultimate origin of time is unnecessary to derive these quantities.

Still required: evolve the spatial field from the radiation-fed branch with declared initial/boundary conditions, include energy and momentum at the same operational level, show why sources/detectors occupy their proposed environments, and compare the resulting forces and frequency/timing predictions to observations. The periodic profile does not answer any of those questions. It must not be relabeled a self-consistent physical solution or a fitted void law.

Verification: 18 cases integrate rays, phase and arrival Jacobians, compare the latter with an analytic integral and compare measured event stretching against independent proper-clock intervals. All checks pass the declared tolerances; all zero-field and uniform controls give zero measured redshift. The local-speed equality is an identity of the postulates, and shared messenger arrival is imposed, not independently discovered. Reproduce with `python research_work/results/inhomogeneous-clock-completion/run.py`. No observation was fitted or held-out data accessed. The nine-goal objective remains active and incomplete.
