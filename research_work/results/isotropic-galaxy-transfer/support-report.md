# Fitted interception profiles need non-isotropic support in many galaxies

The preceding well-strength fit made progress by rejecting an unhelpful transfer modification. This calculation keeps the original attenuated parameters and asks whether their deposited profiles can be supported by isotropic bound particles. It applies an analytic condition to the actual fitted scales rather than another toy-density experiment.

## Known necessary condition

For a stationary spherical bound-particle distribution depending only on relative energy epsilon=Psi-v^2/2, nonnegative f(epsilon) gives the known velocity-space integral

    rho(Psi)=4pi sqrt(2) integral_0^Psi f(epsilon) sqrt(Psi-epsilon) d epsilon,
    d rho/d Psi=2pi sqrt(2) integral_0^Psi f(epsilon)/sqrt(Psi-epsilon) d epsilon >=0.

These are the standard isotropic distribution relations; see the [spherical distribution-function derivation](https://galaxiesbook.org/chapters/I-04.-Equilibria-of-Collisionless-Stellar-Systems_6-Spherical-distribution-functions.html). With positive spherical gravitational mass, Psi decreases outward, so density cannot rise outward in this particular equilibrium. Adding deposited self-gravity does not reverse that necessary sign. This is not a criterion for arbitrary field states or anisotropic particle distributions.

## Exact threshold for this capture profile

Let x=r/a and K=k0 a. Expand the previously derived optical depth near the center:

    tau(x,mu)=K[pi/4+mu x-(3pi/8)(1-mu^2)x^2+O(x^3)].

Angular averaging removes odd terms. Including the local opacity factor gives

    rho_d(x)/rho_d(0)=1+[K^2/6+pi K/4-2]x^2+O(x^4).

Therefore this profile rises away from the center when

    K > K_critical = [-3pi/2+sqrt((3pi/2)^2+48)]/2
                   =1.833274745.

This is an analytic consequence of the project's postulated opacity and known absorption law. It is not a newly discovered universal stability law or a claimed unique formula in the literature. With the retained fit, k0=0.122551 kpc^-1 and a=4.032983 R_disk, the condition corresponds to R_disk>3.709236 kpc. This threshold is conditional on fitted parameters; their uncertainties are not propagated here.

## Evaluation at actual fitted scales

| Sample | Profiles | Isotropic support excluded by the central condition |
|---|---:|---:|
| SPARC | 149 | 36 |
| Milky Way at stipulated scale | 1 | 0 |
| SLACS morphology-transfer sample | 6 | 5 |

The SLACS result inherits the explicitly postulated equivalent-exponential size mapping. These are model profiles evaluated at observed-derived scales, not observations that real galaxies lack stable particle populations. Direct small-radius angular integration verifies the analytic expansion for every listed case. support-results.json preserves all objects and coefficients.

The other 113 SPARC cases and the Milky Way merely pass this central necessary condition. They are not thereby proven to have positive distribution functions, full equilibrium or collective stability. No full Eddington inversion is claimed.

## Consequences for the proposed mechanism

Interception can empty the central deposited component while loading it farther out. Isotropic orbits tend to carry particles through the interior, so they cannot maintain the resulting outward-rising central density in a stationary spherical bound population. A positive Jeans pressure integral would not fix this distribution-function obstruction.

The mechanism therefore needs a specified alternative in the affected systems: tangentially biased orbital support, redistribution after capture, or physically defined field/bound-state stresses. A circular-orbit snapshot can avoid the isotropic restriction, as the earlier toy study discussed, but forming that distribution requires angular momentum and energy accounting, and its stability is separate. None is silently supplied by the current deposited-density formula.

This does not change the reported conditional rotation/lensing predictions, but it prevents describing those profiles as dynamically realized until their support is demonstrated. No cosmic age, photon supply or new gravity response was chosen in this check. All six goals remain open.

Reproduce: `python research_work/results/isotropic-galaxy-transfer/support.py`.
