# CC-1: a shared local cone instead of a photon-speed patch

20 September 2026. [Protocol](protocol.md), [constitutive model and spherical source](model.py), [saved summary](evidence-v1/summary.json), [read-only verifier](verify_archive.py).

The approach has changed: the scalar and directional fields now determine the same local propagation geometry for photons and mediator waves. This removes FM-1's mismatch between a field with fixed wave speed and photons with field-dependent speed. All 600 declared local equation samples and 108 spherical-source checks pass. This is analytic and numerical evidence for a candidate to implement, not a completed nonlinear 3D simulation or an observational solution.

## Equations and shared cone

Let alpha=exp(g phi), beta=alpha b, and b=kappa eta A/sqrt(1+eta^2 |A|^2), where 0<=kappa<1. Let E=sqrt(m^2+p^2). The particle Hamiltonian is

```text
H_particle = alpha E + beta dot p
velocity = alpha p/E + beta
force = -E grad alpha - (grad beta)^T p
```

The three-dimensional field Hamiltonian is

```text
H_field = integral [alpha (Pi^2 + |grad F|^2)/2
                   - beta dot sum_a(Pi_a grad F_a)
                   + omega^2 |F|^2/2] d^3x
F = (phi,A_x,A_y,A_z)
```

The protocol gives its full canonical equations, including derivatives of alpha and beta that act as field self-sources. Omitting those terms would break Hamiltonian energy exchange. Point-source regularization remains unresolved; the tested spherical profile is a regulator for the planned solver.

For frozen coefficients, the principal frequencies are beta dot k +/- alpha |k|. A photon has velocity beta+alpha n; a massive particle lies inside that same cone because |p|/E<1. Equivalently, rays have the local optical metric

```text
ds^2 = -alpha^2 dt^2 + |dx-beta dt|^2.
```

This is an effective geometry on a nonexpanding background, not a derivation of Einstein gravity. The Hamiltonian wave modes share this principal metric. The chosen mass potential is not the generally covariant scalar-field mass term, and A's three components are not a covariant vector theory merely because they transform under spatial rotations.

Both local speeds can differ from one in coordinates; their common cone is the relevant comparison. A light clock measures the local geometry. A coordinate-speed bound is not restored by silently redefining only the photon speed. The construction does not yet establish global causality, Lorentz covariance, or consistency with measured clocks and rulers.

## Why a simple fixed-speed patch is restrictive

Suppose a smooth massless Hamiltonian is H=|p|h(n,F), with ordinary vacuum h(n,0)=1. Euler homogeneity gives n dot velocity=h. If |velocity|<=1 for every small field displacement of either sign, then h<=1 and h is at its maximum at F=0. Consequently every first field derivative of h at vacuum must vanish. A nonzero linear light-field coupling cannot satisfy all those assumptions simultaneously.

This short argument applies only to a smooth, energy-homogeneous ray Hamiltonian saturating a fixed vacuum speed bound, with two-sided field variations. It does not exclude restricted fields, nonsmooth laws, dispersive particles, or changing the mediator cone. CC-1 takes the last route. This is an elementary bound argument, not a claimed new general no-go theorem.

## Positivity and the stability actually checked

Since |beta|/alpha=|b|<kappa, particle energy is at least alpha(E-kappa|p|)>0 for any nonzero-energy particle. For each field component the kinetic/gradient quadratic block has smallest eigenvalue alpha-|beta|>0. The remaining mass potential is nonnegative.

The particle momentum Hessian is alpha(I/E-pp^T/E^3). It is positive definite for m>0 and semidefinite for m=0, with the expected photon radial zero mode. This is a stronger property than checking energy alone.

The source-free zero-field vacuum has four stable linear wave modes with frequency squared |k|^2+omega^2. Freezing nonzero alpha,beta gives (frequency-beta dot k)^2=alpha^2 |k|^2+alpha omega^2. Nonzero homogeneous fields are generally not stationary solutions of the nonlinear equations, so that frozen calculation does not certify their full stability. Inhomogeneous field terms, particle-driven instabilities, collapse and long-lived circulation require actual evolution tests. Positive local kinetic blocks do not establish global nonlinear stability.

At weak fields the source is (g E,kappa eta p). Static linear elimination retains scalar attraction and a current-current term. It does not demonstrate that nonlinear self-coupling preserves FM-1's measured bending enhancement. That measurement must be repeated with CC-1 evolution.

## Results

The 600 samples span both coupling signs, zero couplings, four drift limits, massless and massive particles, random momenta and 3D field directions. They are local mathematical stress samples, not 600 independently fitted theories.

| Test | Largest error or result |
|---|---:|
| Particle velocity vs Hamiltonian finite difference | 4.80e-10 scaled |
| Matter field-source derivative | 5.65e-11 scaled |
| Momentum Hessian vs finite difference | 1.37e-9 scaled |
| Photon local-cone residual | 2.22e-16 |
| Massive local-cone excess | 2.22e-16, roundoff |
| Negative momentum-Hessian eigenvalue residual | 9.99e-16, roundoff |
| Field kinetic block and particle energy | positive in all 600 |
| Frozen Fourier frequency mismatch | 2.66e-15 |
| Fourier eigenvalue real-part residual | 1.11e-15 |
| Rotation covariance residual | 1.78e-15 |
| Photon momentum-doubling velocity change | zero |
| Wrong field-drift sign negative control | detected in all 361 nonzero-drift samples |

The spherical profile is proportional to max(1-r^2/a^2,0)^3. Its exact volume normalization is 64 pi a^3/315 and its mean squared radius is 3a^2/11. Independent radial quadrature reproduced both to the reported numerical precision. The grid profile includes the derivative of its normalization, and all 108 weight, force, deposition and rotation checks pass.

| Physical profile radius | Rotation spread, N=32 | N=48 | N=64 |
|---|---:|---:|---:|
| 0.6 | 0.2406% | 0.01364% | 0.002451% |
| 0.9 | 0.02394% | 0.004760% | 0.0008082% |
| 1.2 | 0.006120% | 0.0007474% | 0.0002967% |

All three finest-grid spreads satisfy the declared 2% ceiling. These are rotations of the source plus sampled field in a 3D interpolation fixture, not a rotated self-gravitating galaxy. At N=64 the reference sampled values are 0.781121, 0.741874 and 0.691454 for the three radii: smoothing-scale sensitivity remains substantial even when rotation error is small. The archive retains each value rather than selecting the most favorable radius. A spherical shape removes the prior separable profile's built-in anisotropy, while lattice errors remain and decrease here with refinement. Finite spatial averaging still does not establish point locality.

## Attribution and observation limits

| Ingredient | Attribution / project status |
|---|---|
| Lapse/drift optical geometry and directional null propagation | Established geometry; [Gibbons, Herdeiro, Warnick and Werner](https://arxiv.org/abs/0811.2877) is a direct precedent. |
| Canonical Hamilton equations, convexity and characteristic analysis | Established mathematical tools, not project inventions. |
| Specific alpha,beta constitutive functions and four-field Hamiltonian | Project candidate; no historical priority claim, Einstein constraints, or observational validation inherited. |
| Spherical compact source averaging | A numerical regularization choice, not a discovered physical emission law. |
| Local common cone and finite-grid rotation results | Demonstrated only within the stated audit. |

The original [GW170817/GRB analysis](https://dcc-lho.ligo.org/LIGO-P1700308/public) constrains light/gravity propagation differences under its emission and source assumptions. A shared local principal cone is relevant but does not by itself reproduce that observation, source radiation, or polarization. We have not inserted expansion or a cosmological-distance fit to claim a pass. No dark matter was introduced and none of the galaxy or cluster scores was improved by this audit.

## Reproduction, first-run error, and next action

Protocol commit b2592f5; numerical source commit c817ce9. The manifest and hashes pin the complete first scientific output. The process exited nonzero only after saving every artifact: its final stdout JSON printer did not handle a NumPy integer. The saved summary uses the separate correct serializer. The one-line printing correction is explicitly checked by verify_archive.py; evidence was neither overwritten nor rerun to hide the error.

Run `python -B research_work/experiments/common_cone/verify_archive.py` for the read-only audit; all 1,823 checks pass. The campaign refuses to overwrite evidence-v1. Its first-run status must not be reported as a clean process exit.

Next implement the full nonlinear 3D canonical field/matter equations with spherical deposition, measure the discrete Hamiltonian gradient, and test source recoil and outward energy flux before a long-lived swirl experiment. Use domain-size and discretization controls. Full coupled derivatives, isolated boundaries, light bundles, clock/polarization physics, source longevity and the joint observational fit all remain incomplete in the goal ledger.
