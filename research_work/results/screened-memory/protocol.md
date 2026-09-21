# SM-1: screened square-root total potential and reciprocal memory

Date: 2026-09-21. Baseline inspected: `28f154242aa6e6ae010a97af8124d01cc63fed3b`.
Status: exploratory candidate, not an adopted or complete theory. This declaration precedes the SM-1 evolution and observational scoring. Algebraic construction, the inverse implementation, and point-Sun illustrative numbers were checked while designing it; those are not claimed as preregistered discoveries.

## Purpose and non-duplication

PM-1 already fitted a local square-root acceleration law. PM-2A already constructed both total-potential and independently sourced auxiliary-field completions. Do not repeat those as new discoveries. This experiment adds a decaying high-acceleration excess to the total-potential completion and an explicit reciprocal time-dependent field/matter energy model. It does not replace the microscopic Phase Junction program or erase its failures.

The supplied alternative-gravity review proposed increasing the power of the auxiliary-field constitutive law at high acceleration. That reduces the fractional correction but leaves the absolute auxiliary force increasing. For an isolated spherical source with fixed boundary conditions and a locally elliptic independently sourced scalar, `gN = mu_chi(gchi/a0)*gchi` has positive derivative with respect to `gchi`, so its inverse cannot have a decreasing absolute excess. This is a restriction on that implementation, not on all companion mechanisms.

## Frozen construction

Let `x = gN/a0`, `z = g/a0`, and freeze

```
a0 = 6.54e-11 m/s^2
h(x) = x + sqrt(x)/(1+x^2)
z = h(x)
mu(z) = h^{-1}(z)/z
```

The numerical a0 is the rounded existing PM-1 fitted scale, not a new fit or microscopic prediction. There is no fitted crossover, exponent, or per-object halo. The constitutive choice itself is phenomenological and is not uniquely derived from the junction algebra.

The field equation is for the TOTAL potential:

```
div[mu(|grad U|/a0) grad U] = 4 pi G rho_b.
```

It is not the review's separately sourced additive chi equation. In spherical symmetry it retains `g = gN + sqrt(a0*gN)` asymptotically at low acceleration and has `g-gN ~ a0^(5/2)/gN^(3/2)` at high acceleration. Define `W(z) = integral_0^z mu(v)*v dv`; use its consistent boundary inverse, relaxation weight, and field functional in nonspherical solves.

For x>0,

```
h'(x) = 1 + (1-3x^2)/(2 sqrt(x) (1+x^2)^2)
h'(x) > b = 1 - (75/128)*(3/5)^(3/4) > 0.6005.
```

Thus mu>0 and `mu+z*mu' = 1/h'(x)>0`. This is convexity/ellipticity of the static scalar gradient functional; it does not establish nonlinear stability of the complete universe. The zero-gradient principal symbol is degenerate, not uniformly elliptic/hyperbolic.

## Reciprocal dynamical extension

One declared preferred-frame effective extension is

```
Utt/cstar^2 + Gamma*Ut/cstar^2 - div[mu(|grad U|/a0) grad U] = -4 pi G rho_b.
m*qddot = -m*grad U.
```

For Gamma=0, field kinetic energy, the gradient functional, particle kinetic energy and the same matter-field interaction form a conservative Hamiltonian. For Gamma>0 record local heat production `Gamma*Ut^2/(4 pi G cstar^2)`; do not call that a derived microscopic bath or a closed bath momentum model. Memory here is the field's retained state and velocity, not a separately assigned halo or a demonstrated long-lived density deposit.

The one-dimensional numerical control uses 4*pi*G=a0=1, cstar=0.7, two mass-0.2 compact extended particles initially at +/-1 and at rest, an initially empty field, domain [-4,4], Dirichlet U=0, source width 0.4, and 161 nodes. Duration 2; dt=0.004, 0.002, 0.001. Compare Gamma=0 and Gamma=1. No initial hidden field energy or source excitation is inserted. The upper local principal speed bound is cstar/sqrt(b), below the reference speed 1 for this numerical choice; a physical common cone is not derived.

## Fixed numerical checks

- Constitutive inverse relative error below 1e-11 over x=1e-14..1e14.
- Gradient-energy derivative relative error below 1e-6 on z=1e-6..1e2, using independent finite differences.
- Full discrete particle/field Hamiltonian derivative residual below 1e-6 on a nonuniform fixture; include source-normalization derivatives.
- Evolving energy error is normalized by the maximum sum of absolute channel energies, not by the initially zero total energy. Require below 1e-4 at the finest dt and decreasing error under timestep refinement. Heat must be nondecreasing. Report trajectories without promoting 1D toy attraction to a galaxy result.
- Preserve numerical failures and diagnostics. Do not adjust equations, couplings, thresholds or selections in response to observational scores.

## Observational reuse and nonspherical extension

Reuse the repository's 149 SPARC names and 89/29/31 split from `companion-extensions/mond-inventory-results.json`, raw `temporal_candidate_audit/data/Rotmod_LTG.zip`, disk mass-to-light 0.5 and bulge 0.7, and signed gas squared speeds. No observed velocity enters the prediction. Report all rows, per-galaxy RMSE and equal-galaxy aggregate statistics. Compare baryons, the unscreened PM-1 law at the SAME rounded a0, and the repository's simple-MOND control at a0=8.563e-11. This is exposed-data algebraic screening, NOT a full disk PDE prediction or a new blind test.

Where executable, reuse PM-2A `fields.Solver` without editing it. The declared synthetic source is an exponential disk with radial scale 1, vertical exponential scale 0.15 and compactness GM/(a0*Rd^2)=0.1 or 10. Correct the unresolved inner mass with the spherical integral, not the inherited thin-cylinder proxy. Compare grids (128,24) and (256,48), radii 0.3,0.5,1,2,4,8,12, inner radius 0.002 and outer radius 80. An outer-radius-160 control must retain at least the same logarithmic radial resolution. Require solver convergence and below 2% force refinement change and 0.5% boundary change before calling these synthetic predictions resolved. Do not substitute algebraic radial mapping for the nonspherical solve.

## Boundaries on claims

A tiny isolated point-Sun correction is NOT a Solar-System pass. Galactic external-field quadrupoles, ephemeris fitting, preferred-frame tests, binary/scalar radiation, lensing coupling, cluster memory, and microscopic derivation remain open. In particular see Desmond, Hees & Famaey (2024), arXiv:2401.04796, on joint RAR/Cassini constraints. AQUAL-type variational gravity and MOND asymptotics have prior art (Bekenstein & Milgrom 1984); do not label them new Phase Junction discoveries. This experiment creates a specific falsifiable construction, not a proof that nature follows it.
