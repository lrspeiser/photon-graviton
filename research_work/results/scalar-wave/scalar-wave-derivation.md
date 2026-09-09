# What the explicit scalar–photon action actually does

This step derives and numerically tests a restricted sector of candidate B from the action register. It is preliminary T03–T05 evidence, not a completed nonexpanding cosmology, a halo fit, or a final choice of law.

## Action and field equations

Suppress charges and absorbers for this propagation test. Use natural units, metric signature (-,+,+,+), and:

```
L = -Z(chi) F_mu_nu F^mu_nu/4
    - partial_mu chi partial^mu chi/2 - V(chi)
Z(chi)=1+g chi
V(chi)=m² chi²/2
```

Here g denotes the dimensionful coefficient g_gamma/Lambda in the candidate action, rescaled into the numerical unit system. We require Z>0. The test uses the decoupled-gravity flat-space limit; it does not claim that a gravitating homogeneous scalar has an exactly flat solution of the full Einstein equations.

Variation gives:

```
partial_mu[Z F^mu nu]=0
Box chi - V'(chi) - Z'(chi) F²/4=0
```

In a local rest frame, the electromagnetic Lagrangian is Z(E²-B²)/2. Thus its permittivity is epsilon=Z and its permeability is mu=1/Z in the chosen units. Their product is one:

```
wave speed² = 1/(epsilon mu) = 1
impedance = sqrt(mu/epsilon) = 1/Z
```

This distinction is why changing Z is not equivalent to imposing a refractive index n=Z. At leading geometric-optics order the principal equation retains the metric-null characteristic. Time derivatives of Z still change amplitude and can exchange energy with the field beyond that leading limit.

## Wave equation and canonical variable

Project the transverse vector potential onto a Fourier mode of wave number k and assume a spatially homogeneous scalar. The projected mode equation is:

```
A_double_dot + (Z_dot/Z) A_dot + k² A = 0
```

The scalar source is spatially averaged in this truncation; other scalar/wave modes are omitted. This is a controlled diagnostic of the projected equations, not a proof that the full nonlinear field evolution remains in this truncation.

Writing u=sqrt(Z) A removes the first derivative:

```
u_double_dot + [k² - (sqrt(Z))_double_dot/sqrt(Z)] u = 0
```

When changes in Z are slow compared with the optical oscillation and the correction term is small relative to k², the leading solution has phase frequency k and amplitude A proportional to Z^-1/2. The field amplitude varies, while the leading traveling-mode energy is approximately constant. This is not the assumed omega=ck/n(t) law. Rapid/background-resonant evolution can produce stronger mode mixing and energy exchange, but that requires its own spectrum and driver calculation; the present tests do not exhaust those regimes.

## An explicit conserved energy test

Use two real quadratures a,b to represent complex mode amplitude A=a+i b. The projected Lagrangian is:

```
L_mode = chi_dot²/2 - m² chi²/2
         + Z/2 [a_dot²+b_dot²-k²(a²+b²)]
```

It gives:

```
chi_double_dot = g/2 [a_dot²+b_dot²-k²(a²+b²)] - m² chi
a_double_dot = -(g chi_dot/Z) a_dot - k² a
b_double_dot = -(g chi_dot/Z) b_dot - k² b
```

The conserved projected energy is:

```
H = chi_dot²/2 + m² chi²/2
    + Z/2 [a_dot²+b_dot²+k²(a²+b²)]
```

In particular:

```
H_EM_dot = Z_dot/2 [k²(a²+b²)-a_dot²-b_dot²]
H_scalar_dot = -H_EM_dot
```

The sign of transfer depends on the state. These equations do not give an independently imposed one-way positive photon loss h U. The homogeneous scalar's energy is a driver/field account; it has not been identified with a bath of traveling companion particles or a deposited halo.

The exact mode equation also conserves W=Z(a b_dot-b a_dot). This Wronskian is distinct from the instantaneous electromagnetic energy. Both were checked, so a nearly conserved large scalar-energy total cannot hide arbitrary drift of the much smaller wave mode.

## Numerical results

Initial values were chi=1, chi_dot=0, a=10^-3, b=0, a_dot=0, b_dot=-k×10^-3, with k=1. Three cases were integrated for 400 time units and repeated with tighter solver tolerances. All maintained Z>0.

| g | m/k | Z range | Instantaneous mode frequency / k | Maximum wave-energy ratio |
|---:|---:|---|---|---:|
| 0 | 0.03 | 1 | 1 within numerical precision | 1 within numerical precision |
| 0.2 | 0.03 | approximately 0.8–1.2 | 0.999906–1.000075 | 1.00000481 |
| 0.2 | 0.3 | approximately 0.8–1.2 | 0.990389–1.007588 | 1.00061207 |

The maximum relative total-energy drift was below 6.3×10^-13; relative Wronskian drift below 2.5×10^-10. The largest frequency difference under tolerance refinement was below 7.1×10^-10. These are numerical consistency checks, not observational constraints.

The slower scalar example changes Z by a substantial amount while leaving the photon-mode frequency close to k. The faster example gives larger oscillatory frequency effects. Neither numerical trajectory demonstrates a monotonic cosmological redshift or irreversible conversion into captured companions. Other potentials, modes, backgrounds or interactions remain untested rather than excluded.

## What operator would change the propagation law?

A useful algebraic comparison adds a distinct electric response relative to a normalized timelike field u:

```
L_EM = -Z F²/4 + Y (u^mu F_mu nu)(u^rho F_rho^nu)/2
u_mu u^mu=-1
```

In the local rest frame of u this becomes:

```
L_EM = (Z+Y) E²/2 - Z B²/2
speed² = Z/(Z+Y)
n² = (Z+Y)/Z
```

Positive Z and Z+Y are necessary photon-sector kinetic/gradient conditions in that frame. They do not prove health of the u/driver/gravity system. Constant or slowly varying coefficients are assumed for the local characteristic calculation.

This operator shows which extra structure can distinguish electric and magnetic responses and hence alter the light cone relative to matter. It is **not an adopted extension**. A complete action would need dynamics/constraints for u, finite functions or parameters for Y and Z, a matter coupling, an energy source for changes, companion propagation and stability checks. If u is constructed from a scalar gradient, that construction and its domain must also be specified. A prescribed preferred frame alone does not close the conservation equations.

## Atomic and cavity consequences cannot be skipped

For the original Z-only candidate with fixed charges and masses, leading Coulomb binding has effective strength proportional to 1/Z. In a simple hydrogenic, nonrelativistic comparison:

```
atomic length proportional to Z
optical transition frequency proportional to Z^-2
```

If a cavity spacer tracked that atomic length, its frequency would scale as Z^-1 because the leading wave speed remains constant. The cavity/optical ratio would therefore scale as Z. Actual nuclei, solids, relativistic corrections and measured responses require further derivation; these leading scalings are not a complete laboratory prediction. They demonstrate why even an amplitude/normalization interaction affects the source/detector side of observed redshift.

A new electric-response term changes these relations again. It cannot inherit the previous clock or brightness results without recalculation.

## Next work and status

The result identifies a specific limitation of a specific candidate sector: Z F² alone does not derive the desired homogeneous refractive-index law in its leading adiabatic regime. The candidate remains open for other conversion/background behavior. The distinct electric-response operator is recorded as an unselected extension requiring a complete action.

The pending user question about smooth wave/clock mechanisms versus tiny transfer events remains relevant. No final propagation or matter law was chosen here. T03, T04 and T05 remain incomplete, with their original broader requirements unchanged.

Evidence: `scalar-wave-checks.json` contains the three cases, conservation/convergence diagnostics and a sampled slow-background trajectory. `check_scalar_wave.py` reproduces the calculation. These outputs verify the projected equations derived above, not the full interacting field theory or its stability in all modes.
