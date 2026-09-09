# Stronger gravitational response must predict motion and lensing together

This is a preliminary T09/T13 derivation for one explicit comparison action. It tests the meaning of the earlier response multiplier eta. It does not select a final field law or close the stronger-gravity path.

## Specified comparison action

Use natural units and signature (-,+,+,+). Introduce a canonical scalar phi and two pressureless matter sectors: ordinary matter b and an assumed deposited sector d. Let:

```
S = integral sqrt(-g) [M_P² R/2 - (partial phi)²/2 - m_phi² phi²/2] d4x
    + S_b[exp(2 beta_b phi/M_P) g, matter_b, photons]
    + S_d[exp(2 beta_d phi/M_P) g, matter_d]
G_* = 1/(8 pi M_P²)
```

The two beta values are constants. The scalar has a positive canonical kinetic term and nonnegative squared mass in the comparison vacuum. These are limited vacuum properties, not a proof of stable deposited configurations or viable cosmology. Ordinary light belongs to the b sector and is minimally electromagnetic in its material metric. There is no added nonconformal photon interaction in this example.

The d sector is assumed to exist for the response calculation; the action does not specify its photon-fed production or permanent capture. If the deposits are instead ordinary matter excitations, treating them as a separately coupled sector would require justification. No new source energy is created by either conformal coupling.

## Linear, static field solution

Expand around phi=0, flat spacetime, and a static nonrelativistic deposited source of mass M_d. Work at weak fields, small conformal perturbations and negligible scalar-stress backreaction. Variation gives, at this order:

```
(nabla²-m_phi²) phi = beta_d rho_d/M_P
nabla² Phi_E = 4 pi G_* rho_d
Psi_E = Phi_E
```

For a point source outside its finite physical radius:

```
phi = -beta_d M_d exp(-m_phi r)/(4 pi M_P r)
Phi_E = Psi_E = -G_* M_d/r
```

The point-source idealization is not used to claim finite field energy at r=0. A finite source and its scalar-field energy are required before computing a complete budget.

Write the ordinary-matter metric as ds_b²=-(1+2 Phi_b)dt²+(1-2 Psi_b)dx². Expanding its conformal factor yields:

```
Phi_b = Phi_E + beta_b phi/M_P
Psi_b = Psi_E - beta_b phi/M_P
Phi_b = -G_* M_d [1+2 beta_b beta_d exp(-m_phi r)]/r
Psi_b = -G_* M_d [1-2 beta_b beta_d exp(-m_phi r)]/r
```

Ordinary slow matter responds to -grad Phi_b. Its radial attraction is:

```
g_b(r) = G_* M_d/r² × [1+2 beta_b beta_d (1+m_phi r) exp(-m_phi r)]
```

Equal-sign couplings strengthen attraction in this regime. Finite scalar range changes the radial shape as well as the amplitude. A single constant eta cannot represent every radius in that case. For extended sources, the Yukawa contribution must be integrated over the source rather than applying a point-mass enclosed-mass formula indiscriminately.

## Light bending from the same potentials

At the stated order, the ordinary-light deflection depends on the transverse gradient of Phi_b+Psi_b. The conformal contributions cancel:

```
Phi_b+Psi_b = 2 Phi_E
deflection = integral grad_perpendicular(Phi_b+Psi_b) dz
point-source magnitude = 4 G_* M_d/b
```

Here b is the impact parameter, and c=1; restore a factor 1/c² in the point-source formula in SI units. This is the standard weak-deflection geometry with asymptotic reference frames, not a full cosmological lens equation with source-distance factors.

Consequently this particular scalar can strengthen motion without supplying the matching direct light-bending enhancement. Scalar stress can itself source the metric at higher order or in other regimes, but its energy and backreaction then have to be computed. They cannot be neglected for the derivation and simultaneously invoked to fix lensing.

The conformal-lensing distinction and investigation of disformal extensions also appear in the primary research article [Lensing with Generalized Symmetrons](https://www.mdpi.com/2674-0346/2/2/9). That is a comparison source, not evidence that an extension already solves our energy, capture or galaxy requirements. The equations here follow directly from the stated linear action and assumptions.

## The laboratory value of G is part of the comparison

G_* is the action parameter, not automatically the measured Newton constant. In the same unscreened theory, for ordinary laboratory sources at distances much shorter than the scalar range:

```
G_lab = G_* (1+2 beta_b²)
eta_motion(r) = [1+2 beta_b beta_d (1+m_phi r) exp(-m_phi r)]
                / (1+2 beta_b²)
eta_lensing = 1/(1+2 beta_b²)
```

These compare the deposited source's contributions with a baseline using G_lab. If the calibration scale is not short compared with the scalar range, the corresponding Yukawa factor belongs in the calibration denominator too. Screening, altered laboratory physics or other interactions would require a new calculation; none is inserted here.

This prevents the scalar enhancement relative to an unmeasured G_* from being mistaken for the same enhancement relative to the G used in the archived galaxy mass estimates.

## Couplings also control deposit self-gravity

In the massless unscreened limit define excess force strengths relative to G_*:

```
A_bb = 2 beta_b²
A_bd = 2 beta_b beta_d
A_dd = 2 beta_d²
A_bd² = A_bb A_dd
```

For example, (beta_b,beta_d)=(0.01,100) gives A_bd=2, A_bb=0.0002 and A_dd=20,000. This is a conditional coupling example, not an observationally allowed halo solution. It shows that suppressing an ordinary-matter force while retaining a cross-sector force can imply a strong force among deposits, which changes their equilibrium, collapse and merger behavior.

The very large response factors in the earlier energy-shortfall tables are not substituted here. Doing so without checking field amplitudes, scalar energy, source structure, local calibration and stability would exceed what the linear calculation establishes. It is possible to investigate other couplings and nonlinear mechanisms; each needs those same accounting checks.

## Numerical checks

The companion script checks four force-gradient cases by differentiating the derived potential numerically, eight line-of-sight deflection integrals using both potentials, and three coupling-product identities. All passed. The largest relative force-gradient error was 1.19×10^-10.

The check parameters are mathematical examples with modest cross-force enhancements. They do not establish galaxy fits, observational bounds, a stable source, or validity of arbitrary large beta values. The formulas are checked only within their declared approximation.

## Consequences for the research queue

1. Keep this conformal scalar as an explicit motion/lensing comparison branch. Do not use its motion enhancement as an automatic solution of the photon-energy and lensing budgets.
2. For a candidate that aims to enhance lensing too, derive an interaction affecting the metric sum seen by light or a fully accounted additional stress source. An extra lensing multiplier cannot be added independently.
3. Derive deposited-sector equilibrium and perturbations using its own self-coupling, not just the force it exerts on stars.
4. Calculate the measured local G and source/detector calibration within the same law before comparing with the archived inferred masses.
5. Join this gravity sector to a derived production, transport and capture model before treating it as a unified candidate.

No physical branch was adopted or rejected generally. T09 and T13 remain incomplete: this is one weak-field comparison, without complete nonlinear field solutions, stability, lensing data or a photon-fed deposited sector.

Files: `gravity-response-checks.json` and `check_gravity_response.py` contain the numerical evidence and reproducible checks.
