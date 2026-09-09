# What companion self-binding would require

## Question and scope

The user suggested graviton attraction analogous to atomic bonding. We separate ordinary massless gravitational waves from a deliberately conditional massive-wave proxy. The latter tests whether a gravitationally supported wave cloud can exist and what capture would cost. It is not adopted as the identity of the companions, a new independent matter population, or an ordinary-graviton calculation.

The distinction is supported by the theoretical literature: [Aoki et al.](https://arxiv.org/html/1710.05606v2) study massive spin-2 geons in a specified bigravity theory, with stability and decay qualifications. Their tensor states are not the scalar state calculated here. We do not import their cosmological motivations or gravity action. A primary scalar-wave comparison, [Silveira and de Sousa](https://www.sciencedirect.com/science/article/pii/037026939501327X), gives a Newtonian spherical binding coefficient against which our independently computed value can be checked.

## 1. Weak gravity is insufficient to balance the pressure of a free massless cloud

A localized massless radiation cloud has stress of the same order as its energy density. Its integrated outward stress is therefore of order its total energy E. In a weak self-generated gravitational potential, binding corrections scale as

```
|W|/E ~ G E/(R c^4) = G M/(R c^2),   M=E/c^2.
```

For gravity characterized by nonrelativistic circular speeds, this compactness is very small. An order-one pressure term cannot be balanced by a correction smaller by many orders of magnitude. Relativistic pressure changes order-one coefficients in this estimate; it does not supply the missing orders of magnitude. This is a weak-field scaling argument, not a solved relativistic geon or an exact collapse threshold.

The earlier [deposit-support derivation](../deposit-support/derivation.md) likewise distinguishes radiation pressure from slow orbital support. Ordinary massless self-confinement must address strong-field geometry, leakage and lifetime. This argument does not exclude compact strongly gravitating objects moving in a weak galactic potential, external confinement, a different force or a different dispersion law. Those are different models requiring their own formation and stability calculations.

## 2. An explicit massive-wave proxy

For a nonrelativistic complex scalar with constituent mass m and conserved occupation N, use the Schrödinger-Poisson Hamiltonian

```
N = integral |psi|^2 d^3x,        M=N m,
H_binding = integral hbar^2 |grad psi|^2/(2m) d^3x
          - G m^2/2 double-integral |psi(x)|^2 |psi(y)|^2/|x-y| d^3x d^3y.
```

Variation gives

```
i hbar psi_t = -hbar^2 Laplacian(psi)/(2m) + m Phi psi,
Laplacian(Phi) = 4 pi G m |psi|^2.
```

The positive gradient term supplies wave pressure; gravity supplies attraction. Rest energy M c^2 is added when comparing total energy with the photon budget. This proxy assumes an occupation already exists. Photon production of that occupation, applicable charge conservation and the actual photon-binding interaction are not derived here.

An analytic normalized Gaussian gives

```
T = 3 hbar^2 M/(4 m^2 R^2),
W = -G M^2/(sqrt(2 pi) R),
R_equilibrium = (3 sqrt(2 pi)/2) hbar^2/(G M m^2).
```

Its radius is a minimum of T+W within the Gaussian family: contraction increases wave pressure faster than attraction. This is a collective cloud, not a collection of assumed atomic bonds.

The code goes beyond that trial shape and solves the coupled spherical boundary-value equations in units hbar=m=G=M=1. With u=r psi,

```
u'' = 2(Phi-mu)u,
M_enclosed' = 4 pi u^2,
Phi' = M_enclosed/r^2.
```

Regularity at the origin, unit total mass, a decaying exterior and Phi(infinity)=0 select the nodeless stationary branch. Independent energy integrals test normalization, the virial relation 2T+W=0, and mu=T+2W. The numerical solution has E_binding about -0.0542564, mu about -0.1627692, and half-mass radius x_half about 3.92510. This is a stationary solution; its full perturbation spectrum, formation and astrophysical lifetime are not calculated.

## 3. Physical scales are predictions of the selected branch

Restoring units gives length scale ell=hbar^2/(G M m^2), velocity scale V=G M m/hbar, and binding energy E0 M V^2, where E0 is the numerical dimensionless energy. At the half-mass radius,

```
R_half = x_half ell,
v_c^2(R_half) = G M/(2 R_half),
m = hbar sqrt(x_half/2)/(R_half v_c).
```

A fixed m therefore predicts a fixed product R_half*v_c for this isolated single-state family. Arbitrarily adjusting m for each object would discard that prediction. Baryonic forces, excited states, multiple components and interactions would change the model and need separate calculations.

The illustrative 10 kpc, 200 km/s cloud requires m c^2 about 1.34e-24 eV. Its ordinary Newtonian mass is about 1.86e11 solar masses, with rest energy about 3.32e58 J. These are consequences of stipulated target scales, not measured masses or a fit. All such energy would still need an accounted source in the fictional universe. The binding correction is only about 1.90e-7 of rest energy and reduces the total; it does not amplify the available energy into additional mass.

Even this massive example is not ordinary two-particle atomic bonding. In the same Newtonian proxy, two equal constituents have

```
alpha_G = G m^2/(hbar c),
a_pair = 2 hbar/(m c alpha_G),
E_pair = -m c^2 alpha_G^2/4.
```

At the illustrative mass above, a_pair is about 2.43e121 metres and |E_pair| about 4.91e-233 eV. These formal weak-coupling scales show that the macroscopic equilibrium relies on a vast collective occupation, not compact pair bonds. They are not calculations of the interactions of ordinary massless gravitons.

## 4. A bound state is not automatically populated

At fixed m, the equilibrium binding energy scales as N^3. Its nonrelativistic chemical potential is mu=dE_binding/dN=3E_binding/N. Adding a slow constituent with kinetic energy K_in requires an outgoing energy channel

```
Q_release approximately K_in - mu = K_in + |mu|.
```

The exact finite-addition calculation uses the difference E_binding(N+DeltaN)-E_binding(N); twelve ledgers verify it. In the illustrative 10 kpc case, the incremental release for an initially slow constituent is about 5.69e-7 of its rest energy. The released energy must appear in radiation, ejected excitations, recoil or another explicit channel. It cannot remain silently counted in both the bound cloud and a second reservoir.

An isolated conservative wave theory can send energy and occupation into escaping waves, but this equilibrium calculation supplies no capture rate or channel efficiency. A classical field solution does not explain how photon-produced companions acquire its coherence. The Schrödinger-Poisson occupation constraint also cannot be assumed to survive an unspecified photon-production interaction.

## 5. Does a bound wave keep stretching light?

A stationary complex bound state has psi(r,t)=psi0(r) exp(-i mu t/hbar), with a time-independent density. An optical coupling to |psi|^2 is therefore static and gives fixed delay, not continuing cumulative stretching under the current ray law. A phase-sensitive coupling could oscillate. It may also break the occupation symmetry used above; its backreaction and charge/energy balance must be derived rather than inferred from the stationary solution.

There is a useful general constraint for any prescribed periodic optical factor n(x,t+P)=n(x,t), with matched stationary clock standards and transparent, single-path propagation. Let F(t) be arrival time for a signal emitted at t. Uniqueness of the ray equation gives

```
F(t+P) = F(t)+P,
J(t) = F'(t)>0,
(1/P) integral over a period J(t) dt = 1.
```

Thus such a repeating field cannot stretch every emission phase by a common factor greater than one. Stretching phases must be compensated by compression phases. For a steady source of equal-energy photons with emission times uniformly covering the cycle, the existing frequency/event relation gives E_out/E_in=1/J. Convexity implies

```
mean_emission_phase(1/J) >= 1/mean_emission_phase(J) = 1.
```

Nontrivial phase variation therefore produces a mean photon-energy gain in this prescribed transparent model, not a uniform energy drain into companions. Maintaining that profile under finite illumination would require an accounted source of work or the field would change. This is not a claim of energy creation or a closed coupled binding calculation. Absorption, phase-selected sources, different clock responses, nonperiodic evolution or a different propagation law require separate treatment.

The optical test uses an explicitly declared cosine profile to check this theorem; it is not fitted to the computed bound-state wavefunction and is not a joined photon-binding simulation. A newly growing occupation or another nonperiodic field state could evade the periodic premise, but its capture, energy supply and time evolution must be calculated.

## What this establishes

The massive proxy demonstrates a conditional collective bound equilibrium and quantifies its binding-energy budget. It does not establish ordinary-graviton bonds, capture, long-term stability, a galaxy population, or sufficient photon supply. Static or merely repeating optical response does not supply the desired persistent redshift under the tested arrival law. The next step must join a specified photon/companion interaction to binding and nonperiodic occupation evolution, including applicable conservation laws, instead of adopting attraction alone as the solution.
