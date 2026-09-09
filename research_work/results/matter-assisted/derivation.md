# An explicit matter-assisted interaction and its predictions

This pass supplies a concrete local interaction involving an ordinary material target, a surviving photon and a scalar companion. It moves beyond assigning a conversion probability by hand. The coupling and target abundance are still unspecified, so it predicts a rate's shape and dependence, not the observed absolute rate. The scalar is a comparison companion, **not** an ordinary spin-two graviton. No target population may be inserted as independent dark matter.

## 1. A target makes energy and momentum balance possible

Consider gamma + T -> gamma' + companion + T', with the target initially at rest, unchanged rest mass M, and a massless companion. Use c=1. Let initial photon energy be E, final photon energy E', initial/final photon directions n,n', and companion direction m. Define:

```
Delta = E-E'
B = E n-E' n'
omega = [M Delta-E E'(1-n.n')] / [M+Delta-B.m]
P_target_final = B-omega m
K_target = sqrt(M^2+|P_target_final|^2)-M
Delta = omega+K_target
```

The formula follows by imposing the final target mass shell, (M+Delta-omega)^2-|B-omega m|^2=M^2. A valid branch requires positive omega and future-directed final target energy. This is exact kinematics, not an interaction probability. Of 300 chosen target scenarios with randomized outgoing directions, 295 have positive companion energy and satisfy the conservation checks; five are disallowed by the specified energies and directions rather than counted as valid reactions.

When the photon stays straight, B=Delta n and:

```
omega = M Delta / [M+Delta(1-n.m)]
```

A companion can leave in a different direction while the target receives recoil. Thus partial conversion does not kinematically require photon-image blurring. Whether an actual interaction favors straight photons is a separate question. Unlike the earlier empty-space scalar vertex, this reaction has an additional momentum receiver.

## 2. A local polarizability interaction

Use signature (+---), hbar=c=1, standard Maxwell photons and a canonical massless scalar chi. For each existing polarizable material target, parameterize its trajectory X(tau), four-velocity u=dX/dtau, u.u=1, and rest-frame electric four-vector E_mu=F_mu_nu u^nu. One possible effective action is:

```
S_bulk = integral d^4x [-F_mu_nu F^mu_nu/4 + (partial chi)^2/2]
S_target = -M integral d tau
S_polarizability = -1/2 integral d tau [a_E + beta chi(X)] E_mu E^mu
```

In the target rest frame the last integrand is +(a_E+beta chi)|electric_field|^2/2. a_E is the ordinary electric polarizability coefficient; beta describes its linear coupling to the companion. beta has energy dimension -4. This is a **hypothesized effective interaction**, not a measured property of the recovered galaxies. The target must be specified as permitted matter with a justified response and abundance. Its motion is dynamical, rather than an external momentum sink with missing energy.

The expression is electromagnetic-gauge invariant because it uses F. It is autonomous and translation invariant with the dynamical target included, so the full equations carry energy and momentum conservation. These properties do not establish a stable ultraviolet completion or unlimited-range validity. The point-target, low-frequency expansion requires photon/companion wavelengths larger than target structure scales, energies below the effective cutoff, and small perturbations. In particular a_E+beta chi cannot be extrapolated to arbitrary scalar amplitudes as a proven stable constitutive law.

Expanding the beta term gives a contact process with two photon legs and one scalar leg on the target. In the leading heavy-target rest-frame approximation, define its matching normalization by:

```
amplitude = 2 M beta E E' (polarization_initial . polarization_final*)
```

This normalization is the convention used below. Covariantly each electric leg is k_mu(epsilon.u)-epsilon_mu(k.u), so replacing a polarization by its own momentum gives zero. The numerical checks verify both Ward identities. Summing final polarizations and averaging the two initial polarizations gives A(theta)=(1+cos^2(theta))/2. This is nonzero at zero photon deflection, but is not sharply concentrated there.

The scalar example is not the only matter-assisted possibility. A primary calculation of charged-scalar scattering with photons and an ordinary graviton exists in [Ahmadiniaz et al., Compton-like scattering of a scalar particle with N photons and one graviton](https://arxiv.org/abs/1908.03425). That is an ordinary-tensor comparison branch, not the source of the polarizability formulas derived here and not proof of an adequate astrophysical conversion rate. Its amplitude, infrared treatment and target population require their own calculation.

## 3. Rate and emitted-energy distribution

For E much smaller than M, neglect recoil only in the leading phase-space expression: omega=E-E'. Use x=E'/E and the relativistic scattering normalization flux=4ME. Integrating the final target momentum and energy delta functions gives:

```
d sigma / (dx dOmega_photon dOmega_companion)
  = beta^2 E^6 x^3(1-x) A(theta) / [8(2pi)^5]

sigma_total(E) = beta^2 E^6 / (480 pi^3)
p(x | event) = 20 x^3(1-x),  0<x<1
mean(1-x) = 1/3
mean((1-x)^2) = 1/7
```

The companion angle is isotropic in this leading contact approximation. The photon-angle density in mu=cos(theta) is 3(1+mu^2)/8, whose mean mu is zero. The runner independently integrates the energy phase space and sums explicit polarization vectors to check the coefficient, moments and E^6 scaling. Finite-recoil corrections are not supplied by calling the leading amplitude exact; only the kinematics in section 1 is exact at finite M.

For dilute uncorrelated stationary targets of number density n_T, and negligible inverse-channel occupation:

```
Gamma(E) = n_T sigma_total(E)
alpha_initial(E) = Gamma(E)/3 = n_T beta^2 E^6 / (1440 pi^3)
```

Neither beta nor n_T is obtained from the redshift data here. Their product could be adjusted to match one reference energy's initial drift, but that would be calibration, not a first-principles prediction. The same constants then make the initial fractional loss 64 times larger when photon energy doubles. Target form factors, dispersion, resonances, ordinary elastic scattering and inverse reactions have been omitted; adding them requires a new derived rate rather than retaining this simple formula.

## 4. Propagating a photon population with the derived kernel

Large losses per event mean the mean energy does not obey a closed exponential ODE. In this approximation:

```
d< E >/ds = -[n_T beta^2/(1440 pi^3)] <E^7>
```

Replacing <E^7> by <E>^7 would be uncontrolled once the spectrum broadens. Instead the runner simulates 20,000 independent histories for each of four initial energies, drawing waiting paths from the energy-dependent rate and the retained fraction x from the derived Beta(4,2) distribution. It credits every lost unit to companion energy in the heavy-target limit; finite recoil belongs in the target ledger when that approximation is relaxed.

For illustration only, normalize the initial fractional drift at E_ref to the preceding empirical alpha_ref and propagate alpha_ref R=0.1. This corresponds to about 402 Mpc, outside the preceding nearby fit range. The results are:

| Initial E/E_ref | All-direction centroid-energy shift | Added fractional energy RMS | Fraction with no event |
|---|---:|---:|---:|
| 0.5 | 0.00135 | 0.02297 | 0.99555 |
| 1 | 0.09837 | 0.19335 | 0.73970 |
| 2 | 1.26540 | 0.29877 | 0 |
| 4 | 3.51989 | 0.29910 | 0 |

These are stochastic simulation estimates, with standard errors reported in the JSON, not observed galaxy redshifts or an observational fit. In particular the centroid counts photons scattered into **all directions**. It is not the spectrum a narrow telescope beam necessarily receives. Unscattered photons retain their original line frequency; scattered photons have a broad distribution and substantial angular redistribution. A ray/acceptance calculation is required for an observed spectrum. The analytic no-event fraction is exp[-Gamma(E_initial)R] and is checked against simulation; it is not a loss of total photon number.

## 5. What this candidate establishes and what it fails to supply

- It supplies an explicit local interaction with a nonzero inelastic amplitude and a target that can receive recoil without special time/void behavior.
- It supplies a conditional cross section, energy kernel and photon-angle distribution from the same leading approximation. It is more specific than a chosen exponential energy-loss equation.
- In this unsuppressed contact limit, it produces strongly color-dependent loss, large energy jumps and angular scattering. It therefore does not supply the desired achromatic narrow-line propagation law.
- It does not identify ordinary gravitons, create a stable deposited reservoir, derive gravity/lensing, or stretch whole events. Inverse channels and target populations remain to be treated.

These are restrictions on this candidate in its stated regime. They are not a rejection of every interaction with matter or every companion mechanism. Potential revisions include a derived collective response that favors small energy and momentum transfers, or an ordinary-tensor emission calculation. A tuned beta(E) proportional to E^-3 would cancel the explicit E^6 rate factor in this formula, but neither its origin nor the resulting spectrum would thereby be explained; it is a target response, not an adopted fix.

The next calculation should connect a permitted target's frequency/momentum response to the joint transfer kernel. The requirements for small fractional loss, narrow angles, reverse events, energy supply and event timing must remain linked. Full goal completion remains unproved.

See [numerical results](matter-assisted-results.json) and [reproducible runner](check_matter_assisted.py). This is preliminary T03/T04/T05 work; the original 20 tasks and 32 observation areas remain in scope.
