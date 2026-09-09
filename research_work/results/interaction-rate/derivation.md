# From an energy-loss formula to a physical interaction rate

This advances the conversion-first goal beyond its single resonant mode example. It derives necessary restrictions for particular interactions and compares explicit transport kernels. No complete photon-graviton theory, astronomical validation, special time-stretching mechanism, or independent dark-matter source is adopted.

## 1. What an interaction has to predict

Use physical path length s. Let W(E,epsilon) be the differential probability per unit path for a photon of energy E to give a fraction epsilon of its energy to a receiving sector, while leaving one lower-energy photon. Then the conditional first and second jump moments are:

```
alpha(E) = integral epsilon W(E,epsilon) d epsilon
d<Delta E>/ds = -E alpha(E)
d<Delta E^2>/ds = E^2 integral epsilon^2 W(E,epsilon) d epsilon
```

These are forward-channel moments. Inverse channels add a positive energy drift and additional variance. The second jump moment is not itself the time derivative of the total ensemble variance; the full distribution equation supplies that derivative.

The fitted alpha in the preceding pass is approximately 0.0002488993 per Mpc. For the achromatic exponential law to apply, the *net* fractional first moment must be approximately energy independent across the photon bands being tested. This is a target for an interaction calculation, not a measured microscopic cross section.

Writing a differential loss-energy density R(E,Delta), with units events per path per energy, one sufficient scale-free form is:

```
R(E,Delta) = Gamma/E * p(Delta/E),    0 < Delta < E
integral_0^1 p(epsilon) d epsilon = 1
alpha = Gamma * mean(epsilon)
```

Here Gamma is an event rate per length, not the capture rate used elsewhere. Gamma and p must come from a physical coupling, available receiving modes, occupations, and momentum constraints. Many distributions have the same first moment but different line profiles and direction changes.

## 2. A fixed companion-energy quantum is not automatically achromatic

Suppose each transfer has fixed energy q and the event rate is independent of photon energy, with negligible probability of exhausting a photon. The mean loss per path is Gamma*q, an absolute loss. It gives E_received approximately E_emitted-Gamma*q*R, not exponential fractional loss.

Calibrate this model to give z=0.1 for a 1 eV reference photon. The mean transferred energy is 0.090909 eV. At the same path, our calculation gives:

| Initial photon energy | Fixed-quantum model shift | Fractional-loss model shift |
|---|---:|---:|
| 0.5 eV | 0.22222 | 0.10000 |
| 1 eV | 0.10000 | 0.10000 |
| 2 eV | 0.04762 | 0.10000 |
| 4 eV | 0.02326 | 0.10000 |

This table compares centroid-energy shift factors, not a fitted observed spectrum. Its approximately 383 Mpc path is an illustrative extrapolation beyond the 10–93 Mpc distance range used for the preceding fit. A fixed-q model could compensate with an energy-dependent event rate; then that dependence must be derived and its noise tested. This result does not exclude all interactions with preferred energy scales.

## 3. Does the simple scalar-photon action give vacuum partial conversion?

Consider the comparison operator L_int proportional to phi F_mu_nu F^mu_nu, in a locally flat region with standard photon dispersion and a stable scalar of nonnegative mass squared. A scalar is not an ordinary graviton; this operator is a test case already related to the repository's scalar-field branch. Scalar/photon interactions in specified magnetic media are discussed by [Ganguly et al.](https://arxiv.org/abs/1907.00565); that source is background context, not evidence for the fictional transfer law. The following vacuum result is derived directly here.

Use signature (+---), c=hbar=1. For gamma(k) -> gamma(k') + phi(q), with real future-directed photons:

```
q = k-k'
q^2 = -2 k.k' = -2 E E' (1-cos(theta)) <= 0
```

A massive outgoing scalar would require q^2=m_phi^2>0, which is unavailable in this empty-space channel. A massless scalar requires collinear photon momenta. The scalar two-photon vertex, up to its normalization and an overall phase, is:

```
V = (k.k') (epsilon.epsilon') - (k.epsilon') (k'.epsilon)
```

For real collinear transverse photons each term vanishes. Replacing either polarization by its photon momentum also makes V vanish, as required by electromagnetic gauge invariance. One hundred numerical kinematic cases check these identities and the spacelike noncollinear q. This is a tree-level on-shell statement for this operator and these dispersion assumptions, not a theorem about all photon-companion interactions.

It explains why a successful two-state Hamiltonian toy does not automatically yield a real vacuum conversion rate. An interaction with matter, an appropriate dynamical environment, modified dispersion, or a different completed field theory requires a new amplitude and momentum/energy calculation. Stationary linear photon-mode mixing in a prescribed static field remains a same-frequency conversion comparison; it cannot be borrowed as an inelastic redshift law.

## 4. Reverse processes and spectral noise

The Hermitian conjugate of an emission vertex permits the inverse reaction. For a specified bosonic receiving mode, creation/annihilation matrix elements have occupation factors n+1 and n. As a deliberately limited kinetic example, let a dilute photon jump on an equally spaced energy ladder by q with energy-independent rates:

```
downward rate = A (n+1)
upward rate   = A n
mean(E(s)) = E(0) - A q s
Var(E(s)) = Var(E(0)) + A q^2 (2n+1) s
```

The moment formulas hold away from ladder boundaries for independent constant-rate jumps. This is not a universal formula for photons in a thermal bath: the matrix elements, state densities, photon occupation factors and energy dependence generally change the rates. A is a stipulated rate per length, and n is a fixed occupation in the reservoir approximation, not a newly assumed primordial companion population.

We solved a 2,001-state master equation by sparse matrix exponentiation. With q=0.001 eV, E(0)=1 eV, A*s=100, and n=0,1,10, all cases end at mean E=0.9 eV. Their variances are respectively 0.0001, 0.0003 and 0.0021 eV squared. Normalization, boundary probabilities and analytic moments pass the specified tolerances. Every upward transfer debits the environment and every downward transfer credits it; net environment gain is 0.1 eV per photon in each case.

This open-system calculation does not prove a finite environment can keep its occupation fixed indefinitely. A completed model must evolve or replenish that environment and include its energy. It also demonstrates why reporting only net lost photon energy can hide large spectral changes. The chosen q and occupations are illustrative, not observationally constrained values.

## 5. Turning spectral precision into a rate requirement

For the earlier independent equal-fraction Poisson model:

```
relative energy RMS^2 = S^epsilon - 1
epsilon <= ln(1+delta^2)/ln(S)
Gamma >= alpha * ln(S)/ln(1+delta^2)
```

S is the centroid energy-shift factor and delta an allowed added fractional RMS. At S=2 and an illustrative width delta=(10 km/s)/c, our fitted alpha requires epsilon at most 1.6052e-9 and Gamma at least 155,057 events per Mpc, corresponding to a mean spacing at most about 6.45 parsecs. The long path accumulates about 4.32e8 events. This does not establish that such events exist or that 10 km/s is the observational bound; a real line-profile analysis must supply the bound and instrumental/intrinsic widths.

These requirements apply to independent equal-fraction events. Coherent transfer, correlated events, different energy kernels or inverse reactions need their own spectral calculation. A mean spacing alone gives no scattering angle or image-blur prediction.

## 6. Consequences and next computation

The next candidate must specify an interaction and a receiving environment that jointly predict:

1. A nonzero allowed inelastic amplitude, with recoil and ordinary or explicitly modified dispersion.
2. Approximately energy-independent net fractional loss, or a quantified, testable departure from it.
3. Small enough line-profile distortion, including reverse events and companion occupations.
4. Photon direction, polarization and event-arrival transfer from that same process.
5. Companion trajectories, capture and a supported stored state supplied by the accounted photon losses.

We have not selected an environment or a new fundamental law. Matter-assisted and collective-field interactions remain permissible research branches. The present scalar vacuum example and fixed-q example are constrained under their assumptions; they do not close the overall fictional hypothesis. The timing and motion/lensing gaps remain unchanged.

Results and verification are in [rate-constraints.json](rate-constraints.json); [check_interaction_rate.py](check_interaction_rate.py) reproduces them. This is preliminary T03/T04/T05 evidence within the original 20-task, 32-area scope.
