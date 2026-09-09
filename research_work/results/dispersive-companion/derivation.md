# A dispersion scale limits companion energy and momentum together

## Why this is the next test

The linear slow-wave model tied small energy packets to large photon deflections. That tradeoff was specific to its scale-free dispersion. Here a positive spatial-gradient term introduces a fixed momentum scale and limits the band of companion waves a photon can emit. This can make energy transfers and angular kicks small together. The important remaining question is whether all photon colors then lose the same fraction of energy.

Without checking the rate as well as the packet size, we could mistake an improvement in image preservation for a complete redshift mechanism. This is an optional effective-field comparison, not a newly identified material, ordinary graviton, or adopted universe law.

## 1. One action supplies the scale

In units hbar=c_photon=1 and the same preferred frame as the preceding pass, use

    L = (1+g phi) E_field^2/2 - B^2/2
        + phi_t^2/2 - v^2 |grad phi|^2/2 - (Laplacian phi)^2/(2 M^2),
    0<v<1, M>0, 1+g phi>0.

There is one scalar degree of freedom, electromagnetic gauge invariance and translation/rotation symmetry. Only spatial derivatives are added, so the time equation remains second order. The new scalar Hamiltonian contribution +(Laplacian phi)^2/(2M^2) is positive. The scalar dispersion is

    omega(Q)^2 = v^2 Q^2 + Q^4/M^2.

M is a declared model scale, not a value fitted to observations or an independently derived material property. The coupling and illuminated-background limitations of the preceding scalar/electric model remain. Radiative corrections, matter clocks and a microscopic source of the preferred frame are not supplied.

This spatially dispersive model is not a demonstrated Lorentz-invariant or causal ultraviolet completion. Group speeds may exceed the photon speed, and the formal polynomial dispersion is not to be extended to arbitrarily high momentum. Introduce an independent validity cutoff Lambda_UV. The calculation using the entire emission band requires photon energies and the band edge to remain below that cutoff, with perturbative rates. The hierarchy Q_c much less than E much less than Lambda_UV is an explicit condition, not an established hierarchy in nature. Cutoff sensitivity of scalar Cherenkov calculations is also emphasized, for a different coupling, by [Dalang, Fleury and Lombriser](https://arxiv.org/abs/2109.10812); their rates are not used here.

The equation phi_tt-v^2 Laplacian(phi)+Laplacian^2(phi)/M^2=(g/2) E_field^2 and the scalar energy flux

    S_phi = -v^2 phi_t grad phi
            + [phi_t grad(Laplacian phi) - (Laplacian phi) grad phi_t]/M^2

give the same local transfer +(g/2) phi_t E_field^2 into scalar energy and its negative in the electromagnetic ledger. Translation symmetry accounts for companion momentum. No unrecorded target supplies energy or recoil.

## 2. The emission band follows from kinematics

For gamma(k)->gamma(k')+phi(Q), put E'=E-omega(Q) and let mu be the angle cosine between the two photons. Energy/momentum conservation gives

    mu = 1 - [Q^2-omega(Q)^2]/[2E(E-omega(Q))].

Emission requires both

    omega(Q)<=Q,
    Q+omega(Q)<=2E.

The first yields a fixed upper band edge

    Q_c = M sqrt(1-v^2).

When E>=Q_c the whole interval 0<Q<Q_c is available. For E<Q_c the upper limit is instead the root Q+omega(Q)=2E. The finite band is a consequence of the action's dispersion, not an independently chosen low-pass filter. The endpoints are understood as limits of the phase-space integral.

## 3. Rate from the same electric vertex

The polarization-averaged vertex remains

    mean |M_vertex|^2 = g^2 E^2(E-omega)^2 P(mu),
    P(mu)=(1+mu^2)/2.

Integrating the energy delta function in the angle between the incident photon and Q gives

    dGamma/dQ = g^2 (E-omega)^2 Q P(mu)/(16 pi omega),
    alpha(E) = (1/E) integral omega dGamma
             = g^2/(16 pi E) integral Q(E-omega)^2 P(mu) dQ.

This is the spontaneous per-photon fractional energy-loss rate, with initially empty companion modes, not a full populated-bath transport solution. It includes all outgoing directions. At low E/M it reproduces the preceding linear-dispersion result and its E^3 color scaling.

For E much greater than Q_c, omega/E and the deflection are small throughout the band. Consequently

    alpha(E) = g^2 E Q_c^2/(32 pi) [1+O(Q_c/E)].

The added scale reduces the leading color dependence from cubic to linear in this regime. It does not remove it.

There is a stronger statement throughout the full-band regime E>=Q_c. Every Q in the integral obeys Q<=E, so the outgoing photon is in the forward hemisphere: mu>=0. As E increases, both (1-omega/E)^2 and P(mu) increase at fixed Q. Therefore

    alpha(E)/E is nondecreasing,
    alpha(2E)/alpha(E) >= 2, for E>=Q_c.

This establishes failure of exact achromatic fractional drift in the regime used to obtain small packets and deflections. It does not claim that every possible transition-band interval, dispersion, vertex or finite-temperature response obeys this bound.

## 4. Improvement in packet size and direction is real but incomplete

At fixed v and M in the high-energy regime, typical companion energy stays of order Q_c, so fractional packet size decreases as Q_c/E. Photon deflections also decrease as Q_c/E. This avoids the earlier linear-dispersion tradeoff. It does not make fractional loss independent of photon energy.

For the illustrative v=0.5 parameter, not fitted to data:

| E/Q_c | Mean fractional loss per event | RMS deflection per event | alpha(2E)/alpha(E) |
| ---: | ---: | ---: | ---: |
| 0.01 | 0.20236 | 38.08 degrees | 7.995 |
| 1 | 0.16171 | 17.08 degrees | 4.312 |
| 10 | 0.03138 | 1.739 degrees | 2.118 |
| 100 | 0.003273 | 0.1724 degrees | 2.011 |
| 10,000 | 0.00003288 | 0.001722 degrees | 2.00011 |

These are single-event all-direction quantities, not simulated astronomical line widths, images or arrival times. Many events accumulate both fluctuations and angular redistribution; no observational image success is inferred from a small single-event angle. Absolute conversion strength still depends on g and the allowed validity range.

## 5. What this changes in the research plan

The last tradeoff is not universal: a dispersion scale can make both transfers and deflections small. This is useful positive structural evidence. The same calculation also shows why this particular electric vertex still cannot give the common fractional redshift law in its high-energy, small-packet regime.

Retain this as a scoped comparison. Do not set g proportional to E^-1/2 by hand to cancel the remaining color dependence; that is a new frequency-dependent interaction requiring its own derivation, causality, fluctuations and energy accounting. A different protected coupling or a derived dynamical response is needed before fitting astronomical redshifts. The same-action clock/timing, capture, supported deposition, gravity and source-history requirements remain open.

## Checks and scope

The script compares independent momentum and frequency integrations, random energy/momentum configurations, the low-energy linear-dispersion limit and the high-energy finite-band asymptote. It samples three v values and seven E/Q_c values without fitting data. Positive free-mode energies follow analytically from the action. None of these checks proves full background stability, quantum consistency, a real material or observational agreement.
