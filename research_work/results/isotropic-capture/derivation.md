# All-direction illumination is not automatically balanced capture

## Why this calculation matters

The last calculation showed that perfectly balanced delivered momentum reduces the radial shove from capture. But a moving receiver does not sample an ambient radiation bath equally from every direction. We now calculate that sampling instead of assuming cancellation. Otherwise, a supposedly balanced store could still slow down and move inward.

We also ask whether the receiver's preference for different companion energies changes the answer. This is an optional comparison with specified absorption cross sections, not a selected material or a new final law. No dark-matter population, expansion or cosmic chronology is imported.

The [radiation-force framework of Bini, Jantzen and Stella](https://arxiv.org/html/0808.1083) provides related context for coupling radiation stress to motion. Our equations below are independently derived for full absorption that increases receiver mass; their constant-mass scattering problem is not substituted for it.

## 1. Count the energy and momentum actually intercepted

Use a local lab frame with isotropic monochromatic companion energy density u. Let a small spherical receiver move at speed beta c along the z axis. Its rest-frame absorption cross section sigma(epsilon') depends on the energy it sees. For a lab ray with directional cosine mu,

    epsilon' = gamma epsilon (1-beta mu).

The lab encounter rate has the relative-flux factor 1-beta mu. This can be obtained by transforming each ray's number density to the receiver frame and converting proper time back to lab time. Therefore

    P = c u/2 integral[-1,1] sigma(epsilon') (1-beta mu) dmu
    F_z = u/2 integral[-1,1] mu sigma(epsilon') (1-beta mu) dmu.

P is lab energy entering the receiver per time; F_z is its lab momentum gain per time. The receiver has total energy E=gamma M c^2 and momentum p=gamma M v. Its radiation-induced velocity change is

    dv/dt = [F_z - v P/c^2]/(gamma M).

These integrals assume complete absorption, a rest-frame orientation-independent cross section, a prescribed bath, and no outgoing radiation. Actual target geometry, shelving radiation, saturation and reverse processes require additional terms.

## 2. Energy-independent absorption gives a definite slowing

For constant sigma, angular integration gives exactly

    P = c sigma u
    F_z = -sigma u beta/3
    dv/dt = -(4/3) P v/(gamma M c^2)
    dM/dt = gamma P/c^2 (1+beta^2/3).

Thus the incoming field has zero net momentum in the lab, but the absorbed subset does not. The receiver preferentially encounters opposing rays. This differs from imposing zero delivered momentum by hand.

An independent frame check transforms the radiation stress to the receiver's momentary rest frame:

    u' = gamma^2 u(1+beta^2/3)
    F'_z/c = -(4/3) gamma^2 u beta.

Absorbing that rest-frame energy and momentum, transforming back, and converting proper time to lab time reproduces P and F_z. Numerical angular integration checks both routes at four speeds.

For a freely moving receiver in a constant bath, with c=1 and time units chosen so P=1, the exact evolution is

    E(t)=E0+t
    p(t)=p0 [E(t)/E0]^(-1/3)
    v(t)=v0 [E(t)/E0]^(-4/3)
    M(t)=sqrt(E(t)^2-p(t)^2).

The code compares integrated energy, momentum and mass growth to this solution at three initial speeds. The incoming energy increases total receiver energy; its kinetic energy can decrease while rest energy grows. This is a local prescribed-bath test, not a closed-universe energy supply.

## 3. Conditional orbital implication

For slow motion and slowly varying nearly circular orbits in a fixed dominant potential, the actual bath force gives

    dlnL/dlnM approximately -1/3,

where L is the receiver's total angular momentum. The factor arises from F_z, while the increasing mass separately changes its angular momentum per unit mass. Consequently,

    Fixed flat circular speed: r proportional to M^(-4/3)
    Fixed Kepler potential:    r proportional to M^(-8/3).

A mass doubling gives radius ratios about 0.397 and 0.157, respectively. The earlier perfectly balanced delivered-momentum comparison gave 0.5 and 0.25. Neither comparison is a full orbit simulation with growing self-gravity. Their difference shows why the angular absorption law matters; these numbers are not predictions for observed galaxies.

## 4. Frequency selection can change the result

Take the comparison sigma(epsilon')=sigma0(epsilon'/epsilon0)^s, with epsilon0 equal to the monochromatic lab energy. This law is stipulated only over the Doppler-sampled band, not as a physical cross section valid to zero or infinite energy. Defining P0=c sigma0 u gives

    P/P0 = gamma^s/2 integral (1-beta mu)^(s+1) dmu
    c F_z/P0 = gamma^s/2 integral mu(1-beta mu)^(s+1) dmu.

At small speed,

    c F_z/P approximately -(s+1) beta/3
    dv/dt approximately -(s+4) P v/(3 M c^2).

| Cross-section energy power s | Small-speed deceleration coefficient multiplying P v/(M c^2) |
|---|---:|
| -6 | -2/3: acceleration along existing motion in this driven comparison |
| -4 | 0 |
| -1 | 1 |
| 0 | 4/3 |
| 2 | 2 |
| 6 | 10/3 |

For s=-4 the ideal angular integrals give P=P0 and c F_z/P0=beta exactly, so F_z=vP/c^2 and the free receiver's velocity does not change. Incoming radiation supplies the momentum needed as the receiver's inertia grows. This cancellation is verified at three speeds; it is not proof of stable galaxy formation, finite capacity, or a material with that response.

The negative-damping comparison is a driven, nonthermal, monochromatic absorption model consuming radiation energy, not perpetual motion or an equilibrium heat bath. Thermal populations and reverse transitions cannot be omitted when assessing an actual receiver. A local cancellation at one band likewise does not ensure cancellation over a broad companion spectrum or after occupation changes.

This is a useful conditional escape from the constant-cross-section slowing result. It identifies a response property to derive and test, not a permission to choose an arbitrary function solely to preserve an orbit. In particular, the earlier oscillator's response-amplitude tail is not automatically this absorption cross section: vertex factors, phase space, spectral thresholds and all exit channels intervene.

## 5. Absorption also makes an external isotropic bath locally directional

The recovered external-capture comparison has uniform incident intensity around a sphere of radius R and constant absorption coefficient kappa. For a stationary receiver at r, a ray with radial direction cosine mu has travelled

    ell = r mu + sqrt(R^2-r^2+r^2 mu^2)
    I(r,mu)=I0 exp(-kappa ell).

The delivered radial momentum per absorbed energy divided by c is

    xi(r) = integral mu I(r,mu) dmu / integral I(r,mu) dmu.

Numerical integration shows xi=0 at the center and xi<0 away from it: attenuation leaves a net inward-directed local field. At the outer boundary, with tau=kappa R,

    A0=1+(1-exp(-2tau))/(2tau)
    A1=-1/2+[1-(1+2tau)exp(-2tau)]/(4tau^2)
    xi(R)=A1/A0.

These analytic boundary expressions independently check the angular calculation. As tau grows, xi(R) approaches -1/2. This adds an inward force to the earlier outer-heavy deposition profile; it does not by itself transport or support the stored reservoir. We have not combined this anisotropic attenuated field with a moving, frequency-selective receiver in a self-consistent evolution.

## What follows

An all-direction incoming bath cannot simply be labeled momentum-balanced at every receiver. Motion, attenuation and frequency-dependent absorption determine the delivered force. The constant-cross-section comparison slows receivers, while a particular stipulated spectral response cancels the velocity change in the ideal free-receiver calculation.

The next physical candidate must jointly derive the receiver's spectrum and bandwidth, actual illumination, all outgoing/reverse channels, finite storage capacity and evolving orbits. Its energy source and gravitational response must then predict both dynamics and lensing. No response law or receiver is adopted here, and the original 20 tasks and 32 observational requirements remain unfinished.

Run `python -X utf8 research_work/results/isotropic-capture/check_isotropic_capture.py`; saved evidence is in `isotropic-capture-results.json`. These are mathematical and transport diagnostics, not astronomical validation.
