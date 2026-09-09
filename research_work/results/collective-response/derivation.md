# Does spatial coherence repair the matter-assisted candidate?

This pass resolves the spatial extent of the preceding scalar/polarizability target. It tests a specific collective modification rather than assuming that the word "coherent" makes an interaction achromatic or lossless. No ordinary-graviton identity, material target population, special clock behavior or final field law is adopted.

## 1. A spatially extended target

Use the preceding local scalar/polarizability vertex, but distribute its coupling over a normalized Gaussian target profile of width L. Work in the coherent channel that leaves the target's internal state unchanged, with negligible center-of-mass recoil energy in the heavy-target limit. Contributions from its different locations add as amplitudes:

```
f(r) = exp[-r^2/(2L^2)] / (2pi L^2)^(3/2)
F(Q) = integral f(r) exp(i Q.r) d^3r = exp[-L^2 Q^2/2]
|F(Q)|^2 = exp[-L^2 Q^2]
Q = E n - E' n' - omega m,       omega = E-E'
```

Here hbar=c=1, Q is the momentum delivered to the target, and n,n',m are the incoming photon, outgoing photon and companion directions. At finite target mass its recoil energy must be restored as in the previous pass. The Gaussian shape and unchanged internal state are specified assumptions, not a derivation of a rigid galactic medium. Binding, finite response times, constituent matter, ordinary scattering, occupation and target abundance must all be supplied for a physical candidate.

Keep the **total** polarizability coefficient beta fixed when comparing the form factor to the point target. Increasing the size of a real object can also increase its number of constituents and total beta. The calculations below isolate the geometrical suppression; they do not assert that every larger physical target has a smaller total cross section.

## 2. Straight paths do not require small energy transfer

Let x=E'/E and y=1-x. If both outgoing particles travel in the incident direction, then:

```
Q = E n - x E n - (1-x) E n = 0
```

This is true for any split 0<x<1. The spatial form factor equals one even when the companion carries a large fraction of the initial energy. Spatial coherence can favor nearly collinear particles, but supplies no small-energy-transfer condition by itself. A response that suppresses large companion frequencies is a distinct physical ingredient.

## 3. Integrate the actual angular distribution

The preceding point-target differential cross section was proportional to:

```
beta^2 E^6 x^3(1-x) A(mu) dx dOmega_photon dOmega_companion
A(mu)=(1+mu^2)/2,    mu=n.n'
```

Define b=E L (or E L/(hbar c) in conventional units), and v=|n-x n'|. Integrating the companion direction exactly and changing the photon-angle variable from mu to v gives the angular factor:

```
I_b(x) = 2pi^2/[b^2 x(1-x)] * integral_(1-x)^(1+x)
          A[(1+x^2-v^2)/(2x)]
          {exp[-b^2(v-1+x)^2]-exp[-b^2(v+1-x)^2]} dv
```

For computation use t=b[v-(1-x)] so that the narrow forward peak remains resolved as b increases. The numerical runner also uses expm1 to avoid cancellation between the exponentials. The point limit is I_0=32pi^2/3. Independent 256- and 512-node energy quadratures agree to the recorded tolerances, and the Gaussian Fourier transform is checked by direct integration.

## 4. Large-target limit

For b much larger than one, away from shrinking endpoint regions, the angular integral is dominated by v approximately 1-x:

```
I_b(x) approximately pi^(5/2) / [b^3 x(1-x)]
p(x | event) -> 3 x^2
mean(1-x) -> 1/4
mean((1-x)^2) -> 1/10
sigma / sigma_point -> 5 sqrt(pi)/(8 b^3)
sigma -> beta^2 E^3 / [768 pi^(5/2) L^3]
mean(photon_angle^2) -> 1/[b sqrt(pi)]
```

The endpoint regions do not dominate these integrated leading moments; the finite-b quadrature tests their approach to the limit. The normalization uses the same amplitude convention as the preceding matter-assisted derivation.

Two results follow together: angles become smaller, while typical energy transfers remain large. The rate's energy dependence weakens from E^6 to approximately E^3 for a fixed target and fixed beta. It does not become energy independent. These are conditional low-energy/large-spatial-extent limits; a completed target response can have additional frequency dependence.

## 5. Numerical results

| b = E L/(hbar c) | Rate / point-target rate | Mean fractional energy loss per event | Photon-angle RMS, radians |
|---:|---:|---:|---:|
| 0, point limit | 1 | 0.33333 | 1.74525 |
| 1 | 0.30453 | 0.33605 | 1.24683 |
| 10 | 9.8586e-4 | 0.26382 | 0.25299 |
| 100 | 1.0953e-6 | 0.25141 | 0.07551 |
| 1,000 | 1.1065e-9 | 0.25014 | 0.02376 |

Holding the target and coupling fixed, doubling photon energy increases the initial fractional loss rate by about 8.023 at reference b=100, approaching the E^3 factor of eight. The point-target factor was 64. This is a conditional improvement in color dependence, not agreement with the desired same-fraction shift across colors.

For illustration, the large-b expression gives L about 4.74 km for a one-eV photon and a **per-event** photon-angle RMS of one arcsecond. That is neither a measured image limit nor a proposed astrophysical target. Many events need a complete angular-transport calculation. At that size the single-event energy transfers still approach 25%, rather than the tiny transfers required by the earlier independent-event narrow-line benchmark. The target must also exist, have the stipulated spatial response, and meet material, energy and local-test constraints.

## 6. Meaning for the full theory

Spatial coherence provides a calculable way to favor straighter propagation. It does not by itself fix spectral sharpness, color dependence, inverse reactions, companion identity, timing, capture or gravitational response. The completed theory must predict these together rather than borrow one favorable property from each incompatible model.

A next candidate can add a **derived dynamical response** of the material or collective state that favors small energy transfers, or evaluate the ordinary spin-two matter-assisted emission amplitude. A material's frequency response concerns its interaction dynamics; it does not automatically mean clocks run differently in voids. A stationary response also does not automatically stretch complete transient light curves: the previously derived event-timing requirement remains.

No phenomenological frequency filter is adopted as a solution here. Its coupling, internal modes, causal response, energy receiver and fluctuations would need calculation. No unseen material reservoir is introduced to supply the desired gravity. This pass remains preliminary T03/T04/T05 evidence within the full 20-task, 32-area goal.

[Results and convergence checks](collective-response-results.json) and [runner](check_collective_response.py) are available alongside this derivation.
