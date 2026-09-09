# Can electric and magnetic couplings cancel the color dependence?

## Why this test is needed

The dispersive companion model improved energy-packet size and direction but retained color-dependent fractional energy loss. Light has electric and magnetic fields, so a natural next comparison lets the companion couple to both. Different terms can interfere in an individual polarization amplitude. We must calculate their combined effect before assuming cancellation works across a spectrum.

The question is not merely whether a rate can be flat at one chosen color. It must remain sufficiently uniform over different colors using the same parameters. Otherwise calibration at one spectral line conceals incorrect shifts of other lines.

This pass covers constant, local, rotationally invariant couplings linear in one scalar and quadratic in electromagnetic fields, without additional derivatives in the vertex. It retains the preceding finite companion emission band and an initially empty companion bath. It does not cover arbitrary nonlocal response, extra operators, background occupation, or polarization-selected transport. No candidate is adopted.

## 1. Specify the enlarged action

Keep the scalar kinetic/spatial terms from the dispersive-companion pass and replace the electromagnetic part by

    L_EM = (1+g_E phi) E_field^2/2 - (1+g_B phi) B^2/2
           + g_O phi E_field dot B.

The real constant parameters g_E, g_B and g_O have inverse-energy units. The last term is parity odd if phi is a scalar. The familiar phi F tilde(F) interaction is one representation of the electric-magnetic term, reviewed in the [PDG axion review](https://pdg.lbl.gov/2023/reviews/rpp2022-rev-axions.pdf). That operator comparison does not identify this companion as an axion, import an axion population, or borrow axion cosmological constraints.

With epsilon=1+g_E phi and b=1+g_B phi, the canonical displacement is D=epsilon E_field+g_O phi B. For epsilon>0 and b>0,

    H_EM = |D-g_O phi B|^2/(2 epsilon) + b B^2/2.

The scalar free energy remains positive, including (Laplacian phi)^2/(2M^2). These are local Hamiltonian conditions, not a complete stability, radiative-naturalness or matter-clock proof. The preferred frame, finite EFT validity cutoff and required hierarchy Q_c << E << Lambda_UV remain explicit assumptions for the high-energy comparisons.

The three electromagnetic structures form the local isotropic quadratic basis in E_field and B under these restrictions. Frequency-dependent coefficients or extra derivatives are outside this bounded family and require a separate action and validity analysis.

## 2. Polarization determines what can cancel

For an incoming photon of energy E, outgoing energy E'=E-omega, and deflection cosine mu, the parity-even linear-polarization amplitudes, divided by E E', are

    g_E mu - g_B,    g_E - g_B mu.

The parity-odd term supplies off-diagonal amplitudes g_O(1-mu), up to polarization-basis signs. Averaging over initial polarizations and summing over final ones gives

    mean |M|^2/(E^2 E'^2)
       = [A(1+mu)^2+B(1-mu)^2]/4,
    A=(g_E-g_B)^2 >= 0,
    B=(g_E+g_B)^2+4 g_O^2 >= 0.

The script verifies this decomposition from electric and magnetic polarization vectors for random angles and all three real couplings. Interference has been included before the polarization average. The two remaining contributions add with nonnegative weights; treating them as freely signed rates would invent an unphysical cancellation. A polarized beam needs the full matrix rather than this averaged reduction.

Equal electric and magnetic coefficients eliminate the forward channel. They do not eliminate all emission into the subluminal dispersive scalar. Adding the parity-odd term changes the nonnegative weight of the same backward-angle channel; it supplies no new unpolarized energy-shape basis in this model.

## 3. Exact rate functions over the full emission band

Retain omega^2=v^2 Q^2+Q^4/M^2 and Q_c=M sqrt(1-v^2). For E>=Q_c, all 0<Q<Q_c is available. Define D_Q=Q^2-omega^2, which is nonnegative in this band. The fractional loss rate becomes

    alpha(E)=A F_plus(E)+B F_minus(E),
    F_plus(E)=1/(64 pi E) integral_0^Qc Q [2(E-omega)-D_Q/(2E)]^2 dQ,
    F_minus(E)=1/(256 pi E^3) integral_0^Qc Q D_Q^2 dQ.

Thus the minus contribution is exactly proportional to E^-3 in the full-band regime. The plus contribution grows as E at high energy. Pure electric coupling reproduces the preceding result with A=B=g_E^2; the numerical checks verify that comparison.

For a transparent algebraic check, define

    I0=integral Q dQ, I1=integral Q omega dQ, I2=integral Q omega^2 dQ,
    IA=integral Q D_Q dQ, IAw=integral Q D_Q omega dQ,
    IA2=integral Q D_Q^2 dQ.

Then

    F_plus(E)=[E I0-2 I1+(I2-IA/2)/E+IAw/(2E^2)+IA2/(16E^3)]/(16 pi),
    F_minus(E)=IA2/(256 pi E^3).

Here I2-IA/2=Q_c^4(1+v^2)/8>0, so the coefficients of the inverse powers in F_plus are positive. Its constant term is negative, but its original integral is positive. Both basis functions are strictly convex on this full-band interval.

No nonzero constant alpha can hold on an open full-band energy interval with fixed A and B. If A>0, the analytic rate contains a positive coefficient of E which no other term can cancel identically. If A=0, only E^-3 remains, or the rate is zero. This is an exact statement about this bounded unpolarized operator family, not every possible photon-companion interaction.

## 4. How close can a mixture get across twofold photon energies?

Opposite slopes can cancel the derivative at one energy. To test a finite range, allow both the mixing ratio and an overall normalization to vary, and minimize the largest fractional departure from a constant over E in [E0,2E0]. This is a mathematical flatness diagnostic, not a fit to astronomical spectra or a measured tolerance.

At high E0/Q_c, write y=E/E0 and the normalized shapes as y and y^-3. The optimal positive mixture has equal endpoint rates:

    f(y)=y+(8/7)y^-3,
    f(1)=f(2)=15/7,
    y_min=(24/7)^(1/4)=1.3607498666...

The maximum/minimum rate ratio is 1.18107148. After optimally choosing the common constant normalization, the smallest maximum fractional error is

    (f_max-f_min)/(f_max+f_min) = 0.08301951...

So even the best high-energy mixture has about an 8.30% departure from a constant somewhere across a factor-of-two photon-energy interval. Equivalently, the rate maximum is about 18.11% above its minimum. These are two descriptions of the same shape, not two observational bounds.

For finite E0/Q_c the script uses the exact moment functions, equates endpoint rates and finds the unique interior minimum. Strict convexity places the maximum at an endpoint. The ratio F_minus/F_plus decreases with energy, so shifting the mixture away from endpoint balance worsens the controlling endpoint-to-minimum ratio. An independent dense-grid mixing optimization checks the result.

Illustrative v=0.5 results:

| E0/Q_c | Best maximum fractional departure from a constant |
| ---: | ---: |
| 1 | 14.68% |
| 10 | 8.78% |
| 100 | 8.35% |
| 10,000 | 8.3024% |

The high-energy optimum also needs a large B/A ratio because the minus basis is suppressed as E^-3. The required ratio is recorded, not assumed radiatively stable or fixed by a symmetry. Multiplying all couplings by a small number changes the total rate but not the relative color variation.

## 5. Decision and remaining routes

This calculation closes a specific tuning route: constant electric, magnetic and parity-odd scalar couplings cannot turn this full-band spontaneous kernel into nonzero exactly achromatic fractional drift. Cancellation at one color is insufficient; adding the parity-odd term does not add a third independent unpolarized spectral shape.

The result does not exclude a narrow transition-band approximation, polarization-dependent transport, a populated/coherent medium, frequency-dependent response, extra derivative operators or other fields. Those are different physical models. A next extension must specify and derive such a change, including its noise, energy destinations, matter response and validity range. Do not continue searching constant coefficients in this already bounded family as though that could give exact broadband cancellation.

No fit of gravity profiles, no-loss companion storage rule or arbitrary energy-dependent coupling repairs this source-kernel result. The complete timing, image, source supply, capture, supported deposit and common motion/lensing requirements remain open.

## Verification scope

Checks cover the polarization decomposition, direct versus moment-integrated rates, recovery of the prior electric-only result, exact inverse-cubic scaling, and independent optimization of finite-band flatness. No astronomical data are used in this pass. Small finite-distance redshift errors, line widths and image changes would still require a full propagation calculation; the rate-flatness percentages are not those observables.
