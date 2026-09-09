# Companion accumulation: reverse energy transfer and noise

## Why this test is necessary

The preceding emission calculations started with empty companion modes. Long-lived companions can accumulate, so the same interaction must allow photons to absorb them as well as emit more of them. Counting only emission would overestimate net redshift and ignore fluctuations from energy going back and forth.

This pass treats an isotropic, incoherent companion occupation n(Q), held fixed during a local calculation. Photon occupation is dilute, so final-state photon stimulation is neglected here. Thermal distributions are illustrative choices, not an adopted companion temperature or a Big-Bang premise. The actual source, depletion and capture history of the bath remains to be derived.

## 1. Use the same forward and inverse vertex

Retain the dispersive scalar action and constant electric/magnetic couplings from the preceding passes. Write

    omega^2=v^2 Q^2+Q^4/M^2, Q_c=M sqrt(1-v^2),
    D_Q=Q^2-omega^2,
    A=(g_E-g_B)^2, B=(g_E+g_B)^2+4g_O^2.

For E>=Q_c both emission and absorption can access the full band. An emitted companion leaves photon energy E-omega; an absorbed companion leaves E+omega. Define the bare per-photon rates per momentum interval

    r_down(Q,E)=Q/(64 pi omega) {A[2(E-omega)-D_Q/(2E)]^2+B D_Q^2/(4E^2)},
    r_up(Q,E)=Q/(64 pi omega) {A[2(E+omega)-D_Q/(2E)]^2+B D_Q^2/(4E^2)}.

The physical downward rate is (1+n) r_down and the upward rate is n r_up. The factor 1+n is bosonic enhancement; n is required to absorb an already present companion. These factors are the usual quantum collision bookkeeping described in [Tong's kinetic-theory notes](https://www.damtp.cam.ac.uk/user/tong/kintheory/two.pdf). The amplitudes, normalization and kinematic band here come from the project's specified action, not from that reference.

The photon and companion energy gains are opposite in every event. A bath that heats photons loses the corresponding energy. Treating n as fixed is a local reservoir approximation, not evidence of unlimited free energy or an indefinitely stationary universe.

## 2. Net loss differs from gross emission

Define positive alpha to mean net photon energy loss:

    alpha(E)=1/E integral omega [(1+n)r_down-n r_up] dQ.

With F_plus and F_minus from the previous pass and

    J0=integral_0^Qc Q n(Q) omega dQ,
    J2=integral_0^Qc Q n(Q) D_Q omega dQ,

the result is exactly

    alpha(E)=A[F_plus(E)-J0/(4 pi)+J2/(16 pi E^2)]+B F_minus(E).

This identity assumes n is independent of the incident photon energy, the whole emission band is available and the displayed moments converge. It does not require n to be thermal.

In the B channel the upward and downward bare rates are equal at the same E. Bosonic enhancement and absorption cancel in the net loss, leaving the spontaneous E^-3 drift exactly unchanged. They do not cancel the number or spread of energy exchanges.

In the A channel the bath correction is negative for nonzero occupation: D_Q<4E^2 in the full band, so J2/(16 pi E^2)<J0/(4 pi). The bath reduces net photon energy loss and can make it negative, meaning photons gain energy. A larger gross emission rate is therefore not proof of more energy available for gravitational deposits.

## 3. The remaining color obstruction

The bath changes only the constant and E^-2 terms of the previous full-band rate. It leaves the positive leading E term when A>0. If A=0, the rate remains proportional to E^-3. Therefore fixed, finite isotropic bath moments cannot make a nonzero exactly constant fractional loss over an open full-band interval within this dilute-photon, constant-coupling model.

This extends the previous obstruction to a populated but incoherent, fixed isotropic bath. It does not establish the same numerical 8.30% flatness bound for all bath choices; that number belonged to the empty-bath calculation. Coherent fields, anisotropic distributions, evolving baths, finite photon occupation and other interactions are outside this result.

## 4. Energy noise can grow while the mean loss stays fixed

The local second energy-jump moment is

    D_E(E)=integral omega^2 [(1+n)r_down+n r_up] dQ.

Both directions add, so D_E increases with n even when their contributions to net drift cancel. This is the infinitesimal second jump moment, not an already integrated line width. A finite-distance spectrum requires the complete energy-dependent transport equation, including feedback and the relevant soft-mode treatment.

For an illustrative v=0.5, E=10Q_c and thermal n=1/(exp(omega/T)-1), the electric-only case gives:

| T/Q_c | Net fractional loss relative to empty-bath loss | Local second jump moment relative to empty bath |
| ---: | ---: | ---: |
| 0 | 1 | 1 |
| 1 | 0.662 | 4.16 |
| 10 | -3.36, net photon heating | 41.73 |

In the B-only case at the same scales, net loss is unchanged, while the second moment rises by factors about 3.69 and 35.92 at T/Q_c=1 and 10. These temperatures and parameters are not observations or proposed universe values. Numerical A/B values fix relative weights; an overall sufficiently small coupling must be restored before interpreting physical rates or claiming the perturbative validity conditions are met. The reported ratios do not depend on that common normalization.

## 5. Detailed balance and the soft-mode limitation

Per-photon rates include final photon phase space. Equilibrium must also include the photon-state density proportional to E^2. For a thermal companion bath, the dilute photon detailed-balance identity is

    E^2 exp(-E/T) (1+n) r_down(Q,E)
       = (E-omega)^2 exp[-(E-omega)/T] n r_up(Q,E-omega).

The script checks this identity directly. A common temperature does not imply one-way photon redshift. Restoring full bosonic photon factors gives Bose-Einstein equilibrium; since this interaction preserves photon number, it does not by itself force the photon chemical potential to zero or guarantee a Planck spectrum. Other photon-number-changing processes would have to be supplied.

There is also an infrared issue that an energy ledger alone would miss. As Q approaches zero in a thermal bath, n approaches T/(vQ). In the A channel the downward event-rate integrand then behaves as

    (1+n)r_down ~ A E^2 T/(16 pi v^2 Q).

The total event count has a logarithmic soft-momentum divergence. The energy-transfer and second-jump moments calculated here remain finite because they carry powers of omega and, for net transfer, forward/reverse cancellation. This does not justify an unregulated finite-rate Poisson simulation. A finite-volume/infrared cutoff or an appropriate inclusive/resummed treatment is required. The script checks the logarithmic cutoff dependence rather than hiding it. The B channel is more strongly suppressed at small Q.

## Decision and next requirements

Accumulating companions cannot be treated as a passive permanent energy ledger while the same reversible coupling remains active. In this model they can return energy, amplify fluctuations and require careful treatment of soft modes. A fixed isotropic bath also does not repair exact achromaticity of the tested full-band kernel.

The next physical response must specify whether companion modes are incoherent, coherent, directional, captured or otherwise modified, and derive that state from the source and receiver dynamics. It must retain reverse transfer and energy/momentum accounting. No companion lifetime, permanent storage, atomic exemption or new gravitational multiplier is adopted here.

## Verification scope

Checks compare direct forward-minus-reverse integrals to the analytic bath correction, verify noise enhancement and B-channel drift cancellation, test photon-state-density detailed balance, and check soft event-count scaling. The suite does not establish a complete bath evolution, an observed spectrum, supported deposits or gravitational predictions.
