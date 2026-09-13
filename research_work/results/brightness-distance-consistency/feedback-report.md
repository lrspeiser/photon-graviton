# Companion inhibition: fit nearby, predict farther

**Outcome:** this specified feedback law fits the lower-redshift brightness group but worsens the frozen farther-group prediction. It leaves the nearby-galaxy redshift results almost unchanged. Do not adopt it as a demonstrated repair.

## Hypothesis and derivation

**New phenomenological postulate for this project comparison; novelty unverified:** attach a companion-energy account to a packet and let its accumulated energy inhibit further conversion:

    E_gamma+E_c=E0, initially E_c=0
    -d ln(E_gamma)/dD = alpha0/[1+chi E_c/E_gamma], chi>=0.

The energy ledger conserves the stipulated packet total. It is not a local stress-energy or transport derivation. In particular, E_c/E_gamma refers to energy associated with that packet, not automatically a locally measurable radiation-density ratio: companions may have departed or photons may have different arrival histories. A local implementation requires that distinction to be resolved.

Define A=ln(E0/E_gamma)=ln(1+z). Then E_c/E_gamma=exp(A)-1, so ordinary integration gives

    alpha0 D = A + chi[exp(A)-1-A]
    D(z) = {ln(1+z)+chi[z-ln(1+z)]}/alpha0.

The integration and exponential/logarithmic identities are known mathematics, not new fundamental laws. For chi>=0 the expression is monotone and can be inverted at a fixed published distance. At chi=0 it recovers the previous constant fractional-loss law. At small D it preserves the original local slope alpha0.

For the brightness comparison we additionally **assume**, without deriving, conserved photon number and event stretch S_t=S_E. The same static flux law gives F=L/[4 pi D(z)^2 S_E^2]. The feedback energy equation does not itself produce that arrival map. Thus an observational pass would still leave the principal timing/interaction gap open; the calculation does not repeat the earlier mistake of assigning a clock-induced shift automatically to a companion current.

## Declared calibration and actual transfer result

Keep alpha0=0.0002488993286382367/Mpc fixed. Calibrate luminosity on the same 77 Cepheid-host rows, with the candidate attenuation at their stipulated geometric distances. Fit only chi, bounded to [0,10], on 466 standardized supernovae at 0.1<=zHD<0.3. Freeze it for 494 farther objects and the 164 nearby galaxy-distance rows. The inherited conditional zHD/zHEL convention and full released covariance remain unchanged. All data were previously exposed; this is not a blind test.

The fitted chi is **0.444683**, not at a boundary.

| Sample | Baseline chi-square | Feedback chi-square |
|---|---:|---:|
| 466 fitted lower-redshift supernovae | 468.030 | 446.192 |
| 494 frozen farther supernovae | 424.482 | 450.456 |

These point-prediction scores use the same covariance. They do not include fitted-parameter prediction uncertainty and are not a complete model-selection likelihood.

| Redshift group | Baseline GLS mean residual, mag | Feedback GLS mean residual, mag |
|---|---:|---:|
| 0.1-0.3 | +0.1249 | +0.0336 |
| 0.3-0.6 | +0.1869 | +0.0210 |
| 0.6-1.0 | +0.2048 | -0.0869 |
| 1.0-3.0 | +0.3497 | -0.1488 |

Negative residual means the model predicts too much dimming. These correlated bin summaries are descriptive, not separately fitted or independent tests. As in the preceding comparisons, the standardized B-band products, covariance, source-population corrections and frame corrections retain published assumptions. This is not raw bolometric photometry.

The galaxy RMS residuals in c times redshift units change from 457.577 to 457.525 km/s in the 104 historical training rows, 437.065 to 434.379 in the 35 validation rows, and 415.414 to 412.246 in the 25 test rows. These small changes do not rescue the farther brightness prediction and do not identify companions as the cause. All three galaxy partitions are reused.

## Why the previous removal model looked similar during fitting

The feedback brightness change relative to chi=0, aside from the small common calibrator adjustment, is

    Delta m = (5/ln10) ln{1+chi[exp(A)-1-A]/A}.

Its small-A expansion is

    Delta m = (5/ln10){chi A/2 + [chi/6-chi^2/8] A^2 + O(A^3)}.

The earlier constant photon-removal model gives Delta m=(2.5/ln10) epsilon A. Thus the two models agree to first order when epsilon=chi, despite different energy accounts and photon counts. Their similar fitted strengths (0.4812 and 0.4447) are not independent evidence for either cause. For the fitted feedback strength, the quadratic term is positive, increasing dimming farther away and contributing to the poorer transfer.

This equivalence is a conditional algebraic deduction, not a new physical law or a claim that the models are identical at all redshifts. The more distant observations are precisely where their differences matter.

**Decision:** retain the failed transfer result. Do not tune another exponent on the same farther sample and call that a prediction. Future candidates need an independently motivated change to the rate's distance/environment dependence or the source/flux model, plus a derived arrival map. The original six objectives remain open. feedback.py records per-object brightness and galaxy predictions, inversion checks, and input hashes; the protocol and prior results are preserved.
