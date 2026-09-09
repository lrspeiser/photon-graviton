# Nonrotating temporal-transport candidates: first consistency and data tests

8 September 2026. Scope: an unscreened, nonexpanding fictional universe confronted with real measurements. No expansion or particulate dark matter is assumed as the proposed explanation. The tests below retain unsuccessful candidates and identify the particular assumptions responsible for their limitations.

## Finding

No tested candidate has yet met all three requirements: protected material clock ratios, shared photon/tensor propagation with observable redshift, and a stable nonexpanding background with the required distances. The calculations sharpen the requirements rather than confirm a completed theory. A universal static metric removes the observable temporal redshift; the earlier separate matter–light interaction has unacceptable leading clock behavior; and the negative curvature favored by the exploratory brightness fit cannot be supported by ordinary positive-enthalpy matter plus a canonical scalar under the Einstein equations.

This report combines exact restricted-model arguments, comparisons with published clock summary measurements, and a rerun of the existing supernova catalog/covariance calculation. It does not contain new raw clock data, a new blind supernova test, a full atomic-structure calculation, or a completed modified-gravity action. Derivative interactions without specified operators and coefficients cannot yet be assigned predictions or data-fit scores.

## 1. Clock protection: how small must the differential response be?

The fixed temporal coefficient implies a present local rate gamma = 7.7315 × 10⁻¹¹ per year in the affine-index realization n(t) = 1 + gamma t, normalized at the present epoch. Parameterize a clock-ratio response by d ln(nu_i/nu_j)/dt = q_ij gamma. The q values refer to observable frequency ratios, not to an assumed interpretation in terms of a varying fine-structure constant.

Lange et al. report the following direct secular frequency-ratio slopes [1]. The intervals below divide each reported Gaussian 95% slope interval by the fixed gamma. They are separate marginal intervals, not a joint likelihood; the two measurements need not be independent.

| Ratio | Reported fractional drift per year | Required residual response q, marginal 95% interval |
|---|---:|---:|
| Yb+ E3 / Yb+ E2, optical / optical | (−6.8 ± 7.5) × 10⁻¹⁸ | [−2.78 × 10⁻⁷, +1.02 × 10⁻⁷] |
| Yb+ E3 / Cs, optical / hyperfine | (−3.1 ± 3.4) × 10⁻¹⁷ | [−1.26 × 10⁻⁶, +4.61 × 10⁻⁷] |

Zero drift lies within one reported standard error of each measurement. This supports requiring protected ratios, but does not identify a mechanism that protects them. Compared with an order-one response to n, the differential response must be suppressed to approximately one part in a million or better, with the optical pair tighter still.

In the earlier constitutive model, epsilon_r = Z n and mu_r = n/Z. Choosing Z = 1/n gives epsilon_r = 1 and mu_r = n². For fixed masses, charges, magnetic moments and leading wavefunctions, optical binding frequencies are constant while hyperfine frequencies scale as n². The predicted optical/hyperfine slope is then −2 gamma = −1.5463 × 10⁻¹⁰ per year. Its magnitude is about 4.55 million times the quoted measurement standard error. This is a scale comparison for the leading toy model, not an exact multi-million-sigma exclusion: actual Yb/Cs magnetic, relativistic and nuclear responses have not been derived for this action.

A compensation exercise makes the required precision explicit. Multiply the leading hyperfine scaling by n^b; then q = −(2+b). The Cs comparison requires b between −2.000000461 and −1.999998737 at the stated marginal 95% level. Setting b = −2 cancels the leading drift. It does not establish a symmetry, protect every other transition, or show that the required magnetic response can arise in one consistent quantum theory. This fitted cancellation is a constraint for model construction, not a successful physical completion.

## 2. Cavity comparison: a useful measurement, but the published analysis removes our leading signal

For a fixed material cavity length and fixed atomic frequency, the ideal nondispersive cavity frequency is proportional to 1/n. Its present fractional drift relative to the atom is −gamma. Over the 2,826,942-second hydrogen-maser/cavity record described by Kennedy et al. [2], the leading expected fractional change would be −6.93 × 10⁻¹². This is a forecast for the fixed-length branch; a solid-state calculation must determine whether the spacer length and mirror phase actually remain fixed.

The published study searches for oscillating signals and explicitly subtracts the overall linear drift of the cavity or maser. Consequently its published oscillation bounds cannot be used as a bound on our leading secular drift. A synthetic projection check in the accompanying code subtracts a fitted intercept and slope from the predicted linear signal; only numerical roundoff remains. That check demonstrates the loss of identifiability, not a real-data detection or null result. The exact affine law has higher-order temporal curvature, so the statement applies to its leading linear term rather than every term for arbitrary duration.

A meaningful secular test needs undetrended frequency records and an independently constrained model of cavity aging, temperature, stress, mirror response and maser drift. Multiple material cavities and atomic species would help distinguish those responses. No such raw-record fit has been performed here, and no communication requesting data has been sent.

## 3. Explicit universal-metric action: clocks and shared waves, but no temporal redshift

Consider the concrete restricted action, in units with c0 = hbar = 1,

S = integral sqrt(−q) [M² R[q]/2 − Lambda − q^(mu nu) partial_mu(phi) partial_nu(phi)/2 − V(phi)] d⁴x + S_SM[q, Psi, A].

Here all Standard Model matter and light couple minimally to q with constant masses and dimensionless couplings. Einstein tensor waves have the same leading characteristic null cone as photons. Local atomic and material-cavity ratios are unchanged by a homogeneous lapse reparameterization; ordinary tidal and environmental corrections are not a secular effect of that lapse.

Use the homogeneous nonrotating metric

ds_q² = −dt²/n(t)² + a0² dSigma_k².

The change of time coordinate dTau = dt/n(t) makes this metric stationary with constant spatial scale. For comoving sources and receivers there is no cosmological frequency shift measured by material clocks. Expressed in t, both the received-wave frequency and atomic frequency carry the same lapse factor, which cancels from the ratio. This conclusion does not depend on the chosen n(t) history. It precedes the question of whether particular matter content solves the background equations.

As a descriptive data check, applying zero total redshift to the original 25 final-test galaxy groups gives RMSE 3,208.4 km/s, versus 412.4 km/s for the original temporal formula. These previously inspected objects are not a new blind test. The zero-redshift number includes no peculiar-velocity model and is not a formal likelihood for arbitrary moving galaxies; it illustrates that removing the systematic redshift term does not reproduce the existing catalog relationship.

Giving q an evolving spatial scale b(Tau) instead yields the usual comoving frequency ratio b_o/b_e. That restores redshift but makes source separations evolve relative to material rulers coupled to q. Calling a different auxiliary metric static does not meet the intended nonexpansion requirement. A successful alternative therefore needs a physical distinction beyond a common lapse, while protecting clock ratios through a demonstrated interaction.

## 4. Brightness and angular distance derived together at the kinematic level

For the separate fixed-material-clock branch, assume a constant-curvature material space, isotropic rays, conserved photon number, no absorption, and both photon-energy reduction and pulse stretching by S = exp(kappa R). If the geometrical beam area at radius R is 4 pi S_k(R)², then

F = L / [4 pi S_k(R)² S²],

D_L = S S_k(R), and D_A = S_k(R).

The angular-distance expression uses transverse physical source lengths measured with fixed material rulers. Thus D_L = (1+z) D_A in this specified transport construction. It is not an assumed universal distance-duality theorem and does not apply automatically to the universal-metric action above. For the ideal isotropic calculation, z denotes the transport redshift; the catalog fit retains the existing separate heliocentric energy/time factor and corrected distance redshift convention.

For negative curvature, S_k(R) = Rc sinh(R/Rc); for positive curvature it is Rc sin(R/Rc). The latter fit is restricted before the first conjugate point. These paired distance laws are derived from the same kinematic assumptions, but still lack a common successful matter–gravity action.

Rerunning the same 960 Pantheon+ supernovae with the released full statistical/systematic covariance, Cepheid-only calibration and original kappa held fixed gives [3]:

| Geometry | Chi-squared | Fitted curvature radius |
|---|---:|---:|
| Flat | 898.475 | Infinite |
| Positive curvature | 898.475 | Flat boundary |
| Negative curvature | 873.041 | 4.547 Gpc |

The negative-curvature radius has a conditional approximate 95% profile interval of 3.847–5.834 Gpc, using Delta chi-squared = 3.841 with the fixed calibration treatment, covariance and functional form. The original four-bin flat-model residual trend has chi-squared 36.743 for three degrees of freedom, corresponding to p = 5.21 × 10⁻⁸ under that specified Gaussian bin model. This is a conditional diagnostic of the trend, not a model-independent discovery significance. The full flat-model chi-squared remains below 960, so it should not be described as a catastrophic global goodness-of-fit rejection.

The fixed expanding reference retained from the previous calculation gives chi-squared 838.483. It is a comparison curve, not a newly optimized cosmology. The curvature fit reuses data already inspected; it is exploratory rather than independent confirmation. Published standardization, selection, host, dust and distance-calibration assumptions remain present.

The paired negative-curvature distance predictions, in the ideal common redshift frame, include D_A = 2.919 Gpc and D_L = 5.839 Gpc at z = 1, and D_A = 5.055 Gpc and D_L = 15.164 Gpc at z = 2. The complete small prediction table is included in the archive. No new angular-size dataset has been fitted, and a model-dependent BAO or CMB ruler has not been imported as if it were a direct distance measurement.

## 5. Does the supernova-favored curvature have a supported static background?

In Einstein gravity, a homogeneous constant-scale background obeys

rho_energy + p = k c⁴ / (4 pi G a0²),

where k = +1, 0, −1 and a0 is the curvature radius for nonzero k. At the fitted negative-curvature radius, this requires

rho_energy + p = −4.89 × 10⁻¹⁰ joules per cubic metre.

Equivalently, dividing by c² gives −5.44 × 10⁻²⁷ kg/m³. This is a required total effective enthalpy, not a measurement of a new material substance. Ordinary matter and radiation have positive rho+p; a homogeneous canonical scalar adds positive kinetic energy to rho+p, while a scalar potential or cosmological constant contributes zero to that sum. Those ingredients therefore cannot support this negative-curvature static solution, whatever potential is chosen. Flat static space has the corresponding zero-enthalpy obstruction for these positive components.

Positive curvature avoids this particular sign obstruction but must still pass stability and distance tests. Classical Einstein-static stability is mode-dependent; acceptable behavior for some perturbations does not establish full stability [4]. No compatible positive-curvature stable solution is demonstrated by the present fit. A scalar potential alone is not a demonstrated repair.

The sign argument applies to Einstein gravity with the stated stress-energy. A modified gravitational action, derivative coupling or additional interacting field changes the effective equations. It may supply the required effective terms without introducing particulate dark matter, but the resulting modes, energy exchange and photon/tensor characteristics must be derived and tested. This report does not claim an impossibility theorem for all such theories.

## Decision and next concrete calculation

The universal lapse-only construction is retained as a mathematically informative branch that fails observable redshift. The earlier fixed-atom constitutive construction is retained as a branch whose leading optical/hyperfine prediction misses the data by an enormous scale. Exact compensation remains a requirement without a derived protective symmetry. The negative-curvature brightness fit remains a useful phenomenological target, now accompanied by angular-distance predictions and a quantitative background stress requirement.

The next viable action must simultaneously provide a physical difference between propagation and material standards, protect both optical/optical and optical/hyperfine responses at the measured levels, align photon and tensor characteristics over the path, and solve stable background equations. Operator choices and symmetries must be specified before assigning these properties. Arbitrarily tuning a magnetic exponent, imposing a shared wave speed and choosing a curvature radius are not substitutes for that derivation. Until that action exists, confirmation of the combined theory is unavailable; the present work establishes constraints and failures of explicit restricted versions.

## Sources and reproduction

1. Lange et al., Improved limits for violations of local position invariance from atomic clock comparisons, Physical Review Letters 126, 011102 (2021). Direct measured ratio drifts in the paper: https://arxiv.org/abs/2010.06620 .
2. Kennedy et al., Precision Metrology Meets Cosmology: Improved Constraints on Ultralight Dark Matter from Atom-Cavity Frequency Comparisons, Physical Review Letters 125, 201302 (2020). The analysis subtracts overall linear drift: https://arxiv.org/abs/2008.08773 .
3. Pantheon+SH0ES DataRelease, commit c447f0fea703fcd0fff57de5000947b5ca81286b, distance table and STAT+SYS covariance: https://github.com/PantheonPlusSH0ES/DataRelease . Brout et al.: https://arxiv.org/abs/2202.04077 . The existing nearby-group CSV is retained from the original paper supplement.
4. Barrow et al., On the Stability of the Einstein Static Universe (2003): https://arxiv.org/abs/gr-qc/0302094 .

Extract the archive and run `python nonrotating_action_test/run_tests.py` with NumPy, pandas and SciPy installed. It regenerates the original brightness outputs in the included sibling directory and writes the new results and conditional distance table. Source hashes and the code are included. Published papers are referenced, not bundled; the clock summary values used by the code are printed above. The source-data rerun, summary-data comparison, analytical restrictions and synthetic detrending demonstration are deliberately identified separately.
