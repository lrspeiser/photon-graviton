# Stellar light constrains the freely inferred lens masses

## Completed update using the newer I-band sizes

The subsequent rerun uses the newer I-band effective radii for all 33 training lenses, retaining the old-size results. Under the same representative seeing and boundary assumptions, the ordinary-matter lens-angle RMS becomes **0.1840 arcsec** and the empirical companion value **0.2083 arcsec**. The companion model still does not outperform the ordinary-matter benchmark in this conditional comparison.

For the 32 mass-matched systems, the updated dispersion-inferred to population-mass ratios are:

| Dynamical model | IMF | Energy loss only | Energy loss plus event stretch |
|---|---|---:|---:|
| Ordinary matter | Chabrier | 3.42 | 2.84 |
| Ordinary matter | Salpeter | 1.92 | 1.60 |
| Empirical companion | Chabrier | 2.81 | 2.37 |
| Empirical companion | Salpeter | 1.60 | 1.35 |

Thus the size correction reduces the discrepancy but does not resolve it under fixed population assumptions. The earlier table below is retained as the initial old-size calculation, not the final size-consistent result. A changed angular half-light radius still does not provide a direct deprojection or a full stellar-population fit. Run lens-training-pilot/run.py --updated-profile --refined, then this directory's run.py --updated-profile to reproduce. A separately evaluated coarser updated-profile run confirms angle changes below 0.001 arcsec on refinement; the exact value is in verification.json.

The public multiband SLACS stellar-population tables are now archived. All 33 pilot training lenses have entries; 32 have published stellar masses under both Chabrier and Salpeter stellar initial-mass functions. J0109+1500 has no published population mass and remains missing. The full release covers 85 systems; no validation/test mass comparisons were scored here.

The first normalization check shows a substantial mass discrepancy. Under the branch with photon-energy loss and event stretching by 1+z, the companion pilot requires a median ordinary-matter mass **2.73 times the Chabrier estimate or 1.56 times the Salpeter estimate** after the stated normalization conversion. The corresponding ordinary-matter-only factors are 3.35 and 1.91. Thus the earlier lens-angle agreement cannot establish an admissible stellar mass model simply because its mass was freely adjusted to the measured dispersion.

This is **not a new stellar-population fit or a calibrated discrepancy significance**. The source populations' age, dust, metallicity, star-formation history and initial-mass-function assumptions are held fixed. In addition, updated I-band sizes differ from those used in the earlier dynamical pilot. These limitations must be resolved before rejecting or selecting a gravity model.

## Why the redshift mechanism changes the stellar-mass inference

Suppose N photons leave a source in source time Δt_e. Let their observed energies be reduced by 1+z, let the arrival interval be Δt_o=SΔt_e, and assume photon number is retained with isotropic geometric spreading over area 4πD². **Known energy/flux bookkeeping under these explicit propagation assumptions** gives

\[
F_{\rm bol}=\frac{L_{\rm bol}}{4\pi D^2(1+z)S},
\qquad D_L^2=D^2(1+z)S.
\]

S=1 is the stationary energy-loss-only comparison; S=1+z is the event-stretching comparison. These are hypotheses for the arrival-time behavior, not two established mechanisms. The prior timing audit found a problem with stationary loss; including it here is a sensitivity calculation, not reinstating it as observationally successful.

For a fixed emitted spectral energy distribution, the corresponding specific flux is

\[
F_{\nu,o}(\nu_o)=\frac{L_{\nu,e}((1+z)\nu_o)}{4\pi D^2S}.
\]

The photon-energy reduction and frequency-bin Jacobian cancel in this expression; integrating over observed frequency restores the bolometric equation. There is therefore no extra arbitrary redshift factor to insert separately into a fixed-template mass rescaling.

**Known luminosity normalization scaling, applied conditionally:** if the same stellar population has the same mass-to-light ratio, then

\[
\frac{M_{*,\rm conditional}}{M_{*,\rm published}}
=\frac{D^2(1+z)S}{D_{L,\rm published}^2}.
\]

We use the previously declared D=ln(1+z)/α and its frozen fitted α. The publication's flat Ω_m=0.3, Ω_Λ=0.7, H₀=70 km/s/Mpc distance is reconstructed **only to undo its luminosity normalization**, not to adopt expansion or dark matter in the fictional universe. This transformation does not remove cosmological age priors, refit stellar spectra, or establish a new posterior. The mathematical formulas above are not claimed as unique to this project; their application to the hypothesized propagation branches is conditional modeling.

In ordinary language: if the same observed brightness represents photons that lost energy and arrived more slowly, the original galaxy must have emitted more light. More intrinsic light generally means more stars for the same stellar population. That connects the timing problem directly to our ordinary-matter baseline.

## Executed training comparison

Ratios below are the dispersion-inferred ordinary-matter mass divided by the normalization-adjusted population mass. They compare 32 systems using the previous pilot's assumed 1.5-arcsec seeing and 20-effective-radius companion boundary. No mass, coupling, distance or timing parameter is fitted in this comparison.

| Dynamical model | Stellar population IMF | Energy loss only, S=1 | Energy loss plus event stretch, S=1+z |
|---|---|---:|---:|
| Ordinary matter | Chabrier | 4.15 | 3.35 |
| Ordinary matter | Salpeter | 2.35 | 1.91 |
| Empirical companion | Chabrier | 3.36 | 2.73 |
| Empirical companion | Salpeter | 1.91 | 1.56 |

These are medians of paired log-mass ratios, exponentiated. Chabrier and Salpeter are alternative assumptions about the distribution of stellar birth masses, not two independent measurements. The published mass-error columns are retained but not treated as covering our propagation, profile or population-prior changes.

## New photometry and a profile inconsistency to resolve

Among the 33 training systems, the tables provide B-band values for 11, V for 32, I for 33 and H for 22. The V filter flag distinguishes F555W from F606W; they must not be treated as the same passband in a spectral fit. Missing bands are retained as null. The downloaded tables do not provide per-band photometric error columns, so they cannot by themselves supply the complete multiband likelihood.

For the 32 mass-matched lenses, the median updated I-band effective radius is **0.829 times** the radius used in the earlier pilot; ratios range from 0.464 to 1.000. The old dynamical masses have not been silently recomputed using those new sizes. A smaller, more concentrated profile can change the mass required by an aperture dispersion. The profile and photometry must be made consistent before interpreting the remaining mass ratios as gravitational evidence.

Published physical Einstein radii, lens masses, stellar fractions and luminosities evolved to redshift zero stay in the raw source archive. They are not substituted for direct constraints in the fictional universe. The stellar-population paper uses informative priors and explicitly studies an age-limited alternative; its cosmic-time evolution is not carried over as a fact of this hypothesis.

## Next calculation and safeguards against a circular fit

First rerun the dynamical prediction with the updated angular light profile consistently, preserving the previous result. Recover photometric error/calibration information and passband responses for a stellar-population likelihood. Refit the population with explicit age, metallicity, dust, star-formation and IMF alternatives compatible with the stated fictional assumptions; do not merely rename a rescaled published posterior as independent mass evidence. Hold the propagation/arrival-time prescription common to photometry and timing.

Only then use the photometric mass constraint jointly with aperture dispersion and lens imaging. If the frozen companion model still fails, test a declared physical revision. Do not conceal the discrepancy with separate per-galaxy light-loss factors, lensing multipliers or companion amplitudes. Validation/test scores remain unopened.

## Reproduction and provenance

Run run.py followed by verify.py. The latter verifies training membership, mass normalization identities and the publication's reference distance against an independent numerical integral. There are 256 normalization scenarios from 32 systems × two dynamical models × two IMFs × two arrival-time branches. No full stellar-population fit was executed.

Primary source: [Auger et al. (2009), SLACS IX](https://arxiv.org/html/0911.2471). Electronic tables and byte descriptions: [CDS/VizieR J/ApJ/705/1099](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJ/705/1099). Retrieved 2026-09-10. Original table hashes, training-only photometry, all normalization rows and numerical checks are archived beside this report. The earlier lens-training-pilot contains the dynamical operator and its original assumptions.
