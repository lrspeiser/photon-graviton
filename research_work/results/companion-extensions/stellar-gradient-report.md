# A bounded stellar mass-to-light gradient at the lens-required mass

13 September 2026. Conditional stellar-model diagnostic with reference companions fixed.

## Outcome

Allowing a radial stellar mass-to-light gradient improves the exact-lens-constrained motion fit, but does not reconcile it. The six-system inner chi-squared decreases from 611.48/613.61 to 471.64/476.53 for Chabrier/Salpeter proxies, still far above 45.81/45.27 for the original unconstrained inner-motion fits. Conditional outer residual sums remain 129.03/127.53, versus 55.01/53.42 for those free fits.

Two systems per population require the imposed maximum central/outer mass-to-light contrast of ten, and two hit the existing anisotropy upper bound. These are fitted boundary values, not measured stellar properties. A follow-up at contrasts 100, 1000 and the limiting concentrated profile gives only modest further improvement in the two gradient-bound cases. No stellar-gradient prescription is adopted into the theory.

## Hypothesis and mass construction

Let nu(r) be the existing spherical deprojection of the observed light. Use the exploratory stellar mass-to-light profile

\[
\Upsilon(r)=\Upsilon_{\rm out}\left[1+\frac{h}{1+(r/R_e)^2}\right],\qquad -0.8\leq h\leq9.
\]

The measured projected half-light radius Re fixes the gradient scale; no scale is fitted. The central/outer contrast is 1+h, ranging from 0.2 to 10. This smooth bounded form is a modeling choice, not a new law or a spectroscopically established population gradient. Positive stellar density is maintained throughout this range. The light-tracer profile is unchanged; its mass weighting is what varies.

For H(r)=[1+(r/Re)^2]^-1, define the luminosity-weighted mean Hbar, the normalized original cumulative mass fraction F0, and the normalized nu H cumulative fraction FH. The mass fraction becomes

\[
F_h(r)=\frac{F_0(r)+h\overline H F_H(r)}{1+h\overline H}
=(1-w)F_0(r)+wF_H(r),\qquad w=\frac{h\overline H}{1+h\overline H}.
\]

For negative h, w is negative; this is an algebraic basis coefficient, not negative physical mass. The full density is proportional to nu(1+hH)>0. The h=0 limit exactly returns the prior constant-M/L model. At h tending to infinity, Fh tends to FH; this limits the concentration achievable with this particular fixed-scale family.

These are known linear mass/force constructions applied to a chosen nuisance profile. The companion density remains the SPARC-amplitude-adjusted exact-third reference, with A=2 C0 and all capture parameters frozen. No photon-transfer or retention formula is modified by this experiment.

## Lens constraint and fitting

At each gradient, the stellar mass is fixed by the catalogue lens angle, using the same effective spherical lens equation. Stellar deflection per unit mass is (1-w)d0+w dH, so the lens-required mass is determined algebraically after subtracting companion bending. The angle is consumed calibration, not a prediction.

Only h and constant orbital beta are then optimized on the inner stellar Vrms bins and released covariance. Beta retains the previous -2..0.45 bounds. The outermost bin is held out of this nuisance fit and evaluated with the inherited covariance-conditioned residual. Compared with constant M/L, this introduces one extra nuisance degree of freedom per system. There are three underlying stellar parameters (normalization, h, beta), with one fixed by the consumed lens constraint. It is not a parameter-free improvement.

Both population proxies, all six systems, the inherited optical geometry and luminosity mapping remain fixed. This analysis is on previously exposed systems, not independent validation or a full lens-image likelihood.

## Per-system fits

The following values are Chabrier results; Salpeter results show the same boundary pattern and are fully recorded.

| System | Central/outer M/L | beta | Inner chi-squared | Bound reached |
|---|---:|---:|---:|---|
| J0037-0942 | 10.000 | -0.0408 | 151.16 | Gradient maximum |
| J1112+0826 | 0.666 | -0.0515 | 13.01 | None |
| J1204+0358 | 10.000 | -0.7611 | 33.82 | Gradient maximum |
| J1402+6321 | 0.513 | 0.4500 | 267.43 | Anisotropy maximum |
| J1621+3931 | 0.368 | 0.4500 | 3.42 | Anisotropy maximum |
| J1630+4520 | 0.623 | 0.3124 | 2.80 | None |

The two large central enhancements are not a common trend across all lenses: other systems select lower central M/L. This is evidence about the flexibility demanded by this conditional fit, not proof of real population variations. J1402+6321 remains particularly poorly fitted within the family.

## Checking whether the gradient bound caused the failure

After observing h=9 boundaries, hold the same gradient scale and extend only those cases to h=99, h=999, and the exact pure-nu-H mass-profile limit. Refit beta directly for each. In the infinite-contrast limit, J0037-0942 still has inner chi-squared 145.34/149.14 (Chabrier/Salpeter), compared with 151.16/154.86 at contrast ten. J1204+0358 approaches 25.89/26.83, compared with 33.82/34.84 at contrast ten.

These endpoint checks show that simply increasing this contrast far beyond ten does not restore the earlier good fit in those cases. They are not a proof covering every gradient shape, scale, orbital model or intermediate contrast. The extreme profiles are diagnostic limits, not astrophysically supported stellar populations. The two beta-bound cases require separate orbital assessment; this experiment does not establish whether broader orbital freedom is physical or sufficient.

## Interpretation

The result narrows a specific explanation for the lensing/motion tension. A fixed-Re mass-to-light gradient can absorb part of it, but substantial residuals and boundary dependence remain. This does not establish that the companion law is the sole problem: the assumed distance mapping, three-dimensional light/mass structure, orbital distribution and image-model lens angle are still conditional.

Future changes should be connected to independent population or kinematic information. Fitting progressively more flexible per-lens gradients solely to force agreement would not derive the proposed gravity mechanism. Positive densities and Jeans second moments also do not guarantee a nonnegative distribution function or dynamical stability; neither is shown here. Source supply, propagation timing, capture states and joint wider tests remain unresolved.

## Verification and reproduction

Run `python research_work/results/companion-extensions/stellar-gradient.py`; see the [pre-execution protocol](stellar-gradient-protocol.md). Fifteen starts explore h and beta using interpolated orbital coefficients, then the best three are polished using direct coefficients. Final scores use direct calculations. The largest interpolated/direct score difference at those optima is 5.36e-6. The h=0 reconstruction recovers the previous exact-angle inner chi-squared to within 1.45e-12, and its lens-required mass to the checked tolerance. All fitted masses and predicted squared speeds remain positive and within the stated numerical mass bounds; lens constraints are satisfied to 1e-10 fractional tolerance.

These checks establish implementation consistency on the inherited numerical grid, not a refined error bound for every new gradient or physical validation. The result records input hashes, all fitted gradients, boundary flags, exact/constant controls, outer residuals and the large-gradient follow-up. No existing reference result is overwritten or reclassified as successful.
