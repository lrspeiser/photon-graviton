# NFW shape comparisons at fixed added-component lens bending

13 September 2026. Controlled mass-distribution comparison, not a free best-fitting dark-halo analysis.

## Outcome

The retained companion profile gives lower aggregate motion chi-squared than all three tested NFW shape controls. Companion totals are 113.56/111.27 for Chabrier/Salpeter; NFW totals range from 162.62 to 194.57 across the three fixed scale choices. All versions use the same four stellar fitting freedoms, data, necessary orbital conditions and lens calibration.

This is a restricted comparison. Each NFW amplitude is set to reproduce the companion component's bending at the catalogue lens radius, and its scale radius is tied to the companion capture scale by one of three prescribed ratios. These constraints are not a general dark-matter fit or a cosmological halo prediction. The best aggregate NFW result is at the largest tested ratio, so the scan does not establish its optimal scale. No claim of superiority over dark matter as a framework follows.

The result shows that the distribution of the additional density matters even when its lensing contribution at one radius is matched. It does not establish a photon source for that density. The unresolved companion-model motion residuals remain, particularly in J1402+6321.

## Profile and mathematical provenance

The comparison uses the established Navarro-Frenk-White profile:

    rho_NFW(r) = rho_s/[x(1+x)^2],  x=r/rs,
    M_NFW(<r) = 4*pi*rho_s*rs^3 [ln(1+x)-x/(1+x)].

See [Navarro, Frenk and White, A Universal Density Profile from Hierarchical Clustering](https://arxiv.org/abs/astro-ph/9611107). The profile and its density integral are known literature/mathematics, not part of our proposed photon-conversion mechanism. We borrow the radial shape as a control without assuming its cosmological formation scenario as a fact of the hypothetical universe.

Choose rs/a_capture=0.3, 1 or 3. At the catalogue lens impact parameter b, set the positive mass amplitude B=4*pi*rho_s*rs^3 by

    deflection_NFW(b; B, rs) = deflection_companion(b).

The deflection is calculated using the same spherical weak-field mass integral as for the companion profile. It is linear in B, so this matching introduces no optimized amplitude parameter. This is an inherited-normalization experiment: the target additional bending comes from the companion model, not a separate observation of dark matter.

With equal added bending at b, the stellar mass required by the catalogue lens angle is the same for a given stellar gradient in both profiles. Away from that radius, their enclosed three-dimensional mass and bending need not agree. Stellar motions therefore distinguish their distributions even though one bending value is identical. This does not imply that matching one catalogue angle matches all observed lens images.

The ideal NFW profile is untruncated for this shape diagnostic. Its total mass diverges logarithmically, while its finite-radius acceleration, relative potential and deflection remain finite. It is not treated as a finite companion inventory, a prediction of a universe's size, or an energy-conserving formation model. A finite physical halo implementation would require an outer boundary and provenance.

## Data, freedoms and constraints

Each scale run uses the same six previously exposed systems and forty motion bins under each of two population-proxy assumptions. All bins are fitted using the full covariance. The catalogue lens angle calibrates stellar mass. Both populations describe the same observed systems; they are not independent samples.

The four fitted stellar parameters are the mass-to-light gradient h, beta0, beta_infinity and orbital transition radius ra. Their bounds and the fixed stellar-gradient scale Re are unchanged. Enforce the same central cusp condition and sampled gamma>=2 beta condition as the preceding reports, with an 8193-point final check. These are necessary conditions for the specified orbital construction, not positive distribution-function or stability proofs.

No additional halo parameter is optimized on motion data within a single run. Testing three scale choices is nevertheless additional model selection, and the inherited companion normalization contains its own prior fitted parameters. The comparison is not presented as an equal-total-complexity likelihood ranking.

## Aggregate motion scores

Lower chi-squared means a better descriptive fit under these stipulated choices.

| Extra component | Chabrier total | Salpeter total |
|---|---:|---:|
| Stellar-only control | 175.088 | 175.088 |
| Reference companions | 113.565 | 111.274 |
| NFW rs/a=0.3 | 193.599 | 194.574 |
| NFW rs/a=1 | 177.024 | 177.067 |
| NFW rs/a=3 | 163.227 | 162.625 |

| Galaxy | Population | Companion | NFW 0.3 | NFW 1 | NFW 3 |
|---|---|---:|---:|---:|---:|
| J0037-0942 | Chabrier | 10.444 | 39.171 | 32.735 | 27.072 |
| J0037-0942 | Salpeter | 9.739 | 39.554 | 32.752 | 26.803 |
| J1112+0826 | Chabrier | 18.668 | 18.579 | 18.732 | 19.607 |
| J1112+0826 | Salpeter | 18.700 | 18.602 | 18.762 | 19.677 |
| J1204+0358 | Chabrier | 14.333 | 17.355 | 17.144 | 17.120 |
| J1204+0358 | Salpeter | 14.424 | 17.559 | 17.355 | 17.337 |
| J1402+6321 | Chabrier | 64.974 | 113.768 | 103.626 | 94.461 |
| J1402+6321 | Salpeter | 63.243 | 114.133 | 103.408 | 93.829 |
| J1621+3931 | Chabrier | 3.409 | 3.008 | 3.051 | 3.096 |
| J1621+3931 | Salpeter | 3.426 | 3.005 | 3.049 | 3.097 |
| J1630+4520 | Chabrier | 1.737 | 1.718 | 1.737 | 1.869 |
| J1630+4520 | Salpeter | 1.742 | 1.721 | 1.741 | 1.881 |

The aggregate result is not uniform improvement in every galaxy. J1621+3931 has slightly lower residuals under these NFW controls, and some other small residual differences depend on the scale. The large differences are concentrated in the difficult systems, especially J0037-0942 and J1402+6321. All per-system outcomes are retained.

At the largest tested NFW scale, J1402+6321 still has total chi-squared 94.46/93.83 and conditional outer residual 6.50/6.48, compared with companion total 64.97/63.24 and outer residual 5.47/5.40. Neither is a demonstrated successful joint description. The conditional residuals depend on covariance and assumptions; they are not model-independent exclusion levels.

## Execution and verification

Each scale/system/population run uses 28 optimizer starts, including the previous companion optimum evaluated with the substituted force. Twenty-seven or twenty-eight converge per case. Every selected objective is no worse than that feasible starting configuration. Multistart optimization does not prove a global optimum over unrestricted halos or orbit models.

The NFW enclosed-mass expression is checked against numerical integration of x/(1+x)^2 at x=1e-6, 1e-3, 1 and 1e3; the maximum relative difference is 4.70e-13. A series expansion protects the small-x subtraction. Tightening the deflection normalization quadrature changes it by at most 6.67e-16 in the executed runs; the added-component bending matches the reference to within 7.78e-16 fractionally. These are numerical consistency checks, not observational precision claims.

Explicit reference-force evaluations reproduce the recorded companion variable-radius objectives exactly at recorded precision. The total motion score decomposes into inner and squared conditional-outer residual contributions to within 7.11e-14. The final lens equation and positivity checks remain active, as do the central and refined slope conditions. None establishes full orbital distribution-function positivity or stability.

## Interpretation and next step

A simple replacement of the present added-density shape by these lens-matched NFW shapes does not improve the aggregate fit. That is useful information about where added gravity would have to reside under the current stellar and optical assumptions. It does not independently measure a dark halo or derive where companions deposit.

A genuinely broader NFW comparison would release or independently constrain its normalization and scale while retaining the same stellar freedoms and accounting for the additional parameters. The decreasing aggregate score toward the largest tested scale makes that limitation especially important. Separately, a physically derived redistribution law would need to predict the useful profile from companion transport and support, rather than selecting it after fitting motions. These are different research tasks and should not be conflated.

## Reproduction

Run gradient-orbits.py with `--slope-constrained --all-motion-bins --free-orbit-radius --nfw-shape 0.3`, then the same command with `--nfw-shape 1` and `--nfw-shape 3`. The script is in research_work/results/companion-extensions/. Each run writes its own nfw-shape-*-results.json and preserves historical outputs. Results contain fitted stars, every motion bin, reference scores, NFW scale/density/amplitude, matching and integration checks, optimizer outcomes and input hashes. Choices are recorded in nfw-shape-control-protocol.md.
