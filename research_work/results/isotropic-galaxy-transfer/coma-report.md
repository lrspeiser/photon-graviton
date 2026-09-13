# Coma shear shapes transfer, but do not distinguish interception

The previous goal turn made progress by deriving capture energy requirements. This calculation compares two specified deposition shapes with actual reconstructed Coma shear points, fitting three inner bins and predicting three outer bins. It does not perform an absolute-mass cluster fit or claim a complete cluster solution.

## What is fitted

The input is the six-point vector-figure reconstruction of [Kubo et al.](https://arxiv.org/abs/0709.0506), with the plotted uncertainties preserved. Use only ratios of plotted radii; no h-dependent distances or NFW mass estimates are adopted. Unknown source geometry is absorbed into a positive amplitude. Unknown physical cluster scale is fitted as a radial ratio. These choices remove any interpretation of the fit as an absolute capture-exposure prediction.

The same postulated opacity shape produces rho(x) proportional to (1+x^2)^-2 J(x). Two optical-depth cases were declared before results: transparent J=1 and strong interception with central chord depth T=100. This is not a transfer of SPARC's physical k0 to Coma. Known spherical projection gives

    Sigma(R)=2 integral_0^infinity rho(sqrt(R^2+z^2)) dz,
    DeltaSigma(R)=2/R^2 integral_0^R Sigma(u)u du-Sigma(R),
    gamma_t(R)=B DeltaSigma(R/a).

The arbitrary unit normalization of DeltaSigma is absorbed into B. The projection and shear relation are known lensing mathematics; capture is a postulate, not a new fundamental derivation. Ordinary stars/gas are not separately subtracted, so this asks whether a capture-like shape could dominate the total observed trend. It cannot determine how much of that trend is actually due to deposits. Reduced-shear corrections, exact annular averaging and full covariance are unavailable here.

## Frozen outer-bin comparison

| Shape | Inner diagonal residual-square sum | Outer diagonal residual-square sum |
|---|---:|---:|
| Transparent capture | 1.162 | 2.702 |
| Strong interception | 1.308 | 2.502 |
| Point-mass reference | 5.037 | 2.326 |

Each capture shape fits an amplitude and radius using bins 1-3 only. The point-mass reference fits one amplitude. There are only three outer bins, and these are exposed data. Scores use plotted diagonal errors rather than a full likelihood; no significance or confidence claim is warranted. In particular the 0.20 outer-score difference between capture shapes is not evidence of a mechanism.

| Outer bin | Transparent predicted shear | Strong-interception predicted shear | Point-mass predicted shear |
|---|---:|---:|---:|
| 4 | 0.001209 | 0.001107 | 0.000357 |
| 5 | 0.000783 | 0.000707 | 0.000217 |
| 6 | 0.000549 | 0.000491 | 0.000146 |

The capture scales are 0.80770 and 0.18271 times the first plotted radius, respectively, with neither at a declared bound. These are shape-fit scales, not physical cluster radii. The amplitudes are likewise not masses or source energies. All six observed shears, including negative values, errors and model predictions are retained in coma-results.json.

Both capture shapes match the inner trend better than the point-mass reference, which has fewer parameters. Their outer predictions remain weakly distinguished by the noisy data, and the reference has a slightly smaller outer score. A low residual is therefore not identification of companions or strong interception.

## Consequence for the full objective

The shape is not immediately incompatible with this limited Coma profile, but no absolute source normalization, ordinary-matter decomposition, velocity prediction, substructure or collision behavior has been tested. Those requirements cannot be replaced by this radial shape comparison. A future physically normalized Coma model needs ordinary mass and source geometry and must preserve the observed shear rather than import a fitted dark-halo mass as the target. All six goals remain open.

Reproduce with `python research_work/results/isotropic-galaxy-transfer/coma.py`. An initial line-of-sight quadrature failed the exact transparent-profile check near its center and was replaced before reporting results by z=sqrt(1+R^2)tan(theta), which resolves the core. The final transparent projection agrees with the analytic Sigma=pi/[2(1+R^2)^(3/2)] within the declared 0.1% check. No observational values or fit choices were changed by that numerical correction.
