# Coma shear recovered from the published figure

The preceding goal turn made progress on receiver consistency. This turn acquires an observational comparison product. All six research objectives remain open.

## Source and extraction

[Kubo et al., Coma weak lensing in SDSS](https://lss.fnal.gov/archive/2007/pub/fermilab-pub-07-453-a-cd.pdf), Figure 2 on PDF page 15, displays six tangential-shear points, one-sigma bars and a cross-component control. The paper reports zero-shear chi-square values 23.33 and 5.65, respectively. These shape measurements are distinct from its fitted NFW curve. The published horizontal axis is in h^-1 Mpc.

extract-kubo.py downloads and pins the PDF by SHA-256, extracts its vector drawing coordinates, identifies the six square centers and error-bar endpoints, and preserves the cross component. The fitted curve is excluded. The original rendered figure was visually inspected to verify axis labels, linear scaling, markers and error-bar grouping. Numeric output is in kubo-figure-data.json. The local source and rendered page are cached under research_work/generated/coma-kubo2007, outside git.

## Approximate reconstructed measurements

| Published radial coordinate (h^-1 Mpc) | Tangential shear | Plotted one-sigma error |
|---:|---:|---:|
| 1.209 | 0.007847 | 0.002534 |
| 2.753 | 0.003093 | 0.001500 |
| 4.465 | 0.003051 | 0.001191 |
| 6.187 | 0.000373 | 0.001021 |
| 7.930 | 0.001415 | 0.000898 |
| 9.657 | -0.000449 | 0.000805 |

These are figure reconstructions, not an author-supplied numeric catalog. We retain more digits in the machine-readable record for reproducibility, not as additional observational precision. A working coordinate allowance of 0.25 PDF point corresponds to about 1.0e-5 in shear; it is an extraction allowance, not a calibrated systematic-error bound.

## Independent checks

Using the recovered values and plotted errors, sum (shear/error)^2 is 23.3359 for tangential shear and 5.6453 for the cross component, agreeing with both published summaries to less than 0.1%. This is an ingestion check, not a test of our companion model. Vector error-line origins agree with independently extracted filled-square centers within 0.21 PDF point. Upper/lower error lengths also agree to this scale. The negative outer shear point and all cross-component measurements remain intact.

## How this can be used

This is an exploratory observational target, not an untouched holdout. It establishes a recoverable measured shear profile rather than requiring us to fit an inferred halo mass. It does not yet justify a final likelihood: bin boundaries and within-bin source distributions, full covariance, shape-calibration uncertainties, and source-distance geometry are not recovered here.

The radial numbers preserve the paper's original coordinate convention. They are not stipulated distances in our fictional universe, and we have not used an expansion law to convert them. Ratios of those coordinates preserve relative angular separations under the paper's single lens-distance scaling, but absolute angular bins need recovery before a physical prediction. In particular, these are representative bin positions, not a license to replace bin-averaged model predictions with point evaluations.

Next, recover the angular binning and weighting or acquire the calibrated source catalog; construct Coma's ordinary-matter and emitter model; then compare bin-averaged reduced-shear predictions with a supported deposit model. No fit parameters were optimized and no final holdouts were opened. No claim of matching Coma has been made.
