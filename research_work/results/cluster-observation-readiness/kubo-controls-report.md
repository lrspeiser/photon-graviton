# Coma source archive and blank-field control audit

The preceding goal turn made progress by recovering Coma's shear figure. This continuation audits the author archive and incorporates the published control measurement. All six objectives remain open.

The [arXiv source archive](https://arxiv.org/src/0709.0506) contains ms.tex and three EPS figures. We inspected its complete member list and the manuscript's binning discussion. No separate data table, analysis pipeline, full covariance, or exact intermediate angular bin boundaries were recovered. This describes this archive, not all possible public or author-held data. kubo-controls-results.json records archive and member hashes; cached source files remain outside git. The archive is read without extracting executable files.

Figure 3 of the [source paper](https://lss.fnal.gov/archive/2007/pub/fermilab-pub-07-453-a-cd.pdf) shows the mean control signal from six blank regions, using matching annuli. We visually inspected PDF page 16 and reconstructed its six tangential and cross-component means and error bars. The original six individual field profiles are not present in the recovered archive, so their covariance cannot be reconstructed from the mean alone.

The reconstructed blank-field null chi-square values are 1.8223 and 12.5664, versus the paper's 1.84 and 12.58. Agreement is within 1%, consistent with the limited plot precision. An initial attempt incorrectly reused Figure 2's vector-item ordering; the consistency check failed. Inspection located Figure 3's zero baseline before its triangles. The final code detects that baseline geometrically and verifies positive error bars, shared center coordinates and axis scales. No failed extraction was accepted or used in a result.

## Conditional offset sensitivity

To test a specific nuisance assumption, subtract each blank-field mean from the corresponding Coma bin and propagate both uncertainties:

    g_adjusted,i = g_Coma,i - g_blank,i,
    variance_adjusted,i = variance_Coma,i + variance_blank,i.

These are known mean-subtraction and independent-error propagation formulas. The assumptions are a common additive offset by bin, independent Coma/control estimates, and diagonal errors between bins. This is an optional sensitivity analysis, not a demonstrated correction or complete covariance model. Matching is by bin identity; the different plotted representative radii are not treated as identical individual source positions.

The resulting zero-shear chi-square is 18.3104 for six degrees of freedom, giving a conditional tail probability 0.00550. Thus this particular subtraction weakens but does not erase the departure from zero. It is not evidence favoring companions over another gravitational source, and the probability is conditional on the error assumptions. No parameters of our gravity model were fitted. We do not subtract the cross component from tangential shear or treat it as proof that all systematics vanish.

## Remaining requirements

The actual intermediate angular bin edges, within-bin source distribution, source geometry, full covariance and calibration uncertainties still need acquisition or reconstruction from a source catalog. Until then, the recovered figures support exploratory checks, not the final joint rotation/lensing likelihood. A study of covariance sensitivity cannot manufacture the missing measured covariance. The existing Coma spectroscopy acquisition remains available for building the separate ordinary-matter and dynamical model. Final observational holdouts remain unopened.
