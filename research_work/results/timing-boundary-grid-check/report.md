# Free-scatter timing fits on finer event grids

We refit both previously exposed difficult synthetic samples with scatter free over 0-0.6, mean width 5-100 days and exponent -2 to 3. The 641-node arrays contain genuinely recomputed midpoint likelihoods; original 321-node entries, object IDs and redshifts match exactly. Source and array hashes are recorded. This extends the earlier fixed-scatter grid check to the estimator now undergoing recalibration.

| Case | b at 321 nodes | b at 641 nodes | Scatter at 321 | Scatter at 641 |
|---|---:|---:|---:|---:|
| split_gaussian-snr5-b1-seed902 | 1.030900450 | 1.030475160 | 1.68077778e-08 | 0.0199554306 |
| split_gaussian-snr5-b0-seed903 | -0.043127943 | -0.043128004 | 0.0288319325 | 0.0293261739 |

The largest exponent change is 0.0004253. However, seed 902 changes from effectively zero intrinsic scatter to about 0.01996. Successful fine-grid starts near zero have objective about 49.18190, versus 49.17940 for the selected finite-scatter candidate: the preference is shallow (about 0.00250 in negative log likelihood). This is neither a significance claim nor proof of a nonzero intrinsic scatter. Several starts remain at inferior local solutions; all are retained.

Raw centered objective values across grids are not directly comparable because each event curve is recentered on its own maximum. Comparisons of candidate objectives within one grid are meaningful. The width-population model and Gaussian integration are known statistical constructions; this experiment introduces no new physical law.

The timing exponent is comparatively stable in these two cases, while the location of the scatter optimum is not. We therefore cannot treat an apparent zero-scatter fit as a robust physical result. The running 160-case recalibration remains unchanged and uses its declared 321-node grid. Its results will be conditional on that grid; fixed-truth likelihood-ratio sensitivity and wider grid refinement are still needed before confident interval claims. Two exposed cases do not establish population-wide convergence. All six scientific goals remain open.
