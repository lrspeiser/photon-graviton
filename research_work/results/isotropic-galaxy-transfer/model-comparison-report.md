# Matched rotation benchmark favors the tested MOND prescription over our candidate

With the same galaxies, ordinary-matter inputs and fitting objective, the simple MOND interpolation predicts rotation better than the current exact-one-third companion model in both comparison groups, using fewer shared parameters. A restricted shared NFW model does worse, but its imposed halo-to-light mapping reaches a boundary; that is not evidence against dark matter generally. Target-specific NFW fits perform much better and are reported separately because they use extra target information.

## Common conditions

All models use the same 149 SPARC galaxies and 3150 accepted radial points, fixed published distances and inclinations, and disk/bulge mass-to-light ratios 0.5/0.7. Fits minimize the same mean of per-galaxy mean squared log10 speed residuals. Parameters are calibrated on 89 galaxies and frozen for 29 validation/31 test galaxies. These partitions have been exposed repeatedly. The companion predictions are read from the existing exact-one-third run and verified against identical radial observations.

Each body has equal weight in the reported RMS. No distance, inclination or stellar-mass nuisance refitting is allowed here. Published velocity errors are included identically in a separate diagnostic; neither that diagonal diagnostic nor these descriptive RMS values is a complete uncertainty-marginalized likelihood.

## Galaxy-to-galaxy predictions

| Model | Shared fitted parameters | Validation RMS km/s | Test RMS km/s | Validation log RMS | Test log RMS |
|---|---:|---:|---:|---:|---:|
| Ordinary matter | 0 | 58.22 | 47.77 | 0.27036 | 0.25608 |
| Companion, exact 1/3 | 3 | 32.49 | 23.59 | 0.11539 | 0.09106 |
| MOND simple, fitted a0 | 1 | **26.88** | **16.40** | **0.09543** | **0.07834** |
| MOND simple, fixed a0 | 0 | 26.24 | 16.49 | 0.09765 | 0.08733 |
| NFW with imposed shared scaling | 3 | 35.18 | 40.54 | 0.12033 | 0.10168 |

The fit selects a0=8.56334e-11 m/s^2. The fixed benchmark uses 1.2e-10. A lower km/s RMS for a fixed-a0 row does not contradict the trained value: training minimized logarithmic error, not km/s error, and the test samples were not optimized. Neither a0 value is a new determination with systematic uncertainties propagated.

The established simple MOND algebraic law is g=g_b/2+sqrt(g_b^2/4+a0*g_b). It solves mu(g/a0)*g=g_b with mu(x)=x/(1+x), verified numerically. This is a phenomenological disk approximation, not a full nonspherical field calculation or an external-field treatment. See [Hees et al., interpolation and rotation constraints](https://academic.oup.com/mnras/article/455/1/449/985618). It is not our empirical extra-force fit renamed MOND.

The known NFW profile is rho=rho_s/[x*(1+x)^2] with x=r/rs ([Navarro, Frenk and White](https://arxiv.org/abs/astro-ph/9611107)). Our additional benchmark mapping rs=s*R_disk, rho_s=rho0*X^p is explicitly our imposed three-parameter predictor. It is not the general prediction of dark matter or a cosmological halo-population model. It uses no expansion or critical-density normalization as a premise.

## NFW boundary follow-up

Training selects log10 rho0=5.07316, s=100 (upper bound) and p=0.09917. The finite shared model's numerical result is therefore conditional on its stated bounds. To distinguish an arbitrary boundary from its qualitative performance, we subsequently calculated the exact large-scale inner-NFW limit: g_extra=A*R_disk*X^p. Fitting its two shared parameters on training gives validation/test RMS 35.16/39.88 km/s and log RMS 0.11979/0.10083. This limiting version still trails the companion candidate and MOND on both comparisons. It is a post-boundary diagnostic, not a secretly enlarged prior or a generic dark-matter test.

## Dark halos with target-specific inner information

For each of the 60 comparison galaxies, fit two NFW parameters to the inner 60% of its radial points, then predict the remaining outer points. The following table scores ONLY those same outer points for every method:

| Model | Target inner speeds used? | Outer validation RMS km/s | Outer test RMS km/s | Outer validation log RMS | Outer test log RMS |
|---|---|---:|---:|---:|---:|
| Frozen companion | No | 33.23 | 29.20 | 0.11261 | 0.09417 |
| Frozen fitted-a0 MOND | No | 22.14 | 16.79 | 0.06361 | 0.06459 |
| Per-galaxy NFW inner fit | **Yes, two parameters each** | 16.14 | 20.07 | 0.06557 | 0.06003 |

The adaptive halos beat our current model on these outer predictions, but the information advantage must remain explicit. They are not a frozen galaxy-to-galaxy prediction. 41/60 adaptive fits reach at least one declared parameter bound, indicating strong shape/scale degeneracy or pressure on the bounds; no precise halo parameter inference is claimed. Outer measurements were not used in those fits. We do not rank NFW versus MOND by cherry-picking validation km/s versus test log error.

## Meaning for the project

The retention revision improves our earlier model, but **it does not yet rival the tested MOND rotation prescription**. A universal claim against dark matter would be unsupported: our restricted shared halo mapping is weaker than the flexibility normally available to individual halos, and broad cosmological evidence was not tested.

The next useful question is why the companion model misses a tighter connection between ordinary matter and the observed acceleration that the MOND benchmark captures. Any revised companion law must predict that connection from source/transport/retention physics, not insert the measured extra acceleration as its input. Lensing, source funding, retention support and the Sun/Earth/Moon interpretation remain separate open requirements.

All global optimizers converge. The small-x NFW mass term uses its analytic expansion to avoid cancellation; the MOND implicit equation and exact companion-data correspondence are checked. Model complexity is reported rather than hidden in a nominal fit score. Full velocity observations, measurement errors and predictions are in `model-comparison-predictions.json`; all summary/error and adaptive fit information is in `model-comparison-results.json`. Reproduce with `python research_work/results/isotropic-galaxy-transfer/model-comparison.py`. All six goals remain open.
