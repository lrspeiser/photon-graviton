# Move the extra source outward, refit inner stars, predict outer stars and lensing

**Outcome:** the fixed outward redistribution worsens the aggregate inner fit, outer-motion predictions and lensing predictions across the six retained galaxies. Two galaxies improve in outer motion, so this is not a claim of uniform failure. The shared change does not provide a common repair.

## Formula and physical meaning

The earlier outer-shell calculation added mass while holding stellar inference fixed. This test instead redistributes the existing finite extra-source profile. At a fixed stellar-mass normalization, prescribe

    M_c,new(<r)=M_c,old(<r/s)
    rho_c,new(r)=rho_c,old(r/s)/s^3
    g_c,new(r)=g_c,old(r/s)/s^2.

These are known mass-conserving rescaling identities, used as a hypothetical spatial change. The source remains nonnegative and its total mass is unchanged before refitting. The old empirical profile is finite because of its retained outer cutoff. s=1 is the original; s=2 moves each source radius outward by a factor two. This is not spatial expansion of the universe and is not a derived capture or redistribution process.

Keep the observed tracer light, ordinary-matter shape, empirical coupling coefficients, conditional distances and equal-potential lens response fixed. For each scenario, refit stellar mass and constant orbit anisotropy on the **inner bins only**, using the prior four optimizer starts and parameter bounds. Consequently the final fitted normalizations can differ: the rescaling conserves extra mass at fixed stellar mass, not necessarily between the two separately fitted models.

Compute the outer velocity prediction with the existing conditional measurement covariance. Compute lensing from the very same ordinary-plus-redistributed acceleration; do not change its multiplier. These are established Jeans projection, Gaussian conditioning and weak-field lensing calculations. The profile rescaling is the additional phenomenological choice, with no originality claim.

## Actual retained measurements and result

The scale s=2 was fixed in the protocol before running; it was not selected to minimize the outer/lens discrepancies. All six galaxies are previously exposed training systems. Neither their outer bins nor lens angles enter the new stellar parameter fits, but this history is not a fresh blind test.

| Separate diagnostic | Original s=1 | Outward s=2 |
|---|---:|---:|
| Inner-bin summed chi-square | 41.490 | 49.848 |
| Outer conditional standardized residual squared sum | 49.394 | 72.468 |
| Lens-angle fractional RMS | 12.04% | 12.83% |
| Fits at orbit-anisotropy boundary | 0 | 0 |

Do not combine these statistics into a global probability. The outer scores use fitted parameters and conditional measurement errors, without parameter uncertainty. Lens RMS compares spherical predictions with published SIE summaries, without lens-summary uncertainty.

| Galaxy | Outer standardized residual, original → outward | Lens prediction, original → outward, arcsec | Catalog SIE angle |
|---|---:|---:|---:|
| J0037-0942 | +1.164 → +2.675 | 1.568 → 1.531 | 1.530 |
| J1112+0826 | +4.529 → +5.093 | 1.234 → 1.203 | 1.490 |
| J1204+0358 | -1.425 → -0.599 | 1.368 → 1.336 | 1.310 |
| J1402+6321 | +3.959 → +4.818 | 1.406 → 1.395 | 1.350 |
| J1621+3931 | -0.660 → +0.032 | 1.124 → 1.132 | 1.290 |
| J1630+4520 | +3.063 → +3.975 | 1.439 → 1.402 | 1.780 |

Positive outer residuals mean observed motions exceed predicted motions. These standardized numbers are not calibrated physical exclusion significances. J0037's lens summary improves substantially while its outer-motion prediction worsens: matching one observable is insufficient. J1204 and J1621 improve their outer predictions; J1112 and J1630 worsen in both outer motion and lensing.

## Verification and interpretation

redistribution.py applies the analytic force rescaling directly, propagates it through the stellar projection and the lens integral, and uses only successful finite optimizer results. The unchanged s=1 calculation reproduces the existing inner-fit chi-square within 0.01 and lens angles within 0.002 arcsec for every system. The inherited integration grid is used; no new high-precision convergence or full posterior claim is made. All per-galaxy masses, anisotropies, predictions and outcomes are saved.

This does not rule out deposition near galaxy edges. It rejects the idea that a common factor-two outward stretch of this empirical profile is already an adequate correction under the retained assumptions. A more physical capture/support law, better ordinary-matter structure or more general stellar orbits may change the result, but must be specified and tested without assigning arbitrary profiles from each lens discrepancy.

**Decision:** retain the baseline and failed common redistribution as separate models. Do not adopt outward rescaling from the isolated improvements. The next deposition model must predict its shape from the receiving material, incident companions and support dynamics, then supply both observables. The redshift/timing interaction, three-dimensional Milky Way, clusters and final untouched prediction remain unfinished; all six goals are open.
