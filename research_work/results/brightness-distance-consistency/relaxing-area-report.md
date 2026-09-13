# A relaxing optical response does not improve the joint comparison

The one-parameter alternative slightly improves the aggregate farther-supernova score but worsens the highest-redshift residual and the six-system lens error. It is not adopted as a joint improvement. No new lens parameter was fitted. All six research goals remain open.

## Formula and status

The new **postulated optical response**, using known exponential mathematics, is

    f = z/(1+z)
    B = 1 + f exp(-q f)
    D_G = D sqrt(B), D_A = D_G/(1+z)

It replaces the retained postulate B=1+f/(1+q f). These are empirical response functions, not first-principles derivations or established photon-to-companion interactions. The known conditional Jacobi/reciprocity relations propagate each response to lens distances. Both have B(0)=1 and B'(0)=1, the necessary regular-observer conditions under the previously declared metric branch.

The training fit gives **q=1.7952170991**, with alpha=0.0002488993286382367/Mpc and event-stretch exponent b=1 held fixed. The latter remains an assumption, not a derived timing mechanism. Maximum beam-area ratio is 1.204922, with ratio 1.203772 at z=1. Calibration uses the same 77 rows; q uses only the 466 nearer training rows. It is then frozen for 494 farther rows and all six lenses.

## Actual-versus-predicted comparison

Lower scores mean smaller residuals using the same covariance. These are exposed standardized-data comparisons, not blind significance claims or independent confirmation of the mechanism.

| Quantity | Retained rational response | Relaxing response |
|---|---:|---:|
| Training brightness covariance score | 444.543837 | 444.517950 |
| Farther brightness covariance score | 389.177267 | 388.100190 |
| Mean observed-minus-predicted magnitude, z>=1 | +0.115906 | +0.155375 |
| Lens fractional RMS, attenuated capture | 13.11169% | 13.20168% |
| Outer stellar conditional residual-square sum | 63.03850 | 63.00835 |

The no-area farther brightness score is 424.481502. Both area candidates improve that aggregate. The additional gain from relaxing the response is small and does not cure the highest-redshift discrepancy. A positive magnitude residual means these supernovae are dimmer than predicted.

| Lens system | Catalog lens angle (arcsec) | Frozen relaxing prediction (arcsec) |
|---|---:|---:|
| J0037-0942 | 1.53 | 1.799676 |
| J1112+0826 | 1.49 | 1.456965 |
| J1204+0358 | 1.31 | 1.516108 |
| J1402+6321 | 1.35 | 1.630968 |
| J1621+3931 | 1.29 | 1.370812 |
| J1630+4520 | 1.78 | 1.723219 |

Catalog angles are published lens-model summaries, not error-free direct deflections. No angle uncertainty is invented; the aggregate here is descriptive fractional RMS, not a lens likelihood. Ordinary stellar mass and orbital anisotropy are fitted to inner stellar bins only; all shared SPARC capture parameters remain frozen. No stellar anisotropy bound is reached. The transparent-capture control has lens RMS 13.91568%, with inner score 41.69793 and outer score 48.08813; the attenuated branch has inner score 48.90883.

## Verification and limits

The calibration objective is not unimodal over the full allowed q interval. An initial single bounded minimizer landed near q=100 and would incorrectly favor endpoint q=0. Before lens transfer, a 401-point training-only scan with local minimization around interior minima located q=1.795217. Candidate scores are saved in the result. No farther observations were used to choose the numerical minimum. The initial endpoint output was an optimizer failure to find the better training fit, not a separate scientific fit accepted for transfer.

The independent Jacobi initial-value integration agrees with the lens-distance integral within 9.1e-12 in relative ratio. Existing input hashes and covariance checks pass. This verifies these calculations, not a global spacetime solution.

Published supernova standardization and covariance assumptions, the conditional distance interpretation, early-type equivalent-disk mapping, spherical lens calculation, and unresolved physical source/support of deposits all remain. Detailed supernova predictions are in `relaxing-area-predictions.json`; stellar observations, predictions, conditional residuals and lens geometry are in `../isotropic-galaxy-transfer/relaxing-optics-results.json`.

**Decision:** retain the previous rational response as the working comparison. This alternative demonstrates that changing brightness geometry in this way is insufficient to reconcile the lens systems. Further formula changes need a stated physical reason and the same frozen cross-observable comparison; these results do not justify adding object-specific lens corrections.

Reproduce:

    python research_work/results/brightness-distance-consistency/bounded-area.py --relaxing
    python research_work/results/isotropic-galaxy-transfer/lensing.py --relaxing-optics
    python research_work/results/isotropic-galaxy-transfer/regular-optics-check.py --relaxing
