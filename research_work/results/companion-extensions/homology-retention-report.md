# Smaller reservoir with reduced retention: no useful improvement

13 September 2026. One new shared parameter, trained on the existing 89 galaxies
and frozen for the existing 29 validation and 31 test galaxies. All data exposed.

## Result

The tested modification does not merit replacing the reference. A 101-point
training scan selected s=1; polishing the endpoint interval found
s=0.99961406, effectively no change. Independent higher-resolution polishing
gives 0.99961449. The chosen change reduces stored mass and radial scale by only
0.0386%. Training mean squared log-speed loss improves by 2.69e-9, while the
speed RMSE worsens slightly in every split. This is not a substantive predictive
gain or a measurement of the new parameter.

| Variant | Training RMSE km/s | Validation RMSE km/s | Test RMSE km/s |
|---|---:|---:|---:|
| Reference, recomputed | 29.0254 | 32.4950 | 23.5910 |
| Training-selected s=0.999614 | 29.0277 | 32.4959 | 23.5945 |
| Fixed s=0.95 | 29.3419 | 32.6404 | 24.0972 |
| Fixed s=0.90 | 29.7099 | 32.8579 | 24.7245 |
| Fixed s=0.80 | 30.6113 | 33.5294 | 26.3460 |
| Fixed s=0.50 | 34.9714 | 37.7068 | 34.3402 |
| Existing fitted simple MOND | 19.8896 | 26.8761 | 16.3978 |

The actual fitting objective is equal-galaxy mean squared log10 speed, not the
km/s metric. Selected log-RMS values change from 0.1384967495, 0.1154012261,
0.0910754795 to 0.1384967398, 0.1153970099, 0.0910742608. Thus the tiny optimum
improves all three logarithmic summaries while worsening km/s summaries. A small
validation log-score improvement also occurs at s=0.95, but training and test
log scores worsen there. Neither detail justifies promoting a substantial change.
The candidate has four shared fitted constants including the three inherited
reference constants, versus one for the existing fitted MOND comparator.

## Formula and what changed

    rho_s(r) = s^-2 rho_reference(r/s)
    M_s(<r) = s M_reference(<r/s)
    v_s squared(r) = v_b squared(r) + G s M_reference(<r/s)/r

The range was fixed to 0.5<=s<=1. The one-third retention factor and original
opacity-defined profile were retained inside rho_reference; ordinary-matter
inputs were unchanged. This transforms the deposited density, rather than
re-solving photon transport in a newly derived interaction. The integration and
change-of-variable identities are established mathematics; tying size and total
retention through s is the untested project hypothesis under examination.

Total stored mass is s times the original. Relative to the original deposited
energy, fraction 1-s must remain outside this stored component. For this trial it
is labeled unretained traveling energy; no capture/escape rate or receiver law has
been derived. Contraction binding release is also not calculated here. Thus the
trial preserves an explicit bookkeeping obligation, not a complete dynamical
energy-conserving formation model. The current best fit makes that change tiny.

## Why the intuitive correction failed

For the companion component alone, the transformation implies

    v_companion,s squared(r) = v_companion,reference squared(r/s).

This shifts the companion speed curve horizontally. Reducing s lowers the speed
only where the original companion contribution is already declining. The earlier
large-radius argument applied beyond most stored mass, not automatically beyond
three stellar disk scale lengths. Stellar outer radii and reservoir outer radii
are different scales, especially for hollow, extended profiles.

At s=0.90, 758 of 1284 sampled outer points instead gain speed. The mean outer
residual increases from +10.41 to +12.75 km/s. At s=0.80 it reaches +14.62 km/s.
The selected near-identity change leaves the original inner/middle/outer pattern
essentially intact. Mean-sign diagnostics were useful, but did not guarantee
that this particular transformation would move enclosed mass in the needed way.

## Numerical verification and limits

Reference reconstruction uses direct angular attenuation and radial mass
integration, including the corrected 2 C0 normalization. At s=1 it differs from
archived speeds by at most 0.02696 km/s; split RMSE differences are below
0.00032 km/s. Both reference and candidate use the same new quadrature, so the
within-run comparison is matched. Doubling 2048 radial/96 angular nodes to
4096/192 changes any selected prediction by at most 0.00232 km/s. Both complete
grid scans and endpoint polishing confirm the near-identity optimum. Extra
digits distinguish calculations, not physical precision. No observational or
mass-model systematic uncertainty has been marginalized.

At the selected parameter, 46/89 training, 12/29 validation and 18/31 test galaxies
have smaller logarithmic losses, with negligible aggregate improvement. Their
individual curves and the full scans are saved. No galaxies were removed and
no target-specific nuisance parameters were added. The bounds were not expanded
to rescue the result. This is a numerical scan/polish result, not a theorem of a
global optimum for every possible profile family.

No lensing or stability success is claimed: those calculations were not rerun.
The proposed replacement already fails to provide a useful rotation improvement,
so it is retained as a diagnostic branch and the exact-third reference stays in
place. The paper supplement and checkpoints record this outcome.

## Next consequence

Do not continue shrinking and reducing retention in this tied proportion. A
future redistribution law must be judged by its effect on enclosed mass at the
actual observed radii, including the hollow profile, before fitting another
parameter. A more physical candidate would specify where incoming energy is
retained or released and then calculate its radial density. This result rejects
promotion of the tested modification, not every radial-shape change or the
companion concept.

Run `python research_work/results/companion-extensions/homology-retention.py`.
The protocol, output JSON, input hashes, fixed-size diagnostics, numerical
refinement and all per-radius predictions are preserved alongside this report.
