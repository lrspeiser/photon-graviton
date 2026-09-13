# Identical halo density gives identical companion gravity

Frozen numerical substitution, 13 September 2026. Six previously fitted SLACS systems; no optimization or new observations. This tests gravitational equivalence, not whether the capture mechanism produces the target distribution.

## Result

Replacing the NFW halo by the identical density interpreted as deposited companions reproduces all 40 stellar-motion bins to at most **0.000201 km/s**. Total motion chi-squared is **49.713056** for the archived analytic halo and **49.712678** for the numerically integrated companion density. The small difference is integration error, not a better fit. Independent density projection reproduces halo lensing to maximum relative difference **5.71e-14** for nonzero halos. J1204+0358 has zero fitted halo amplitude and reproduces exactly as a zero-source control; its reported zero lens difference is absolute, not a defined fractional error.

| Galaxy | Maximum motion difference (km/s) |
|---|---:|
| J0037-0942 | 0.000138 |
| J1112+0826 | 0.000056 |
| J1204+0358 | 0 |
| J1402+6321 | 0.000200 |
| J1621+3931 | 0.000133 |
| J1630+4520 | 0.000069 |

The production axisymmetric gravity solver used for flattened capture also passes a spherical NFW density control: maximum fractional error in circular speed squared falls from **0.000190** to **0.0000473** when the radial grid doubles. This checks its spherical limit, not every flattened or strongly structured density. Earlier oblate-homeoid checks remain the independent nonspherical control.

## What was fixed and what was independently calculated

All archived stellar masses, stellar gradients, orbital anisotropies, halo scale radii, halo amplitudes and distance geometry remain fixed. No lens normalization is recalibrated. The harness reuses only setup from free-nfw.py, stopping before optimizer code. Stellar motions use the same tracer, projection and covariance calculation, with the halo force replaced by direct integration of the specified density on 8193- and 16385-point log-radius grids. The same stellar and orbital assumptions are essential to this equivalence.

Known NFW mathematics, not a new project formula:

    x = r / rs
    rho(r) = A / [4 pi rs^3 x (1+x)^2]
    M(<r) = A [ln(1+x) - x/(1+x)]
    g(r) = G M(<r) / r^2

The companion calculation obtains M by numerically integrating 4 pi r^2 rho. Lensing is independently obtained by integrating the density in the projected cylinder: a shell outside impact radius b contributes the fraction 1-sqrt(1-b^2/r^2). This is compared with the archived enclosed-mass deflection integral. Newtonian force, weak-field GR bending and spherical Jeans projection are established methods already used in the project. Calling this density companion deposits is a hypothetical interpretation; it does not derive the NFW profile or establish originality.

The same weak-field gravitational coupling is imposed for both sources. Equal density would not alone guarantee equal lensing if a future companion theory introduces different stresses or a different relation between the two metric potentials. This test makes no such modification. The formal NFW tail is used for equivalence; its unbounded total mass is not accepted as a physical companion energy budget. No formation time, energy supply, cutoff, retention or orbital support has been inferred.

## Interpretation

No gravity-solver failure appears in these checks. The old companion capture profile can fail while an imposed halo-shaped companion profile succeeds at exactly the halo's level. The unresolved step is producing and supporting the required distribution with a shared capture/migration law. These inherited fits are not perfect or independent predictions: they have six fitted parameters per galaxy and several boundary solutions. Their chi-squared is not evidence for a newly successful companion formation theory.

Keep the original one-third capture law as the existing reference. Do not silently replace its predicted deposits with the target NFW distribution. A next inverse-transport diagnostic can quantify where its deposits exceed or fall short of that target, at fixed ordinary-matter inputs, before specifying a conservative migration or outer-storage rule. Different stellar/orbital calibrations between earlier models must first be matched; subtracting their separately optimized densities would confound the comparison.

Reproduce with `python research_work/results/companion-extensions/halo-equivalence.py`. Per-system predictions, both numerical resolutions, fit scores, solver controls and source hashes are in halo-equivalence-results.json. An initial diagnostic run encountered a zero-amplitude halo when computing a relative lens error; the harness now handles it explicitly as a zero-source control. No production solver was changed.
