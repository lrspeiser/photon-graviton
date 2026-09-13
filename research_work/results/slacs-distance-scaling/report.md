# A common distance rescaling cannot repair these lensing predictions

Under the current empirical completion, changing the normalization of the
redshift-distance law leaves the stellar and lensing predictions unchanged when
the corresponding mass and acceleration scales are transformed consistently.
This closes one apparent way to repair the discrepancies; it does not establish
that the geometry or the physical theory is correct.

## Derivation and provenance

Our adopted effective prescriptions are

    D(z) = ln(1+z)/alpha,
    a_star = c^2 alpha,
    g_extra = A a_star^(1-p) g_b^p.

The exponential loss relation and power-law force parameterization are known
mathematical forms, used here as empirical postulates. Neither their adoption
nor this scaling argument establishes a microscopic photon interaction.

Take a positive scale s and transform

    alpha -> alpha/s,       D -> s D,
    r -> s r,               M -> s M,
    a_star -> a_star/s.

Keep measured angular profiles, seeing, velocities, redshifts and anisotropy
fixed. Keep A and p fixed. The physical source cutoff remains 20 times the
angularly inferred reference scale and therefore also grows by s.

For the same dimensionless ordinary mass profile, established Newtonian scaling
gives g_b -> g_b/s. The prescribed extra force then transforms as

    g_extra -> A (a_star/s)^(1-p) (g_b/s)^p = g_extra/s.

The stellar Jeans velocity-squared prediction contains acceleration times
length; both factors cancel. Angular aperture and seeing convolution are also
unchanged. Thus M -> s M gives the same velocity predictions at every fixed
anisotropy, not just one fitted velocity.

The weak-field deflection integral contains g times a path length and is
unchanged. All distances in D_ls/D_s scale together, so that ratio is unchanged.
Consequently the angular Einstein-radius prediction is unchanged as well.
These are conditional consequences of the stated equations using established
Jeans and lensing mathematics, not independently new physical laws.

The mass posterior transforms accordingly when the prior bounds scale with M.
Uniform mass and uniform log mass each retain their form after normalization.
Fixed hard physical mass bounds can break the invariance through the prior;
that would not supply observational evidence for a distance correction.

## Executed check

`run.py` rebuilds all six stellar component models with physical radii, seeing,
cutoff, mass and acceleration scale multiplied by the required factors. It
checks s=0.1 and 10 for both baryonic and extra-force models, at three fixed
anisotropies, using the previously inferred posterior median masses. It computes
the lens roots afresh rather than reusing the old angles. This gives 24 model
and scale cases, 72 velocity-profile comparisons and 24 lens-angle comparisons.

Declared numerical gates are relative velocity error below 1e-6 and angle
difference below 1e-5 arcseconds. `results.json` records every case and source
hash. These are symmetry checks, not new fits or observations.

All 24 cases pass. The maximum relative velocity difference is 2.00e-15;
the maximum lens-angle difference is 5.97e-10 arcseconds. These small numerical
differences support the analytic invariance at the tested scales.

## What can and cannot change the outcome

Simply changing the common alpha normalization cannot improve these lensing
residuals while retaining the current a_star=c^2 alpha relation and free dynamical
mass inference. Holding a_star fixed while varying alpha changes a separate
physical assumption and is not a pure distance rescaling.

A different *shape* of D(z) can change D_ls/D_s, so it is not covered by this
invariance. Independent distances or externally constrained stellar masses can
also break it. Those require actual constraints and consistent photometric and
X-ray conversion, not a separate distance chosen for each discrepant lens.

The earlier [force-normalization degeneracy](../theory-priority-review/report.md)
varied alpha while compensating A at fixed galaxy distances. This is a different
invariance: A stays fixed while redshift-derived distances and inferred masses
scale together. Neither permits changing independently measured distances freely.

Full joint image/stellar inference, source geometry and the companion-derived
force remain unresolved. No reserved observations were opened.
