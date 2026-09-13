# Reconstructing spectral aging from actual observation dates

We recovered all 35 published epoch measurements for the 13 distant supernovae
and refitted their spectral-aging slopes. The results closely reproduce the
published rates. Accounting for uncertainty in the first spectrum is important,
but it does not reveal a missing correction in the published slope errors.

## Inputs and their limits

[Blondin et al. Table 2](https://arxiv.org/html/0804.3595v1#S3.T2) supplies Julian
observation dates and spectral-age estimates. The dates are observer times; the
ages remain template-derived. This is a reconstruction of the published table,
not a new analysis of the original pixel spectra.

The paper explicitly notes that template ages include a redshift correction
and that the historical reference objects have z<=0.05. For an interval from
one reference object, undoing that correction rescales it by at most 1.05.
This alone is not a bound on all template-matching errors: different epochs
can select different reference objects. A full rerun is needed to quantify that
effect. It is inaccurate to say the paper simply imposed the distant source's
stretch on its measured observer dates.

## Formula provenance: known Gaussian regression

For each supernova fit the spectral age y_i against its observed date t_i:

    y_i = intercept + rate * (t_i-t_0).

The intercept is free; an offset in the estimated date of maximum does not
automatically change the inferred rate. Assuming independent absolute age
errors sigma_i, ordinary weighted regression gives the fit and covariance.

The equivalent difference formulation has d_i=y_i-y_0 for i>0 and

    Cov(d_i,d_j) = delta_ij sigma_i^2 + sigma_0^2.

The common first-spectrum error appears in every pair. Generalized least
squares with this covariance yields exactly the same slope and slope variance
as the intercept fit. The code checks that identity for every object. These
are established statistical formulas, not a new propagation law.

## Representative reconstructed results

| Supernova | Observed span, days | Reconstructed aging rate | Published rate |
|---|---:|---:|---:|
| 1996bj | 10.05 | 0.5274 +/- 0.3702 | 0.527 +/- 0.369 |
| 1997ex | 30.95 | 0.7450 +/- 0.0759 | 0.745 +/- 0.076 |
| 2001go | 37.47 | 0.6512 +/- 0.0621 | 0.652 +/- 0.062 |
| 2006mk | 32.07 | 0.7507 +/- 0.0599 | 0.753 +/- 0.060 |
| 2006tk | 14.02 | 0.8352 +/- 0.1818 | 0.835 +/- 0.181 |

All 13 results and all epoch dates are retained in `results.json`. Rounded
input ages and errors prevent exact reproduction. The largest rate difference
is about 0.00234, much smaller than its quoted individual error; this comparison
does not establish a complete systematic error budget.

For SN 1997ex, the recorded dates are JD 2450815.08, 2450839.96 and
2450846.03. The spectral ages are -1.6, 17.4 and 21.2 days. We fit those ages
against date differences 0, 24.88 and 30.95 days directly, without dividing
those observed intervals by 1+z before the fit.

If the first spectral age were incorrectly treated as exact, the SN 2006tk
slope error would be about 0.0605 instead of 0.1818. The published error is
0.181, consistent with retaining the origin uncertainty. We therefore do not
use the naive fixed-origin calculation as a new criticism of the publication.

## Consequence for the hypothesis

The earlier shared-normalization sensitivity remains relevant, but simple
origin-error bookkeeping does not remove the published aging-rate pattern.
Any no-stretch explanation must supply additional source or calibration physics.
Any stretching explanation must derive the effect rather than insert it into
the timestamps being compared.

Unknown covariance among absolute spectral ages or across objects remains
outside this reconstruction. The agreement of two algebraically equivalent
fits verifies the implementation, not independence of the templates or a
complete cosmological interpretation. The actual spectral arrays, observing
conditions and historical library still need a consistent reproduction.

Run `run.py`. The source HTML is cached with a content hash; if missing it is
retrieved from the versioned paper URL. An existing recorded HTML hash must
match before a rerun proceeds. No reserved scientific observations were opened.
