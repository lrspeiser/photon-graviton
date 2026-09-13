# Does freeing the spectral-clock normalization remove the timing pattern?

Not for one shared normalization across all 35 published rows. Under the
diagonal-error approximation, a redshift-dependent aging rate still fits better
than a constant rate. However, allowing a different calibration for the distant
group changes the comparison substantially. Consistency of the nearby and
distant source clocks remains an essential physical assumption.

## Observations and provenance

The existing repository table transcribes 22 nearby and 13 distant spectral-aging
estimates from [Blondin et al. (2008)](https://arxiv.org/abs/0804.3595). Spectral
ages are template-based estimates, not direct readings from clocks inside the
supernovae. This calculation reuses those exposed estimates; it does not provide
new spectra or an independent replication of their age determination.

Only object, redshift, observed aging rate and quoted uncertainty are read.
Previously stored model-prediction columns do not enter the fit. No distance,
expansion history or universe age is supplied. The natural sample gap is split
at z=0.1, reproducing the 22/13 groups.

## Effective fit and established statistical derivation

The effective comparison is

    aging_rate(z) = a (1+z)^(-b).

Here a allows a common multiplicative calibration of the spectral clock. b=0
means no redshift-dependent event stretching; b=1 is the matched stretching
relation. This power law is a known empirical form, not a newly derived photon
interaction. If a represents intrinsic source evolution or calibration, its
origin still needs to be specified.

For fixed b, put f_i=(1+z_i)^(-b), w_i=1/sigma_i^2. Established weighted least
squares gives

    a_hat(b) = sum(w_i f_i y_i) / sum(w_i f_i^2),
    chi2_profile(b) = sum[w_i (y_i-a_hat(b) f_i)^2].

The implementation verifies the zero weighted residual derivative and searches
b between -8 and 8. All reported optima and interval endpoints are interior.
It also reproduces the original unit-normalization residual scores before
releasing a. The error model is the same diagonal approximation; shared template
covariance is not supplied by this table and has not been manufactured.

## Results with one normalization within each listed sample

| Sample | Best b | Best a | Formal delta-chi2-one b interval | No stretch chi2, free a | b=1 chi2, free a |
|---|---:|---:|---|---:|---:|
| All 35 | 0.867 | 0.969 | 0.745 to 0.992 | 81.94 | 25.72 |
| Nearby 22 alone | 3.293 | 1.002 | 0.860 to 5.761 | 22.05 | 21.09 |
| Distant 13 alone | 0.811 | 0.951 | 0.196 to 1.430 | 5.12 | 3.47 |

These intervals describe this approximate likelihood only, not calibrated
coverage or a full uncertainty budget. The nearby-only large best b is weakly
constrained because its redshift range is short; it is not evidence for a
separate propagation law. Nor should the all-sample interval endpoint just below
one be interpreted as a significant failure of b=1.

Previously fixing a=1 gave chi2=150.569 for b=0 and 26.949 for b=1. Freeing the
common calibration reduces the difference to 56.221 in chi-square. We do not
convert this difference into a detection significance: templates, selection and
cross-object correlations remain outside the calculation.

## What a group-specific source explanation would require

Under b=0, the fitted nearby normalization is 0.95687, and the distant one is
0.71460. Thus the distant group would need an aging-rate scale about 25.3 percent
lower relative to the nearby group. Giving the groups those two independent
normalizations reduces total chi2 to 27.17. This is an additional fitted freedom,
not evidence that the required source/calibration difference exists.

The two-group result also shows why the distant sample alone does not provide
the same comparison as the anchored combination: an overall lower distant
aging rate can be absorbed into its free a. Nearby source information is doing
real work in the joint constraint. A population or calibration model must
predict that relationship rather than choose it after seeing the distant rates.

## Consequence for the photon-companion hypothesis

A model that only reduces photon energy and leaves event spacing unchanged still
has to explain these aging-rate observations. A single universal clock-calibration
rescaling does not resolve the approximate combined mismatch. A source-evolution
alternative is mathematically possible in this fit but needs independent evidence.

Conversely, choosing b=1 does not derive the physical cause of event stretching.
Our conversion mechanism must predict it along with redshift and brightness,
using a consistent source and detector clock convention.

The ongoing artificial-light-curve calibration answers a different question:
whether the timing estimator recovers known input stretch under its simulation
assumptions. It remains unchanged and is not replaced by this spectral-table
analysis. No reserved observations were opened. Reproduce with `run.py`; numerical
results and hashes are recorded in `results.json`.
