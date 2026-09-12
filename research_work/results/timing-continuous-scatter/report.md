# Continuous source-scatter integration: development revision

## Purpose and status

The frozen 160-case timing calibration has fits at its imposed scatter floor. This separate component permits a continuous normal distribution of intrinsic log widths, including exactly zero scatter. It changes no running code, checkpoint, success threshold or coverage result. It repairs a numerical representation, not the source-physics assumptions, and is not yet a calibrated replacement estimator.

## Formula provenance

These are known normal-distribution integration and interpolation identities, not new physical laws. Let x be log width, mu=a+b ln(1+z), and L_i(x) the event likelihood. Between tabulated nodes we linearly interpolate **likelihood**, not log likelihood. For sigma>0 we calculate

I_i = integral[L_i(x) N(x;mu,sigma) dx] / integral[N(x;mu,sigma) dx]

on the existing finite width support. Each interval's linear function has an analytic Gaussian integral, evaluated with normal CDF differences and first moments. At sigma=0, I_i=L_i(mu) by interpolation. The distribution thus need not snap to grid nodes when its width becomes small. At zero scatter a mean outside support has zero likelihood; positive-scatter cases with unrepresentably small support probability raise an explicit error.

The interpolation itself remains an approximation to the underlying light-curve likelihood. It can create derivative corners at zero scatter, so exact integration does not guarantee smooth optimization or adequate grid resolution.

## Executed checks

35 comparisons with independent adaptive quadrature cover sigma=0, 1e-8, 1e-4, .01, .03, .1 and .6, at interior and boundary means. Maximum absolute difference was 4.45e-16. This establishes integration accuracy for the declared interpolant only.

We then refitted a and b at six fixed scatter values on the first two already-exposed faint boundary samples, using their hash-verified saved event likelihoods. We evaluated both the 321-node grid and its 161-node subsample. No observed supernova flux enters these checks.

| Sample | Original b at floor .03 | Revised b at .03 | Revised b at zero scatter | Zero-scatter b on 161 nodes |
|---|---:|---:|---:|---:|
| Injected stretch, seed 902 | 1.03003 | 1.03000 | 1.03094 | 1.03696 |
| Injected no stretch, seed 903 | -0.04309 | -0.04304 | -0.04805 | -0.05988 |

The .03 fits agree closely with the earlier representation. At zero scatter, halving the grid changes b by about .0060 and .0118. One of three optimizer starts failed for the second zero-scatter 321-node fit; two succeeded. These findings must remain visible, rather than declaring the numerical repair complete. All profile values and successful-start counts are retained in results.json. This sparse fixed-scatter profile does not establish a global best scatter or a confidence interval.

Initial broad-bound exploration encountered negligible probability inside the finite support. The reported development profiles therefore use mean width 5-80 days and b=-.5 to 1.5, keeping means inside support for these data. Those are **narrower computational bounds than the frozen experiment**, not physical exclusions or a replacement coverage protocol. Future full-domain implementation must handle the support constraint explicitly. We do not claim to have repaired all boundary behavior by changing those bounds.

## Remaining work

Refine the underlying event-width grid near the narrow-scatter optimum and establish optimizer reliability there. Then freeze the revised estimator, its domain and uncertainty procedure for a new calibration including zero and small scatter. Boundary estimates must be permitted with appropriately tested intervals; an ordinary interior-only confidence argument is insufficient. Retain all failures from the original experiment. The six observational demonstrations remain incomplete.

Reproduce:

```powershell
python research_work/results/timing-continuous-scatter/verify.py
python research_work/results/timing-continuous-scatter/run.py
```
