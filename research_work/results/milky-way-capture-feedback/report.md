# Capture responds to the gravity of its own deposits

## Result

The previous Milky Way calibration does not survive the proposed gravitational feedback at its fitted normalization. With the same incoming-intensity/age combination, the p=3 training-bin RMS discrepancy rises from8.76 to370.33km/s; p=2 rises from13.78 to193.15km/s. The final p=3 prediction spans about351–809km/s over the existing training bins, whose inferred circular-speed proxies are about223–242km/s.

This tests a missing physical coupling, not a new observational fit. It rejects reuse of the previous fitted amplitudes in this chronological feedback model. It does not rule out every alternative capture law or a recalibration of the coupled equations.

## What changed

The earlier profile used only ordinary-matter well depth. We now let accumulated deposits contribute to the depth that controls subsequent capture. **Project growth postulate, not a known microscopic law:**

\[
\frac{\partial\rho_D}{\partial s}
=C\left[\frac{W_{\rm total}}{W_{\rm total}+W_*}\right]^p,
\quad W_{\rm total}=\max[-\Phi_b-\Phi_D(s),0].
\]

s is dimensionless exposure age from0 to1. C is the same calibrated density scale from the previous100kpc case, Wstar=40000(km/s) squared, and p=2 or3. Start with zero deposits. **Known Newtonian response used conditionally:** solve Laplacian Phi_D=4pi G rho_D with zero potential at infinity. The same potential gives stellar forces.

We integrate the growth history. We do not replace it with rho_D=C*f(W_final), which would incorrectly treat the final capture rate as if it had operated throughout formation. The fixed-depth control reproduces the previous profile, because then the capture rate is constant in exposure time.

Constant isotropic optically thin supply and permanent in-place retention remain assumptions. C combines supply, opacity and age; none is independently measured here. The100kpc cutoff and ordinary-matter geometry are unchanged. A quasistatic Poisson response is used; finite propagation during formation, kinetic support and gravitational work are not closed by this test.

## Why the earlier shape disappears

Deposits deepen the potential throughout the capture region. The saturating depth factor then approaches one over most of the volume. Capture loses much of its original preference for the ordinary-matter shape, and the deposited density becomes much more uniform.

At final exposure, rho_D/C ranges from0.911 to0.965 for p=2 and0.885 to0.951 for p=3. A nearly uniform spherical density produces circular speed increasing roughly with radius, unlike the desired nearly flat contribution. The computed field includes axisymmetric structure rather than assuming a sphere, but this explains its broad trend.

This is strong positive feedback within a bounded capture-rate law, not a demonstrated finite-time mathematical divergence. Since0<=f(W)<=1, rho_D<=C*s over this finite experiment; the code verifies that bound. Permanent deposits would continue growing under continuing supply. A time-independent total mass has not been established.

## Same parameters, before and after feedback

| Variant | Fitted C, solar masses/kpc cubed | No-feedback training RMS | Feedback training RMS | Final deposit mass within100kpc |
|---|---:|---:|---:|---:|
| p=2 | 5.51336e7 | 13.78km/s | 193.15km/s | 2.12e14 solar masses |
| p=3 | 1.27119e8 | 8.76km/s | 370.33km/s | 4.74e14 solar masses |

Those masses are requirements of the unmodified capture exposure, not observed Galactic masses. No assertion is made that photons can supply them. The previous p=3 no-feedback mass was8.31e11 solar masses under the same cutoff. Adding the proposed feedback is therefore consequential even before energy availability is considered.

Examples of p=3 predictions with feedback:

| Training bin | Inferred Jeans speed proxy | Predicted speed |
|---|---:|---:|
| 6–7kpc | 242.0km/s | 351.2km/s |
| 8–9kpc | 235.6km/s | 425.0km/s |
| 11–12kpc | 232.6km/s | 543.1km/s |
| 17–18kpc | 226.2km/s | 808.6km/s |

All12 predictions and the formation histories are in [results.json](results.json). As before, the observed-side values are approximate Jeans inferences with distance, population and equilibrium assumptions. These are unweighted training diagnostics, not significance estimates or newly tested stars.

## Verification and provenance

The [protocol](protocol.md) was declared before calculation. [run.py](run.py) uses a precomputed Legendre Green kernel and explicit midpoint growth updates on a positive real-space density. It reuses only the542 training-star radii and previously exposed12 bins. No reserved stars or vertical force inferences are read in this test. The preceding fitted-result file is hashed.

The frozen-depth control is evaluated with this independent radial integration method and differs from prior predictions by at most0.0171km/s, or0.00474km/s on the fine grid. This isolates the large change as feedback rather than a changed baseline or force solver.

The initial p=3 time-resolution difference exceeded the declared2km/s tolerance. That failure remains archived. Refining to256 exposure steps gives a128-to256-step prediction difference below0.761km/s on the fine grid. At256 steps,400-to800 radial-node predictions differ by less than0.117km/s; the mass difference is below0.015%. All final declared gates pass. p=2 passed at128 steps. The finite density bounds hold throughout.

Reproduce both stages:

```text
python research_work/results/milky-way-capture-feedback/run.py
python research_work/results/milky-way-capture-feedback/run.py --refine
```

The refinement retains the first stage's runs and failed checks. Limiting BLAS threads to one can improve reproducibility of run time; it does not change the model.

## Next required test

Calibrate the exposure normalization using the **coupled growth equations**, and test the entire resulting radial curve, rather than transferring a no-feedback fit. Keep p, depth scale and cutoff explicit; do not choose a separate density at every radius. If a coupled fit survives, check its stability to source history/cutoff, its actual photon-energy budget, and its vertical/lensing consequences with admissible observations.

An alternative law that responds only to ordinary matter would be a distinct physical hypothesis needing justification. It cannot be introduced silently to recover the prior8.76km/s result. Persistent redshift, physical energy/momentum conservation, mechanical retention and genuine unseen-data evaluation remain incomplete; all nine goals stay active.
