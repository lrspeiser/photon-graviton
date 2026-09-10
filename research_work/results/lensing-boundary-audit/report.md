# Outer deposits affect lensing even when inner rotation is unchanged

The frozen extra-gravity templates now produce 1,788 conditional companion-deflection predictions across the 149 retained SPARC galaxies. Moving a trial source boundary from the last measured rotation radius to ten times that radius changes median extra light bending by about 14%, 23% and 59% at impact parameters of one quarter, one half and one times the last measured radius, respectively. The inner circular-speed prediction stays unchanged under the stipulated spherical-source construction.

This is a calculation of predictive ambiguity, **not an observational lensing fit**. No lensing measurements have been compared here, and no reserved stellar likelihood scores were opened. An added metric assumption is necessary to turn the fitted stellar force into light bending.

## Frozen setup and results

We retain A=0.2422960666, p=0.4624587420 and a*=7.249608969×10⁻¹⁰ m/s² from the archived empirical galaxy fit. At every existing radius, the extra radial acceleration remains g_c=A a*(g_b/a*)^p. This is the project's empirical power-law template using known mathematics; its photon interpretation is unproven. We introduce no newly fitted coefficient.

Each galaxy's radial extra-force template is interpolated on the original measured interval. Between the last measured radius R_max and a trial boundary r_t, we assume g_c(r)=g_c(R_max)(r/R_max)^(-2p), corresponding to the finite-baryonic-mass exterior approximation. Beyond r_t, enclosed equivalent source mass remains constant, so g_c(r)=g_c(r_t)(r_t/r)². This is a **conditional spherical geometry and exterior assumption**, not a measured capture profile. R_max is the endpoint of the available rotation measurements, not an observed physical edge. Some ordinary matter can lie beyond it; the exterior approximation does not establish otherwise.

| Closest ray distance from galaxy center | Median extra deflection ratio, r_t=10R_max versus R_max | Range across all 149 templates |
|---|---:|---:|
| 0.25R_max | 1.137 | 1.093–1.249 |
| 0.5R_max | 1.231 | 1.195–1.353 |
| R_max | 1.590 | Same ratio by exterior construction |

These ratios concern the **extra companion contribution only**. They are not percentage changes in total lensing, image separation, Einstein radius or shear. Ordinary matter also bends light and would need to be included in an observational prediction.

For example, the NGC2403 template at b=10.435 kpc gives extra deflections of 0.11634, 0.12991, 0.13828 and 0.14123 arcseconds for trial cutoffs of 1, 2, 5 and 10 times R_max=20.87 kpc. These are predictions under the stated assumptions, not measured deflections of NGC2403.

In plain language, a star inside a spherical outer shell feels no net pull from that shell. A light ray passing through and beyond the system samples a longer path through its field. Extending the outer deposits can therefore leave the inner rotation unchanged while changing the bending. The capture boundary is observationally consequential, not merely a bookkeeping choice.

## Formulas and the extra lensing assumption

**Known weak-field metric and lensing construction:** write

\[
ds^2=-(1+2\Phi/c^2)c^2dt^2+(1-2\Psi/c^2)d\mathbf x^2.
\]

Slow stellar motion constrains −∇Φ. Light deflection depends on the transverse gradient of Φ+Ψ. Our empirical stellar-force fit does not determine Ψ.

For this calculation we explicitly assume **Ψ_c=Φ_c**, a familiar no-slip weak-field choice, not a derived companion field equation. With static spherical sources, straight unperturbed rays and distant endpoints, the extra deflection magnitude is

\[
\hat\alpha_c(b)=\frac{4}{c^2}\int_0^\infty
g_c\!\left(\sqrt{b^2+z^2}\right)\frac{b}{\sqrt{b^2+z^2}}\,dz.
\]

All these are known lensing equations applied to proposed source profiles. No claim of mathematical novelty is made. With a constant Ψ_c=ηΦ_c, the predicted companion deflection instead scales by (1+η)/2 relative to the table. A spatially varying relation requires a new integral. In particular, matching stellar rotation alone cannot select η. The previously studied minimal conformal scalar example has canceling direct scalar lensing contributions at its stated order; the equal-potential assumption here does not repair that example automatically.

**Known projected-source identity, used as an independent check under the same equal-potential assumption:**

\[
\hat\alpha_c(b)=\frac{4G M_{c,\mathrm{cylinder}}(<b)}{c^2b}.
\]

For a spherical shell at r>b, the fraction of its mass within that cylinder is 1−√(1−b²/r²). We independently integrated these projected shells and compared with the ray-force integral for every prediction. The largest relative discrepancy is 2.1×10⁻¹². Point-source and analytic power-law integral limits were also checked. These validate the numerical calculation under its assumptions; they do not validate the assumptions against observations.

The standard weak-field relationship is described in [Bartelmann & Maturi, Weak gravitational lensing](https://arxiv.org/abs/1612.06535). This calculation uses the local isolated-lens deflection, not a cosmological distance-redshift law. Image positions and shear require additional source/lens geometry and an observational likelihood. Those remain uncomputed.

## Three source-eligibility warnings

Unlike the earlier 63-point Milky Way check, extending the spherical-source reconstruction across these 149 different templates exposes three cases with negative equivalent density in the force interpolant: NGC6015, NGC6946 and UGC07125.

- **NGC6015 and UGC07125:** enclosed equivalent mass increases at the original nodes, but the interpolated acceleration produces a negative-density interval. Monotone interpolation of enclosed mass is an available alternative that preserves node values. It has not been silently substituted here.
- **NGC6946:** enclosed equivalent mass already decreases between some original template nodes. Interpolation alone cannot make every node exact and the spherical equivalent density nonnegative. A different geometry, adjusted inputs within uncertainties or a modified response would be needed. This is not evidence that the galaxy has negative physical mass.

All three are retained and flagged; their curves are formal effective-potential calculations, not accepted positive deposited-mass models. Restricting the summary to the 146 templates positive on the refined probe grid gives essentially the same median boundary ratios: 1.138, 1.231 and 1.590. Grid positivity is not a proof over all radii.

Original train/validation/test labels are preserved as provenance. These galaxy predictions and their parent measurements were previously exposed; calculating field eligibility is not a new blind test of observed rotation. No lensing observation is used to choose a cutoff, interpolation or metric response.

## What this changes in the work plan

To make a joint motion-and-lensing prediction, the companion model must specify both **where its deposited source ends or changes character** and **how that source affects the spatial as well as temporal metric**. The current rotation fit supplies neither. A free outer radius and a free lensing multiplier could be tuned to data, but that would be additional empirical modeling, not proof of photon conversion.

The next observational lensing comparison must carry source/lens distances and their assumptions, measured shear or image data rather than only a fitted dark-halo mass, and ordinary-matter uncertainties. The same capture and metric rule must then predict those observations with parameters fixed from the declared training analysis. The independent supernova-timing, propagation, conservation and source-energy requirements remain open where previously recorded.

Run run.py and verify.py to reproduce this audit from the existing archived rotation predictions and coefficients. results.json records input hashes, all 1,788 prediction rows are in predictions.json, and source eligibility plus independent numerical checks are archived alongside them. The three impact-parameter fractions and four cutoff ratios are declared sensitivity choices, not optimized measurements or confidence bounds.
