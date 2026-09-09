# A shared companion-deposition profile fitted to galaxy data

9 September 2026. This is an empirical fit within the fictional companion-wave hypothesis. It identifies a positive extra-source distribution capable of improving galaxy rotation predictions. It does not derive capture, establish an energy supply, or demonstrate a new light–gravity interaction.

## Result

Among three positive-density cored profiles, the validation-selected formula is

\[
\boxed{\rho_X(r)=\frac{V_X^2}{4\pi G(r^2+r_c^2)}}
\]

with shared parameters

\[
\boxed{r_c=1.0723 R_d},
\qquad
\boxed{V_X^2=1.5723\times10^4\,(\mathrm{km/s})^2
\left(\frac{L_{3.6}}{10^{10}L_{\odot,3.6}}\right)^{0.3395}
\left(\frac{R_d}{3\,\mathrm{kpc}}\right)^{0.3792}}.
\]

R_d is the measured stellar disk scale length. L_3.6 is the observed 3.6-micron luminosity, not bolometric luminosity or a measured companion power. The fit supplies four constants shared by all galaxies; no individual halo amplitude or radius is fitted on the comparison objects. Rounded coefficients are for presentation; results.json contains full precision. These coefficients do not have established physical significance or quoted confidence intervals.

The corresponding enclosed source and rotation prediction follow by integration:

\[
M_X(<r)=\frac{V_X^2}{G}\,[r-r_c\arctan(r/r_c)],
\]
\[
\boxed{v_{\rm pred}^2(r)=v_b^2(r)+V_X^2
\left[1-\frac{r_c}{r}\arctan\left(\frac r{r_c}\right)\right]}.
\]

Near the center the density is finite and the added squared speed scales as r². At large radius the density falls as r⁻² and the added squared speed approaches V_X². The source must eventually truncate or steepen; continuing it to infinity would require infinite mass and energy. No outer truncation radius is determined by these fits.

This is a familiar cored isothermal halo shape with a newly fitted shared photometric scaling, not a newly discovered fundamental gravity law. It can represent a deposited field in the fictional model, but rotation data do not identify what the source is made of.

## What “depositing it” requires

If stored companion energy density u_d produces equivalent source density ρ_X=ηu_d/c², the required stored energy is

\[
u_d(r)=\frac{c^2V_X^2}{4\pi G\eta(r^2+r_c^2)}.
\]

For constant local deposition power per volume q_d, zero initial deposit, residence time τ and no spatial redistribution,

\[
\dot u_d=q_d-u_d/\tau,\qquad
u_d(T)=q_d\tau(1-e^{-T/\tau}).
\]

Consequently a candidate deposition rule is

\[
\boxed{q_d(r)=\frac{c^2V_X^2}
{4\pi G\eta T_{\rm eff}(r^2+r_c^2)}},
\quad T_{\rm eff}=\tau(1-e^{-T/\tau}).
\]

Use consistent units, e.g. SI, to obtain q_d in W/m³. This formula specifies how much deposition would produce the fitted gravitational profile. It does not show that photons can supply that power. Only the combination η q_d T_eff is determined. For absorption q_d=α_c v_c u_c, it constrains α_c v_c u_c η T_eff, not each factor separately.

At a chosen outer boundary R, integrate to obtain

\[
P_{\rm required}(<R)=\frac{c^2V_X^2}{G\eta T_{\rm eff}}
[R-r_c\arctan(R/r_c)].
\]

A physical external-collection model must independently satisfy P_required ≤ πR²v_cu_c for capture fraction at most one, or calculate a different effective cross-section from dynamics. The source profile does not resolve the earlier energy-budget shortage at η=1. Amplification η, residence time, and incoming companion energy remain explicit unknowns.

## Data and fitting procedure

Source: [SPARC](https://astroweb.cwru.edu/SPARC/), [Lelli et al. 2016](https://arxiv.org/abs/1606.09251). Starting from 175 galaxies, retain the same 149-galaxy, 3,150-point selection as the prior companion-wave study: quality flag ≤2, inclination ≥30 degrees, positive disk scale length, and at least five valid rotation points. Fixed disk and bulge mass-to-light ratios are 0.5 and 0.7; signed gas contributions are retained. Catalog distances and inclinations are fixed.

The previously used split has 89 training, 29 validation and 31 test galaxies. Fit each candidate on training galaxies only, minimizing the galaxy-equal mean squared log10 speed residual. Choose among the three cored models by validation log error. The test subset is then reported with shared parameters unchanged. These subsets have been exposed in earlier investigations, so this is exploratory reuse, not a new blind test.

Three starts per optimizer were used. All optimizers returned success. No observational-error likelihood or parameter posterior was computed. Galaxy-equal km/s RMSE is a descriptive score; fit selection used log-speed error. Photometry, baryonic mass modeling and rotation-distance calibration share assumptions, so this is not an entirely independent measurement pipeline.

## Comparison

| Model | Shared parameters | Training RMSE (km/s) | Validation RMSE | Reused test RMSE |
|---|---:|---:|---:|---:|
| Baryons only | 0 | 52.57 | 58.22 | 47.77 |
| Area-motivated cored profile | 2 | 26.58 | 32.94 | 23.19 |
| Luminosity-scaled cored profile | 3 | 19.39 | 28.10 | 15.75 |
| Luminosity-and-size cored profile | 4 | 19.42 | **26.00** | **14.76** |
| Acceleration-response comparison | 2 | 20.01 | 27.38 | 17.20 |

The selected cored profile has validation/test log-speed RMSE 0.09046/0.07040 dex. The luminosity-only version gives 0.09849/0.07719 dex. The extra size parameter improves these descriptive scores, but the small difference is not established as statistically significant. The simple three-parameter version is a useful lower-complexity alternative:

\[
V_X^2=1.3450\times10^4(\mathrm{km/s})^2(L_{3.6}/10^{10}L_{\odot,3.6})^{0.4476},
\quad r_c=0.9254R_d.
\]

The area-motivated model uses V_X² proportional to R_d and r_c proportional to R_d, so M_X at a fixed multiple of R_d scales as R_d². It replaces the previous observing-limit R_max with a photometric size. Its best fit is V_X²=2.2234×10⁴(km/s)²(R_d/3 kpc), r_c=2.1204 R_d. This differs from the earlier outer-mass-only area regression; the scores should not be compared as identical tests.

The acceleration comparison uses

\[
g_X=A a_{\rm ref}(g_b/a_{\rm ref})^q,
\quad a_{\rm ref}=10^{-10}\,\mathrm{m/s^2},
\quad A=0.70275,\quad q=0.46246.
\]

It is a force prescription, not automatically a positive deposited-energy profile. Its spherical-equivalent M_X=r²g_X/G decreases over 1 of 3,001 neighboring radial intervals at fixed baryonic inputs. This finite-difference diagnostic can be affected by noisy baryon modeling; it is not a statistical exclusion. The selected cored profile, by construction, has nonnegative density at every radius.

## Cluster transfer: what could actually be tested

The cored galaxy formula uses stellar disk size and luminosity. Clusters are not stellar disks, and the available cluster table does not contain comparable structural measurements. We did **not** replace R_d with R_500 and claim a valid transfer of this formula.

The acceleration comparison can be transferred to the 11 clusters from the preceding study using enclosed baryonic mass and R_500. Use g_b=G M_b/R_500², M_b=M_500(f_gas+0.02), and freeze the galaxy-fitted A,q. Sources: [Ettori et al. Table 1](https://arxiv.org/abs/1805.00035), [Eckert et al. Table 2](https://arxiv.org/abs/1805.00034). The 2% stellar fraction is a stated assumption; hot-gas fractions are published hydrostatic estimates. The previous exclusions and standard distance/HSE assumptions are retained.

The required extra cluster mass exceeds this prediction by a median factor **1.90**, range **1.70–2.42**. This result belongs to the acceleration formula, not to the selected cored galaxy formula. Allowing stronger deposition or retention could address it, but this test does not determine the responsible interaction. No cluster radial profile or lensing map was fitted.

## Interpretation and next discriminating step

The useful result is a shared positive-density formula that places an additional source broadly through galaxies and predicts their rotation curves materially better than baryons alone. It avoids assigning an arbitrary deposit amplitude separately to every comparison galaxy. It is still a fitted halo description, with ordinary astrophysical parameter assumptions held fixed.

A stronger next step is to derive an absorption/redistribution rule whose integrated solution yields this profile and predicts its normalization from an independently estimated radiation supply. Such a rule should use quantities meaningful in both gas-dominated clusters and stellar disks. Until then, the luminosity and size exponents describe correlations; they do not establish photons as the source of extra gravity or solve the local-clock and signal-stretching issues.

The archive contains fit.py, protocol.json, full numerical results, comparison predictions, and source data. Extract preserving both companion_deposition_fit and companion_wave_test directories; run `python companion_deposition_fit/fit.py` with NumPy and SciPy installed. The cluster source transcription is retained in the input CSV and prior cluster_test.py. No per-object fitting occurs for the comparison galaxies.
