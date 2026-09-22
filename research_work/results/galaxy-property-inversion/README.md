# JR-2: what galaxy-property changes let the unchanged formula fit?

**21 September 2026, America/Los_Angeles. Executed inverse inference, not corrected observations.**

## Theory first

Our working hypothesis is that matter and radiation generate an extended companion state with connected stellar-motion and lensing responses. Here the direction of inference changes: hold the exact JR-1 R10 force law and every universal coefficient fixed, and determine which galaxy properties would need to change. Such properties are conditional predictions of the assumed formula, not evidence the published observations are wrong.

Baseline GitHub main was `df6bb0a20062aa3bbcdfecccbecba1cc3a80ae7f`. The original JR-1 package supplies the exact parameters, data and inherited forward calculations. All original input and formula-file hashes are unchanged. Scope is **149 SPARC galaxies and six SLACS galaxies, not clusters**; no cluster forward model is implemented here.

## What changed and what did not

For disks, vary one distance, one inclination, one stellar mass-to-light multiplier and one gas multiplier per galaxy. No independent radial corrections, companion amplitudes or gravity coefficients. At fixed angular profile, distance factor d changes radius and effective radius by d, photometric/gas masses by d^2 and baryonic squared speeds by d. Stellar/gas multipliers propagate through both ordinary gravity and the unchanged companion-source prescription. The R10 disk/spheroid relation is unchanged.

Inclination variations keep the reconstructed projected velocity `v_catalog*sin(i_catalog)` and its uncertainty fixed and compare them to `v_model*sin(i_new)`. This is a projected major-axis proxy, not a new fit of full velocity maps or inclination-dependent image deprojection. The original signed-gas and clipping conventions remain.

The restricted fit allows distance/inclination within two quoted marginal catalog errors, stellar M/L within 0.30 dex and gas within 0.20 dex. It penalizes departures using the catalog errors and working stellar/gas scales of 0.15/0.10 dex. **The stellar/gas allowances are sensitivity assumptions, not individual measured uncertainties.** The original fixed-distance policy is relaxed only in separate counterfactual copies at the user's request.

The wide diagnostic allows distance factors 0.2–5, inclination 5–90 degrees, stellar multipliers 0.1–10 and gas multipliers 0.25–4. These are not plausible error bars. Both the original quoted-error objective and a separately labelled post-primary fractional-RMS objective are preserved. Optimizing the latter does not mean the quoted errors were enlarged or changed.

For lenses, freeze observed angular light profiles, 40 V_rms measurements with full within-object covariance, Einstein-angle summaries and the existing conditional distance geometry. Change stellar mass and constant orbital anisotropy, using every new mass in both ordinary and companion gravity. Mass bounds are +/-0.30 dex, beta in [-1,0.35]. No independent photon-response coefficient. The inherited 5% lens-deflection scale in joint optimization is a working tradeoff, not a measured uncertainty.

## Rotation results: improvement but no single universal correction

| Input treatment | Below 10% rotation RMS | Below 20% rotation RMS | Median fractional RMS |
|---|---:|---:|---:|
| Original R10 properties | 36/149 | 95/149 | 15.05% |
| One common four-parameter property recalibration | 32/149 | 93/149 | 14.98% |
| Separate restricted property fits | 91/149 | 120/149 | 8.34% |
| Separate wide fits, fractional-RMS objective | 122/149 | 142/149 | 5.75% |

These are descriptive thresholds, not a statistical acceptance criterion. All 3,152 points and all 149 galaxies are retained. The original outlier subgroup is defined before adjustment as fractional RMS >20%: 54 galaxies. Restricted changes rescue 26 of those; the wide fractional-RMS calculation rescues 47. One previously sub-20% galaxy crosses the threshold in the restricted uncertainty-weighted fit, leaving 29 total above it.

Restricted fits reach at least one limit in **70/149** cases. Wide fractional-RMS fits reach a limit in **148/149**. The latter are not evidence that small catalog corrections explain the model. Raw rotation chi-square remains 26428.5 in the restricted fit and 28407.3 in the broad fractional-RMS fit, over 3,152 measurements and before accounting for all fitted property parameters. Minimizing fractional RMS does not minimize measurement-error chi-square.

The wide quoted-error joint fit, preserved separately, reaches 136/149 below 20%, median fractional RMS 6.50%, and chi-square 14541.5. Its 146 boundary cases expose substantial nonuniqueness. None of these fits proves a global optimum or validates the microscopic companion mechanism.

## The common disk-property changes

Of the 54 original outliers, **38 are overpredicted** on mean signed fractional residual and **16 underpredicted**. Their restricted joint fits ask for different changes:

| Original class | Count | Median distance change | Median inclination change | Median stellar M/L change | Median gas change |
|---|---:|---:|---:|---:|---:|
| Model predicts speeds too high | 38 | -8.6% | -5.39 degrees | -29.7% | +1.4% |
| Model predicts speeds too low | 16 | +8.6% | +1.64 degrees | +87.2% | -5.6% |

The overpredicted class generally needs less inferred source gravity and/or a more face-on orientation. A smaller inclination raises the intrinsic speed inferred from a fixed projected Doppler signal. The other class needs the opposite correction, especially more stellar mass per unit light. These are conditional directions selected using the residuals, not independent evidence of a catalog bias.

One-property-at-a-time wide diagnostics rescue 27/54 through distance, 25 through inclination, 20 through stellar mass and 11 through gas. Only 8 of the distance-only rescues and 6 of the inclination-only rescues stay within two quoted catalog errors. These alternative groups overlap and their counts must not be added.

### What the difficult objects share

| Original input characteristic | 54 outliers | Other 95 |
|---|---:|---:|
| Median effective 3.6-micron surface brightness [solar luminosities/pc^2] | 30.83 | 178.99 |
| Median input gas fraction Mgas/(Mgas+Mstar) | 0.639 | 0.386 |
| Median morphology code | 9 | 6 |
| Medium-quality Q=2 rotation curves | 28/54, 51.9% | 28/95, 29.5% |

The difficult subset is typically fainter per unit area, more gas-rich and later-type, with a larger medium-quality fraction. Catalog Q=1 and Q=2 are both retained. These are descriptive, post-hoc associations in an already exposed sample. They do not determine whether source reconstruction, viewing geometry, non-circular motion or a missing physical response causes the discrepancy. The residual correlation with inclination itself is weak (Spearman 0.039): not all failures are face-on objects.

### Examples of the actual changes required

| Galaxy | Original -> restricted RMS | Distance [Mpc] | Inclination [deg] | Stellar M/L multiplier | Gas multiplier |
|---|---:|---:|---:|---:|---:|
| UGC07125 | 55.8% -> 3.7% | 19.80 -> 8.59 | 90.0 -> 90.0 | 1.095 | 0.866 |
| UGC07866 | 30.5% -> 5.7% | 4.57 -> 4.48 | 44.0 -> 35.2 | 0.866 | 0.906 |
| UGC07399 | 32.6% -> 2.3% | 8.43 -> 13.49 | 55.0 -> 57.9 | 1.828 | 1.136 |
| UGC09037 | 24.2% -> 6.2% | 83.60 -> 79.82 | 65.0 -> 63.1 | 0.501 | 1.012 |
| ESO444-G084 | 25.5% -> 8.3% | 4.83 -> 5.12 | 32.0 -> 40.8 | 0.890 | 1.366 |
| DDO170 | 22.6% -> 6.7% | 15.40 -> 10.42 | 66.0 -> 64.0 | 0.872 | 1.430 |

UGC07125's 57% distance reduction is about 1.90 quoted flow-distance errors, illustrating why a two-error allowance need not be a small physical change. UGC07399 instead needs a 60% increase and reaches its distance bound. UGC09037 reaches the lower stellar M/L bound. They are competing demands, not one calibration shift.

The common four-parameter calibration chooses distances x0.799, inclination +0.0065 degrees, stellar M/L x1.260 and gas x1.585 (the upper bound). Mean galaxy velocity RMS barely changes, **17.71 -> 17.53 km/s**, and the number below 20% goes **95 -> 93**. It does not explain the success of separate property changes.

Seven curves remain above 20% even in the executed wide fractional-RMS fits: **UGC11557, D631-7, NGC2903, NGC2915, NGC4217, CamB and NGC3741**. This is not a proof that no other source reconstruction could fit them.

## Lenses: extra stellar mass fixes rings but worsens the stars

All six Einstein radii can be matched exactly by solving for stellar mass, without changing the force formula. That is a conditional mass inference from one ring constraint, not an independent prediction.

| Lens | Stellar mass change for exact ring | Original stellar RMS | Stellar RMS after mass change at old orbit | After also refitting orbit |
|---|---:|---:|---:|---:|
| J0037-0942 | -4.7% | 1.24% | 2.46% | 2.46% |
| J1112+0826 | +20.5% | 1.69% | 6.95% | 4.75% |
| J1204+0358 | -4.7% | 1.37% | 2.76% | 2.74% |
| J1402+6321 | -4.7% | 1.95% | 2.98% | 2.81% |
| J1621+3931 | +27.1% | 3.50% | 6.05% | 5.19% |
| J1630+4520 | +25.3% | 1.82% | 8.17% | 6.40% |

The three under-bent systems need about 20–27% more stellar mass than the existing R10 source. Their accurate stellar-motion profiles do not favor that simple repair. Orbital adjustment helps but leaves total chi-square 348.4 over 40 bins in the exact-ring-plus-orbit calculation; two of the other systems hit beta=0.35. The corresponding mass shifts are 0.90, 1.74 and 1.40 published log-mass errors, but only as scale comparisons: R10 masses already contain fitted normalizations, and the published errors do not describe all conditional-geometry/population systematics.

### Required environmental or geometric input

Keep the R10 internal galaxy unchanged. Let F be its fraction of the required deflection at the observed ring. Then required effective external convergence is `kappa=1-F`, and the alternative required multiplier of `Dls/Ds` is `1/F`.

| Under-bent lens | Required extra convergence | Alternative Dls/Ds increase | Foreground redshift |
|---|---:|---:|---:|
| J1112+0826 | 0.1065 | +11.9% | 0.2730 |
| J1621+3931 | 0.1277 | +14.6% | 0.2449 |
| J1630+4520 | 0.1245 | +14.2% | 0.2479 |

This is a required line-of-sight/geometry property, **not a discovered external mass distribution**. A free convergence can algebraically close one ring condition; it is not a full lens-image fit or a same-theory environmental construction. Actual foreground/background matter, its tidal effects and lens-plane geometry must be modeled. The distance-ratio alternative is degenerate at this level and is not a percentage change in galaxy distance itself. Negative required convergence for the other three (-0.028 to -0.031) represents over-focusing relative to the baseline, not literal negative matter.

**A shared feature worth isolating:** all three under-bent lenses have foreground redshift **0.2449–0.2730**, while the other three are at **0.1644–0.2046**. This makes distance/redshift or stellar-population normalization a concrete possible contributor. Six post-hoc objects do not establish a causal redshift trend. Their apparent axis ratios do not divide the two groups similarly cleanly.

## Interpretation

The exercise finds two different kinds of demand. Disk outliers are often low-surface-brightness, gas-rich systems requiring object-specific changes in inferred source strength or viewing geometry. Lenses need more projected focusing than their stellar dynamics imply, which points toward geometry or line-of-sight structure as alternatives to simply adding central stellar mass. Neither conclusion shows the published data are wrong or that R10 is right.

Distance and inclination are nearly degenerate: ordinary squared speeds scale approximately as d, while the companion amplitude scales as d^(2p) with 2p=1.0525. Projecting by sin(i) permits compensations. Extreme boundary estimates are not precise galaxy measurements.

No source-energy, microscopic-conversion, disappearing-time, Solar-system, cluster or merger solution is claimed. The broader companion hypothesis is not identified with this single empirical construction.

## Reproduce and evidence

The exact executed scripts are [inverse_properties.py](inverse_properties.py) and [supplement.py](supplement.py). Use the original `Photon-Graviton-JR1` package directory as `--bundle`, or `original-JR1` inside the new complete package:

```sh
OPENBLAS_NUM_THREADS=1 python inverse_properties.py --bundle original-JR1 --output fresh-JR2
OPENBLAS_NUM_THREADS=1 python supplement.py original-JR1 fresh-JR2
```

The full originating-conversation package preserves all original inputs, all changed-property copies, optimizer/bound flags, observations, predictions, one-at-a-time inversions, metadata, and this report. It does not rely on a running process or temporary Actions retention. Original byte hashes remain unchanged. Baseline velocities reproduce archived R10 to **5.68e-14 km/s**.

A separate readback audit reconstructs **1,192** disk predictions from explicit algebra rather than calling the fitted forward function. Together with stored lens-covariance checks it performs **3,636** array/metric comparisons, maximum projected-speed discrepancy **8.53e-14 km/s**, maximum scaled metric discrepancy **8.94e-15**. This is numerical verification of these calculations, not observational validation or a rerun of every historical suite.

Primary context: Lelli et al. (2016), arXiv:1606.09251, SPARC; Li et al. (2018), arXiv:1803.00022, uncertain distance/inclination/M/L methodology; the pinned SLACS/TDCOSMO inputs in JR-1; Schneider and Sluse (2013), arXiv:1306.0901, lens-model degeneracies. These sources do not independently establish the counterfactual properties inferred here.
