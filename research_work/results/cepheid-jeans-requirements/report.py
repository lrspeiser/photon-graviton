"""Make an explicit conditional report of the training-only Jeans audit."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
s=json.loads((HERE/'results.json').read_text())
rows=json.loads((HERE/'requirements.json').read_text())
tilt=json.loads((HERE/'tilt-training.json').read_text())
v=json.loads((HERE/'verification.json').read_text())
primary=[r for r in rows if r['model']=='completion']
floor_gap=np.array([r['minimum_vc_for_nonincreasing_radial_pressure_kms']-r['model_vc_kms'] for r in primary])
table='\n'.join(f"| {r['R_kpc']:.2f} | {r['minimum_vc_for_nonincreasing_radial_pressure_kms']:.1f} | {r['model_vc_kms']:.1f} | {r['local_density_efold_growth_kpc']:.3f} | {r['required_midplane_d_vRvz_dz_kms2_per_kpc']:.0f} |" for r in primary)
tilttable='\n'.join(f"| {r['R_low_kpc']}–{r['R_high_kpc']} | {r['above']} / {r['below']} | {r['reflection_linear_slope_kms2_per_kpc']:.0f} | {r['bootstrap_slope_percentiles_2_5_16_50_84_97_5'][0]:.0f} to {r['bootstrap_slope_percentiles_2_5_16_50_84_97_5'][-1]:.0f} | {r['required_slope_for_completion_approx']:.0f} |" for r in tilt)
x=np.arange(4);observed=np.array([r['reflection_linear_slope_kms2_per_kpc'] for r in tilt])
interval=np.array([r['bootstrap_slope_percentiles_2_5_16_50_84_97_5'] for r in tilt])
required=np.array([r['required_slope_for_completion_approx'] for r in tilt])
fig,ax=plt.subplots(figsize=(9,5.5),layout='constrained')
ax.errorbar(x-.08,observed,yerr=[observed-interval[:,0],interval[:,-1]-observed],fmt='o',capsize=4,color='#174364',label='Training-star linear estimate; star-bootstrap range')
ax.scatter(x+.08,required,marker='D',color='#bb4b21',label='Approximate value required by current added field')
ax.axhline(0,color='gray',lw=.8)
ax.set_xticks(x,[f"{r['R_low_kpc']}–{r['R_high_kpc']}" for r in tilt])
ax.set(xlabel='Radius from Galactic axis (kpc)',ylabel='Radial–vertical motion slope ((km/s)²/kpc)',title='Does the omitted vertical-motion term supply the needed correction?')
ax.legend(fontsize=8,loc='upper right');ax.grid(axis='y',alpha=.2)
ax.text(.02,.02,'Exploratory finite-height fit: selection, density gradients and position errors are not fully modeled.',transform=ax.transAxes,fontsize=8,bbox=dict(facecolor='white',edgecolor='none',alpha=.85))
fig.savefig(HERE/'tilt-comparison.png',dpi=170);plt.close(fig)
report=f'''# Can orbital assumptions repair the Cepheid discrepancy?

**A different declining radial profile alone cannot reconcile this frozen field with the training-star moments under steady, axisymmetric, zero-tilt assumptions.** All twelve bins fall below the corresponding minimum circular speed. This eliminates one simple proposed repair within that framework; it does not exclude all companion physics or prove that all the inferred moments are unbiased.

The added-field prediction is below the conditional floor by {floor_gap.min():.2f}–{floor_gap.max():.2f} km/s (bin RMS {np.sqrt(np.mean(floor_gap**2)):.2f} km/s). These are model-to-condition gaps, not statistical significances. Ordinary-matter results are also retained in `requirements.json`.

## What was calculated, and why

**Known physics:** the steady axisymmetric radial Jeans equation relates the potential to the stellar population:

`Vc² = Sphi - SR [1 + d ln(n SR)/d ln R] - T`

`T = (R/n) d(n Q)/dz`, with `SR=<vR²>`, `Sphi=<vphi²>`, and `Q=<vR vz>`.

Here `n` is stellar tracer number density, not ordinary-matter mass density or companion density. The term `T` describes how radial and vertical stellar motions vary together with height. It is not itself vertical gravitational acceleration. These equations are known Jeans mechanics, not formulas unique to our hypothesis. [Primary formulation and moment-error subtraction](https://arxiv.org/html/1810.09466v2).

**Conditional algebra, not a new physical law:** if `T=0` and `n SR` does not increase outward, then `Vc² >= Sphi - SR`. This floor allows any nonincreasing radial-pressure profile, so it does not depend on choosing 4 kpc or 27.3 kpc scale lengths. To reach our lower predicted speed with zero tilt instead requires

`d ln(n SR)/d ln R = (Sphi - Vmodel²)/SR - 1`.

The required slopes are positive in every bin, spanning {s['required_completion_log_pressure_slope_range'][0]:.2f}–{s['required_completion_log_pressure_slope_range'][1]:.2f}. If the previously assumed declining radial second moment is retained, the implied local stellar-density growth lengths are {s['required_completion_density_efold_kpc_range'][0]:.3f}–{s['required_completion_density_efold_kpc_range'][1]:.3f} kpc. A growth length is the distance for a local exponential increase by a factor of e; it is not a measured scale length. No such density profile was fitted or adopted, and raw survey counts cannot determine it without selection corrections.

In plain language: changing how gradually the usual outward-declining stellar population thins out is insufficient. Within this simplified equilibrium picture, making the weaker field work requires reversing that trend, a substantial extra orbital term, or a change in the inferred measurements or gravity.

## Required changes by radius

The last column asks how fast the radial–vertical correlation would have to change with height if the original radial profile is retained. It uses `T_required = Vproxy² - Vmodel²`. Near a reflection-symmetric midplane with `Q=a z`, `T(0)=R a`.

| Radius (kpc) | Conditional floor (km/s) | Added-field speed (km/s) | Required density growth length if zero tilt (kpc) | Required a ((km/s)²/kpc) |
|---|---:|---:|---:|---:|
{table}

These are local bin requirements. They do not construct a global equilibrium distribution or prove that such a distribution is possible.

## Check against above/below-disk training motions

We used only the existing 542 training stars. Measurement covariance is propagated into `vR*vz` under the same proper-motion/RV and independent 7-percent distance-error scenario as the preceding reduction. Within four broad radial regions, a simple odd linear model `Q(z)=a z` is fitted. A second fit allows an intercept. No field parameter is changed.

| Radius (kpc) | Stars above / below | Estimated a | 2.5–97.5% star-bootstrap range | Approximate a required |
|---|---:|---:|---:|---:|
{tilttable}

All slope columns have units `(km/s)²/kpc`. Required broad-bin values are star-count-weighted summaries of the finer-bin requirements, not an exact forward prediction at each star's height.

![Tilt diagnostic](tilt-comparison.png)

The exploratory slopes do not supply the required positive corrections. However, the intervals above resample stars only: they are **not a complete confidence interval for the physical midplane derivative**. The sample is selected, radial variations are pooled, the data extend to finite heights, position errors are ignored in the regression, and vertical density gradients away from the plane matter. The median absolute heights range from about 0.043 to 0.223 kpc, and the outer region has only 26 stars. Allowing an intercept changes the outer slope from about 34 to 9, illustrating model sensitivity. We therefore do not turn this comparison into a sigma-level rejection or apply it as a validated correction to the rotation curve.

## Consequences for the larger task

This result makes the next choice more concrete. We should not spend further effort selecting a different positive exponential scale length to repair this particular field. The remaining plausible work is to quantify distance/selection and non-equilibrium effects, constrain ordinary-matter components, and test a physically specified three-dimensional companion response. Any change to that response must also address the earlier excessive pull toward the disk, rather than independently increase the radial force in each region.

The result does not yet compare stars inside the bulge; these Cepheids span 6–18 kpc. The large red-giant sample remains necessary for the user's bulge-above/below-versus-plane comparison. Photon redshift, event timing, conversion and permanent storage remain separate unsolved requirements. The total energy-supply budget is still deferred.

## Verification and provenance

This analysis uses previously exposed training results and introduces no fresh holdout claim. Validation/test outcomes remain unopened. The regression archive includes every leave-one-star-out slope range and maximum individual leverage; these are influence diagnostics, not exclusions of stars.

`run.py` records source hashes and checks training-only identifiers. `verify.py` tests the Jeans signs with nine analytic density/velocity fields using independent finite differences (maximum error {v['max_finite_difference_error_kms2']:.4g} (km/s)²), verifies the inverse requirements, and checks the pressure floor across {v['bound_grid_checks']} slope cases. Cross-covariance diagonals reproduce the previous error propagation, and all regional star counts match. These checks validate the calculation under its assumptions, not those assumptions themselves.

Run `run.py`, `verify.py`, then `report.py` using Python. Stellar rows remain in the ignored cache; the code, input hashes, per-bin requirements and exploratory regression summaries are versioned. No earlier model score, calibration or holdout assignment was changed.
'''
(HERE/'report.md').write_text(report,encoding='utf-8',newline='\n')
print(json.dumps(dict(conditional_floor_gap_min_kms=float(floor_gap.min()),max_kms=float(floor_gap.max()),rms_kms=float(np.sqrt(np.mean(floor_gap**2)))),indent=2))
