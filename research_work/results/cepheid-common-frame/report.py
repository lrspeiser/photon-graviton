"""Present the training-only diagnostic, without opening held-out outcomes."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
s=json.loads((HERE/'results.json').read_text())
rows=json.loads((HERE/'training-bins.json').read_text())
sens=json.loads((HERE/'distance-sensitivity.json').read_text())
verification=json.loads((HERE/'verification.json').read_text())
R=np.array([r['R_mean_kpc'] for r in rows]);obs=np.array([r['jeans_proxy_kms'] for r in rows])
low=np.array([r['jeans_proxy_kms'] for r in sens['0.93']]);high=np.array([r['jeans_proxy_kms'] for r in sens['1.07']])
maxshift=float(np.max(abs(np.r_[low-obs,high-obs])))
maxnoise=float(max(abs(r['jeans_proxy_kms']-r['uncorrected_jeans_proxy_kms']) for r in rows))
fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
ax.plot(R,obs,'o-',color='#153c60',label='Training stars: approximate Jeans estimate')
ax.plot(R,[r['ordinary_vc_kms'] for r in rows],color='#777777',label='Frozen ordinary matter')
ax.plot(R,[r['completion_vc_kms'] for r in rows],color='#ba4b20',label='Frozen additional field')
ax.fill_between(R,np.minimum(low,high),np.maximum(low,high),color='#153c60',alpha=.15,label='Distance scale +/-7%, fixed stars/bins (not confidence band)')
for r in rows:ax.annotate(str(r['n']),(r['R_mean_kpc'],r['jeans_proxy_kms']),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8)
ax.set(xlabel='Distance from Galactic rotation axis (kpc)',ylabel='Circular-speed estimate or prediction (km/s)',
    title='542 training Cepheids: improvement, but a remaining shortfall')
ax.text(.02,.03,'Numbers above points are star counts.\nNo validation or test velocities opened; no new gravity parameters fitted.',transform=ax.transAxes,fontsize=9)
ax.grid(alpha=.18);ax.legend(loc='lower left',bbox_to_anchor=(0,.12),fontsize=8)
fig.savefig(HERE/'comparison.png',dpi=170);plt.close(fig)
table='\n'.join(f"| {r['R_mean_kpc']:.2f} | {r['n']} | {r['jeans_proxy_kms']:.2f} | {r['ordinary_vc_kms']:.2f} | {r['completion_vc_kms']:.2f} |" for r in rows)
report=f'''# Individual Cepheids in a shared coordinate frame

**The frozen additional field improves the training comparison, but remains too weak in every radial bin.** Starting from individual measurements does not remove the shortfall found with the published Cepheid curve. This is a conditional training diagnostic, not a validated reconstruction or a fresh holdout result.

| Model | Bin RMS discrepancy (km/s) | Mean predicted minus inferred speed (km/s) |
|---|---:|---:|
| Ordinary matter | {s['training_proxy_scores']['ordinary']['rms_kms']:.2f} | {s['training_proxy_scores']['ordinary']['bias_kms']:.2f} |
| Frozen additional field | {s['training_proxy_scores']['completion']['rms_kms']:.2f} | {s['training_proxy_scores']['completion']['bias_kms']:.2f} |

The comparison averages squared model circular speeds at the actual stellar radii before taking a square root. It uses equal weights across twelve bins; these RMS values are not uncertainty-weighted significances. The last two bins contain only five and six stars.

![Training Cepheid comparison](comparison.png)

## Distances and formulas

**Known empirical calibration, not our invention:** Ripepi et al. (2023), Table 3, gives absolute Wesenheit magnitude `M_W = a + b log10(P/day)`. Fundamental-mode coefficients are `(a,b)=(-2.744,-3.391)`; first-overtone coefficients are `(-3.224,-3.588)`. Its All Sky calibration uses Gaia parallax information. It does not require a cosmological redshift-distance relation.

**Known photometric definitions:** `w = G - 1.90(BP-RP)` and `d/kpc = 10^((w-M_W-10)/5)`, using intensity-averaged variable-star photometry. The color combination reduces extinction sensitivity under an adopted extinction law; it cannot guarantee zero dust systematics. The table's ABL scatter is not a magnitude scatter or a fractional-distance uncertainty. Coefficient covariance, metallicity, extinction-law variation and calibration selection are not fully modeled here.

**Known approximate Jeans equation:** `Vc^2 = <vphi^2> - <vR^2> [1 - R/Rd - 2R/Rsigma]`, with fixed `Rd=4 kpc` and `Rsigma=27.3 kpc` for this diagnostic. This assumes steady axisymmetric dynamics and prescribed radial profiles, omits the vertical mixed-moment term, and uses error-corrected second moments. It is not a one-star speed prediction. The additional field is the previously archived QUMOND-style conservative completion of our empirical response, not a derived photon-capture equation.

In plain language, distances convert angular movement into speed. The orbit approximation then estimates how strong gravity would have to be to support a population with those motions. Each step has assumptions; matching the final numbers alone cannot establish photon conversion.

## Selection and reserved data

Protocol commit `f04e0fe` preceded calculation of individual-star distances and velocities. A fixed hash of Gaia ID assigns 60/20/20 percent roles. Among measurement-ready stars with a supported single pulsation mode, there are 1,389 training, 433 validation and 443 test candidates. Only training coordinates and velocities were processed. The spatial and vertical-velocity cuts leave **542 training stars: 348 fundamental and 194 first overtone**.

No numerical selection was adjusted to match the authors' 903-star sample. Multimode objects are deferred, and we do not claim exact replication of their sky exclusions. The selection here uses `6 <= R <= 18 kpc`, absolute azimuth within 30 degrees of the Sun direction, `|z| <= 0.5 kpc`, and `|vz| <= 100 km/s`.

The role assignment is a random-star reserve, not a geographically independent survey. Gaia systematics and the previously exposed published curve remain shared. Validation and test candidate counts are availability counts, not final spatially selected counts. A stronger held-out claim requires a completed distance/selection/orbit method and frozen model choices before inspecting those outcomes.

## Numerical and assumption checks

- The common frame now matches the archived giant code exactly: Astropy `galcen_distance=8.2 kpc`, `z_sun=0.0208 kpc`, solar Cartesian velocity `[11.1,248,7.25] km/s`. Positive reported rotation follows the Galactic disk. After the first training run, we corrected the interpretation of 8.2 from cylindrical radius to three-dimensional distance, a **0.0264 pc** convention difference. This implementation amendment is recorded in the protocol; the original is preserved in Git.
- The exact coordinate transform to Feng's convention changes median individual rotation by {s['exact_coordinate_frame_diagnostics']['Feng']['median_delta_vphi_kms']:.3f} km/s; to Eilers' convention by {s['exact_coordinate_frame_diagnostics']['Eilers']['median_delta_vphi_kms']:.3f} km/s. These are individual-velocity shifts, not revised published circular-speed curves. Their listed solar radii are interpreted cylindrically in these comparison transforms.
- A first-order error calculation includes proper-motion correlation, radial-velocity error and an independent 7-percent distance-error scenario. It neglects sky-position uncertainty and distance-motion correlations. The error subtraction changes the bin proxy by at most **{maxnoise:.3f} km/s**.
- Rescaling all selected distances by 0.93 or 1.07, retaining nominal stars and bin membership, moves the inferred-speed proxy by at most **{maxshift:.3f} km/s**. This is a fixed-membership sensitivity, not a full refit or a confidence interval; model radii and sample membership would also need updating for a complete alternative-distance evaluation.
- Analytic-coordinate fixtures recover positions and velocities to below `1e-9` in their stated units. Monte Carlo propagation for 20 training stars with 4,000 draws each agrees with predicted variances to a median {100*verification['median_relative_variance_difference']:.2f}% and maximum {100*verification['max_relative_variance_difference']:.2f}%. This checks the implementation of the assumed error model, not its physical completeness.
- Refined and coarse field predictions differ by at most {max(s['max_coarse_refined_bin_change_kms'].values()):.4f} km/s. Frozen source/data hashes, training-only identifiers, and bin counts pass verification.

## Actual training comparison

The inferred column is an approximate population estimate from measurements, not a directly measured gravitational acceleration. All speeds are km/s.

| Mean radius (kpc) | Stars | Inferred from motions | Ordinary matter | Added-field total |
|---|---:|---:|---:|---:|
{table}

## Implication for the larger theory

This result supports continuing to investigate extra gravitational effects, but it does not specifically support photon-generated companions over other explanations. The current field's radial shortfall and earlier excess pull toward the disk need a common physical response, not independent adjustments for each observable. A fuller treatment of mass-model uncertainty and stellar populations is still needed before attributing all residuals to the field law.

Separately, a viable photon mechanism must produce the observed spectral shifts and event timing, transport companion energy consistently, and explain long-lived storage and lensing. Those unresolved requirements are not addressed by this stellar reduction. The total photon-supply budget remains deferred, not passed.

The next step is to test the distance and population assumptions on the training stars and connect the resulting orbit model to the same three-dimensional field used for bulge/disk comparisons. Do not open the reserved outcomes merely to tune away this discrepancy.

## Reproduce

Run `python -X utf8 run.py`, then `verify.py`, then `report.py` from this directory (or use their repository-relative paths). Raw and derived stellar rows are in the ignored cache; scripts, protocol, summaries and figure are versioned. Run.py only writes training-star velocity products. Reported results are derived from the frozen catalog in the protocol.

Sources: [Ripepi et al. 2023 calibration](https://openaccess.inaf.it/bitstream/20.500.12386/36524/1/aa43990-22-compressed.pdf); [Feng et al. 2026 methodology](https://academic.oup.com/mnras/article/546/2/stag011/8416425); [Eilers et al. coordinate and moment conventions](https://arxiv.org/html/1810.09466v2); [known QUMOND field structure](https://arxiv.org/abs/0911.5464).
'''
(HERE/'report.md').write_text(report,encoding='utf-8',newline='\n')
print(json.dumps(dict(max_distance_proxy_shift_kms=maxshift,max_noise_correction_kms=maxnoise),indent=2))
