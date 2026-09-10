"""Summarize the specified absorption test without promoting it to a data fit."""
from pathlib import Path
import json
import hashlib
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
d=json.loads((HERE/'results.json').read_text())
for p,h in d['source_hashes'].items():assert digest(ROOT/p)==h
for name,fit in d['fits'].items():
    assert digest(ROOT/f'research_work/data-cache/baryonic-shielding/{name}.npz')==fit['cache_sha256']
source=json.loads((HERE.parent/'baryon-attached-deposits/results.json').read_text())
target=np.array([r['required_extra_per_baryon'] for r in source['rows']])
fine=d['fits']['fine'];pred=np.exp(fine['log_prediction'])
residual=np.log(pred/target)
assert abs(np.sqrt(np.mean(residual**2))-fine['rms_log_residual'])<1e-12
records=[]
for i,r in enumerate(source['rows']):
    records.append(dict(R_kpc=r['R_kpc'],z_kpc=r['z_kpc'],phi_rad=r['phi_rad'],
                        required_loading=float(target[i]),predicted_loading=float(pred[i]),
                        predicted_over_required=float(pred[i]/target[i])))
(HERE/'predictions.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8',newline='\n')
worst=int(np.argmax(abs(residual)))
assessment=dict(maximum_mismatch_factor=float(np.exp(abs(residual).max())),worst=records[worst],
                fraction_within_factor_two=float(np.mean(abs(residual)<=np.log(2))),
                inferred_from_observations=False,holdouts_opened=False,
                results_sha256=digest(HERE/'results.json'))
(HERE/'assessment.json').write_text(json.dumps(assessment,indent=2)+'\n',encoding='utf-8',newline='\n')
lines=[
    '# Can baryonic shielding create the required uneven deposit pattern?', '',
    '**The tested shielding rule improves the match to the target source shape but does not reproduce it.** Its exponentiated RMS logarithmic residual falls from 6.01 for constant loading to 4.49 with shielding. This is an equal-weight comparison of model values at 240 grid points, not a measured uncertainty, stellar likelihood or independent observational success.', '',
    'This follows the constant-loading audit by giving the incoming companions a specified spatial transport law. It retains one common absorption coefficient and one common normalization; it does not assign a separate loading to each location.', '',
    '## Hypothesis and formula provenance', '',
    'Assume a stationary, isotropic companion intensity at a spherical boundary, straight propagation, no internal emission or scattering, and absorption coefficient per unit length `beta(x)=kappa*rho_b(x)`. Deposits remain attached to their local ordinary matter, and all locations share the same exposure duration. Ordinary matter is held fixed. These are optional physical assumptions for this test, not consequences of the user\'s general deposit concept.', '',
    'Known absorption transport gives:', '',
    '`Sigma(x,n) = integral_from_x_to_boundary rho_b(x-s*n) ds`',
    '`I(x,n) = I_boundary * exp[-kappa*Sigma(x,n)]`',
    '`q_model(x) = C * <exp[-kappa*Sigma(x,n)]>_directions`.', '',
    'Here Sigma is ordinary-matter column density, kappa is absorption cross section per ordinary mass, and q is the extra equivalent source per ordinary mass. The angular brackets average the incident directions. The exponential and angular transport are known mathematics; identifying the absorbing sector as companions and specifying beta are hypothetical. No novel fundamental equation is claimed.', '',
    'C contains the accumulated exposure, capture normalization and assumed conversion from retained energy to an ordinary Newtonian source. Fitting it independently is a source-shape diagnostic; it does not measure the available photon energy or break the source-to-gravity normalization degeneracy. The constant-loading comparison is the earlier shape model, not physical capture at exactly zero absorption with a finite fixed exposure.', '',
    '## Numerical setup', '',
    'The ordinary density uses the same analytic bar, nuclear, stellar-disk, gas-disk and softened-central components as the preceding audit. Their sum is checked against all 240 archived ordinary-density values. Rays end at a 30-kpc sphere. A 60-kpc boundary is an explicit sensitivity test; these radii are modeling choices, not measured halo boundaries.', '',
    'Angular quadrature uses 128 directions for the coarse run and 512 for the finer runs. Radial quadrature uses 64 and 128 nodes on each segment, splitting the path at the disk midplane before integration. The zero-optical-depth angular normalization is verified. The search scans kappa from 1e-12 to 1e-5 kpc squared per solar mass, then refines around the lowest scanned interval. Reported optima lie inside that scan; this is not a proof of a global optimum over every possible capture law.', '',
    'For each kappa the common log-normalization is solved analytically by the mean log target minus mean log attenuation. The loss is the mean squared log residual across the selected grid. It does not use observational error weights, volume weights, a selection function or held-out stars.', '',
    '| Calculation | Best kappa (kpc²/Msun) | RMS log residual | exp(RMS log residual) |',
    '|---|---:|---:|---:|',
    f"| Constant loading | — | {d['constant_loading_rms_log']:.4f} | {np.exp(d['constant_loading_rms_log']):.4f} |",
]
for name,v in d['fits'].items():
    lines.append(f"| {name} | {v['kappa_kpc2_per_Msun']:.6g} | {v['rms_log_residual']:.4f} | {v['geometric_scatter_factor']:.4f} |")
lines += ['',
    f"Only {100*assessment['fraction_within_factor_two']:.1f}% of sampled loadings fall within a factor of two of the target. The largest mismatch is a factor of {assessment['maximum_mismatch_factor']:.1f}. Those summaries depend on the declared grid and loss function; they are not rejection probabilities.", '',
    '## Resolution and boundary checks', '',
    'Holding the coarse fitted parameters fixed, finer quadrature changes predicted log loadings by at most 0.01830, roughly 1.85% multiplicatively. Refitting barely changes kappa or the overall residual. Moving the boundary to 60 kpc leaves the fitted residual factor at approximately 4.489. Thus the tested resolution and boundary changes are much smaller than the remaining factor-level mismatch. More refinement would be necessary for precision capture-rate inference; no such precision claim is made.', '',
    'The target effective source itself has up to 5.4% fine/finer density changes in the previous audit. Ordinary-matter uncertainties, more general outer illumination and ray bending are not covered by these numerical checks.', '',
    '## What remains physically unresolved', '',
    'Shielding can reduce the flux per unit ordinary matter in the interior, giving the desired direction of uneven loading. In this specified geometry, a single absorption coefficient does not supply the full required spatial pattern. Merely adjusting the global exposure cannot fix that shape mismatch.', '',
    'A next physical candidate would need an independently specified change, such as capture dependence on the gravitational environment, internally generated companions with a modeled stellar-emission history, a supported reservoir separate from ordinary matter, or a different gravitational response. Giving each region an arbitrary cross section would make agreement a definition rather than a prediction.', '',
    'This stationary snapshot also omits deposition-induced motion, recoil, drag and changing baryonic mass. It cannot serve as a self-consistent history when accumulated loading is substantial. It does not derive the boundary bath from photon redshift, test total photon supply, or show that ordinary clocks and event timing work. The broader theory remains unvalidated; no holdout has been opened.', '',
    '## Reproduction', '',
    'Run `run.py`, then `report.py`. Input hashes, both quadrature levels, the boundary sensitivity, fitted parameters and frozen-parameter comparisons are retained. `predictions.json` records the required and predicted loading at every probe. Large column-density arrays remain in the ignored cache. The existing target potential, source and stellar catalogs are unchanged.', '',
]
(HERE/'report.md').write_text('\n'.join(lines),encoding='utf-8',newline='\n')
print(json.dumps(assessment,indent=2))
