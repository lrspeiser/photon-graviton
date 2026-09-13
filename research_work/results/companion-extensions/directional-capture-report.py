"""Summarize completed direction-dependent capture calculations without refitting."""
from pathlib import Path
import json,hashlib
import numpy as np
P=Path(__file__).resolve().parent
d=json.loads((P/'directional-capture-refined.json').read_text())
coarse=json.loads((P/'directional-capture-results.json').read_text())
check=json.loads((P/'directional-capture-check-results.json').read_text())
saved=json.loads((P.parent/'isotropic-galaxy-transfer/model-comparison-predictions.json').read_text())
ref={r['galaxy']:r for r in saved if r['model']=='companion_third'}
baryons={r['galaxy']:r for r in saved if r['model']=='baryons'}
summary=[];rng=np.random.default_rng(20260913)
for b in d['branches']:
 name=b['label'].rsplit('_raw',1)[0] if b['label'].endswith('_raw') else b['label'].removesuffix('_fixed_mass')
 raw=b['label'].endswith('_raw');ratios=[(1-b['mixture'])+b['mixture']*g['variants'][name]['mass_ratio'] if raw else 1. for g in d['geometry']]
 delta={}
 for split in ['train','validation','test']:
  changes=[]
  for row in b['predictions']:
   if row['split']!=split:continue
   y=np.array(row['observed_kms']);v=np.array(row['predicted_kms']);v0=np.array(ref[row['galaxy']]['predicted_kms'])
   changes.append(float(np.mean(np.log10(v/y)**2)-np.mean(np.log10(v0/y)**2)))
  values=np.array(changes);boot=values[rng.integers(0,len(values),(1000,len(values)))].mean(axis=1)
  delta[split]=dict(galaxies=len(changes),improved=int(np.sum(values<0)),mean_log_loss_difference=float(values.mean()),descriptive_bootstrap95=np.quantile(boot,[.025,.975]).tolist())
 matched={}
 if raw:
  for split in ['train','validation','test']:
   errors=[];logs=[]
   for row,ratio in zip(b['predictions'],ratios):
    if row['split']!=split:continue
    name0=row['galaxy'];vb=np.array(baryons[name0]['predicted_kms']);v0=np.array(ref[name0]['predicted_kms']);y=np.array(row['observed_kms'])
    v=np.sqrt(vb*vb+ratio*(v0*v0-vb*vb));errors.append(float(np.mean((v-y)**2)));logs.append(float(np.mean(np.log10(v/y)**2)))
   matched[split]=dict(RMSE_kms=float(np.sqrt(np.mean(errors))),log_RMS=float(np.sqrt(np.mean(logs))))
 summary.append(dict(label=b['label'],mixture=b['mixture'],scores=b['scores'],radial=b['radial'],stored_mass_fraction_range=[min(ratios),max(ratios)],paired_log_loss=delta,same_inventory_spherical_control=matched,max_prediction_refinement_change_kms=b['max_prediction_refinement_change_kms']))
assert all(x['max_prediction_refinement_change_kms']<.25 for x in summary),[(x['label'],x['max_prediction_refinement_change_kms']) for x in summary]
out=dict(scope='Frozen coarse-trained mixtures evaluated at refined resolution; no new fits; multiple exposed-data comparisons.',branches=summary,
    refinement_max_kms=max(x['max_prediction_refinement_change_kms'] for x in summary),source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [P/'directional-capture-results.json',P/'directional-capture-refined.json',P/'directional-capture-check-results.json']})
(P/'directional-capture-summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
lines=['# Directional capture: disk geometry changes the errors, but does not resolve both signs','',
'13 September 2026. Four geometries, each with two normalization choices. One shared mixture per branch, fit only to the original 89 training galaxies and frozen for 29 validation and 31 test galaxies. All samples exposed; no new observational validation.','',
'## Actual three-dimensional calculation','',
'The incoming far-field bath is isotropic. The capture region is oblate with m squared=R squared+z squared/q squared and opacity kappa=k0 W(n_z)/(1+m squared/a squared)^2. Top and bottom are symmetric, but paths through the disk plane differ from vertical paths. Unlike the earlier spherical calculation, deposited density depends on both radius and height, and its equatorial gravitational force includes nonspherical multipoles. This is an idealized capture geometry tied to the observed stellar disk scale, not a measured disk/bulge absorber or mapped nearby/distant stellar radiation sky.','',
'The four endpoint choices are q=0.5 with W=1 (flat isotropic capture); q=0.5 with W=3 n_z squared (vertical-sensitive capture); q=0.5 with W=1.5(1-n_z squared) (side-sensitive capture); and q=0.25 with W=1 (thinner isotropic capture). Each W has angular mean one. W represents an unproved directional interaction response, not a claim that more photons are observed arriving from those directions. The same W enters both local absorption and path attenuation.','',
'    rho_q(R,z) = 2 C0 eta(X) [1+m squared/a squared]^-2 <W exp(-tau)>','    tau = integral_along_incoming_path kappa dl','    rho_mix = (1-f_mix) rho_reference + f_mix rho_endpoint','',
'The exact-third eta and all original reference constants remain unchanged. Ray integration, exponential attenuation, angular quadrature, linear source superposition and the Legendre solution of Poisson gravity are known mathematics. The ellipsoidal opacity, W choices and storage interpretation are proposed hypotheses. Standard nonspherical potential methods are described in [Tremaine, potential theory lectures](https://www.ias.edu/sns/tremaine/lectures/ast513/potential).','',
'Two normalization choices separate effects. Raw branches retain the calculated capture-derived density amplitude. Fixed-mass branches rescale the endpoint density to the spherical inventory over the same numerical domain: this isolates redistribution, but is not derived from the transport law. That deterministic per-galaxy rescaling uses no measured rotation residual, yet it inherits the previously fitted reference inventory. The global mixture is an additional fitted parameter, not a derived capture fraction.','',
'## Frozen-mixture rotation scores','',
'Lower RMSE is better. Values are equal-galaxy km/s errors. Fitting minimizes mean squared log10 speed, so a km/s gain need not improve the actual objective. Eight searched families and one additional parameter must be acknowledged.','',
'| Branch | Fitted mixture | Training RMSE | Validation RMSE | Test RMSE |',
'|---|---:|---:|---:|---:|']
rs=d['reference_scores'];lines.append(f"| Spherical reference | 0 | {rs['train']['RMSE_kms']:.3f} | {rs['validation']['RMSE_kms']:.3f} | {rs['test']['RMSE_kms']:.3f} |")
for b in summary:lines.append(f"| {b['label']} | {b['mixture']:.5f} | {b['scores']['train']['RMSE_kms']:.3f} | {b['scores']['validation']['RMSE_kms']:.3f} | {b['scores']['test']['RMSE_kms']:.3f} |")
lines+=['| Existing fitted simple MOND | Not this mixture | 19.890 | 26.876 | 16.398 |','',
'## Inner versus outer signed error','',
'Negative means too slow; positive means too fast. Bins are below one, one to below three, and at least three stellar disk scale lengths. These are not necessarily the outer radii of the companion reservoir. Each galaxy is weighted equally within a bin.','',
'| Branch | Inner mean km/s | Middle mean km/s | Outer mean km/s |','|---|---:|---:|---:|',
'| Spherical reference | -6.726 | -7.486 | +10.412 |']
for b in summary:lines.append('| '+b['label']+' | '+' | '.join(f"{r['bias_kms']:+.3f}" for r in b['radial'])+' |')
lines+=['','The fixed-mass flattened branches can raise inner speeds, but also raise outer speeds that were already too high. The raw thinner branch lowers the outer excess while making the inner/middle shortage slightly worse. Thus direction-dependent capture can alter the ratio, but these simple variants do not repair both errors together. The trained mixtures remain small; the fit does not support replacing the entire spherical component by these flattened endpoints.','',
'The summary JSON also reports signed paired logarithmic-loss changes, per-split improvement counts and descriptive galaxy-bootstrap intervals. These do not include full distance, inclination, stellar-mass or selection uncertainty, and are not blind significance or model-selection evidence. Numerical score improvements are not proof of a physical capture mechanism.','',
'## Inventory and energy','',
'| Endpoint | Raw total / reference total, across galaxies |','|---|---:|']
for name in d['geometry'][0]['variants']:
 ratios=[g['variants'][name]['mass_ratio'] for g in d['geometry']];lines.append(f'| {name} | {min(ratios):.4f} to {max(ratios):.4f} |')
lines+=['','The raw endpoints store less, rather than creating or deleting energy. Relative to the spherical capture assumption, the uncaptured fraction remains in traveling companions in the bookkeeping interpretation; a common self-consistent source history and collision law are not derived. Mixture inventories are (1-f_mix)+f_mix times the raw ratio and are recorded separately. Fixed-mass branches preserve the numerical inventory by construction but require redistribution/support that has not been calculated. No cooling spectrum, stability, vertical stellar-motion or lensing success follows from this comparison.','',
'## Matched inventory-only control','',
'As a post-fit diagnostic, reduce the spherical density by exactly the same per-galaxy total-inventory ratio as each raw mixture, with no further fitting. This separates lower retention from directional shape. It is not an additional trained competitor or an untouched test. The summary JSON retains all split speed/log errors.','',
'| Raw branch | Directional validation / test RMSE | Same-inventory spherical validation / test RMSE |','|---|---:|---:|']
for b in summary:
 if b['same_inventory_spherical_control']:
  c=b['same_inventory_spherical_control'];lines.append(f"| {b['label']} | {b['scores']['validation']['RMSE_kms']:.3f} / {b['scores']['test']['RMSE_kms']:.3f} | {c['validation']['RMSE_kms']:.3f} / {c['test']['RMSE_kms']:.3f} |")
lines+=['',
'For the thin raw branch, the directional test RMSE is 23.334 km/s versus 23.243 for its same-inventory spherical control. The directional branch has slightly better logarithmic scores than that control, so the verdict depends on the stated metric; it does not establish a general advantage. Against the original reference, its validation logarithmic score slightly worsens, and descriptive paired log-loss intervals span zero in all three splits. The speed gain therefore cannot be cleanly attributed to a successful directional mechanism.','',
'## Verification and limitations','',
f"Coarse resolution was {coarse['resolution']}; refined resolution was {d['resolution']}. Mixtures were not refitted on comparison data or at higher resolution. The largest frozen-mixture speed change on refinement is {out['refinement_max_kms']:.6f} km/s. The maximum reference force-reconstruction speed discrepancy is {max(g['reference_max_difference_kms'] for g in d['geometry']):.6f} km/s; the mixture baseline itself uses the exact archived spherical predictions.", '',
f"Independent direct numerical ray integrals agree with the analytic expression within {check['max_ray_integral_relative_error']:.3g} relative error over {check['ray_cases']} cases. Smooth oblate density fields checked against a separate homoeoid integral agree to {max(x['relative_error'] for x in check['ellipsoidal_force_cases'] if x['refined']):.3g} relative error at refined resolution. These validate the numerical machinery, not the proposed physics.",'',
f"The grid extends from 1e-6 to 3000 reference capture radii. Positive exterior density bounds give omitted mass at most {100*max(v['tail_upper_fraction_of_grid_mass'] for g in d['geometry'] for v in g['variants'].values()):.4f}% of the grid inventory. This finite domain is a numerical approximation, not an assumed age or universe size. Fixed-mass normalization is exact only within this common finite domain.",'',
'The model is axisymmetric and reflection symmetric. It does not model the Milky Way bar, independent bulge opacity, individual stellar source positions or environmental asymmetry. Direction-sensitive cross sections and a long-lived deposited state are assumptions. Rotation was computed in the equatorial plane; predicting motion above/below the disk and light bending requires additional calculations. No outward-force term was clipped to simulate a positive circular speed.','',
'## Consequence','',
'Keep the original shared law as the reference. Assess the small raw thinner-capture improvement against the inventory-only control, rather than attributing any gain automatically to direction. None of these fitted mixtures repairs both radial signs. Further work should specify a physically motivated disk/bulge capture relation or use vertical-motion constraints to distinguish geometries, rather than enlarge this sweep solely to find a lower residual.','',
'The user also proposed making capture respond to stronger gravity deeper in the galaxy. This is distinct from the present prescribed opacity profile, although that profile already has greater local opacity toward the center. A gravitational-depth rule must distinguish potential depth from local acceleration, calculate the incoming energy surviving interception along each path, and include feedback from deposits. A deep central well can have small local acceleration by symmetry. Greater central deposition still contributes to outer gravity. Earlier binding-feedback branches should be reviewed before constructing a directional version, rather than labeling all depth feedback untested or repeating an old prescription. No depth-feedback experiment was performed in this run.','',
'## Reproduction','',
'Run directional-capture.py, directional-capture.py --refine, directional-capture-check.py and directional-capture-report.py in that order. Protocol, source hashes, density-derived force arrays, inventories, scans, frozen predictions, refinement and summary are saved alongside this report. The paper supplement indexes this result; PDF v1.5 predates it.']
(P/'directional-capture-report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([{k:b[k] for k in ['label','mixture','scores','stored_mass_fraction_range','max_prediction_refinement_change_kms']} for b in summary],indent=2))
