from pathlib import Path
import json,hashlib
import numpy as np
H=Path(__file__).resolve().parent;R=H.parent
protocol=json.loads((H/'protocol.json').read_text())
for file,digest in protocol['frozen_inputs'].items():assert hashlib.sha256((R/file).read_bytes()).hexdigest()==digest,file
rows=json.loads((H/'predictions.json').read_text());fine=json.loads((H/'predictions-refined.json').read_text())
prep=json.loads((H/'preparation.json').read_text());expected=set(prep['eligible_names'])
assert len(rows)==3*len(expected)==36 and all(r['role']=='test' for r in rows)
by={(r['Name'],r['model']):r for r in rows};models=['baryons','shared_stellar_mass','stationary_wave']
assert all({r['Name'] for r in rows if r['model']==m}==expected for m in models)
ds=max(abs(r['sigma_pred_km_s']-by[r['Name'],r['model']]['sigma_pred_km_s']) for r in fine)
da=max(abs(r['theta_pred_arcsec']-by[r['Name'],r['model']]['theta_pred_arcsec']) for r in fine)
projection=max(r['projected_mass_lensing_identity_relative_error'] for r in rows)
lam=protocol['parameters']['common_stellar_mass_multiplier']
scale=max(abs(by[n,'shared_stellar_mass']['sigma_pred_km_s']/by[n,'baryons']['sigma_pred_km_s']/np.sqrt(lam)-1) for n in expected)
assert ds<.01 and da<1e-5 and projection<1e-6 and scale<1e-12
individual={}
for n in sorted(expected):
    individual[n]={m:float(.5*(np.log(by[n,m]['sigma_pred_km_s']/by[n,m]['sigma_observed_km_s'])**2+
                               np.log(by[n,m]['theta_pred_arcsec']/by[n,m]['theta_SIE_arcsec'])**2)) for m in models}
scores={m:float(np.mean([v[m] for v in individual.values()])) for m in models}
# Post-test influence diagnostic only: keep every object in the reported primary score.
influence=[{'omitted_for_diagnostic_only':n,'wave_minus_mass_score':float(np.mean([
    v['stationary_wave']-v['shared_stellar_mass'] for k,v in individual.items() if k!=n]))} for n in sorted(expected)]
paired={}
for predicted,observed,label in [('sigma_pred_km_s','sigma_observed_km_s','dispersion'),('theta_pred_arcsec','theta_SIE_arcsec','lensing')]:
    paired[label]=sum(abs(by[n,'stationary_wave'][predicted]-by[n,'stationary_wave'][observed])<
                      abs(by[n,'shared_stellar_mass'][predicted]-by[n,'shared_stellar_mass'][observed]) for n in expected)
summary={'frozen_inputs_unchanged':True,'test_systems':len(expected),'excluded':prep['excluded'],
 'joint_scores':scores,'wave_minus_mass_primary_score':scores['stationary_wave']-scores['shared_stellar_mass'],
 'paired_wave_smaller_absolute_error_counts_vs_mass':paired,'individual_scores':individual,
 'post_test_leave_one_out_influence':influence,
 'checks':{'max_sigma_resolution_change_km_s':ds,'max_angle_resolution_change_arcsec':da,
           'max_projected_mass_identity_relative_error':projection,'mass_scaling_identity_relative_error':scale},
 'status':'First fixed-protocol test completed; these residuals are now exposed and cannot be reused as a fresh test after development'}
(H/'verification.json').write_text(json.dumps(summary,indent=2)+'\n',newline='\n')
print(json.dumps(summary,indent=2))
