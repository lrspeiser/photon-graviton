from pathlib import Path
import hashlib,json
import numpy as np
H=Path(__file__).resolve().parent;R=H.parent
result=json.loads((H/'results.json').read_text())
rows=json.loads((H/'predictions.json').read_text())
fine=json.loads((H/'predictions-refined.json').read_text())
selected=json.loads((R/'wave-lens-training/parameter-grid/selected-calibration.json').read_text())
for key in ['field_mass_eV_c2','source_to_stellar_mass_ratio']:
    assert result['protocol'][key]==selected[key]
for file,digest in selected['hashes'].items():
    assert hashlib.sha256((R/'wave-lens-training'/file).read_bytes()).hexdigest()==digest
assert result['input_sha256']['wave-lens-validation\\protocol.json']==hashlib.sha256((H/'protocol.json').read_bytes()).hexdigest()
by={(r['Name'],r['model']):r for r in rows}
names=set(json.loads((H/'preparation.json').read_text())['eligible_names'])
assert len(rows)==14 and {r['Name'] for r in rows}==names and all(r['role']=='validation' for r in rows)
ds=max(abs(r['sigma_pred_km_s']-by[r['Name'],r['model']]['sigma_pred_km_s']) for r in fine)
da=max(abs(r['theta_pred_arcsec']-by[r['Name'],r['model']]['theta_pred_arcsec']) for r in fine)
old=json.loads((R/'lens-training-pilot/predictions-refined-updated-profile-validation.json').read_text())
old={r['Name']:r for r in old if r['model']=='baryons' and r['seeing_fwhm_arcsec']==1.5}
mass_error=max(abs(r['sigma_pred_km_s']-r['sigma_observed_km_s']*np.sqrt(r['stellar_mass_Msun']/old[r['Name']]['mass_from_sigma_Msun']))
               for r in rows if r['model']=='baryons')
projection=max(r['projected_mass_lensing_identity_relative_error'] for r in rows)
assert ds<.01 and da<1e-5 and mass_error<1e-7 and projection<1e-6
scores={}
for model in ['baryons','stationary_wave']:
    ss=[r for r in rows if r['model']==model]
    scores[model]=float(.5*np.mean([np.log(r['sigma_pred_km_s']/r['sigma_observed_km_s'])**2+
                                  np.log(r['theta_pred_arcsec']/r['theta_SIE_arcsec'])**2 for r in ss]))
paired={}
for predicted,observed,label in [('sigma_pred_km_s','sigma_observed_km_s','dispersion'),('theta_pred_arcsec','theta_SIE_arcsec','lensing')]:
    paired[label]=sum(abs(by[n,'stationary_wave'][predicted]-by[n,'stationary_wave'][observed])<
                      abs(by[n,'baryons'][predicted]-by[n,'baryons'][observed]) for n in names)
summary={'frozen_calibration_hashes_unchanged':True,'eligible_validation_systems':len(names),
         'test_rows_scored':0,'three_resolution_check_systems':sorted({r['Name'] for r in fine}),
         'max_sigma_resolution_change_km_s':ds,'max_angle_resolution_change_arcsec':da,
         'max_baryon_scaling_identity_error_km_s':mass_error,'max_projected_mass_identity_relative_error':projection,
         'joint_log_error_scores':scores,'systems_with_smaller_absolute_error':paired,
         'interpretation':'Reused validation sample with frozen training parameters; not blind or full-uncertainty inference'}
(H/'verification.json').write_text(json.dumps(summary,indent=2)+'\n',newline='\n')
print(json.dumps(summary,indent=2))
