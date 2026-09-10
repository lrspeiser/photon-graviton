from pathlib import Path
import json
import numpy as np
H=Path(__file__).resolve().parent
base=json.loads((H/'predictions.json').read_text())
fine=json.loads((H/'predictions-refined.json').read_text())
by={(r['Name'],r['model']):r for r in base}
ds=max(abs(r['sigma_pred_km_s']-by[r['Name'],r['model']]['sigma_pred_km_s']) for r in fine)
da=max(abs(r['theta_pred_arcsec']-by[r['Name'],r['model']]['theta_pred_arcsec']) for r in fine)
old=json.loads((H.parent/'lens-training-pilot/predictions-refined-updated-profile.json').read_text())
old={r['Name']:r for r in old if r['model']=='baryons' and r['seeing_fwhm_arcsec']==1.5}
errors=[]
for r in base:
    if r['model']!='baryons':continue
    expected=r['sigma_observed_km_s']*np.sqrt(r['stellar_mass_Msun']/old[r['Name']]['mass_from_sigma_Msun'])
    errors.append(abs(expected-r['sigma_pred_km_s']))
assert len(base)==64 and all(r['role']=='training' for r in base)
projection=max(r['projected_mass_lensing_identity_relative_error'] for r in base)
assert projection<1e-6
assert ds<.01 and da<1e-5 and max(errors)<1e-7
result={'three_resolution_check_systems':sorted({r['Name'] for r in fine}),
        'max_sigma_resolution_change_km_s':ds,'max_angle_resolution_change_arcsec':da,
        'max_baryon_mass_scaling_identity_error_km_s':max(errors),
        'max_projected_mass_lensing_identity_relative_error':projection,
        'scope':'Numerical checks, not observational validation'}
(H/'verification.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))
