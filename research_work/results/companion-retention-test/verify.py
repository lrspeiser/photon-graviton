from pathlib import Path
import json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;r=json.loads((HERE/'results.json').read_text())
assert hashlib.sha256((HERE/'run.py').read_bytes()).hexdigest()==r['script_sha256']
assert len(r['runs'])==8
for x in r['runs']+[r['double_box_control']]:
 np.testing.assert_allclose(x['final_wave_energy']/x['initial_wave_energy'],x['wave_energy_retained_fraction'],rtol=1e-12)
 loss=x['initial_wave_energy']-x['final_wave_energy']
 assert abs(loss-x['background_energy_gain'])/x['initial_wave_energy']<1e-6
 assert x['energy_error_over_initial_wave']<1e-6
 assert abs(x['integrated_background_work']-x['background_energy_gain'])/x['initial_wave_energy']<1e-5
 np.testing.assert_allclose(x['history']['wave_energy_fraction'][-1],x['wave_energy_retained_fraction'])
 np.testing.assert_allclose(x['homogeneous_adiabatic_retention_prediction'],1/x['mean_field_final'],rtol=1e-10)
 assert x['minimum_field']>0
for x in r['refinements']:assert max(x['absolute_changes'].values())<1e-5
assert r['double_box_retention_change']<1e-4
v={'script_hash':True,'energy_retention_ratios':True,'wave_loss_equals_background_gain':True,'work_integral_independent_of_energy_subtraction':True,'error_scaled_to_initial_wave_energy':True,'grid_and_box_controls':True,'lossless_retention_established':False,'astronomical_validation':False}
(HERE/'verification.json').write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(v,indent=2))
