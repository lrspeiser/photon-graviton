from pathlib import Path
import json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
p=json.loads((HERE/'protocol.json').read_text());r=json.loads((HERE/'results.json').read_text())
assert r['protocol_sha256']==hashlib.sha256((HERE/'protocol.json').read_bytes()).hexdigest()
for path,digest in p['input_hashes'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
assert {x['target_name'] for x in r['rows']}=={x['target_name'] for x in p['targets']} and len(r['rows'])==8
assert sum(x['below_original_calibration_distance_range'] for x in r['rows'])==7
for x in r['rows']:
 assert x['source_position_separation_arcsec']<30 and x['haynes_class']==1
 np.testing.assert_allclose(x['observed_cmb_z'],(1+x['heliocentric_optical_velocity_kms']/p['c_kms'])*x['solar_to_cmb_factor']-1,rtol=1e-12)
 np.testing.assert_allclose(x['exponential_predicted_z'],np.expm1(p['alpha_per_mpc']*x['distance_mpc']),rtol=1e-12)
 np.testing.assert_allclose(x['linear_predicted_z'],p['linear_alpha_per_mpc']*x['distance_mpc'],rtol=1e-12)
for model in ['exponential','linear']:
 residual=np.array([p['c_kms']*(x[model+'_predicted_z']-x['observed_cmb_z']) for x in r['rows']])
 np.testing.assert_allclose(np.sqrt(np.mean(residual**2)),r['scores'][model]['rms_kms'])
 np.testing.assert_allclose(np.mean(residual),r['scores'][model]['bias_kms'])
 np.testing.assert_allclose(np.mean(abs(residual)),r['scores'][model]['mae_kms'])
v={'frozen_protocol_and_input_hashes':True,'exact_eight_target_membership':True,'position_and_quality_checks':True,'frame_and_prediction_reconstruction':True,'metric_reconstruction':True,'refitted_parameters':False,'certified_independent_validation':False,'outcomes_exposed':True}
(HERE/'verification.json').write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(v,indent=2))
