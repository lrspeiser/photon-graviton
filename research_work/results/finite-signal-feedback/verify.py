from pathlib import Path
import json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
r=json.loads((HERE/'results.json').read_text())
assert len(r['runs'])==10
for x in r['runs']+[r['color_control'],r['analytic_homogeneous_control']]:
 assert x['maximum_relative_energy_error']<1e-8 and x['maximum_relative_momentum_error']<1e-8
 assert x['sampled_global_optical_lower_bound']>0
 assert abs(x['driver_energy_lost']+x['signal_energy_lost']-x['field_energy_gain'])<1e-9
 assert abs(x['field_energy_gain']-(x['final_field_energy']-x['initial_field_energy']))<1e-12
 for i,(a,b) in enumerate(zip(x['packets'][:-1],x['packets'][1:])):
  expected=(b['detector_crossing_time']-a['detector_crossing_time'])/(b['source_crossing_time']-a['source_crossing_time'])
  np.testing.assert_allclose(expected,x['intervals'][i]['reference_duration_stretch'],rtol=1e-12)
  expected_q=(b['detector_q1_clock']-a['detector_q1_clock'])/(b['source_q1_clock']-a['source_q1_clock'])
  np.testing.assert_allclose(expected_q,x['intervals'][i]['q1_duration_stretch'],rtol=1e-12)
 for p in x['packets']:
  np.testing.assert_allclose(p['reference_frequency_stretch'],p['source_frequency']/p['detector_frequency'],rtol=1e-12)
  assert abs(p['independent_frozen_field_jacobian']-p['reference_frequency_stretch'])<2e-8
  if 'feedback_stretch_increment' in p:
   np.testing.assert_allclose(p['reference_frequency_stretch']-p['weak_background_stretch_at_same_source_time'],p['feedback_stretch_increment'],atol=1e-12)
for x in r['mode_refinement']:assert x['max_stretch_change']<1e-6
assert r['fixed_energy_color_stretch_difference']<1e-8
for p in r['analytic_homogeneous_control']['packets']:
 assert abs(p['reference_frequency_stretch']-np.exp(.1))<1e-8 and abs(p['q1_frequency_stretch']-1)<1e-8
for name,value in r['source_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==value
checks={'energy_and_momentum':True,'global_spatial_positivity_at_sampled_times':True,'intervals_and_carrier_ratios_reconstructed':True,'matched_emission_time_controls':True,'fixed_energy_color_control':True,'mode_refinement':True,'analytic_homogeneous_clock_control':True,'source_hashes':True,'astronomical_validation':False}
(HERE/'verification.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(checks,indent=2))
