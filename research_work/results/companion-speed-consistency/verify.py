from pathlib import Path
import hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
r=json.loads((HERE/'results.json').read_text());m=json.loads((HERE/'matched-speed-check.json').read_text())
assert len(r['runs'])==12 and len(m['runs'])==6
assert hashlib.sha256((HERE/'matched_speed.py').read_bytes()).hexdigest()==m['script_sha256']
for x in r['runs']:
 np.testing.assert_allclose(np.array(x['sampled_optical_factor_on_signal_path'])*x['wave_speed'],x['sampled_companion_to_photon_speed_ratio'],rtol=1e-12)
 assert x['max_relative_energy_error']<1e-8 and x['max_relative_momentum_error']<1e-8
for x in r['refinements']:assert max(x['mode_changes'].values())<1e-6
assert m['finite_difference_energy_gradient_max_error']<1e-7
for x in m['runs']:
 assert abs(x['frequency_relative_error'])<1e-5 and x['relative_energy_error']<1e-5
 np.testing.assert_allclose(x['continuum_photon_speed'],1/x['background_n'])
for N in [1.,1.2,2.]:
 pair=[x for x in m['runs'] if x['background_n']==N]
 errors=[abs(x['measured_phase_speed']/(1/N)-1) for x in pair]
 assert errors[1]<errors[0]/3.9
for name,value in r['source_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==value
v={'old_half_speed_reference_reproduced_by_driver':True,'speed_ratios_reconstructed':True,'fixed_speed_conservation_and_mode_checks':True,'proposed_energy_gradient_check':True,'linear_wave_speed_and_grid_dispersion':True,'source_hashes':True,'new_field_full_photon_interaction_tested':False,'lossless_companion_transport_proved':False,'astronomical_validation':False}
(HERE/'verification.json').write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(v,indent=2))
