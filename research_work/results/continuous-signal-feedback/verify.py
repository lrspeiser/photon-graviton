from pathlib import Path
import hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent;r=json.loads((HERE/'results.json').read_text())
assert hashlib.sha256((HERE/'run.py').read_bytes()).hexdigest()==r['run_sha256']
for run in r['runs']+[r['coincident_split_control'],r['field_mode_control']]:
 assert run['max_relative_energy_error']<1e-8 and run['max_relative_momentum_error']<1e-8
 assert run['sampled_global_optical_lower_bound']>0
for run in [x for x in r['runs'] if x['quadrature_nodes']==32]:
 c=r['resolved_curves'][str(run['signal_energy'])]
 w=np.array(c['normalized_count_weights']);a=np.array(c['source_times']);b=np.array(c['detector_times']);fe=np.array(c['source_frequency']);fo=np.array(c['detector_frequency'])
 assert len(w)==32 and (w>0).all();np.testing.assert_allclose(w.sum(),1,atol=1e-14)
 assert (b>a).all() and (fe>0).all() and (fo>0).all()
 # Initial coordinates increase from back to front; crossing times decrease.
 assert (np.diff(a)<0).all() and (np.diff(b)<0).all()
 def width(t,ww):
  ww=ww/ww.sum();return np.sqrt(ww@((t-ww@t)**2))
 np.testing.assert_allclose(width(b,w)/width(a,w),run['count_weighted_duration_ratio'],rtol=1e-12)
 np.testing.assert_allclose(width(b,w*fo)/width(a,w*fe),run['energy_weighted_duration_ratio'],rtol=1e-12)
 np.testing.assert_allclose(w@(fe/fo),run['mean_carrier_stretch'],rtol=1e-12)
 A=np.column_stack([np.ones(len(a)),a]);beta=np.linalg.lstsq(A*np.sqrt(w)[:,None],b*np.sqrt(w),rcond=None)[0]
 np.testing.assert_allclose(beta[1],run['best_affine_slope'],rtol=1e-12)
 resid=np.sqrt(w@((b-A@beta)**2))/width(b,w)
 np.testing.assert_allclose(resid,run['affine_residual_rms_over_received_width'],atol=1e-12)
base=next(x for x in r['runs'] if x['signal_energy']==.01 and x['quadrature_nodes']==32)
for ctl in [r['coincident_split_control'],r['field_mode_control']]:
 for key in ['mean_carrier_stretch','energy_weighted_duration_ratio','field_energy_gain']:assert abs(ctl[key]-base[key])<1e-8
for c in r['quadrature_checks']:assert max(c['16_to_32_absolute_changes'].values())<1e-6
v={'source_hash':True,'conservation_checks':True,'crossing_order':True,'weighted_widths_and_carrier_means_reconstructed':True,'independent_weighted_linear_fit':True,'quadrature_refinement':True,'coincident_split_invariance':True,'field_mode_control':True,'astronomical_validation':False}
(HERE/'verification.json').write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(v,indent=2))
