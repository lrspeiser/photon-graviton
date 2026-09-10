"""Independently recompute report metrics from prediction rows and check artifacts."""
from pathlib import Path
import hashlib,json
import numpy as np
H=Path(__file__).resolve().parent; ROOT=H.parents[2]
load=lambda f:json.loads((H/f).read_text())
r=load('results.json'); red=load('redshift-predictions.json'); sp=load('galaxy-rotation-predictions.json'); mw=load('milky-way-predictions.json'); repair=load('repair-predictions.json')
for path,digest in load('input-manifest.json').items(): assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
assert len(red)==492 and len(sp)==6300 and len(mw)==648 and len(repair)==162
for model in ['fixed_alpha','refit_alpha','quadratic_history']:
    for split in ['train','validation','test']:
        rows=[x for x in red if x['model']==model and x['split']==split]
        value=np.sqrt(np.mean([(x['predicted_z']-x['observed_z'])**2 for x in rows]))*299792.458
        assert abs(value-r['redshift'][model]['scores'][split]['RMSE_kms'])<1e-8
for model in ['baryons','power']:
    for split in ['train','validation','test']:
        rows=[x for x in sp if x['model']==model and x['split']==split]
        names={x['galaxy'] for x in rows}
        value=np.sqrt(np.mean([np.mean([(x['predicted_kms']-x['observed_kms'])**2 for x in rows if x['galaxy']==n]) for n in names]))
        assert abs(value-r['sparc']['scores'][model][split]['galaxy_weighted_RMSE_kms'])<1e-9
for prof in load('capture-profiles.json'):
    assert np.all(np.asarray(prof['rho_Msun_kpc3'])>=0)
    assert np.all(np.diff(prof['mass_enclosed_Msun'])>=0)
assert r['gas_numerics']<.005
assert r['checks']['capture_energy_relative_error']<5e-5
assert all(v['optimizer_success'] for v in load('repair-results.json').values())
bulge=load('bulge-inputs.json'); assert len(bulge['fields'])==57
assert len({(x['survey'],x['field']) for x in bulge['fields']})==57
assert all(x['N']>0 and x['dispersion_kms']>0 and x['dispersion_error_kms']>0 for x in bulge['fields'])
assert len(load('bulge-response-predictions.json'))==3615
checks=load('bulge-results.json')['checks']
assert all(x['center_depth']<0 and x['center_acceleration_magnitude']<1e-9 for x in checks.values())
assert checks['equatorial_deposits']['center_vertical_curvature']>0
assert checks['upper_lower_deposits']['center_vertical_curvature']<0
html=(H/'comparison.html').read_text(encoding='utf-8'); assert 'const records=DATA;' not in html and 'Not derived' in html
assert all((H/f).is_file() for f in ['report.md','bulge-local-well.md','summary.png','bulge-summary.png'])
metadata={'status':'passed','scope':'Metric reproduction, data counts/hashes, nonnegative capture density and enclosed mass, gas/capture numerical checks, bulge symmetry/curvature and output presence',
          'input_prediction_rows':len(red)+len(sp)+len(mw)+len(repair),'bulge_field_rows':57,'bulge_synthetic_positions':3615,
          'code_and_protocol_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in H.iterdir() if p.suffix=='.py' or p.name in ['protocol.json','repair-amendment.md']}}
(H/'verification.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf-8',newline='\n')
print('Verified metrics, input hashes, numerical limits, and bulge response checks.')
