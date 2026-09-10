"""Check frozen provenance, reserved-role separation and metric reconstruction."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
p=json.loads((HERE/'protocol.json').read_text());s=json.loads((HERE/'results.json').read_text())
bins=json.loads((HERE/'validation-bins.json').read_text())
for name,digest in p['frozen_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
assert hashlib.sha256((HERE/'protocol.json').read_bytes()).hexdigest()==s['protocol_sha256']
cache=ROOT/'research_work/data-cache/cepheid-stars'
d=pd.read_parquet(cache/'cepheid-validation-selected.parquet')
assert hashlib.sha256((cache/'cepheid-validation-selected.parquet').read_bytes()).hexdigest()==s['selected_catalog_sha256']
split=pd.read_parquet(cache/'cepheid-common-frame-split.parquet')
assert d.source_id.is_unique and len(d)==s['selected_validation_stars']
assert set(d.source_id)<=set(split.loc[split.role.eq('validation'),'source_id'])
assert not set(d.source_id)&set(split.loc[split.role.isin(['train','test']),'source_id'])
assert set(d.role)=={'validation'}
assert sum(r['n'] for r in bins)==len(d)
assert d.R_kpc.between(6,18).all() and d.z_kpc.abs().le(.5).all() and d.phi_deg.abs().le(30).all() and d.vz_kms.abs().le(100).all()
for r in bins:
    assert r['n']==int(d.bin.eq(r['bin']).sum())
    if r['n']<5:assert 'jeans_proxy_kms' not in r
    elif r['moment_valid']:
        independently=r['corrected_azimuthal_second_moment']-r['corrected_radial_second_moment']*(1-r['R_mean_kpc']/4-2*r['R_mean_kpc']/27.3)
        assert abs(independently-r['jeans_proxy_kms']**2)<1e-8
errors=[]
for model,score in s['scores'].items():
    residual=np.array([r[model+'_vc_kms']-r['jeans_proxy_kms'] for r in bins if r.get('moment_valid')])
    errors.append(abs(np.sqrt(np.mean(residual**2))-score['rms_kms']))
    assert errors[-1]<1e-10
    assert int((residual<0).sum())==score['underpredicted_bins']
training=json.loads((HERE.parent/'cepheid-common-frame/training-bins.json').read_text())
mass=json.loads((HERE.parent/'baryon-component-response/predictions.json').read_text())
matched={}
valid=set(s['evaluated_bins'])
for model,key in [('original_ordinary','ordinary_vc_kms'),('original_completion','completion_vc_kms')]:
    residual=np.array([r[key]-r['jeans_proxy_kms'] for r in training if r['bin'] in valid])
    matched[model]=dict(bins=len(residual),rms_kms=float(np.sqrt(np.mean(residual**2))))
mass_rotation=[r for r in mass if r['observable']=='rotation_training']
assert len(mass_rotation)==12
residual=np.array([r['balanced']-r['observed'] for i,r in enumerate(mass_rotation) if i in valid])
matched['balanced_completion']=dict(bins=len(residual),rms_kms=float(np.sqrt(np.mean(residual**2))))
out=dict(frozen_hashes_unchanged=True,validation_ids_disjoint_from_train_and_test=True,
    spatial_and_bin_count_checks=True,minimum_count_rule_obeyed=True,
    max_metric_reconstruction_error=max(errors),training_metrics_on_same_bin_indices=matched,
    qualification='Matched-bin training values use their own stars and mean radii; they are descriptive comparisons, not identical data or independent systematic-error realizations.')
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out,indent=2))
