"""Check table integrity, frozen inputs, scores, and spatial refinement."""
from pathlib import Path
import importlib.util
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
protocol=json.loads((HERE/'protocol.json').read_text())
for path,expected in protocol['hashes'].items():
    assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected
obs=json.loads((HERE/'observations.json').read_text())['rows']
pred=json.loads((HERE/'predictions.json').read_text())
scores=json.loads((HERE/'results.json').read_text())['scores']
assert len(obs)==12 and len(pred)==36
for model in protocol['models']:
    rows=[r for r in pred if r['model']==model['name']]
    assert [(r['R_kpc'],r['vc_kms'],r['error_kms']) for r in rows]==[(r['R_kpc'],r['vc_kms'],r['error_kms']) for r in obs]
    residual=np.array([r['predicted_kms']-r['vc_kms'] for r in rows])
    assert abs(np.sqrt(np.mean(residual**2))-scores[model['name']]['rmse_kms'])<1e-10
source=ROOT/'research_work/results/conservative-field-completion/run.py'
spec=importlib.util.spec_from_file_location('frozen_field_check',source)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
b=mod.Baryons(False);c=mod.Completion(b,False,200)
R=np.array([r['R_kpc'] for r in obs]);z=np.zeros(len(R))
ab=mod.force(b,R,z)[:,0];ac=mod.force(c,R,z)[:,0]
changes={}
for model in protocol['models']:
    lam=model['lambda_mass'];a=lam*ab
    if model['name']!='ordinary':a+=lam**mod.P*ac
    low=np.sqrt(-R*a)
    high=np.array([r['predicted_kms'] for r in pred if r['model']==model['name']])
    changes[model['name']]=float(np.max(abs(low-high)))
    assert changes[model['name']]<.1
out=dict(all_protocol_hashes_unchanged=True,all_12_rows_retained_for_all_three_models=True,
    scores_recomputed=True,max_coarse_refined_prediction_change_kms=changes,
    refinement_is_numerical_check_not_refit=True)
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n')
print(json.dumps(out,indent=2))
