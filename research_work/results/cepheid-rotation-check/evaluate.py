"""Evaluate only the predeclared models against the new catalog table."""
from pathlib import Path
import importlib.util
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
protocol=json.loads((HERE/'protocol.json').read_text())
for relative,expected in protocol['hashes'].items():
    assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==expected,relative
source=ROOT/'research_work/results/conservative-field-completion/run.py'
spec=importlib.util.spec_from_file_location('frozen_completion',source)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
observations=json.loads((HERE/'observations.json').read_text())
rows=observations['rows'];R=np.array([r['R_kpc'] for r in rows]);z=np.zeros(len(R))
print('Building frozen field for Cepheid table radii',flush=True)
baryons=mod.Baryons(protocol['prediction_resolution']['refined'])
completion=mod.Completion(baryons,True,protocol['prediction_resolution']['response_outer_radius_kpc'])
ab=mod.force(baryons,R,z)[:,0];ac=mod.force(completion,R,z)[:,0]
assert np.all(ab<0) and np.all(ac<0)
observed=np.array([r['vc_kms'] for r in rows]);error=np.array([r['error_kms'] for r in rows])
predictions=[];scores={}
for model in protocol['models']:
    lam=model['lambda_mass']
    a=lam*ab
    if model['name']!='ordinary':a=a+lam**mod.P*ac
    pred=np.sqrt(-R*a);residual=pred-observed
    scores[model['name']]=dict(rows=len(rows),rmse_kms=float(np.sqrt(np.mean(residual**2))),
        bias_kms=float(np.mean(residual)),median_prediction_observation_ratio=float(np.median(pred/observed)),
        rms_residual_over_published_bootstrap_error=float(np.sqrt(np.mean((residual/error)**2))),
        within_one_published_error=int(np.sum(abs(residual)<=error)),
        underpredicted_rows=int(np.sum(residual<0)))
    for i,row in enumerate(rows):
        predictions.append(dict(model=model['name'],lambda_mass=lam,**row,
            predicted_kms=float(pred[i]),residual_kms=float(residual[i]),
            residual_over_published_bootstrap_error=float(residual[i]/error[i])))
result=dict(protocol_commit='f7b2980',protocol_sha256=hashlib.sha256((HERE/'protocol.json').read_bytes()).hexdigest(),
    observations_sha256=hashlib.sha256((HERE/'observations.json').read_bytes()).hexdigest(),
    all_frozen_hashes_match=True,parameters_refitted=False,scores=scores,
    scope='First frozen-model evaluation of this Cepheid table; shared Gaia systematics and model assumptions prevent a fully independent raw-data claim',
    now_exposed=True)
for name,obj in [('results.json',result),('predictions.json',predictions)]:
    (HERE/name).write_text(json.dumps(obj,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2),flush=True)
