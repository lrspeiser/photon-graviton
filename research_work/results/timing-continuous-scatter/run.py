"""Development profiles on the two previously exposed boundary cases."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.optimize import minimize
from integration import log_likelihood
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
protocol=json.loads((HERE.parent/'timing-boundary-diagnosis/protocol.json').read_text(encoding='utf-8'))
results=[]
# Development-only domain keeps all population means inside the finite grid.
# It is narrower than the original experiment and does not replace its fit.
for case in protocol['cases']:
 path=ROOT/'research_work/generated/timing-coverage-calibration'/(case['label']+'.npz')
 assert hashlib.sha256(path.read_bytes()).hexdigest()==case['array_sha256']
 with np.load(path) as data: x=data['log_width'];logs=data['event_log_likelihoods'];z=data['redshift']
 logs=logs-logs.max(axis=1)[:,None]
 profiles=[]
 for step in [1,2]:
  assert (len(x)-1)%step==0
  for sigma in [0.,.005,.01,.02,.03,.1]:
   objective=lambda q:-log_likelihood(x[::step],logs[:,::step],z,q[0],q[1],sigma)
   starts=[[case['fit']['a'],case['fit']['b']],[np.log(30),0],[np.log(30),1]]
   fits=[minimize(objective,q,method='L-BFGS-B',bounds=[np.log([5,80]),[-.5,1.5]],options={'ftol':1e-11,'maxiter':600}) for q in starts]
   good=[f for f in fits if f.success and np.isfinite(f.fun)]
   assert good
   best=min(good,key=lambda f:f.fun)
   profiles.append(dict(width_nodes=len(x[::step]),sigma=sigma,a=float(best.x[0]),b=float(best.x[1]),negative_centered_log_likelihood=float(best.fun),successful_starts=len(good)))
 results.append(dict(label=case['label'],input_sha256=case['array_sha256'],original_b=case['fit']['b'],profiles=profiles))
 print(json.dumps(dict(label=case['label'],profiles=profiles)),flush=True)
(HERE/'results.json').write_text(json.dumps(dict(scope='Selected exposed cases, separate continuous-interpolant revision; no replacement coverage claim',cases=results),indent=2)+'\n',encoding='utf-8',newline='\n')
