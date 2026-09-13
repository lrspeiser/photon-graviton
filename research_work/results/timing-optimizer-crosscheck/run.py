"""Independent derivative-free checks on the exposed boundary-case development fits."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
from scipy.optimize import minimize
H=Path(__file__).resolve().parent;R=H.parents[2]
sys.path.insert(0,str(H.parent/'timing-boundary-capable'))
from run import outside,averaged_likelihood
source=H.parent/'timing-boundary-capable/results.json';prior=json.loads(source.read_text(encoding='utf-8'))
rows=[]
for c in prior['cases']:
 p=R/'research_work/generated/timing-coverage-calibration'/(c['label']+'.npz')
 assert hashlib.sha256(p.read_bytes()).hexdigest()==c['array_sha256']
 with np.load(p) as d:x=d['log_width'];logs=d['event_log_likelihoods'];z=d['redshift']
 L=np.exp(logs-logs.max(axis=1)[:,None]);lz=np.log1p(z)
 def objective(q):
  mu=q[0]+q[1]*lz
  if np.any(outside(x,mu,q[2])>1e-6):return 1e12
  try:v=averaged_likelihood(x,L,mu,q[2])
  except ArithmeticError:return 1e12
  if np.any(v<=0):return 1e12
  return float(-np.log(v).sum())
 old=c['revised']['best'];trials=[]
 for fixed_zero in [False,True]:
  bounds=[np.log([5,100]),[-2,3]]+([] if fixed_zero else [[0,.6]])
  q=np.array(old['parameters'][:2] if fixed_zero else old['parameters'])
  f=minimize((lambda a:objective([*a,0])) if fixed_zero else objective,q,method='Powell',bounds=bounds,options={'ftol':1e-12,'xtol':1e-10,'maxiter':1000})
  params=[*f.x,0.] if fixed_zero else f.x.tolist()
  trials.append(dict(method='Powell',fixed_zero=fixed_zero,parameters=params,objective=float(f.fun),success=bool(f.success),message=str(f.message)))
 valid=[v for v in trials if v['success'] and v['objective']<1e11]
 best=min([old]+valid,key=lambda v:v['objective'])
 rows.append(dict(label=c['label'],array_sha256=c['array_sha256'],previous_best=old,trials=trials,best=best,objective_improvement=old['objective']-best['objective'],b_change=best['parameters'][1]-old['parameters'][1]))
 print(json.dumps({k:rows[-1][k] for k in ['label','objective_improvement','b_change']}),flush=True)
out=dict(scope='Optimization development only, exposed synthetic cases; no recalibrated coverage',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),cases=rows)
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
