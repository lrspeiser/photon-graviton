"""Numerical profile-interval sensitivity on two exposed synthetic samples."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
from scipy.optimize import minimize,brentq
from scipy.stats import chi2
H=Path(__file__).resolve().parent;R=H.parents[2]
sys.path.insert(0,str(H.parent/'timing-boundary-capable'))
from run import outside,averaged_likelihood
source=H.parent/'timing-boundary-grid-check/results.json';data=json.loads(source.read_text(encoding='utf-8'))
rows=[];cut=float(chi2.ppf(.95,1))
for c in data['cases']:
 for kind,folder in [('coarse','timing-coverage-calibration'),('fine','timing-scatter-grid-refinement')]:
  p=R/'research_work/generated'/folder/(c['label']+'.npz')
  assert hashlib.sha256(p.read_bytes()).hexdigest()==c[kind+'_sha256']
  with np.load(p) as d:x=d['log_width'];logs=d['event_log_likelihoods'];z=d['redshift']
  L=np.exp(logs-logs.max(axis=1)[:,None]);lz=np.log1p(z)
  best=c[kind+'_fit']['best'];a,b,s=best['parameters'];evaluations=[]
  def profile(t):
   def objective(q):
    try:v=averaged_likelihood(x,L,q[0]+t*lz,q[1])
    except ArithmeticError:return 1e12
    return float(-np.log(v).sum()) if np.all(v>0) else 1e12
   trials=[]
   for q in [[a,s],[a,0],[np.log(30),.1]]:
    f=minimize(objective,q,method='SLSQP',bounds=[np.log([5,100]),[0,.6]],constraints={'type':'ineq','fun':lambda q:1e-6-outside(x,q[0]+t*lz,q[1])},options={'ftol':1e-10,'maxiter':600})
    trials.append(dict(parameters=f.x.tolist(),objective=float(f.fun),success=bool(f.success),outside=float(outside(x,f.x[0]+t*lz,f.x[1]).max())))
   good=[v for v in trials if v['success'] and v['outside']<=1.00001e-6 and v['objective']<1e11]
   if not good:raise RuntimeError('No profile fit')
   lr=2*(min(v['objective'] for v in good)-best['objective'])
   if lr < -2e-6:raise RuntimeError('Profile beats unrestricted fit')
   evaluations.append(dict(b=float(t),lr=float(lr),trials=trials))
   return max(0.,lr)-cut
  endpoints=[]
  for direction in [-1,1]:
   for step in [.2,.4,.8,1.6,3.2]:
    endpoint=float(np.clip(b+direction*step,-2,3))
    if profile(endpoint)>0:break
   else:raise RuntimeError('Interval not bracketed')
   root=brentq(profile,min(b,endpoint),max(b,endpoint),xtol=1e-7)
   endpoints.append(float(root))
  rows.append(dict(label=c['label'],grid=kind,nodes=len(x),array_sha256=c[kind+'_sha256'],best=best,nominal_95_endpoints=endpoints,evaluations=evaluations))
  print(json.dumps(dict(label=c['label'],grid=kind,endpoints=endpoints)),flush=True)
paths=[Path(__file__),source,H.parent/'timing-boundary-capable/run.py',H.parent/'timing-continuous-scatter/integration.py']
out=dict(scope='Numerical nominal interval sensitivity only; coverage uncalibrated, two exposed synthetic samples',threshold=cut,hashes={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},cases=rows)
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
