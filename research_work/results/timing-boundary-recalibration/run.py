"""Recalibrate on all exposed frozen simulations; no observational inference."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2
H=Path(__file__).resolve().parent;R=H.parents[2]
sys.path.insert(0,str(H.parent/'timing-boundary-capable'))
from run import fit,outside,averaged_likelihood
sources=[H/'run.py',H/'protocol.json',H.parent/'timing-boundary-capable/run.py',H.parent/'timing-continuous-scatter/integration.py',H.parent/'timing-coverage-calibration/results.json']
hashes={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
cache=R/'research_work/generated/timing-boundary-recalibration';cache.mkdir(exist_ok=True)
cp=cache/'checkpoint.json'
state=json.loads(cp.read_text(encoding='utf-8')) if cp.exists() else dict(hashes=hashes,cases=[])
assert state['hashes']==hashes
base=json.loads(sources[-1].read_text(encoding='utf-8'))
for c in base['cases']:
 if c['label'] in {v['label'] for v in state['cases']}:continue
 p=R/'research_work/generated/timing-coverage-calibration'/(c['label']+'.npz')
 assert hashlib.sha256(p.read_bytes()).hexdigest()==c['array_sha256']
 with np.load(p) as d:x=d['log_width'];logs=d['event_log_likelihoods'];z=d['redshift']
 fitted=fit(x,logs,z,c['fit']);best=fitted['best'];a,b,s=best['parameters']
 L=np.exp(logs-logs.max(axis=1)[:,None]);lz=np.log1p(z);truth=c['truth_b']
 def objective(q):
  try:v=averaged_likelihood(x,L,q[0]+truth*lz,q[1])
  except ArithmeticError:return 1e12
  return float(-np.log(v).sum()) if np.all(v>0) else 1e12
 def constraint(q):return 1e-6-outside(x,q[0]+truth*lz,q[1])
 trials=[]
 for q in [[a,s],[a,0],[np.log(30),.1]]:
  f=minimize(objective,q,method='SLSQP',bounds=[np.log([5,100]),[0,.6]],constraints={'type':'ineq','fun':constraint},options={'ftol':1e-10,'maxiter':600})
  trials.append(dict(parameters=f.x.tolist(),objective=float(f.fun),success=bool(f.success),outside=float(outside(x,f.x[0]+truth*lz,f.x[1]).max())))
 good=[v for v in trials if v['success'] and v['outside']<=1.00001e-6 and v['objective']<1e11]
 raw_lr=2*(min(v['objective'] for v in good)-best['objective']) if good else None
 valid=bool(good) and raw_lr>=-2e-6 and min(a-np.log(5),np.log(100)-a,b+2,3-b,.6-s)>1e-5
 lr=max(0,raw_lr) if raw_lr is not None and raw_lr>=-2e-6 else None
 row=dict(label=c['label'],family=c['family'],snr=c['snr'],truth_b=truth,array_sha256=c['array_sha256'],fit=fitted,truth_profiles=trials,raw_lr=raw_lr,lr=lr,valid=bool(valid),accepted_95=bool(valid and lr<=chi2.ppf(.95,1)))
 state['cases'].append(row)
 tmp=cp.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n',encoding='utf-8',newline='\n');tmp.replace(cp)
 print(json.dumps(dict(completed=len(state['cases']),label=c['label'],b=b,sigma=s,valid=bool(valid),lr=lr)),flush=True)
assert len(state['cases'])==160
(H/'results.json').write_text(json.dumps(state,indent=2)+'\n',encoding='utf-8',newline='\n')
