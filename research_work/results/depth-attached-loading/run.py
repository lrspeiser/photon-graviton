"""Monotone well-depth attachment bound for existing 3D source target."""
from pathlib import Path
import json,hashlib,importlib.util
import numpy as np
from scipy.optimize import linprog
H=Path(__file__).resolve().parent;R=H.parents[2]
p=H.parent/'baryon-attached-deposits/results.json';d=json.loads(p.read_text(encoding='utf-8'))
for name,value in d['source_hashes'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==value
code=H.parent/'full-bar-orbits/run.py';spec=importlib.util.spec_from_file_location('orbit_depth',code);orbit=importlib.util.module_from_spec(spec);spec.loader.exec_module(orbit)
rows=d['rows'];xyz=np.array([[r['R_kpc']*np.cos(r['phi_rad']),r['R_kpc']*np.sin(r['phi_rad']),r['z_kpc']] for r in rows]);q=np.array([r['required_extra_per_baryon'] for r in rows]);outputs=[];paths=[Path(__file__),p,code]
for field_kind in ['ordinary','full']:
 field=orbit.Field(field_kind);paths+=field.paths
 W=-np.concatenate([field.evaluate(part)[0] for part in np.array_split(xyz,16)])
 assert np.all(np.isfinite(W))
 # Only potential ordering is used: additive reference does not affect monotonicity.
 for name,mask in [('full',np.ones(len(q),bool)),('inner',np.array([r['R_kpc']<=8 and abs(r['z_kpc'])<=1 for r in rows]))]:
  ix=np.flatnonzero(mask);ix=ix[np.argsort(W[ix])];v=q[ix];n=len(v)
  # If W_i <= W_j then monotonic eta demands eta_i <= eta_j.
  mismatch=(v[:,None]-v[None,:])/(v[:,None]+v[None,:]);mismatch[np.tril_indices(n)]=0
  i,j=np.unravel_index(np.argmax(mismatch),mismatch.shape);bound=float(mismatch[i,j])
  A=[];rhs=[]
  for k in range(n):
   a=np.zeros(n+1);a[k]=1;a[-1]=-v[k];A.append(a);rhs.append(v[k])
   a=np.zeros(n+1);a[k]=-1;a[-1]=-v[k];A.append(a);rhs.append(-v[k])
  for k in range(n-1):
   a=np.zeros(n+1);a[k]=1;a[k+1]=-1;A.append(a);rhs.append(0.)
  cost=np.zeros(n+1);cost[-1]=1
  fit=linprog(cost,A_ub=np.array(A),b_ub=rhs,bounds=[(0,None)]*n+[(0,1)],method='highs')
  assert fit.success and abs(fit.fun-bound)<1e-7
  outputs.append(dict(field=field_kind,subset=name,n=n,minimum_worst_relative_error=bound,linear_program_error=float(fit.fun),critical_shallower=dict(location=rows[int(ix[i])],negative_potential_kms2=float(W[ix[i]])),critical_deeper=dict(location=rows[int(ix[j])],negative_potential_kms2=float(W[ix[j]])),ordered_rows=[dict(index=int(k),negative_potential_kms2=float(W[k]),required_loading=float(q[k]),fitted_loading=float(fit.x[t])) for t,k in enumerate(ix)]))
  print(json.dumps({k:outputs[-1][k] for k in ['field','subset','minimum_worst_relative_error']}),flush=True)
out=dict(scope='Conditional test of monotone depth-only attached loading under common exposure; no observation or general capture exclusion',hashes={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},cases=outputs)
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
