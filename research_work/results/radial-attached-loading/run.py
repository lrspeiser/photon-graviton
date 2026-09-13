"""Conditional minimax bound for arbitrary radial loading of ordinary matter."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.optimize import minimize_scalar
H=Path(__file__).resolve().parent;R=H.parents[2]
p=H.parent/'baryon-attached-deposits/results.json';d=json.loads(p.read_text(encoding='utf-8'))
for name,value in d['source_hashes'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==value
rows=d['rows'];out=[]
for name,subset in [('full',rows),('inner_R_le_8_z_le_1',[r for r in rows if r['R_kpc']<=8 and abs(r['z_kpc'])<=1]),('midplane',[r for r in rows if abs(r['z_kpc'])==0])]:
 assert subset
 groups=[]
 for radius in sorted(set(r['R_kpc'] for r in subset)):
  group=[r for r in subset if r['R_kpc']==radius];q=np.array([r['required_extra_per_baryon'] for r in group]);lo=float(q.min());hi=float(q.max())
  eta=2*lo*hi/(lo+hi);err=(hi-lo)/(hi+lo)
  opt=minimize_scalar(lambda e:np.max(abs(e/q-1)),bounds=(lo,hi),method='bounded',options={'xatol':1e-11})
  assert opt.success and abs(opt.fun-err)<1e-7
  groups.append(dict(R_kpc=radius,n=len(group),q_min=lo,q_max=hi,eta_minimax=eta,minimum_worst_fractional_error=err,minimum_worst_multiplicative_factor=float(np.sqrt(hi/lo)),minimum_location=group[int(q.argmin())],maximum_location=group[int(q.argmax())]))
 out.append(dict(subset=name,n=len(subset),radial_groups=groups,minimum_worst_fractional_error=max(g['minimum_worst_fractional_error'] for g in groups),minimum_worst_multiplicative_factor=max(g['minimum_worst_multiplicative_factor'] for g in groups)))
result=dict(scope='Conditional shape test of arbitrary eta(R) attached loading; no new observed fit',source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),subsets=out)
(H/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([{k:v for k,v in s.items() if k!='radial_groups'} for s in out]))
