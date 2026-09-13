"""Finite formation time for the frozen threshold-mixture galaxy model."""
from pathlib import Path
import json,math,hashlib
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
P=Path(__file__).resolve().parent;OLD=P.parent/'isotropic-galaxy-transfer'
old=json.loads((P/'threshold-galaxies-results.json').read_text())
for name,digest in old['input_sha256'].items():assert hashlib.sha256((P.parents[2]/'temporal_candidate_audit/data'/name).read_bytes()).hexdigest()==digest
base=[r for r in old['rows'] if r['half_width_decades']==0]
baryons={r['galaxy']:r for r in json.loads((OLD/'model-comparison-predictions.json').read_text()) if r['model']=='baryons'}
L=6*math.log(10)
def w(y):return math.sin(math.pi/3)/(2*math.pi*(math.cosh(y/3)+.5))
norm=quad(w,-L,L)[0]
def integrate(fn):return quad(lambda y:fn(math.exp(y))*w(y),-L,L,epsabs=1e-11)[0]/norm
def occupied(X,u):return integrate(lambda t:X/(X+t)*(1 if u is None else -math.expm1(-(X+t)*u)))
out=dict(scope='All exposed original galaxies, fixed 12-decade thresholds, empty initial storage and common constant-source dimensionless duration; no fitting or physical ages',rows=[],scores=[],energy_checks=[])
for row in base:
 name=row['galaxy'];b=baryons[name];v0=np.array(row['predicted_kms']);vb=np.array(b['predicted_kms'])**2;y=np.array(b['observed_kms']);X=row['X'];eq=occupied(X,None)
 ref=next(r for r in old['rows'] if r['galaxy']==name and r['half_width_decades']==6)
 assert abs(eq-ref['retention'])<1e-10
 u90=brentq(lambda u:occupied(X,u)-.9*eq,0,math.log(10)/X,xtol=1e-10)
 for u in [.1,1.,10.,100.,1000.,None]:
  f=occupied(X,u);v=np.sqrt(vb+(v0*v0-vb)*f/row['retention'])
  out['rows'].append(dict(galaxy=name,split=row['split'],X=X,duration=u,occupied_fraction=f,equilibrium_fraction=eq,fraction_of_equilibrium=f/eq,duration_to_90_percent=u90,RMSE_kms=float(np.sqrt(np.mean((v-y)**2))),log_mse=float(np.mean(np.log10(v/y)**2))))
for u in [.1,1.,10.,100.,1000.,None]:
 for split in ['train','validation','test']:
  rr=[r for r in out['rows'] if r['duration']==u and r['split']==split]
  score=dict(duration=u,split=split,RMSE_kms=float(np.sqrt(np.mean([r['RMSE_kms']**2 for r in rr]))),log_RMS=float(np.sqrt(np.mean([r['log_mse'] for r in rr]))))
  if u is None:
   ref=next(r for r in old['scores'] if r['half_width_decades']==6 and r['split']==split)
   assert abs(score['RMSE_kms']-ref['RMSE_kms'])<1e-8
  out['scores'].append(score)
for X in [min(r['X'] for r in base),1.,max(r['X'] for r in base)]:
 for u in [1.,100.]:
  def capture(t):
   k=X+t;feq=X/k
   return X*((1-feq)*u+feq*(-math.expm1(-k*u))/k)
  def release(t):
   k=X+t;feq=X/k
   return t*feq*(u+math.expm1(-k*u)/k)
  C=integrate(capture);R=integrate(release);U=occupied(X,u)
  assert abs(C-R-U)<1e-8 and R>=0
  out['energy_checks'].append(dict(X=X,duration=u,captured_energy_per_capacity=C,released_energy_per_capacity=R,stored_energy_per_capacity=U))
(P/'threshold-history-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(out['scores']);print('u90 range',min(r['duration_to_90_percent'] for r in out['rows']),max(r['duration_to_90_percent'] for r in out['rows']))
