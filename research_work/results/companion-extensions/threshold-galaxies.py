"""Frozen SPARC response to finite threshold ranges; no fitting."""
from pathlib import Path
import json,math,io,zipfile,hashlib
import numpy as np
from scipy.integrate import quad
from scipy.special import expit
P=Path(__file__).resolve().parent;OLD=P.parent/'isotropic-galaxy-transfer';BASE=P.parents[2]/'temporal_candidate_audit/data'
reference=json.loads((OLD/'third-radiation-retention-results.json').read_text());saved=json.loads((OLD/'third-radiation-retention-predictions.json').read_text())
for name,digest in reference['input_sha256'].items():assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==digest
catalog={}
for line in (BASE/'SPARC_Lelli2016c.mrt').read_text().splitlines():
 f=line.split()
 if len(f)==19:
  try:catalog[f[0]]=(float(f[11]),float(f[7]))
  except ValueError:pass
def weight(y):return math.sin(math.pi/3)/(2*math.pi*(math.cosh(y/3)+.5))
def eta(X,W):
 if W==0:return X**(1/3)/(1+X**(1/3))
 limit=W*math.log(10);norm=quad(weight,-limit,limit,epsabs=1e-12)[0]
 return quad(lambda y:expit(math.log(X)-y)*weight(y),-limit,limit,epsabs=1e-12)[0]/norm
out=dict(scope='Exposed original SPARC splits, frozen companion shape and amplitude except explicit eta substitution; no optimization',input_sha256=reference['input_sha256'],rows=[],scores=[])
with zipfile.ZipFile(BASE/'Rotmod_LTG.zip') as archive:
 for row in saved:
  name=row['galaxy'];arr=np.atleast_2d(np.loadtxt(io.BytesIO(archive.read(name+'_rotmod.dat'))))
  vb=arr[:,3]*abs(arr[:,3])+.5*arr[:,4]*abs(arr[:,4])+.7*arr[:,5]*abs(arr[:,5])
  good=np.isfinite(arr).all(axis=1)&(arr[:,0]>0)&(arr[:,1]>0)&(arr[:,2]>0)&(vb>0)
  arr,vb=arr[good],vb[good]
  assert arr[:,0].tolist()==row['R_kpc'] and arr[:,1].tolist()==row['observed_kms']
  pred=np.array(row['predicted_kms']);extra=pred**2-vb;assert extra.min()>-1e-8
  rd,L=catalog[name];X=L/rd**2;e0=eta(X,0)
  for W in [0,3,6,9]:
   e=eta(X,W);v=np.sqrt(vb+extra*e/e0)
   if W==0:assert np.max(abs(v-pred))<1e-10
   out['rows'].append(dict(galaxy=name,split=row['split'],half_width_decades=W,X=X,retention=e,reference_retention=e0,RMSE_kms=float(np.sqrt(np.mean((v-arr[:,1])**2))),log_mse=float(np.mean(np.log10(v/arr[:,1])**2)),max_prediction_change_kms=float(np.max(abs(v-pred))),predicted_kms=v.tolist()))
for W in [0,3,6,9]:
 for split in ['train','validation','test']:
  rr=[r for r in out['rows'] if r['half_width_decades']==W and r['split']==split]
  score=dict(half_width_decades=W,split=split,galaxies=len(rr),RMSE_kms=float(np.sqrt(np.mean([r['RMSE_kms']**2 for r in rr]))),log_RMS=float(np.sqrt(np.mean([r['log_mse'] for r in rr]))),max_prediction_change_kms=max(r['max_prediction_change_kms'] for r in rr))
  if W==0:
   known=reference['models']['attenuated'];assert min(abs(score['RMSE_kms']-known[k][split]['RMSE_kms']) for k in ['scores','finer_scores'])<1e-8
  out['scores'].append(score)
assert len(saved)==149 and sum(len(r['R_kpc']) for r in saved)==3150
out['input_range']=[min(r['X'] for r in out['rows']),max(r['X'] for r in out['rows'])]
(P/'threshold-galaxies-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(out['input_range']);print(out['scores'])
