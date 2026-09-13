"""Transfer fixed pair-calibrated opacity to the entire exposed source grid."""
from pathlib import Path
import json,hashlib
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import logsumexp
H=Path(__file__).resolve().parent;R=H.parents[2]
code=H.parent/'attached-shielding-pair/run.py';ns={'__file__':str(code)}
exec(compile(code.read_text(encoding='utf-8').split('\nrows=[]\n')[0],str(code),'exec'),ns)
density=ns['density'];data=ns['d'];locations=data['rows'];assert len(locations)==240
p=H.parent/'attached-shielding-pair/results.json';prior=json.loads(p.read_text(encoding='utf-8'))
for name,v in prior['hashes'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==v
kappa=next(v['cross_section_kpc2_per_Msun'] for v in prior['cases'] if v['boundary_kpc']==30 and v['angular_mu_nodes']==64)
rows=[]
for n in [16,32]:
 mu,w=leggauss(n);phi=2*np.pi*(np.arange(2*n)+.5)/(2*n)
 dirs=np.array([[np.sqrt(1-m*m)*np.cos(p),np.sqrt(1-m*m)*np.sin(p),m] for m in mu for p in phi]);weights=np.repeat(w/(4*n),2*n);t,wt=leggauss(n//2)
 for index,loc in enumerate(locations):
  point=np.array([loc['R_kpc']*np.cos(loc['phi_rad']),loc['R_kpc']*np.sin(loc['phi_rad']),loc['z_kpc']]);dot=dirs@point;length=-dot+np.sqrt(dot*dot+30**2-point@point);column=np.zeros(len(dirs))
  for lo,hi in zip([0,.01,.03,.1,.3,1,3,10,30],[.01,.03,.1,.3,1,3,10,30,60]):
   a=np.minimum(lo,length);b=np.minimum(hi,length);dist=a[:,None]+(b-a)[:,None]*(t+1)/2
   pos=point[None,None,:]+dirs[:,None,:]*dist[:,:,None];values=density(pos.reshape(-1,3)).reshape(dist.shape)
   column+=(values*wt).sum(axis=1)*(b-a)/2
  transmission=float(np.exp(logsumexp(np.log(weights)-kappa*column)))
  assert 0<transmission<=1
  rows.append(dict(index=index,mu_nodes=n,location=loc,transmission=transmission))
  if (index+1)%40==0:print(json.dumps(dict(mu_nodes=n,completed=index+1,total=240)),flush=True)
q=np.array([v['required_extra_per_baryon'] for v in locations]);summaries=[]
for name,mask in [('full',np.ones(240,bool)),('inner',np.array([v['R_kpc']<=8 and abs(v['z_kpc'])<=1 for v in locations])),('midplane',np.array([v['z_kpc']==0 for v in locations]))]:
 for n in [16,32]:
  T=np.array([v['transmission'] for v in rows if v['mu_nodes']==n]);a=(q/T)[mask];lo=float(a.min());hi=float(a.max());normal=2*lo*hi/(lo+hi)
  summaries.append(dict(subset=name,mu_nodes=n,minimum_worst_relative_source_error=(hi-lo)/(hi+lo),optimal_common_normalization=normal,required_normalization_range=[lo,hi],n=int(mask.sum())))
T0=np.array([v['transmission'] for v in rows if v['mu_nodes']==16]);T1=np.array([v['transmission'] for v in rows if v['mu_nodes']==32])
out=dict(scope='Transfer across exposed model grid, not held-out observations; no per-point fitting',kappa_kpc2_per_Msun=kappa,boundary_kpc=30,hashes={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),code,p]},summaries=summaries,maximum_relative_transmission_refinement=float(np.max(abs(T1/T0-1))),rows=rows)
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(summaries),flush=True)
