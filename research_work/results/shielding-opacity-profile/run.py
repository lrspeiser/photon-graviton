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
rows=[];all_columns={}
for n in [16,32]:
 columns_at_n=[]
 mu,w=leggauss(n);phi=2*np.pi*(np.arange(2*n)+.5)/(2*n)
 dirs=np.array([[np.sqrt(1-m*m)*np.cos(p),np.sqrt(1-m*m)*np.sin(p),m] for m in mu for p in phi]);weights=np.repeat(w/(4*n),2*n);t,wt=leggauss(n//2)
 for index,loc in enumerate(locations):
  point=np.array([loc['R_kpc']*np.cos(loc['phi_rad']),loc['R_kpc']*np.sin(loc['phi_rad']),loc['z_kpc']]);dot=dirs@point;length=-dot+np.sqrt(dot*dot+30**2-point@point);column=np.zeros(len(dirs))
  for lo,hi in zip([0,.01,.03,.1,.3,1,3,10,30],[.01,.03,.1,.3,1,3,10,30,60]):
   a=np.minimum(lo,length);b=np.minimum(hi,length);dist=a[:,None]+(b-a)[:,None]*(t+1)/2
   pos=point[None,None,:]+dirs[:,None,:]*dist[:,:,None];values=density(pos.reshape(-1,3)).reshape(dist.shape)
   column+=(values*wt).sum(axis=1)*(b-a)/2
  columns_at_n.append(column.copy())
  transmission=float(np.exp(logsumexp(np.log(weights)-kappa*column)))
  assert 0<transmission<=1
  rows.append(dict(index=index,mu_nodes=n,location=loc,transmission=transmission))
  if (index+1)%40==0:print(json.dumps(dict(mu_nodes=n,completed=index+1,total=240)),flush=True)
 all_columns[n]=(np.array(columns_at_n),weights.copy())

from scipy.optimize import minimize_scalar
cache=R/'research_work/generated/shielding-opacity-profile';cache.mkdir(exist_ok=True)
q=np.array([v['required_extra_per_baryon'] for v in locations]);outputs=[]
for n,(columns,weights) in all_columns.items():
 archive=cache/f'columns-{n}.npz';np.savez_compressed(archive,columns=columns,weights=weights)
 # Reproduce prior transmissions without changing opacity or geometry.
 expected=json.loads((H.parent/'attached-shielding-grid/results.json').read_text(encoding='utf-8'))
 reference=np.array([v['transmission'] for v in expected['rows'] if v['mu_nodes']==n])
 np.testing.assert_allclose(np.exp(logsumexp(np.log(weights)[None,:]-kappa*columns,axis=1)),reference,rtol=1e-12)
 for name,mask in [('full',np.ones(240,bool)),('inner',np.array([v['R_kpc']<=8 and abs(v['z_kpc'])<=1 for v in locations]))]:
  def assess(logk):
   lt=logsumexp(np.log(weights)[None,:]-10**logk*columns[mask],axis=1)
   logq=np.log(q[mask])-lt;lo=float(logq.min());hi=float(logq.max())
   return dict(log10_kappa=float(logk),log_required_amplitude_min=lo,log_required_amplitude_max=hi,log_amplitude_span=hi-lo,minimum_worst_relative_error=float(np.tanh((hi-lo)/2)),log_optimal_normalization=float(np.log(2)+lo+hi-np.logaddexp(lo,hi)))
  grid=np.linspace(-13,-7,121);scan=[assess(v) for v in grid];candidates=[scan[0],scan[-1]];trials=[]
  for j in range(1,len(grid)-1):
   if scan[j]['log_amplitude_span']<=min(scan[j-1]['log_amplitude_span'],scan[j+1]['log_amplitude_span']):
    f=minimize_scalar(lambda a:assess(a)['log_amplitude_span'],bounds=(grid[j-1],grid[j+1]),method='bounded',options={'xatol':1e-10})
    trials.append(dict(success=bool(f.success),logk=float(f.x),objective=float(f.fun)))
    if f.success:candidates.append(assess(f.x))
  best=min(candidates,key=lambda a:a['log_amplitude_span'])
  outputs.append(dict(subset=name,mu_nodes=n,column_file=str(archive.relative_to(R)),column_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),best=best,scan=scan,optimizer_trials=trials))
  print(json.dumps(dict(subset=name,mu_nodes=n,best=best)),flush=True)
paths=[Path(__file__),code,p,H.parent/'attached-shielding-grid/results.json']
out=dict(scope='Exposed model calibration of one global mass opacity and amplitude; not data validation or all-opacity proof',log10_kappa_search=[-13,-7],boundary_kpc=30,hashes={str(v.relative_to(R)):hashlib.sha256(v.read_bytes()).hexdigest() for v in paths},cases=outputs)
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
