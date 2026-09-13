"""Isotropic straight-ray shielding in the fixed 3D baryonic model."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import logsumexp
from scipy.optimize import brentq
H=Path(__file__).resolve().parent;R=H.parents[2]
sys.path.insert(0,str(H.parent/'bar-field-foundation'))
from field import bar,nuclei
source=H.parent/'baryon-attached-deposits/results.json';d=json.loads(source.read_text(encoding='utf-8'))
for name,v in d['source_hashes'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==v
pairpath=H.parent/'depth-attached-loading/results.json';pair=json.loads(pairpath.read_text(encoding='utf-8'))['cases'][1]
locations=[pair[k]['location'] for k in ['critical_shallower','critical_deeper']]
xyz=np.array([[r['R_kpc']*np.cos(r['phi_rad']),r['R_kpc']*np.sin(r['phi_rad']),r['z_kpc']] for r in locations]);target=locations[0]['required_extra_per_baryon']/locations[1]['required_extra_per_baryon']
def density(x):
 rr=np.maximum(np.hypot(x[:,0],x[:,1]),1e-30);z=x[:,2];v=bar(x)+nuclei(x)
 for sigma,rd,h,hole,kind in d['disk_parameters']:
  surf=sigma*np.exp(-hole/rr-rr/rd)
  vert=np.exp(-abs(z)/h)/(2*h) if kind=='exp' else np.exp(-2*np.logaddexp(z/(2*h),-z/(2*h))+2*np.log(2))/(4*h)
  v+=surf*vert
 v+=3*4.1e6*.001**2/(4*np.pi*(np.sum(x*x,axis=1)+.001**2)**2.5)
 return v
np.testing.assert_allclose(density(xyz),[v['declared_baryon_density_Msun_kpc3'] for v in locations],rtol=1e-12)
rows=[]
for radius in [30.,60.]:
 for n in [16,32,64]:
  mu,w=leggauss(n);phi=2*np.pi*(np.arange(2*n)+.5)/(2*n)
  dirs=np.array([[np.sqrt(1-m*m)*np.cos(p),np.sqrt(1-m*m)*np.sin(p),m] for m in mu for p in phi]);weights=np.repeat(w/(4*n),2*n);assert abs(weights.sum()-1)<1e-14
  t,wt=leggauss(n//2);columns=[]
  for point in xyz:
   dot=dirs@point;length=-dot+np.sqrt(dot*dot+radius*radius-point@point)
   column=np.zeros(len(dirs))
   for lo,hi in zip([0,.01,.03,.1,.3,1,3,10,30,60],[.01,.03,.1,.3,1,3,10,30,60,120]):
    a=np.minimum(lo,length);b=np.minimum(hi,length);dist=a[:,None]+(b-a)[:,None]*(t+1)/2
    pos=point[None,None,:]+dirs[:,None,:]*dist[:,:,None]
    values=density(pos.reshape(-1,3)).reshape(dist.shape)
    column+=((values*wt).sum(axis=1)*(b-a)/2)
   columns.append(column)
  columns=np.array(columns);assert np.all(columns>0)
  def logs(logk):return logsumexp(np.log(weights)[None,:]-10**logk*columns,axis=1)
  root=brentq(lambda l:float(np.diff(logs(l))[0]*-1-np.log(target)),-13,-7,xtol=1e-10)
  values=logs(root)
  rows.append(dict(boundary_kpc=radius,angular_mu_nodes=n,phi_nodes=2*n,segment_nodes=n//2,required_exposure_ratio=target,cross_section_kpc2_per_Msun=10**root,relative_intensities=np.exp(values).tolist(),achieved_ratio=float(np.exp(values[0]-values[1])),mean_columns_Msun_kpc2=(columns@weights).tolist()))
  print(json.dumps(rows[-1]),flush=True)
paths=[Path(__file__),source,pairpath,H.parent/'bar-field-foundation/field.py',H.parent/'bar-field-foundation/published_bar.py']
out=dict(scope='Necessary pairwise shielding feasibility under postulated mass-proportional opacity; not a whole-galaxy fit or supply calculation',hashes={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},locations=locations,cases=rows)
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
