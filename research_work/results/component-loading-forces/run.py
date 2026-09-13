"""Forward gravity of the fitted attached component source; no new fitting."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib,importlib.util
import numpy as np
from scipy.special import roots_legendre
H=Path(__file__).resolve().parent;R=H.parents[2]
source=H.parent/'component-attached-loading/results.json';d=json.loads(source.read_text(encoding='utf-8'))
for name,value in d['hashes'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==value
c=next(v for v in d['cases'] if v['resolution']=='finer' and v['fit_subset']=='inner');eta=c['coefficient_by_component']
code=H.parent/'full-bar-orbits/run.py';spec=importlib.util.spec_from_file_location('orbit_force',code);orbit=importlib.util.module_from_spec(spec);spec.loader.exec_module(orbit);old=orbit.full.old
p=H.parent/'baryon-attached-deposits/results.json';inputs=json.loads(p.read_text(encoding='utf-8'));rows=inputs['rows'];xyz=np.array([[v['R_kpc']*np.cos(v['phi_rad']),v['R_kpc']*np.sin(v['phi_rad']),v['z_kpc']] for v in rows]);rad=np.linalg.norm(xyz,axis=1);cyl=np.hypot(xyz[:,0],xyz[:,1])
ordinary=orbit.Field('ordinary');full=orbit.Field('full')
def acc(field):return np.concatenate([field.evaluate(v)[1] for v in np.array_split(xyz,16)])
ab=acc(ordinary);target=acc(full)-ab
bar=orbit.full.load_field(R/'research_work/data-cache/bar-field/bar-L64.npz')
abar=acc(bar);assert eta['nuclei']==eta['softened_centre']==0
predictions=[]
for refined in [False,True]:
 r=np.geomspace(1e-5,500,3072 if refined else 1536);ell=np.arange(0,257 if refined else 129,2);mu,w=roots_legendre(1024 if refined else 512)
 leg,_=old.basis(ell,mu);angular=leg*(w[:,None]*(2*ell+1)[None,:]/2);coeff=np.empty((len(r),len(ell)))
 for start in range(0,len(r),32):
  rr=r[start:start+32,None];rc=rr*np.sqrt(1-mu**2);z=rr*mu;rho=np.zeros_like(rc)
  for j,(sigma,rd,h,hole,kind) in enumerate(inputs['disk_parameters']):
   surf=sigma*np.exp(-hole/rc-rc/rd);vert=np.exp(-abs(z)/h)/(2*h) if kind=='exp' else np.exp(-2*np.logaddexp(z/(2*h),-z/(2*h))+2*np.log(2))/(4*h)
   rho+=eta['disk_'+str(j)]*surf*vert
  coeff[start:start+len(rr)]=rho@angular
 inn,out=old.scaled_integrals(r,coeff*r[:,None],coeff*r[:,None],ell);fac=4*np.pi*old.G/(2*ell+1)
 field=old.Expansion(r,ell,-fac*(inn+out),fac*((ell+1)*inn-ell*out)/r[:,None]);f=old.force(field,cyl,xyz[:,2])
 pred=eta['bar']*abar+np.c_[f[:,0]*xyz[:,0]/cyl,f[:,0]*xyz[:,1]/cyl,f[:,1]]
 predictions.append(pred);print(json.dumps(dict(refined=refined,completed_points=len(rows))),flush=True)
low,high=predictions;norm=np.linalg.norm(target,axis=1);assert np.all(norm>0)
subsets={'inner':np.array([v['R_kpc']<=8 and abs(v['z_kpc'])<=1 for v in rows]),'full':np.ones(len(rows),bool),'midplane':xyz[:,2]==0};summaries=[]
for name,mask in subsets.items():
 err=np.linalg.norm(high-target,axis=1)/norm;ref=np.linalg.norm(high-low,axis=1)/norm
 summaries.append(dict(subset=name,n=int(mask.sum()),median_relative_vector_error=float(np.median(err[mask])),maximum_relative_vector_error=float(err[mask].max()),rms_absolute_vector_error_kms2_per_kpc=float(np.sqrt(np.mean(np.sum((high[mask]-target[mask])**2,axis=1)))),maximum_refinement_over_target=float(ref[mask].max())))
paths=[Path(__file__),source,code,p,H.parent/'conservative-field-completion/run.py',H.parent/'full-bar-completion/run.py',*ordinary.paths,*full.paths]
result=dict(scope='Forward force comparison of fixed inner component-density fit against empirical target; no observed orbit or lensing validation',hashes={str(v.relative_to(R)):hashlib.sha256(v.read_bytes()).hexdigest() for v in paths},coefficients=eta,summaries=summaries,rows=[dict(location=rows[i],ordinary_acceleration=ab[i].tolist(),target_extra_acceleration=target[i].tolist(),coarse_extra_acceleration=low[i].tolist(),fine_extra_acceleration=high[i].tolist()) for i in range(len(rows))])
(H/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(summaries))
