#!/usr/bin/env python3
"""R10: explicitly post-transfer stellar-population follow-up (AMENDMENT-3).

No observations change. R9 remains the original frozen transfer record.
This extension adds one common spheroidal stellar normalization, not a
per-object companion amplitude or a separate photon response.
"""
from __future__ import annotations
import argparse,json,hashlib,io,zipfile,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import run as J
BASE_SOURCE=J.source_parameters
BASE_LENS=J.Lens
J.DEFAULTS['logusph']=0.;J.BOUNDS['logusph']=(-.3,.4)
SPEC=dict(J.SPECS[9]);SPEC['name']='R10_spheroid_stellar_population';SPEC['free']=SPEC['free']+['logusph']

def source(Mstar,Mgas,Re,p,saturation=False,spheroid_fraction=0.):
 f=np.asarray(spheroid_fraction);u=10**p.get('logusph',0.);factor=1+(u-1)*f
 return BASE_SOURCE(np.asarray(Mstar)*factor,Mgas,Re,p,saturation,u*f/factor)

class Lens(BASE_LENS):
 def local_parameters(self,p,spec):
  local=super().local_parameters(p,spec).copy()
  local['logu']+=p.get('logusph',0.);local['logusph']=0.
  return local

J.source_parameters=source
J.Lens=Lens

class Experiment(J.Experiment):
 def __init__(self,*args,**kwargs):
  super().__init__(*args,**kwargs);bul=[]
  with zipfile.ZipFile(J.ROOT/J.INPUT_NAMES[0]) as archive:
   for g in self.galaxies:
    rows=np.atleast_2d(np.loadtxt(io.BytesIO(archive.read(g['name']+'_rotmod.dat'))));rows=rows[rows[:,0]>0];bul.append(.7*rows[:,5]**2)
  self.bulge_v2=np.concatenate(bul)
 def sparc_predictions(self,p,spec,baryons_only=False):
  u=10**p.get('logusph',0.);v2=np.maximum(self.vg2+10**p['logu']*(self.vs2+(u-1)*self.bulge_v2),0.)
  if not baryons_only:
   fs=self.Msph/self.Ms;updated_fs=u*fs/(1+(u-1)*fs)
   A,rc,rt=source(self.Ms,self.Mg,self.Re,p,spec['saturation'],fs)
   q=p['q']+(p['q_sph']-p['q'])*updated_fs if spec.get('separate_q',False) else np.full_like(self.Ms,p['q'])
   v2+=self.r*J.companion_force(self.r,A[self.ids],rc[self.ids],rt[self.ids],q[self.ids])
  return np.sqrt(v2)
 def residual(self,x,spec):
  r=super().residual(x,spec);p=J.DEFAULTS|dict(zip(spec['free'],x))
  return np.r_[r,p['logusph']/.15]

def main():
 ap=argparse.ArgumentParser(__doc__);ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args();out=args.output_dir;target=out/(SPEC['name']+'.json')
 if target.exists():raise FileExistsError('Preserve existing follow-up')
 old=json.loads((out/'R9_spheroid_shape_stellar_priors.json').read_text());E=Experiment();free=SPEC['free'];lo=np.array([J.BOUNDS[k][0] for k in free]);hi=np.array([J.BOUNDS[k][1] for k in free]);base=J.DEFAULTS|old['parameters'];x0=np.clip(np.array([base[k] for k in free]),lo+1e-7,hi-1e-7)
 rng=np.random.default_rng(20260931);starts=[x0]+[np.clip(x0+rng.normal(0,.08,len(free))*(hi-lo),lo+1e-6,hi-1e-6) for _ in range(3)];fits=[];runs=[]
 for n,x in enumerate(starts):
  t=time.monotonic();f=least_squares(E.residual,x,args=(SPEC,),bounds=(lo,hi),max_nfev=300,ftol=2e-8,xtol=2e-8,gtol=2e-7,x_scale='jac');fits.append(f);runs.append(dict(start=x.tolist(),objective=float(f.fun@f.fun),success=bool(f.success),nfev=f.nfev,parameters=f.x.tolist(),seconds=time.monotonic()-t));print('START',n,f.fun@f.fun,f.success,flush=True)
 best=min(fits,key=lambda f:float(f.fun@f.fun));p=J.DEFAULTS|dict(zip(free,best.x.tolist()));record=dict(specification=SPEC,parameters=p,free_parameter_count=len(free),objective=float(best.fun@best.fun),optimizer_starts=runs,scope='Post-R9-transfer exploratory development; no fresh independent confirmation',bounds_touched=[k for k in free if min(abs(p[k]-np.array(J.BOUNDS[k])))<1e-3],source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 # Freeze parameters before calculating follow-up comparison scores.
 (out/'R10-parameters-before-scoring.json').write_text(json.dumps(record,indent=2)+'\n')
 metrics=E.metrics(p,SPEC,True);record['metrics']=metrics
 fine=Experiment(n=3201,order=128,deproj_order=256);checks=[]
 for a,b in zip(E.lenses,fine.lenses):
  va,_,_=a.prediction(p,SPEC);vb,_,_=b.prediction(p,SPEC);aa=a.angle(p,SPEC);ab=b.angle(p,SPEC)
  checks.append(dict(name=a.name,max_vrms_refinement_relative=float(np.max(abs(vb/va-1))),angle_refinement_relative=ab/aa-1))
 record['numerical_checks']=checks
 target.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');print('R10',json.dumps(dict(parameters=p,objective=record['objective'],sparc=metrics['sparc'],lenses=metrics['lenses'],numerical_checks=checks)),flush=True)
if __name__=='__main__':main()
