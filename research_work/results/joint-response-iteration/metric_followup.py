#!/usr/bin/env python3
"""R11: one universal companion temporal/spatial ratio (AMENDMENT-4).

Post-transfer development. The same chi field supplies both potentials,
but U_s=Phi_b+gamma_chi*chi is now a fitted effective relation, not a
microscopic relativistic derivation or a per-object lensing adjustment.
"""
from __future__ import annotations
import argparse,json,hashlib,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares,brentq
import population_followup as P
J=P.J
J.DEFAULTS['gamma_chi']=1.;J.BOUNDS['gamma_chi']=(0.,3.)
SPEC=dict(P.SPEC);SPEC['name']='R11_universal_metric_response';SPEC['free']=SPEC['free']+['gamma_chi']

class Lens(P.Lens):
 def prediction(self,p,spec,baryons_only=False):
  v,_,(A,rc,rt)=super().prediction(p,spec,baryons_only);local=self.local_parameters(p,spec)
  q=local['q_sph'] if spec.get('separate_q',False) else local['q'];w=(1+p.get('gamma_chi',1.))/2
  ray=10**local['logu']*self.lens_gb+w*J.companion_force(self.lens_r,A,rc,rt,q)
  alpha=4/J.CLIGHT**2*np.dot(self.lens_w,ray*self.lens_r);f=self.ratio*alpha/(self.theta/J.ARCSEC)-1
  return v,float(f),(A,rc,rt)
 def angle(self,p,spec,baryons_only=False):
  local=self.local_parameters(p,spec);q=local['q_sph'] if spec.get('separate_q',False) else local['q'];w=(1+p.get('gamma_chi',1.))/2
  if baryons_only:A,rc,rt=0.,self.Re,20*self.Re
  else:A,rc,rt=P.source(self.Mstar,0.,self.Re,local,spec['saturation'],1.)
  def residual(b):
   r=b/np.cos(self.lens_t);g=10**local['logu']*self.baryon_force(r)+w*J.companion_force(r,A,rc,rt,q)
   return self.ratio*4/J.CLIGHT**2*np.dot(self.lens_w,g*r)-b/self.Dl
  grid=self.Re*np.geomspace(1e-6,1e4,200);vals=np.array([residual(b) for b in grid]);idx=np.flatnonzero((vals[:-1]>0)&(vals[1:]<=0))
  if not len(idx):return 0.
  i=idx[-1];return float(brentq(residual,grid[i],grid[i+1],xtol=1e-10)/self.Dl*J.ARCSEC)
J.Lens=Lens


def main():
 ap=argparse.ArgumentParser(__doc__);ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args();out=args.output_dir;target=out/(SPEC['name']+'.json')
 if target.exists():raise FileExistsError('Preserve existing follow-up')
 old=json.loads((out/'R10_spheroid_stellar_population.json').read_text());E=P.Experiment();free=SPEC['free'];lo=np.array([J.BOUNDS[k][0] for k in free]);hi=np.array([J.BOUNDS[k][1] for k in free]);base=J.DEFAULTS|old['parameters'];x0=np.clip(np.array([base[k] for k in free]),lo+1e-7,hi-1e-7)
 rng=np.random.default_rng(20260932);starts=[x0]+[np.clip(x0+rng.normal(0,.08,len(free))*(hi-lo),lo+1e-6,hi-1e-6) for _ in range(3)];fits=[];runs=[]
 for n,x in enumerate(starts):
  t=time.monotonic();f=least_squares(E.residual,x,args=(SPEC,),bounds=(lo,hi),max_nfev=300,ftol=2e-8,xtol=2e-8,gtol=2e-7,x_scale='jac');fits.append(f);runs.append(dict(start=x.tolist(),objective=float(f.fun@f.fun),success=bool(f.success),nfev=f.nfev,parameters=f.x.tolist(),seconds=time.monotonic()-t));print('START',n,f.fun@f.fun,f.success,flush=True)
 best=min(fits,key=lambda f:float(f.fun@f.fun));p=J.DEFAULTS|dict(zip(free,best.x.tolist()));record=dict(specification=SPEC,parameters=p,free_parameter_count=len(free),objective=float(best.fun@best.fun),optimizer_starts=runs,scope='Post-R9-transfer effective metric-coupling development; not fresh independent confirmation',bounds_touched=[k for k in free if min(abs(p[k]-np.array(J.BOUNDS[k])))<1e-3],source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (out/'R11-parameters-before-scoring.json').write_text(json.dumps(record,indent=2)+'\n');metrics=E.metrics(p,SPEC,True);record['metrics']=metrics
 fine=P.Experiment(n=3201,order=128,deproj_order=256);checks=[]
 for a,b in zip(E.lenses,fine.lenses):
  va,_,_=a.prediction(p,SPEC);vb,_,_=b.prediction(p,SPEC);aa=a.angle(p,SPEC);ab=b.angle(p,SPEC)
  checks.append(dict(name=a.name,max_vrms_refinement_relative=float(np.max(abs(vb/va-1))),angle_refinement_relative=ab/aa-1))
 record['numerical_checks']=checks;target.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');print('R11',json.dumps(dict(parameters=p,objective=record['objective'],sparc=metrics['sparc'],lenses=metrics['lenses'],numerical_checks=checks)),flush=True)
if __name__=='__main__':main()
