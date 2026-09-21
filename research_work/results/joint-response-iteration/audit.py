#!/usr/bin/env python3
"""JR-1 frozen-parameter numerical audit and nuisance-matched baryon control.

This does not refit the selected companion model or select a revision from
transfer observations. The additional baseline is descriptive development,
not a complete evidence-ratio comparison.
"""
from __future__ import annotations
import argparse, json, hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares
import run as J


def main():
 ap=argparse.ArgumentParser(__doc__);ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args();out=args.output_dir
 target=out/'numerical-audit.json'
 if target.exists():raise FileExistsError('Audit already exists; preserve it')
 summary=json.loads((out/'summary.json').read_text());name=summary['selection']['selected_model'];record=json.loads((out/(name+'.json')).read_text());p=J.DEFAULTS|record['parameters'];spec=record['specification']
 coarse=J.Experiment();fine=J.Experiment(n=3201,order=128,deproj_order=256);rows=[]
 for a,b in zip(coarse.lenses,fine.lenses):
  v0,f0,_=a.prediction(p,spec);v1,f1,(A,rc,rt)=b.prediction(p,spec);angle0=a.angle(p,spec);angle1=b.angle(p,spec)
  local=b.local_parameters(p,spec);q=local['q_sph'] if spec.get('separate_q',False) else local['q'];bend=[]
  for impact in [b.b,angle1/J.ARCSEC*b.Dl]:
   fun=lambda t:float((10**local['logu']*b.baryon_force(impact/np.cos(t))+J.companion_force(impact/np.cos(t),A,rc,rt,q))*impact/np.cos(t))
   independent,err=quad(fun,0,np.pi/2,epsabs=1e-7,epsrel=1e-9,limit=300)
   rr=impact/np.cos(b.lens_t);direct=float(np.dot(b.lens_w,(10**local['logu']*b.baryon_force(rr)+J.companion_force(rr,A,rc,rt,q))*rr))
   bend.append(dict(impact_kpc=impact,relative_quadrature_difference=direct/independent-1,adaptive_estimated_relative_error=err/abs(independent)))
  rows.append(dict(name=a.name,maximum_Vrms_refinement_relative=float(np.max(abs(v1/v0-1))),theta_refinement_relative=angle1/angle0-1,deflection_at_observed_change=f1-f0,lensing_quadrature=bend,light_mass_normalization_error=b.total_light_check-1))
 print('NUMERICAL',json.dumps(rows),flush=True)
 # Independent source integrals use adaptive quadrature of the same photometry.
 import zipfile,io
 source=[]
 catalog={}
 for line in (J.ROOT/J.INPUT_NAMES[1]).read_text().splitlines():
  f=line.split()
  if len(f)==19:
   try:catalog[f[0]]=float(f[11])
   except ValueError:pass
 with zipfile.ZipFile(J.ROOT/J.INPUT_NAMES[0]) as archive:
  for name in ['NGC2403','NGC5055','UGC06667']:
   d=np.atleast_2d(np.loadtxt(io.BytesIO(archive.read(name+'_rotmod.dat'))));rd=catalog[name]
   mask=(d[:,0]>0)&(d[:,6]>0);r=d[mask,0];s=.5e6*d[mask,6]
   adaptive=np.pi*r[0]**2*s[0]+2*np.pi*s[-1]*rd*(r[-1]+rd)
   for l in range(len(r)-1):
    def integrand(x):return 2*np.pi*x*np.exp(np.log(s[l])+(x-r[l])/(r[l+1]-r[l])*(np.log(s[l+1])-np.log(s[l])))
    adaptive+=quad(integrand,r[l],r[l+1],epsabs=1e-3,epsrel=1e-11)[0]
   direct=J.stellar_disk_mass(d,rd);source.append(dict(name=name,relative_mass_difference=direct/adaptive-1))
 # Analytic positive-density identity at the selected parameters, every source.
 full=json.loads((out/'selected-full-results.json').read_text());density=[]
 for d in full['sparc_rows']:
  gal=next(g for g in coarse.galaxies if g['name']==d['name']);fs=gal['Mbulge']/gal['Mstar'];q=p['q']+(p['q_sph']-p['q'])*fs if spec.get('separate_q',False) else p['q']
  r=np.geomspace(d['rc_kpc']*1e-4,d['rt_kpc']*1e4,1000)
  slope=1+q/(1+(r/d['rc_kpc'])**q)-r*r/(r*r+d['rt_kpc']**2)
  density.append(float(slope.min()))
 # An explicit clock-rate scan: common rate rescaling changes relaxation,
 # not stationary occupation. This is not evidence time is absent.
 a,b=.3,.4;P=2.;N=2*P/(a+np.sqrt(a*a+8*b*P));clock=[]
 for k in np.geomspace(1e-12,1e12,13):
  nn=2*k*P/(k*a+np.sqrt((k*a)**2+8*k*b*k*P));tau=1/(k*a+4*k*b*nn)
  clock.append(dict(common_rate_multiplier=float(k),stationary_occupation=float(nn),relaxation_time_in_reference_units=float(tau),relative_occupation_difference=float(nn/N-1)))
 # Extra baseline: fit only ordinary stellar normalization/orbital nuisances,
 # retaining exactly the same training data blocks and priors as R9.
 keys=['logu']+[J.nuisance_key(k,n) for k in ('dm','b') for n in J.TRAIN_LENSES]
 lo=np.array([J.BOUNDS[k][0] for k in keys]);hi=np.array([J.BOUNDS[k][1] for k in keys]);bspec=dict(name='matched_nuisance_baryon_control',free=keys,saturation=False,radial_beta=False,stellar_nuisance=True)
 def residual(x):
  pp=J.DEFAULTS|dict(zip(keys,x));v=coarse.sparc_predictions(pp,bspec,True);parts=[((v-coarse.y)*coarse.sweight)[coarse.train]]
  for l in coarse.lenses:
   if not l.training:continue
   v,f,_=l.prediction(pp,bspec,True);parts.extend([solve_triangular(l.chol,v-l.y,lower=True)/np.sqrt(len(l.y)*4),np.array([f/.05/2]),np.array([pp[J.nuisance_key('dm',l.name)]/l.mass_log_error/2,pp[J.nuisance_key('b',l.name)]/.3/2])])
  return np.concatenate(parts)
 starts=[np.zeros(len(keys)),np.r_[.25,np.zeros(len(keys)-1)],np.r_[.39,np.full(4,.12),np.zeros(4)]];fits=[]
 for x in starts:
  f=least_squares(residual,np.clip(x,lo+1e-7,hi-1e-7),bounds=(lo,hi),max_nfev=200,ftol=2e-8,xtol=2e-8,gtol=2e-7,x_scale='jac');fits.append(f)
 best=min(fits,key=lambda z:float(z.fun@z.fun));bp=J.DEFAULTS|dict(zip(keys,best.x.tolist()));baseline=coarse.metrics(bp,bspec,True,True)
 (out/'matched-nuisance-baryon-full-results.json').write_text(json.dumps(baseline,indent=2,allow_nan=False)+'\n')
 baseline_record=dict(specification=bspec,parameters=bp,objective=float(best.fun@best.fun),optimizer_starts=[dict(success=bool(f.success),nfev=f.nfev,objective=float(f.fun@f.fun)) for f in fits],bounds_touched=[k for k in keys if min(abs(bp[k]-np.array(J.BOUNDS[k])))<1e-3],metrics={k:baseline[k] for k in ['sparc','lenses']})
 result=dict(scope='Frozen selected-parameter refinement; no selected-parameter refitting after transfer exposure',selected_model=spec['name'],lens_checks=rows,source_mass_checks=source,minimum_selected_enclosed_mass_log_slope=min(density),clock_rate_scan=clock,matched_nuisance_baryon_control=baseline_record,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),passed_numerics=bool(max(r['maximum_Vrms_refinement_relative'] for r in rows)<.002 and max(abs(r['theta_refinement_relative']) for r in rows)<.002 and max(abs(c['relative_quadrature_difference']) for r in rows for c in r['lensing_quadrature'])<.002))
 target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print('AUDIT',json.dumps(result),flush=True)
if __name__=='__main__':main()
