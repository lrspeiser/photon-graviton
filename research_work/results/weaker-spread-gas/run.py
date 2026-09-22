#!/usr/bin/env python3
"""JR-4: fixed-source environmental weakness and mass-preserving spreading.
Run: python run.py --jr1 /path/Photon-Graviton-JR1 --output /fresh/result-dir
The original JR-1 package is supplied with the reproducibility bundle.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time, platform
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, minimize_scalar, brentq
from scipy.integrate import quad, solve_ivp
from scipy.linalg import solve_triangular
from scipy.stats import spearmanr

SPECS=[('baseline',[],False),('constant_total',['a'],True),
 ('weaker_companion',['a','logSigma'],False),('spread_only',['b','logSigma'],False),
 ('weaker_spread_companion',['a','b','logSigma'],False),
 ('weaker_total',['a','logSigma'],True),('weaker_total_spread',['a','b','logSigma'],True)]
BOUNDS={'a':(0.,np.log(10.)), 'b':(0.,np.log(10.)), 'logSigma':(-1.,3.)}
DEFAULT={'a':0.,'b':0.,'logSigma':1.}

def dump(path,data):
 path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')

def modifiers(Sigma,pars,name):
 w=np.ones_like(np.asarray(Sigma,float)) if name=='constant_total' else 1/(1+np.asarray(Sigma)/10**pars['logSigma'])
 return np.exp(-pars['a']*w),np.exp(pars['b']*w)

def score(y,v,e):
 f=(v-y)/y
 return dict(fractional_RMS=float(np.sqrt(np.mean(f*f))),RMSE_kms=float(np.sqrt(np.mean((v-y)**2))),
  chi2=float(np.sum(((v-y)/e)**2)),mean_fractional=float(f.mean()))

def summarize(rows):
 return dict(galaxies=len(rows),points=sum(r['points'] for r in rows),
  mean_galaxy_RMSE_kms=float(np.mean([r['RMSE_kms'] for r in rows])),
  median_fractional_RMS=float(np.median([r['fractional_RMS'] for r in rows])),
  mean_fractional_square=float(np.mean([r['fractional_RMS']**2 for r in rows])),
  raw_chi2=float(sum(r['chi2'] for r in rows)),
  below10=sum(r['fractional_RMS']<.1 for r in rows),below20=sum(r['fractional_RMS']<.2 for r in rows))

class Study:
 def __init__(self,jr1):
  self.jr1=jr1;src=jr1/'repository/research_work/results/joint-response-iteration'
  sys.path.insert(0,str(src));import population_followup as P
  self.P=P;self.J=P.J;self.saved=json.loads((jr1/'evidence/R10_spheroid_stellar_population.json').read_text())
  self.p=self.saved['parameters'];self.spec=self.saved['specification'];self.E=P.Experiment();E=self.E;p=self.p
  u=10**p['logusph'];fs=E.Msph/E.Ms;updated=u*fs/(1+(u-1)*fs)
  self.A,self.rc,self.rt=P.source(E.Ms,E.Mg,E.Re,p,False,fs)
  self.q=p['q']+(p['q_sph']-p['q'])*updated
  self.mb=10**p['logu']*(E.Ms+(u-1)*E.Msph)+E.Mg
  self.Sigma=self.mb/(2*np.pi*E.Re**2*1e6)
  self.fgas=E.Mg/self.mb
  self.vb2=np.maximum(E.vg2+10**p['logu']*(E.vs2+(u-1)*E.bulge_v2),0.)
  self.base=self.prediction(DEFAULT,'baseline',False)
  self.base_rows=self.rows(self.base,DEFAULT,'baseline')
  self.original_over=np.array([r['fractional_RMS']>.2 and r['mean_fractional']>0 for r in self.base_rows])
  self.original_under=np.array([r['fractional_RMS']>.2 and r['mean_fractional']<0 for r in self.base_rows])
  sizes=np.bincount(E.ids);self.sizes=sizes
  self.hashes={f:hashlib.sha256((self.J.ROOT/f).read_bytes()).hexdigest() for f in self.J.INPUT_NAMES}
 def prediction(self,p,name,total):
  E=self.E;s,lam=modifiers(self.Sigma,p,name);ids=E.ids
  gc=self.J.companion_force(E.r,self.A[ids]*s[ids]/lam[ids],self.rc[ids]*lam[ids],self.rt[ids]*lam[ids],self.q[ids])
  return np.sqrt((s[ids] if total else 1.)*self.vb2+E.r*gc)
 def rows(self,v,pars,name):
  E=self.E;s,lam=modifiers(self.Sigma,pars,name);rr=[]
  for i,g in enumerate(E.galaxies):
   ix=E.ids==i;sc=score(E.y[ix],v[ix],E.error[ix]);rr.append(dict(name=g['name'],split=g['split'],points=int(ix.sum()),
    Sigma_proxy_Msun_pc2=float(self.Sigma[i]),fgas=float(self.fgas[i]),strength=float(s[i]),spread=float(lam[i]),
    radius_kpc=E.r[ix].tolist(),observed_kms=E.y[ix].tolist(),errors_kms=E.error[ix].tolist(),predicted_kms=v[ix].tolist(),**sc))
  return rr
 def stats(self,rows):
  out={sp:summarize([r for r in rows if r['split']==sp]) for sp in ('train','validation','test')};out['all']=summarize(rows)
  for label,mask in [('original_overpredicted_outliers',self.original_over),('original_underpredicted_outliers',self.original_under)]:
   out[label]=summarize([r for r,m in zip(rows,mask) if m])
  return out
 def lens(self,l,pars,name,total):
  local=l.local_parameters(self.p,self.spec);J=self.J
  A,rc,rt=self.P.source(l.Mstar,0.,l.Re,local,False,1.);q=local['q_sph'];mass=l.Mstar*10**local['logu']
  Sigma=mass/(2*np.pi*l.Re**2*1e6);s,lam=map(float,modifiers(Sigma,pars,name));bf=s if total else 1.
  force=lambda r:bf*10**local['logu']*l.baryon_force(r)+J.companion_force(r,s*A/lam,lam*rc,lam*rt,q)
  v=l.moments(force(l.r),local,False)
  def ring(b):
   r=b/np.cos(l.lens_t);return l.ratio*4/J.CLIGHT**2*np.dot(l.lens_w,force(r)*r)-b/l.Dl
  grid=l.Re*np.geomspace(1e-6,1e4,140);val=np.array([ring(b) for b in grid]);idx=np.flatnonzero((val[:-1]>0)&(val[1:]<=0))
  angle=0. if not len(idx) else brentq(ring,grid[idx[-1]],grid[idx[-1]+1],xtol=1e-10)/l.Dl*J.ARCSEC
  e=v-l.y;wh=solve_triangular(l.chol,e,lower=True)
  # Independent adaptive quadrature at measured Einstein impact.
  b=l.b;integral,err=quad(lambda t:float(force(b/np.cos(t))*b/np.cos(t)),0,np.pi/2,epsabs=1e-7,epsrel=1e-9,limit=250)
  fixed=np.dot(l.lens_w,force(l.lens_r)*l.lens_r)
  return dict(name=l.name,role='prior_R10_fit' if l.training else 'prior_R10_transfer',Sigma_proxy_Msun_pc2=Sigma,strength=s,spread=lam,
   observed_vrms=l.y.tolist(),predicted_vrms=v.tolist(),covariance=l.cov.tolist(),
   stellar_fractional_RMS=float(np.sqrt(np.mean((e/l.y)**2))),stellar_chi2=float(wh@wh),
   theta_observed=l.theta,theta_predicted=float(angle),theta_fractional_error=float(angle/l.theta-1),
   deflection_quadrature_relative=float(fixed/integral-1))

def fit_forms(S,objective):
 E=S.E;train=E.split[E.ids]=='train';den=(E.y if objective=='fractional' else E.error)*np.sqrt(S.sizes[E.ids]*np.sum(E.split=='train'))
 records=[]
 for name,keys,total in SPECS:
  def res(x):
   p=DEFAULT|dict(zip(keys,x));return ((S.prediction(p,name,total)-E.y)/den)[train]
  if keys:
   lo=np.array([BOUNDS[k][0] for k in keys]);hi=np.array([BOUNDS[k][1] for k in keys]);starts=[]
   for logSig in (-.5,1.,2.5):
    for val in (.1,1.5):starts.append(np.array([logSig if k=='logSigma' else val for k in keys]))
   if name=='constant_total':starts=starts[:2]
   fits=[least_squares(res,x,bounds=(lo,hi),max_nfev=250,ftol=1e-11,xtol=1e-11,gtol=1e-10) for x in starts]
   f=min(fits,key=lambda z:z.fun@z.fun);p=DEFAULT|dict(zip(keys,f.x.tolist()))
   optimization=[dict(success=bool(f.success),nfev=f.nfev,objective=float(f.fun@f.fun),parameters=dict(zip(keys,f.x.tolist()))) for f in fits]
  else:p=DEFAULT.copy();optimization=[]
  v=S.prediction(p,name,total);rows=S.rows(v,p,name)
  # Only training and validation metrics used here; comparison values released below after selection.
  record=dict(name=name,total_response=total,free=keys,parameters=p,optimizer=optimization,fit_objective=objective,
   train=summarize([r for r in rows if r['split']=='train']),validation=summarize([r for r in rows if r['split']=='validation']),
   bounds_touched=[k for k in keys if min(abs(p[k]-np.array(BOUNDS[k])))<1e-4])
  print('FIT',objective,name,json.dumps({k:record[k] for k in ('parameters','train','validation','bounds_touched')}),flush=True)
  records.append(record)
 return records

def inverse_diagnostics(S):
 E=S.E;out=[]
 for i,g in enumerate(E.galaxies):
  ix=E.ids==i;r=E.r[ix];vb2=S.vb2[ix];y=E.y[ix]
  models={}
  for mode in ('total_strength','companion_strength','spreading','strength_and_spread'):
   keys=2 if mode=='strength_and_spread' else 1
   def calc(x):
    a=x[0] if mode!='spreading' else 0.;b=x[1] if keys==2 else x[0] if mode=='spreading' else 0.
    s=np.exp(-a);lam=np.exp(b);gc=S.J.companion_force(r,S.A[i]*s/lam,S.rc[i]*lam,S.rt[i]*lam,S.q[i])
    return np.sqrt((s if mode=='total_strength' else 1.)*vb2+r*gc)
   starts=[np.full(keys,j) for j in (.01,.6,1.8)]
   fits=[least_squares(lambda x:(calc(x)-y)/y,x,bounds=(np.zeros(keys),np.full(keys,np.log(10.))),max_nfev=150,gtol=1e-10,ftol=1e-10,xtol=1e-10) for x in starts]
   f=min(fits,key=lambda z:z.fun@z.fun);a=f.x[0] if mode!='spreading' else 0.;b=f.x[1] if keys==2 else f.x[0] if mode=='spreading' else 0.
   models[mode]=dict(strength=float(np.exp(-a)),spread=float(np.exp(b)),bound=bool(np.any(f.x>np.log(10.)-1e-4)),success=bool(f.success),**score(y,calc(f.x),E.error[ix]))
  out.append(dict(name=g['name'],split=g['split'],Sigma=float(S.Sigma[i]),fgas=float(S.fgas[i]),original=S.base_rows[i]['fractional_RMS'],over=bool(S.original_over[i]),under=bool(S.original_under[i]),models=models))
 groups={}
 for mode in out[0]['models']:
  groups[mode]={}
  for group,mask in [('all',np.ones(149,bool)),('over',S.original_over),('under',S.original_under)]:
   rr=[z for z,m in zip(out,mask) if m];vals=[z['models'][mode] for z in rr]
   groups[mode][group]=dict(n=len(rr),below20=sum(z['fractional_RMS']<.2 for z in vals),median_RMS=float(np.median([z['fractional_RMS'] for z in vals])),median_strength=float(np.median([z['strength'] for z in vals])),median_spread=float(np.median([z['spread'] for z in vals])),upper_bounds=sum(z['bound'] for z in vals))
 strengths=np.array([z['models']['total_strength']['strength'] for z in out]);corr=spearmanr(strengths,S.fgas)
 return dict(rows=out,groups=groups,strength_gas_fraction_spearman=float(corr.statistic),scope='Per-object feasibility fits, not a universal law or causal evidence')

def gas_checks():
 # Uniform cold sphere shell collapse, GM=R0=1. Stop before r=0 singularity.
 stop_r=.001;rows=[]
 for s in (1.,.5,.25,.1):
  def event(t,y):return y[0]-stop_r
  event.terminal=True;event.direction=-1
  f=solve_ivp(lambda t,y:[y[1],-s/y[0]**2],(0,20),[1.,0.],events=event,rtol=2e-10,atol=1e-12,max_step=.01)
  measured=float(f.t_events[0][0]);analytic=(np.arccos(np.sqrt(stop_r))+np.sqrt(stop_r*(1-stop_r)))/np.sqrt(2*s)
  rows.append(dict(s=s,collapse_time=measured,analytic=analytic,relative_error=measured/analytic-1,tff_ratio=s**-.5,Jeans_mass_ratio=s**-1.5,
   escape_speed_ratio=np.sqrt(s),Q_background_only=np.sqrt(s),Q_background_and_cloud=s**-.5))
 return dict(rows=rows,assumptions='Fixed density cold spherical cloud and separately a fixed-shape scaled galactic potential. No measured SFR or cloud coupling was fitted.')

def main():
 ap=argparse.ArgumentParser(__doc__);ap.add_argument('--jr1',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=args.output
 if (out/'summary.json').exists():raise FileExistsError('Preserve completed results')
 out.mkdir(parents=True,exist_ok=True);start=time.monotonic();S=Study(args.jr1)
 saved={r['name']:np.array(r['predicted_kms']) for r in S.saved['metrics']['sparc_rows']}
 baseline_error=max(float(np.max(abs(np.array(r['predicted_kms'])-saved[r['name']]))) for r in S.base_rows)
 print('BASELINE',baseline_error,flush=True)
 dump(out/'manifest.json',dict(baseline_commit='7386427ad7228b8d57319b1be78bbea3c47f31e8',inputs=S.hashes,base_parameters=S.p,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),specifications=SPECS,bounds=BOUNDS,Python=platform.python_version()))
 primary=fit_forms(S,'fractional');sensitivity=fit_forms(S,'quoted_errors')
 selected=min(primary,key=lambda z:z['validation']['mean_fractional_square'])
 dump(out/'selection-before-comparison.json',dict(primary=primary,sensitivity=sensitivity,selected=selected['name'],criterion='Primary fractional fit: minimum validation mean equal-galaxy fractional square; original baseline eligible.'))
 full={}
 for rec in primary:
  name=rec['name'];p=rec['parameters'];v=S.prediction(p,name,rec['total_response']);rows=S.rows(v,p,name)
  lenses=[S.lens(l,p,name,rec['total_response']) for l in S.E.lenses]
  full[name]=dict(parameters=p,total_response=rec['total_response'],stats=S.stats(rows),galaxies=rows,lenses=lenses)
 dump(out/'all-model-predictions.json',full)
 inverse=inverse_diagnostics(S);dump(out/'individual-diagnostics.json',inverse)
 gas=gas_checks();dump(out/'gas-physics-controls.json',gas)
 # Selected lens refinement only; immutable parameters.
 fine=S.P.Experiment(n=3201,order=128,deproj_order=256);refinement=[]
 for l,coarse in zip(fine.lenses,full[selected['name']]['lenses']):
  ff=S.lens(l,selected['parameters'],selected['name'],rec['total_response'] if False else selected['total_response'])
  refinement.append(dict(name=l.name,stellar_max_relative=float(np.max(abs(np.array(ff['predicted_vrms'])/coarse['predicted_vrms']-1))),angle_relative=float(ff['theta_predicted']/coarse['theta_predicted']-1)))
 # Exact preservation of total companion normalization under pure dilation.
 lam=3.;mass_error=float(np.max(abs((S.A/lam)*(S.rt*lam)/(S.A*S.rt)-1)))
 # Change in companion potential depth at Re under selected law.
 s,ls=modifiers(S.Sigma,selected['parameters'],selected['name']);depth=[]
 for i,g in enumerate(S.E.galaxies):
  re=g['Re'];a,rc,rt,q=S.A[i],S.rc[i],S.rt[i],S.q[i]
  old=quad(lambda t:float(S.J.companion_force(re/np.cos(t),a,rc,rt,q)*re*np.sin(t)/np.cos(t)**2),0,np.pi/2,epsabs=1e-5,epsrel=1e-9)[0]
  new=quad(lambda t:float(S.J.companion_force(re/np.cos(t),s[i]*a/ls[i],rc*ls[i],rt*ls[i],q)*re*np.sin(t)/np.cos(t)**2),0,np.pi/2,epsabs=1e-5,epsrel=1e-9)[0]
  depth.append(dict(name=g['name'],companion_escape_speed_ratio=float(np.sqrt(new/old)),strength=float(s[i]),spread=float(ls[i])))
 dump(out/'selected-potential-depth.json',depth)
 hash_ok=all(hashlib.sha256((S.J.ROOT/k).read_bytes()).hexdigest()==h for k,h in S.hashes.items())
 mini={}
 for name,x in full.items():
  lr=x['lenses'];mini[name]=dict(parameters=x['parameters'],total_response=x['total_response'],stats=x['stats'],
   lens_stellar_mean_fractional=float(np.mean([r['stellar_fractional_RMS'] for r in lr])),lens_stellar_chi2=float(sum(r['stellar_chi2'] for r in lr)),
   lens_angle_fractional_RMS=float(np.sqrt(np.mean([r['theta_fractional_error']**2 for r in lr]))),lenses=lr)
 allsel=full[selected['name']]['galaxies'];groups={}
 for label,mask in [('all',np.ones(149,bool)),('over_outliers',S.original_over),('under_outliers',S.original_under),('gas_rich',S.fgas>.5),('gas_poor',S.fgas<=.5)]:
  rr=[r for r,m in zip(allsel,mask) if m];groups[label]=dict(n=len(rr),median_strength=float(np.median([r['strength'] for r in rr])),median_spread=float(np.median([r['spread'] for r in rr])),median_fgas=float(np.median([r['fgas'] for r in rr])),stats=summarize(rr))
 summary=dict(experiment='JR-4',selected=selected['name'],primary=mini,fit_records=primary,sensitivity=sensitivity,
  baseline_replay_error_kms=baseline_error,original_over_outliers=int(S.original_over.sum()),original_under_outliers=int(S.original_under.sum()),
  selected_groups=groups,inverse_summary=inverse['groups'],inverse_gas_spearman=inverse['strength_gas_fraction_spearman'],
  lens_refinement=refinement,mass_preservation_error=mass_error,input_hashes_unchanged=hash_ok,gas_controls=gas,
  elapsed_seconds=time.monotonic()-start,limitations=['All observations previously exposed; development not blind confirmation','No actual SFR, gas evolution, molecular-cloud perturbation response, or cluster calculation','Effective static environmental normalization and redistribution are not microscopic gravity derivations','Original R10 stellar populations and conditional lens geometry unchanged'])
 dump(out/'summary.json',summary);print('DONE',selected['name'],json.dumps({k:v['stats'] for k,v in mini.items()}),flush=True)
if __name__=='__main__':main()
