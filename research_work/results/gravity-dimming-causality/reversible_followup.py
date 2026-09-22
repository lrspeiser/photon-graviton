#!/usr/bin/env python3
"""A two-parameter, energy-conserving reversible channel, not free galaxy opacities."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.linalg import expm
from test_dimming import brief,dump

def main():
 ap=argparse.ArgumentParser(__doc__);ap.add_argument('--jr2',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 target=args.output/'reversible-shared.json'
 if target.exists():raise FileExistsError('Preserve previous results')
 bundle=args.jr2/'original-JR1';sys.path.insert(0,str(bundle/'repository/research_work/results/joint-response-iteration'))
 import population_followup as P
 J=P.J;old=json.loads((bundle/'evidence/R10_spheroid_stellar_population.json').read_text());p=old['parameters'];E=P.Experiment()
 ids=E.ids;fs=E.Msph/E.Ms;us=10**p['logusph'];qn=p['q']+(p['q_sph']-p['q'])*us*fs/(1+(us-1)*fs)
 vb=10**p['logu']*(E.vs2+(us-1)*E.bulge_v2)
 vre=np.array([np.interp(g['Re'],g['r'],vb[ids==i]) for i,g in enumerate(E.galaxies)])
 vge=np.array([np.interp(g['Re'],g['r'],g['gas_v2']) for g in E.galaxies])
 size=np.bincount(ids);mask=E.train;calls=0
 def state(lu):
  A,rc,rt=P.source(E.Ms*np.exp(lu),E.Mg,E.Re,p,False,fs)
  vc=E.Re*J.companion_force(E.Re,A,rc,rt,qn)
  H=vc/(vc+np.maximum(vre*np.exp(lu)+vge,0))
  return H,A,rc,rt
 def transport(lu,tau,R):
  H=state(lu)[0];X=R*H;T=1+.5*(X-1)*(-np.expm1(-tau*H))
  return T,X,H
 def solve(tau,R):
  lo=np.full(149,-np.log1p(R));hi=np.full(149,np.log(2.))
  for _ in range(48):
   mid=(lo+hi)/2;T,_,_=transport(mid,tau,R);v=mid+np.log(T)
   lo=np.where(v<0,mid,lo);hi=np.where(v>=0,mid,hi)
  lu=(hi+lo)/2;T,X,H=transport(lu,tau,R);_,A,rc,rt=state(lu)
  v=np.sqrt(np.maximum(E.vg2+vb*np.exp(lu)[ids],0)+E.r*J.companion_force(E.r,A[ids],rc[ids],rt[ids],qn[ids]))
  return lu,T,X,H,v
 def residual(x):
  nonlocal calls
  calls+=1;v=solve(*x)[-1]
  return ((v-E.y)/E.y/np.sqrt(size[ids]*89))[mask]
 starts=[(0.,0.),(.5,.5),(2.,2.),(5.,1.),(5.,5.),(9.,8.)];fits=[];records=[]
 for i,x in enumerate(starts):
  f=least_squares(residual,np.clip(x,1e-7,None),bounds=([0,0],[10,9]),max_nfev=150,ftol=1e-10,xtol=1e-10,gtol=1e-9,diff_step=1e-4)
  fits.append(f);records.append(dict(start=x,parameters=f.x.tolist(),cost=float(f.fun@f.fun),nfev=f.nfev,success=bool(f.success)))
  print('SHARED REVERSE START',i,records[-1],flush=True)
 best=min(fits,key=lambda f:float(f.fun@f.fun));tau,R=map(float,best.x);lu,T,X,H,v=solve(tau,R)
 rows=[]
 for i,g in enumerate(E.galaxies):
  ix=ids==i;r=v[ix]-g['observed'];rows.append(dict(name=g['name'],split=g['split'],fractional_RMS=float(np.sqrt(np.mean((r/g['observed'])**2))),
   RMSE_kms=float(np.sqrt(np.mean(r*r))),chi2=float(np.sum((r/g['error'])**2)),
   intrinsic_light_factor=float(np.exp(lu[i])),net_dimming_fraction=float(1-T[i]),
   incoming_companion_flux_per_primary_photon_flux=float(X[i]),outgoing_companion_flux_per_primary_photon_flux=float(1+X[i]-T[i]),
   companion_fraction_at_Re=float(H[i]),predicted_kms=v[ix].tolist()))
 # Global coarse root scan on final parameters, all galaxies.
 grid=np.linspace(-np.log1p(R),np.log(2.),130);tt=transport(grid[:,None],tau,R)[0];rr=grid[:,None]+np.log(tt)
 nroots=np.sum((rr[:-1]<=0)&(rr[1:]>=0),axis=0)
 exact=[]
 for i in range(149):
  op=.5*tau*H[i]*np.array([[-1.,1.],[1.,-1.]])
  ex=expm(op)@np.array([1.,X[i]])
  exact.append(float(np.max(abs(ex-[T[i],1+X[i]-T[i]]))))
 result=dict(experiment='JR-3 reversible follow-up',post_primary_development=True,
  parameters=dict(tau=tau,incoming_companion_coefficient=R),optimizer_starts=records,forward_evaluations=calls,
  all=brief(rows),splits={s:brief([r for r in rows if r['split']==s]) for s in ('train','validation','test')},
  count_net_dimming=int(np.sum(T<1-1e-5)),count_net_brightening=int(np.sum(T>1+1e-5)),
  median_photon_flux_multiplier=float(np.median(T)),p10_photon_flux_multiplier=float(np.percentile(T,10)),p90_photon_flux_multiplier=float(np.percentile(T,90)),
  checks=dict(max_selfconsistency_error=float(np.max(abs(np.exp(lu)*T-1))),min_roots=int(nroots.min()),max_roots=int(nroots.max()),
    energy_ledger_max_error=float(np.max(abs(T+(1+X-T)-(1+X)))),matrix_exponential_max_error=max(exact)),
  rows=rows,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  limits=['Incoming companion flux X=R*H is hypothesized, not independently observed or derived from reservoir transport.',
   'Equal-rate two-channel forward solution conserves energy but does not close galactic source replenishment.',
   'Uniform band flux changes at fixed populations; no multiband or dynamical formation model.',
   'Existing validation/comparison data were previously exposed. A better score is not proof of causality.'])
 dump(target,result);print('REVERSE FINAL',json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
if __name__=='__main__':main()
