"""Positive finite Coma inverse mass envelopes, no cosmological mass conversion.

python coma-inverse.py [--output-dir DIR] [--canonical]
Regenerates into a fresh directory and compares with the archived
coma-inverse-results.json; only --canonical overwrites the archive. The
archived file was produced by the pre-entry-point version at 3884b4f.
"""
from pathlib import Path
import json,hashlib,sys
import numpy as np
from scipy.optimize import linprog,minimize,minimize_scalar
from scipy.integrate import quad
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
import evidence_io  # noqa: E402
paths=[P.parent/'cluster-observation-readiness/kubo-figure-data.json',P/'coma-inverse-protocol.md',Path(__file__)]

def compute():
 raw=json.loads(paths[0].read_text());rows=raw['rows'];r=np.array([v['published_radius_h_inverse_Mpc'] for v in rows]);r/=r[0]
 y=np.array([v['shear_t'] for v in rows]);err=np.array([v['plotted_sigma_t'] for v in rows]);out=dict(scope='Total gravitating positive finite spherical distributions fitting exposed Coma weak-shear bins. Not baryon-subtracted companions, physical mass, reduced-shear likelihood or derived transport.',radius_ratios=r.tolist(),observed_shear=y.tolist(),plotted_errors=err.tolist(),cases=[])
 checks=[]
 for R in [.1,1,10]:
  # Unit Plummer with total mass pi has rho=3/4*(1+r^2)^(-5/2).
  surf=2*quad(lambda z:.75*(1+R*R+z*z)**-2.5,0,np.inf,epsabs=1e-12)[0]
  checks.append(abs(surf/(1+R*R)**-2-1))
 assert max(checks)<1e-8
 for extent in [10,30,100]:
  a=np.geomspace(.1,extent*r[-1],25);A=a[None,:];R=r[:,None]
  shear=A*A*R*R/(R*R+A*A)**2;kappa=A**4/(R*R+A*A)**2
  B=shear/err[:,None];t=y/err
  # Positive weights have natural convergence units; objectives scaled for LP conditioning.
  lhs=np.vstack([shear,-shear,kappa]);rhs=np.r_[y+2*err,-y+2*err,np.full(len(r),.05)]
  def lp(obj):
   res=linprog(obj/max(abs(obj)),A_ub=lhs,b_ub=rhs,bounds=(0,None),method='highs')
   assert res.success,res.message
   assert max(lhs@res.x-rhs)<1e-7
   return res.x
  initial=lp(a*a)
  res=minimize(lambda w:float(np.sum((B@w-t)**2)),initial,jac=lambda w:2*B.T@(B@w-t),method='SLSQP',bounds=[(0,None)]*len(a),constraints=[{'type':'ineq','fun':lambda w:.05-kappa@w,'jac':lambda w:-kappa}],options={'maxiter':2000,'ftol':1e-11})
  assert res.success,res.message
  assert min(res.x)>-1e-8 and max(kappa@res.x)<.0500001
  def describe(w):
   gam=shear@w;kap=kappa@w;enc=a*a*r[-1]**3/(r[-1]**2+a*a)**1.5
   return dict(weights=w.tolist(),predicted_shear=gam.tolist(),predicted_kappa=kap.tolist(),chi2=float(sum(((gam-y)/err)**2)),total_mass_units=float(sum(w*a*a)),mass_inside_last_units=float(sum(w*enc)),max_reduced_vs_weak_difference_sigma=float(max(abs(gam/(1-kap)-gam)/err)))
  enc=a*a*r[-1]**3/(r[-1]**2+a*a)**1.5
  case=dict(max_scale_over_last_radius=extent,scales=a.tolist(),best=describe(res.x),minimum_total=describe(lp(a*a)),maximum_total=describe(lp(-a*a)),minimum_enclosed=describe(lp(enc)),maximum_enclosed=describe(lp(-enc)))
  out['cases'].append(case)
  print(extent,'chi2',case['best']['chi2'],'total',[case[k]['total_mass_units'] for k in ['minimum_total','maximum_total']],'enclosed',[case[k]['mass_inside_last_units'] for k in ['minimum_enclosed','maximum_enclosed']],flush=True)
 out['max_projection_check_relative_error']=max(checks)
 # Post-envelope exploratory single-component fit. Known n=5 polytrope gives a
 # self-gravitating Plummer equilibrium; this does not derive its photon supply.
 def single(loga):
  a=np.exp(loga);s=a*a*r*r/(r*r+a*a)**2;k=a**4/(r*r+a*a)**2
  w=np.clip(np.sum(s*y/err**2)/np.sum((s/err)**2),0,.05/max(k))
  pred=w*s
  return float(sum(((pred-y)/err)**2)),float(w),pred
 grid=np.linspace(np.log(.1),np.log(100*r[-1]),401);losses=[single(v)[0] for v in grid]
 i=int(np.argmin(losses));opt=minimize_scalar(lambda z:single(z)[0],bounds=(grid[max(0,i-1)],grid[min(400,i+1)]),method='bounded')
 best=min([grid[0],grid[-1],float(opt.x)],key=lambda z:single(z)[0]);chi,w,pred=single(best);a=np.exp(best)
 out['exploratory_single_plummer']=dict(chi2=chi,scale_in_first_radius_units=float(a),central_kappa=w,total_mass_units=float(w*a*a),predicted_shear=pred.tolist(),note='Post-envelope exploratory two-parameter cluster fit, not a frozen galaxy-to-cluster prediction or a test of one universal polytropic constant.')
 print('single Plummer',out['exploratory_single_plummer'])
 rr=np.geomspace(1e-3,1e3,16001);rho=3/(4*np.pi)*(1+rr*rr)**-2.5
 pressure=rho/(6*np.sqrt(1+rr*rr));gravity=rr/(1+rr*rr)**1.5
 derivative=np.gradient(pressure,np.log(rr))/rr
 equilibrium_error=float(max(abs(derivative[10:-10]/(-rho[10:-10]*gravity[10:-10])-1)))
 total=quad(lambda x:3*x*x/(1+x*x)**2.5,0,np.inf,epsabs=1e-11)[0]
 assert equilibrium_error<1e-4 and abs(total-1)<1e-9
 out['stationary_support_check']=dict(unit_G_M_a=1,relative_force_balance_error=equilibrium_error,total_mass=total,pressure_law='P=K rho^(6/5), K=(4*pi/3)^(1/5)/6 for G=M=a=1',scope='Known self-gravitating Plummer stationary support, not a source-driven evolution or stability proof; ordinary matter omitted in this analytic control.')
 out['source_sha256']={str(p.relative_to(P.parents[2])).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
 return out

def main():
 args=evidence_io.parse(__doc__)
 out=compute()
 # The flexible-mixture SLSQP fits have a flat objective and non-unique weights: under different BLAS
 # threading they stop at slightly different points (chi2 within ~1e-5, derived masses within ~0.1%).
 # LP extrema, the single-Plummer fit and all checks are still compared exactly.
 rules=[(r'/cases\[\d+\]/best/chi2$',0.,5e-5),(r'/cases\[\d+\]/best/weights',None,None),(r'/cases\[\d+\]/best/',5e-3,1e-9)]
 return evidence_io.finish(args,'coma-inverse',json.dumps(out,indent=2)+'\n',P/'coma-inverse-results.json',
                           ignore=['/source_sha256/research_work/results/companion-extensions/coma-inverse.py'],rules=rules)

if __name__=='__main__':
 raise SystemExit(main())
