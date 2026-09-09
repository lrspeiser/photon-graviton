"""Exploratory one-parameter history fit; already-inspected Pantheon+ data.
Run from extracted bundle with NumPy, SciPy and pandas installed.
"""
from pathlib import Path
import sys,contextlib,io,json
import numpy as np
from scipy.optimize import minimize_scalar,brentq
from scipy.linalg import cho_factor,cho_solve
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'shared_interaction_test'))
with contextlib.redirect_stdout(io.StringIO()):import brightness as b
x=np.log1p(b.z);N=len(x);one=np.ones(N)
def history(p):
 y=p*x
 return np.ones_like(x) if abs(p)<1e-12 else np.expm1(y)/y
def curvature(a):
 return np.ones_like(x) if a<1e-12 else np.sinh(a*x)/(a*x)
def evaluate(f,offset=False,sel=None):
 rr=b.r-5*np.log10(f)
 if sel is None: fac=b.vf;oo=one
 else:rr=rr[sel];fac=cho_factor(b.V[np.ix_(sel,sel)]);oo=np.ones(len(sel))
 v1=cho_solve(fac,oo);off=float(v1@rr/(oo@v1)) if offset else 0.
 rr=rr-off
 return float(rr@cho_solve(fac,rr)),off

def fit(fun,bounds,offset=False,sel=None):
 obj=lambda t:evaluate(fun(t),offset,sel)[0]
 s=minimize_scalar(obj,bounds=bounds,method='bounded',options={'xatol':1e-10})
 v=float(s.x);cc,off=evaluate(fun(v),offset,sel)
 grid=np.linspace(*bounds,401);assert cc<=min(obj(t) for t in grid)+1e-6
 ints={}
 for d in [1.,3.841458820694124]:
  target=lambda t:obj(t)-cc-d
  ints[str(d)]=[float(brentq(target,bounds[0],v)) if target(bounds[0])>0 else None,float(brentq(target,v,bounds[1])) if target(bounds[1])>0 else None]
 return dict(parameter=v,chi_squared=cc,offset_mag=off,conditional_profile_intervals=ints)
W=np.array(b.binweights);bc=W@b.V@W.T;bf=cho_factor(bc);bo=np.ones(4);bw=cho_solve(bf,bo);bw/=sum(bw)
def diagnostics(f,off=0):
 rr=b.r-5*np.log10(f)-off;means=W@rr;mean=float(bw@means);d=means-mean
 return dict(bin_mean_residuals_mag=means.tolist(),bin_constant_residual_chi_squared=float(d@cho_solve(bf,d)),
  note='Post-fit diagnostic on four correlated bins; no naive chi-square significance assigned after fitting.',
  last_minus_first_bin_mag=float(means[-1]-means[0]),last_minus_first_standard_error_before_parameter_fit=float(np.sqrt(bc[-1,-1]+bc[0,0]-2*bc[-1,0])),rmse_mag=float(np.sqrt(np.mean(rr**2))))
base=dict(parameter=None,chi_squared=evaluate(np.ones(N))[0],offset_mag=0.)
cal=dict(parameter=None,chi_squared=evaluate(np.ones(N),True)[0],offset_mag=evaluate(np.ones(N),True)[1])
h=fit(history,(-3,3));curv=fit(curvature,(0,3));hc=fit(history,(-3,3),True)
models={'fixed_exponential_flat':(base,np.ones(N),0),'offset_only':(cal,np.ones(N),1),'history_only':(h,history(h['parameter']),1),'negative_curvature_only':(curv,curvature(curv['parameter']),1),'history_plus_offset_diagnostic':(hc,history(hc['parameter']),2)}
for name,(r,f,k) in models.items():
 r.update(diagnostics(f,r['offset_mag']));r['added_fit_parameters']=k;r['delta_AIC_vs_fixed']=r['chi_squared']+2*k-base['chi_squared']
sens={}
for name,mask in [('z_below_0_6',b.z<.6),('z_at_least_0_6',b.z>=.6),('z_below_1',b.z<1)]:
 ids=np.flatnonzero(mask);sens[name]=dict(N=len(ids),**fit(history,(-3,3),sel=ids))
# Check distance integration independently of closed expression.
from scipy.integrate import quad
checks=[]
for pp in [0.,h['parameter'],-.5,.5]:
 for zz in [.1,1.,2.]:
  numerical=quad(lambda n:n**(-pp-1),1/(1+zz),1)[0]
  closed=np.log1p(zz) if pp==0 else np.expm1(pp*np.log1p(zz))/pp
  checks.append(abs(numerical-closed))
assert max(checks)<1e-10
out=dict(status='Exploratory reuse, not independent validation',N=N,kappa_per_Mly=b.cfg['kappa_per_Mly'],M=b.M,
 law='dn/dt=gamma*n^p, n_o=1; R=expm1(p*ln(1+zHD))/(kappa*p); DL=(1+zHEL)*R; DA=R in ideal rest-frame flat model.',
 models={n:v[0] for n,v in models.items()},subset_sensitivity_not_holdouts=sens,
 maximum_integral_check_error=max(checks),
 caveats=['Full released STAT+SYS covariance and Cepheid calibration propagated as in previous test.', 'Standardized published magnitudes and reference-frame corrections are conditional on their original processing.', 'Kappa fixed without its uncertainty; profile intervals conditional on this model and covariance.', 'AIC is descriptive within this likelihood, not validation or a probability of physical truth.', 'No new atomic interaction, angular-distance dataset or sustainable dynamical background tested.'])
(P/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
