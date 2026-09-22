"""Post-primary JR-2 diagnostic: fractional-RMS fit and exact-ring orbital adjustment.
No formula coefficient or observation changes. Minimizing fractional RMS instead
of quoted-error chi2 is a diagnostic objective, not altered measurement errors.
"""
import sys,json,io,zipfile
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares,minimize_scalar
from scipy.linalg import solve_triangular
from inverse_properties import Disk,catalog,dump
root=Path(sys.argv[1]);out=Path(sys.argv[2]);sys.path.insert(0,str(root/'repository/research_work/results/joint-response-iteration'))
import population_followup as P
J=P.J;saved=json.loads((root/'evidence/R10_spheroid_stellar_population.json').read_text());p=saved['parameters'];spec=saved['specification'];E=P.Experiment();cat=catalog(J.ROOT/J.INPUT_NAMES[1]);old=json.loads((out/'disk-results.json').read_text());a=[]
with zipfile.ZipFile(J.ROOT/J.INPUT_NAMES[0]) as z:
 for g in E.galaxies:
  raw=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(g['name']+'_rotmod.dat'))));raw=raw[raw[:,0]>0];d=Disk(g,cat[g['name']],raw,J,P,p)
  oldr=next(r for r in old if r['name']==g['name']);lo=np.array([np.log10(.2),5.,-1.,np.log10(.25)]);hi=np.array([np.log10(5.),90.,1.,np.log10(4.)])
  def res(x):return (d.forward(x)[1]-d.proj)/d.proj
  starts=[d.base,np.array(oldr['cases']['restricted_joint']['x']),np.array(oldr['cases']['wide_joint']['x'])]
  best=min((least_squares(res,np.clip(x,lo+1e-8,hi-1e-8),bounds=(lo,hi),max_nfev=300,
      ftol=1e-11,xtol=1e-11,gtol=1e-9,x_scale='jac') for x in starts),key=lambda f:f.fun@f.fun)
  a.append(dict(name=g['name'],original_outlier=oldr['original_outlier'],**d.metrics(best.x),success=bool(best.success),
    x=best.x.tolist(),bound_indices=[i for i in range(4) if min(best.x[i]-lo[i],hi[i]-best.x[i])<1e-5]))
lens_old=json.loads((out/'lens-results.json').read_text());lr=[]
for l in E.lenses:
 o=next(r for r in lens_old if r['name']==l.name);base=l.local_parameters(p,spec);pd=p.copy();pd['logu']=base['logu']+o['cases']['mass_for_ring']['log_stellar_mass_shift_dex'];pd['logusph']=0.
 ss=dict(spec,stellar_nuisance=False)
 def r(beta):
  pd['beta']=beta;v,f,pars=l.prediction(pd,ss);w=solve_triangular(l.chol,v-l.y,lower=True);return w@w
 fit=minimize_scalar(r,bounds=(-1,.35),method='bounded',options={'xatol':1e-10});pd['beta']=fit.x;v,f,pars=l.prediction(pd,ss)
 lr.append(dict(name=l.name,stellar_mass_factor=o['cases']['mass_for_ring']['stellar_mass_factor'],beta=float(fit.x),
   stellar_chi2=float(fit.fun),stellar_fractional_RMS=float(np.sqrt(np.mean(((v-l.y)/l.y)**2))),
   theta_predicted_arcsec=float(l.angle(pd,ss)),theta_observed_arcsec=l.theta,beta_on_bound=bool(fit.x>.34999 or fit.x<-.99999)))
result=dict(scope='Post-primary diagnostic; fractional-RMS objective is not a measurement-error likelihood; no gravity-law changes',
 disk_rows=a,lens_exact_ring_orbit_adjustment=lr,
 summary=dict(galaxies=len(a),within10pct=sum(x['fractional_RMS']<.1 for x in a),within20pct=sum(x['fractional_RMS']<.2 for x in a),
  original_outliers_rescued_below20pct=sum(x['original_outlier'] and x['fractional_RMS']<.2 for x in a),
  median_fractional_RMS=float(np.median([x['fractional_RMS'] for x in a])),
  still_above20pct=[x['name'] for x in a if x['fractional_RMS']>=.2],
  bounds_touched=sum(bool(x['bound_indices']) for x in a),
  raw_chi2=sum(x['chi2'] for x in a)))
dump(out/'supplement.json',result);print(json.dumps(result['summary']));print(json.dumps(lr))
