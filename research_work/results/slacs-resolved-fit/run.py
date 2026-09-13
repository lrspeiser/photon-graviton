from pathlib import Path
import json,hashlib,importlib.util
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import solve_triangular
from model import AnnularModel
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];PILOT=HERE.parent/'slacs-motion-lensing-pilot'
spec=importlib.util.spec_from_file_location('pilot',PILOT/'model.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
paths=[HERE.parent/'slacs-resolved-input-audit/results.json',PILOT/'results.json',HERE.parent/'lensing-data-readiness/conditional-geometry.json',HERE/'protocol.json']
data=json.loads(paths[0].read_text(encoding='utf-8'));base=json.loads(paths[1].read_text(encoding='utf-8'));geo={r['Name']:r for r in json.loads(paths[2].read_text(encoding='utf-8'))};protocol=json.loads(paths[3].read_text(encoding='utf-8'))
par=base['parameters'];A,p,astar=par['A'],par['p'],par['a_star_m_s2']*3.085677581491367e19/1e6
pilot={r['Name']:r for r in base['rows'] if r['model']=='empirical_extra' and r['cutoff_in_a']==20}
rows=[];excluded=[]
for item in data['systems']:
 if float(item['release_use_flag'])!=1:
  excluded.append(dict(Name=item['Name'],reason='Release use flag is zero; underlying reason unresolved'));continue
 name=item['Name'];a=pilot[name]['scale_a_kpc'];Dl=geo[name]['conditional_Dl_Mpc']*1000;ratio=geo[name]['conditional_Dls_over_Ds'];edges=np.r_[item['inner_arcsec'],item['outer_arcsec'][-1]]*Dl/old.ARCSEC
 psf=item['psf_fwhm_arcsec']*Dl/old.ARCSEC/np.sqrt(8*np.log(2));model=AnnularModel(a,edges,psf,A,p,astar,20)
 y=np.array(item['vrms_kms']);cov=np.array(item['covariance_kms_squared']);chol=np.linalg.cholesky(cov)
 const=len(y)*np.log(2*np.pi)+2*np.log(np.diag(chol)).sum()
 fitted=[]
 for label,extra in [('baryons',False),('empirical_extra',True)]:
  def objective(q):
   residual=y-model.predict(np.exp(q[0]),q[1],extra);v=solve_triangular(chol,residual,lower=True);return float(v@v)
  starts=[[np.log(pilot[name]['inferred_stellar_mass_Msun']),b] for b in protocol['starts_beta']]
  fits=[minimize(objective,start,bounds=[np.log(protocol['mass_Msun_bounds']),protocol['constant_beta_bounds']],method='L-BFGS-B',options={'ftol':1e-11,'maxiter':1000}) for start in starts]
  good=[f for f in fits if f.success and np.isfinite(f.fun)];assert good
  best=min(good,key=lambda f:f.fun);mass=np.exp(best.x[0]);beta=best.x[1];prediction=model.predict(mass,beta,extra)
  theta=old.angle(mass,a,Dl,ratio,A if extra else 0,p,astar,20)
  row=dict(Name=name,model=label,n_bins=len(y),mass_Msun=float(mass),beta=float(beta),beta_at_boundary=bool(min(abs(beta-np.array(protocol['constant_beta_bounds'])))<1e-5),chi2=float(best.fun),nominal_bins_minus_parameters=len(y)-2,normalized_log_likelihood=float(-.5*(best.fun+const)),predicted_vrms_kms=prediction.tolist(),observed_vrms_kms=y.tolist(),predicted_angle_arcsec=theta,catalog_angle_arcsec=pilot[name]['catalog_SIE_arcsec'],angle_ratio=theta/pilot[name]['catalog_SIE_arcsec'] if theta else None,starts=[dict(success=bool(f.success),chi2=float(f.fun),parameters=f.x.tolist()) for f in fits])
  fitted.append((row,extra))
 fine=AnnularModel(a,edges,psf,A,p,astar,20,n=8001,order=256)
 for row,extra in fitted:
  prediction=fine.predict(row['mass_Msun'],row['beta'],extra);coarse=np.array(row['predicted_vrms_kms']);change=prediction-coarse
  row['maximum_relative_vrms_refinement']=float(np.max(abs(change/coarse)));row['refinement_whitened_norm']=float(np.linalg.norm(solve_triangular(chol,change,lower=True)));rows.append(row)
 print(json.dumps([{k:v for k,v in row.items() if k in ['Name','model','beta','beta_at_boundary','chi2','angle_ratio','maximum_relative_vrms_refinement']} for row,_ in fitted]),flush=True)
(HERE/'results.json').write_text(json.dumps(dict(scope=protocol['scope'],rows=rows,excluded=excluded,hashes={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths+[HERE/'run.py',HERE/'model.py',PILOT/'model.py']}),indent=2)+'\n',encoding='utf-8',newline='\n')
