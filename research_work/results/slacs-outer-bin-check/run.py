from pathlib import Path
import json,hashlib,importlib.util
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import solve_triangular
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'slacs-component-refit'))
from model import ComponentModel
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];PILOT=HERE.parent/'slacs-motion-lensing-pilot'
spec=importlib.util.spec_from_file_location('pilot',PILOT/'model.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
paths=[HERE.parent/'slacs-resolved-input-audit/results.json',PILOT/'results.json',HERE.parent/'lensing-data-readiness/conditional-geometry.json',HERE/'protocol.json']
data=json.loads(paths[0].read_text(encoding='utf-8'));base=json.loads(paths[1].read_text(encoding='utf-8'));geo={r['Name']:r for r in json.loads(paths[2].read_text(encoding='utf-8'))};protocol=json.loads(paths[3].read_text(encoding='utf-8'))
par=base['parameters'];A,p,astar=par['A'],par['p'],par['a_star_m_s2']*3.085677581491367e19/1e6
pilot={r['Name']:r for r in base['rows'] if r['model']=='empirical_extra' and r['cutoff_in_a']==20}
profiles_path=HERE.parent/'slacs-light-profile-audit/results.json'
profiles={r['Name']:r for r in json.loads(profiles_path.read_text(encoding='utf-8'))['rows']}
paths.append(profiles_path)
rows=[];excluded=[]
for item in data['systems']:
 if float(item['release_use_flag'])!=1:
  excluded.append(dict(Name=item['Name'],reason='Release use flag is zero; underlying reason unresolved'));continue
 name=item['Name']
 if profiles.get(name,{}).get('status')!='conditional_components_available':
  excluded.append(dict(Name=name,reason='Published component parameters unavailable'));continue
 a=pilot[name]['scale_a_kpc'];Dl=geo[name]['conditional_Dl_Mpc']*1000;ratio=geo[name]['conditional_Dls_over_Ds'];edges=np.r_[item['inner_arcsec'],item['outer_arcsec'][-1]]*Dl/old.ARCSEC
 psf=item['psf_fwhm_arcsec']*Dl/old.ARCSEC/np.sqrt(8*np.log(2));components=[dict(R=c['R_arcsec']*Dl/old.ARCSEC,n=c['n'],amp=c['amp_at_R'],bn=c['bn']) for c in profiles[name]['components']]
 model=ComponentModel(a,edges,psf,A,p,astar,20,components)
 y=np.array(item['vrms_kms']);cov=np.array(item['covariance_kms_squared']);full_chol=np.linalg.cholesky(cov);chol=np.linalg.cholesky(cov[:-1,:-1])
 const=(len(y)-1)*np.log(2*np.pi)+2*np.log(np.diag(chol)).sum()
 fitted=[]
 for label,extra in [('baryons',False),('empirical_extra',True)]:
  def objective(q):
   residual=y[:-1]-model.predict(np.exp(q[0]),q[1],extra)[:-1];v=solve_triangular(chol,residual,lower=True);return float(v@v)
  starts=[[np.log(1e11),b] for b in protocol['starts_beta']]
  fits=[minimize(objective,start,bounds=[np.log(protocol['mass_Msun_bounds']),protocol['constant_beta_bounds']],method='L-BFGS-B',options={'ftol':1e-11,'maxiter':1000}) for start in starts]
  good=[f for f in fits if f.success and np.isfinite(f.fun)];assert good
  best=min(good,key=lambda f:f.fun);mass=np.exp(best.x[0]);beta=best.x[1];prediction=model.predict(mass,beta,extra)
  theta=model.angle(mass,Dl,ratio,extra)
  cross=cov[-1,:-1];weight=np.linalg.solve(cov[:-1,:-1],cross)
  conditional_mean=prediction[-1]+weight@(y[:-1]-prediction[:-1])
  conditional_var=cov[-1,-1]-cross@weight;assert conditional_var>0
  conditional_residual=y[-1]-conditional_mean
  full_residual=solve_triangular(full_chol,y-prediction,lower=True)
  identity_error=abs(float(full_residual@full_residual)-best.fun-conditional_residual**2/conditional_var);assert identity_error<1e-7
  row=dict(Name=name,model=label,n_bins=len(y),mass_Msun=float(mass),beta=float(beta),beta_at_boundary=bool(min(abs(beta-np.array(protocol['constant_beta_bounds'])))<1e-5),chi2=float(best.fun),n_fitted_bins=len(y)-1,nominal_bins_minus_parameters=len(y)-3,outer_raw_prediction_kms=float(prediction[-1]),outer_observed_kms=float(y[-1]),outer_conditional_mean_kms=float(conditional_mean),outer_conditional_error_sigma_kms=float(np.sqrt(conditional_var)),outer_conditional_standardized_residual=float(conditional_residual/np.sqrt(conditional_var)),gaussian_factorization_error=float(identity_error),normalized_log_likelihood=float(-.5*(best.fun+const)),predicted_vrms_kms=prediction.tolist(),observed_vrms_kms=y.tolist(),predicted_angle_arcsec=theta,catalog_angle_arcsec=pilot[name]['catalog_SIE_arcsec'],angle_ratio=theta/pilot[name]['catalog_SIE_arcsec'] if theta else None,starts=[dict(success=bool(f.success),chi2=float(f.fun),parameters=f.x.tolist()) for f in fits])
  fitted.append((row,extra))
 fine=ComponentModel(a,edges,psf,A,p,astar,20,components,n=8001,order=256,deproj_order=512)
 for row,extra in fitted:
  prediction=fine.predict(row['mass_Msun'],row['beta'],extra);coarse=np.array(row['predicted_vrms_kms']);change=prediction-coarse
  row['coarse_total_light_ratio']=model.total_mass_check;row['fine_total_light_ratio']=fine.total_mass_check
  row['maximum_relative_vrms_refinement']=float(np.max(abs(change/coarse)));row['refinement_whitened_norm']=float(np.linalg.norm(solve_triangular(full_chol,change,lower=True)));rows.append(row)
 print(json.dumps([{k:v for k,v in row.items() if k in ['Name','model','beta','beta_at_boundary','chi2','outer_conditional_standardized_residual','maximum_relative_vrms_refinement']} for row,_ in fitted]),flush=True)
(HERE/'results.json').write_text(json.dumps(dict(scope=protocol['scope'],rows=rows,excluded=excluded,hashes={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths+[HERE/'run.py',HERE.parent/'slacs-component-refit/model.py',PILOT/'model.py',HERE.parent/'slacs-resolved-fit/model.py']}),indent=2)+'\n',encoding='utf-8',newline='\n')
