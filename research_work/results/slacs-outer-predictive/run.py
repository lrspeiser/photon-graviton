from pathlib import Path
import json,hashlib,sys
import numpy as np
from scipy.special import logsumexp,ndtr
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'slacs-component-refit'))
from model import ComponentModel,ARCSEC
paths=[HERE/'protocol.json',HERE.parent/'slacs-resolved-input-audit/results.json',HERE.parent/'slacs-light-profile-audit/results.json',HERE.parent/'slacs-motion-lensing-pilot/results.json',HERE.parent/'lensing-data-readiness/conditional-geometry.json']
protocol=json.loads(paths[0].read_text(encoding='utf-8'));data=json.loads(paths[1].read_text(encoding='utf-8'))
profiles={r['Name']:r for r in json.loads(paths[2].read_text(encoding='utf-8'))['rows']};base=json.loads(paths[3].read_text(encoding='utf-8'));geo={r['Name']:r for r in json.loads(paths[4].read_text(encoding='utf-8'))}
pilot={r['Name']:r for r in base['rows'] if r['model']=='empirical_extra' and r['cutoff_in_a']==20};par=base['parameters'];A,p,astar=par['A'],par['p'],par['a_star_m_s2']*3.085677581491367e19/1e6
results=[]
for item in data['systems']:
 name=item['Name']
 if float(item['release_use_flag'])!=1 or profiles.get(name,{}).get('status')!='conditional_components_available':continue
 Dl=geo[name]['conditional_Dl_Mpc']*1000;a=pilot[name]['scale_a_kpc'];edges=np.r_[item['inner_arcsec'],item['outer_arcsec'][-1]]*Dl/ARCSEC;psf=item['psf_fwhm_arcsec']*Dl/ARCSEC/np.sqrt(8*np.log(2))
 comps=[dict(R=c['R_arcsec']*Dl/ARCSEC,n=c['n'],amp=c['amp_at_R'],bn=c['bn']) for c in profiles[name]['components']]
 model=ComponentModel(a,edges,psf,A,p,astar,20,comps)
 y=np.array(item['vrms_kms']);cov=np.array(item['covariance_kms_squared']);C=cov[:-1,:-1];inv=np.linalg.inv(C);cross=cov[-1,:-1];weight=inv@cross;noise_var=float(cov[-1,-1]-cross@weight);assert noise_var>0
 fine_mass,fine_beta=protocol['grids'][-1];lm=np.linspace(*np.log(protocol['mass_Msun_bounds']),fine_mass);betas=np.linspace(*protocol['beta_bounds'],fine_beta);m=np.exp(lm)/1e11
 coefficients=np.array([model.coefficients(b) for b in betas]);out=[]
 for label,extra in [('baryons',False),('empirical_extra',True)]:
  loglike=np.empty((fine_beta,fine_mass));means=np.empty_like(loglike)
  for j,(cb,cc) in enumerate(coefficients):
   pred=np.sqrt(m[:,None]*cb[None,:]+(m[:,None]**p*cc[None,:] if extra else 0))
   residual=y[:-1]-pred[:,:-1];loglike[j]=-.5*np.einsum('mi,ij,mj->m',residual,inv,residual)
   means[j]=pred[:,-1]+residual@weight
  for prior in protocol['priors']:
   summaries=[]
   for step in [2,1]:
    massgrid=lm[::step];betagrid=betas[::step];ll=loglike[::step,::step].copy();mu=means[::step,::step]
    if prior=='uniform_mass_and_beta':ll+=massgrid[None,:]
    # Product trapezoidal weights; common spacings cancel upon normalization.
    ll[[0,-1],:]+=np.log(.5);ll[:,[0,-1]]+=np.log(.5)
    prob=np.exp(ll-logsumexp(ll));mean=float(np.sum(prob*mu));variance=float(noise_var+np.sum(prob*(mu-mean)**2));cdf=float(np.sum(prob*ndtr((y[-1]-mu)/np.sqrt(noise_var))))
    summary=dict(mass_nodes=len(massgrid),beta_nodes=len(betagrid),mean_kms=mean,sd_kms=float(np.sqrt(variance)),observed_cdf=cdf,beta_edge_probability=float(prob[(betagrid<betas[0]+.05)|(betagrid>betas[-1]-.05)].sum()),mass_edge_probability=float(prob[:,(massgrid<lm[0]+.1)|(massgrid>lm[-1]-.1)].sum()))
    if step==1:
     keep=prob>1e-14;dropped=float(prob[~keep].sum());assert dropped<1e-6
     w=prob[keep];w=w/w.sum();u=mu[keep];noise=np.sqrt(noise_var)
     def quantile(level):return float(brentq(lambda v:float(np.sum(w*ndtr((v-u)/noise)))-level,float(u.min()-10*noise),float(u.max()+10*noise)))
     summary.update(interval95_kms=[quantile(.025),quantile(.975)],median_kms=quantile(.5),discarded_quantile_weight=dropped)
    summaries.append(summary)
   coarse,fine=summaries;g=protocol['numerical_gates'];passed=abs(coarse['observed_cdf']-fine['observed_cdf'])<=g['maximum_observed_cdf_change'] and abs(coarse['mean_kms']-fine['mean_kms'])<=g['maximum_mean_change_kms'] and abs(coarse['sd_kms']-fine['sd_kms'])<=g['maximum_sd_change_kms']
   out.append(dict(model=label,prior=prior,observed_outer_kms=float(y[-1]),conditional_measurement_sd_kms=float(np.sqrt(noise_var)),grids=summaries,numerical_pass=bool(passed)))
 results.append(dict(Name=name,predictions=out))
 print(json.dumps(dict(Name=name,predictions=[dict(model=r['model'],prior=r['prior'],mean=r['grids'][-1]['mean_kms'],interval=r['grids'][-1]['interval95_kms'],cdf=r['grids'][-1]['observed_cdf'],numerical_pass=r['numerical_pass']) for r in out])),flush=True)
(HERE/'results.json').write_text(json.dumps(dict(scope=protocol['scope'],systems=results,hashes={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths+[HERE/'run.py',HERE.parent/'slacs-component-refit/model.py',HERE.parent/'slacs-resolved-fit/model.py']}),indent=2)+'\n',encoding='utf-8',newline='\n')
