"""Finish slow EM weight fits without changing the data or basis centres."""
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import softmax
from run import HERE,OUT,save

original=json.loads((HERE/'results.json').read_text())
data=pd.read_parquet(OUT/'training-J-diagnostic.parquet')
launch=pd.read_parquet(OUT/'expanded-launches.parquet')
rows=[]
for model,item in original['models'].items():
    y=data[model+'_J_kms2'].to_numpy();s=data[model+'_J_sigma_kms2'].to_numpy()
    lookup=data.set_index('source_id')[model+'_J_kms2']
    centres=lookup.loc[launch.source_id].to_numpy()
    for fit in item['expanded_mixture_fits']:
        if fit['converged']:continue
        broadening=fit['fixed_J_broadening_kms2'];variance=s*s+broadening*broadening
        logK=-.5*((y[:,None]-centres[None,:])**2/variance[:,None]+np.log(2*np.pi*variance[:,None]))
        offsets=logK.max(axis=1);K=np.exp(logK-offsets[:,None])
        def objective(theta):
            w=softmax(np.r_[theta,0.]);mix=K@w
            value=-np.mean(np.log(mix))
            gw=-np.mean(K/mix[:,None],axis=0)
            grad=(w*(gw-np.dot(w,gw)))[:-1]
            return value,grad
        w=np.maximum(fit['weights'],1e-16);theta=np.log(w[:-1]/w[-1]);theta=np.clip(theta,-40,40)
        opt=minimize(objective,theta,method='L-BFGS-B',jac=True,bounds=[(-40,40)]*len(theta),
                     options=dict(maxiter=1500,ftol=1e-13,gtol=1e-8,maxls=40))
        weights=softmax(np.r_[opt.x,0.]);nll=float(opt.fun-np.mean(offsets))
        assert nll<=fit['training_mean_negative_log_density']+1e-8
        row=dict(model=model,fixed_J_broadening_kms2=broadening,success=bool(opt.success),
                 termination=str(opt.message),iterations=int(opt.nit),weights=weights.tolist(),
                 original_EM_mean_negative_log_density=fit['training_mean_negative_log_density'],
                 refined_mean_negative_log_density=nll,
                 maximum_absolute_theta_gradient=float(np.max(abs(opt.jac))))
        rows.append(row);save('mixture-refinement.json',dict(fits=rows,
             interpretation='Convergence refinement of training-only marginal library fits, not a theory comparison'))
        print(model,broadening,opt.success,opt.nit,flush=True)
