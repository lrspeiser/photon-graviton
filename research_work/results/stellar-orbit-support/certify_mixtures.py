"""Use the convex weight problem directly and report an objective-gap bound."""
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from run import HERE,OUT,save

original=json.loads((HERE/'results.json').read_text())
refined=json.loads((HERE/'mixture-refinement.json').read_text())['fits']
data=pd.read_parquet(OUT/'training-J-diagnostic.parquet')
launch=pd.read_parquet(OUT/'expanded-launches.parquet')
rows=[]
for model,item in original['models'].items():
    y=data[model+'_J_kms2'].to_numpy();s=data[model+'_J_sigma_kms2'].to_numpy()
    centres=data.set_index('source_id').loc[launch.source_id,model+'_J_kms2'].to_numpy()
    for fit in item['expanded_mixture_fits']:
        broadening=fit['fixed_J_broadening_kms2'];variance=s*s+broadening*broadening
        logK=-.5*((y[:,None]-centres[None,:])**2/variance[:,None]+np.log(2*np.pi*variance[:,None]))
        offsets=logK.max(axis=1);K=np.exp(logK-offsets[:,None])
        def objective(w):
            mix=K@w
            return -np.mean(np.log(mix)),-np.mean(K/mix[:,None],axis=0)
        earlier=next((r for r in refined if r['model']==model and r['fixed_J_broadening_kms2']==broadening),None)
        start=np.array(earlier['weights'] if earlier else fit['weights'])
        epsilon=1e-10;k=len(start)
        start=epsilon+(1-k*epsilon)*start/start.sum()
        opt=minimize(objective,start,method='SLSQP',jac=True,bounds=[(epsilon,1)]*k,
            constraints=[dict(type='eq',fun=lambda w:w.sum()-1,jac=lambda w:np.ones(k))],
            options=dict(maxiter=500,ftol=1e-12))
        # Keep exactly feasible weights for the convex Frank-Wolfe gap.
        u=np.maximum(opt.x-epsilon,0);u/=u.sum();w=epsilon+(1-k*epsilon)*u
        f,g=objective(w)
        gap=float((1-k*epsilon)*(np.dot(u,g)-g.min()))
        assert gap>=-1e-10
        nll=float(f-np.mean(offsets))
        row=dict(model=model,fixed_J_broadening_kms2=broadening,weights=w.tolist(),
            scipy_success=bool(opt.success),termination=str(opt.message),iterations=int(opt.nit),
            mean_negative_log_density=nll,minimum_weight=epsilon,
            convex_objective_gap_upper_bound_nats_per_star=max(gap,0),
            gap_below_one_millionth=bool(gap<1e-6),
            interpretation='Optimization accuracy for this training marginal mixture, not evidence for a gravity theory')
        rows.append(row);save('mixture-certification.json',dict(fits=rows))
        print(model,broadening,opt.success,opt.nit,gap,flush=True)
