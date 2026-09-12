"""Joint recovery using independent analytic artificial flux injection."""
from pathlib import Path
import json,sys
import numpy as np
from scipy.optimize import minimize
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'joint-light-forward'))
sys.path.insert(0,str(HERE.parent/'joint-light-likelihood'))
from transport import predict,load_filters,MPC_CM,AB0_FNU,C_ANGSTROM_S
from likelihood import log_likelihood
p=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'));filters=load_filters()
with np.load(HERE.parent/'des-calibration-covariance/des5yr-zero-point-covariance.npz',allow_pickle=False) as a:cov=a['cov']
distances=np.array(p['distances_mpc'],dtype=float);times=np.repeat(p['observer_days'],4).astype(float);bands=np.tile(np.array(p['bands']),len(p['observer_days']));indices=np.tile(np.arange(4),len(times)//4*len(distances))
alpha=p['alpha_per_mpc'];W=p['source_fwhm_days'];scaleL=10**(-.4*p['source_absolute_AB_mag'])
lnu=4*np.pi*(1e-5*MPC_CM)**2*AB0_FNU*scaleL
S=np.exp(alpha*distances);offset=np.array([filters[b][3] for b in bands]);results=[]
for truth in p['truth_b']:
    rng=np.random.default_rng(p['seed']);delta=np.linalg.cholesky(cov)@rng.normal(size=4);measurements=[];errors=[];analytic_means=[]
    for D,s in zip(distances,S):
        peak=1e11*scaleL*(1e-5/D)**2*s**(-truth)*10**(-.4*offset)
        mean=peak*np.exp(-np.log(2)*(2*times/(W*s**truth))**2)
        e=peak/p['peak_flux_snr'];shift=-.4*np.log(10)*mean*delta[np.tile(np.arange(4),len(times)//4)]
        measurements.extend(mean+shift+rng.normal(0,e));errors.extend(e);analytic_means.extend(mean)
    y=np.array(measurements);e=np.array(errors);zobs=S-1+rng.normal(0,p['redshift_sigma'],len(S))
    def forward(q):
        rate=alpha*np.exp(q[0]);b=q[1];amplitude=np.exp(q[2]);width=W*np.exp(q[3])
        sed=lambda time,lam:np.where((lam>=1000)&(lam<=20000),amplitude*lnu*C_ANGSTROM_S/lam**2*np.exp(-np.log(2)*(2*time/width)**2),0.)
        return np.concatenate([predict(times,bands,D,sed,alpha_per_mpc=rate,stretch_exponent=b,filters=filters)['native_fluxcal'] for D in distances])
    truth_error=float(np.max(abs(forward([0,truth,0,0])-analytic_means)/np.maximum(np.array(analytic_means),1e-30)))
    assert truth_error<1e-10
    def objective(q):
        model=forward(q);zp=np.expm1(alpha*np.exp(q[0])*distances)
        return -log_likelihood(y,e,model,indices,cov)+.5*np.sum(((zobs-zp)/p['redshift_sigma'])**2+np.log(2*np.pi*p['redshift_sigma']**2))
    bounds=[np.log([.5,1.5]),[-.5,1.5],np.log([.5,2]),np.log([.5,2])]
    trials=[minimize(objective,[np.log(ar),b,0,0],method='L-BFGS-B',bounds=bounds,options={'maxiter':400,'ftol':1e-11,'gtol':1e-5}) for ar,b in [(.9,0),(1.,.5),(1.1,1)]]
    valid=[v for v in trials if v.success]
    if not valid:raise RuntimeError('All recovery fits failed')
    fit=min(valid,key=lambda v:v.fun);q=fit.x
    errors_out={'relative_alpha_error':float(abs(np.expm1(q[0]))),'absolute_b_error':float(abs(q[1]-truth)),'relative_luminosity_error':float(abs(np.expm1(q[2]))),'relative_width_error':float(abs(np.expm1(q[3])))}
    interior=all(lo+1e-5<v<hi-1e-5 for v,(lo,hi) in zip(q,bounds))
    result={'injected_b':truth,'fitted_alpha_per_mpc':float(alpha*np.exp(q[0])),'fitted_b':float(q[1]),'fitted_common_luminosity_scale':float(np.exp(q[2])),'fitted_common_width_days':float(W*np.exp(q[3])),'errors':errors_out,'screen_pass':bool(interior and all(errors_out[k]<=v for k,v in p['screening_limits'].items())),'analytic_vs_forward_max_relative_error':truth_error,'starts':[{'success':bool(v.success),'nll':float(v.fun),'parameters':v.x.tolist()} for v in trials]}
    results.append(result);print(json.dumps(result),flush=True)
(HERE/'results.json').write_text(json.dumps({'scope':p['scope'],'results':results,'observed_data_fitted':False,'source_evolution_allowed':False,'uncertainty_coverage_established':False},indent=2)+'\n',encoding='utf-8',newline='\n')
