"""Independent dense and scalar-integral checks of the joint likelihood."""
from pathlib import Path
import json,sys
import numpy as np
from scipy.integrate import quad
from scipy.stats import multivariate_normal,norm
from likelihood import log_likelihood
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'joint-light-forward'))
from transport import predict,load_filters,MPC_CM,AB0_FNU,C_ANGSTROM_S
with np.load(HERE.parent/'des-calibration-covariance/des5yr-zero-point-covariance.npz',allow_pickle=False) as a:v=a['cov'];labels=a['labels']
assert labels.tolist()==['DES5YR '+b for b in 'griz']
filters=load_filters();rng=np.random.default_rng(881);errors=[]
for n in [4,12,40]:
    index=np.arange(n)%4;t=np.linspace(-40,60,n);bands=np.array(list('griz'))[index]
    lnu=4*np.pi*(1e-5*MPC_CM)**2*AB0_FNU*10**(.4*19)
    sed=lambda time,lam:lnu*C_ANGSTROM_S/lam**2*np.exp(-np.log(2)*(2*time/30)**2)
    predicted=predict(t,bands,100,sed,filters=filters)['native_fluxcal']
    measurement=np.maximum(10,predicted/20)
    y=predicted+rng.normal(0,measurement)
    for scale in [0.,1.,4.]:
        cov=v*scale;jac=np.zeros((n,4));jac[np.arange(n),index]=-.4*np.log(10)*predicted
        dense=np.diag(measurement**2)+jac@cov@jac.T
        expected=float(multivariate_normal.logpdf(y,mean=predicted,cov=dense))
        calculated=log_likelihood(y,measurement,predicted,index,cov)
        delta=abs(expected-calculated);assert delta<1e-9
        errors.append({'n':n,'covariance_scale':scale,'absolute_log_likelihood_error':delta})
# Independently integrate a single shared Gaussian offset rather than forming covariance.
m=np.array([10.,20.,30.]);e=np.array([.2,.3,.4]);y=np.array([10.1,19.9,30.2]);sd=.006
j=-.4*np.log(10)*m
integral=quad(lambda delta:np.exp(np.sum(norm.logpdf(y,loc=m+j*delta,scale=e))+norm.logpdf(delta,scale=sd)),-.08,.08,epsabs=1e-11,epsrel=1e-11)[0]
scalar=log_likelihood(y,e,m,np.zeros(3,dtype=int),np.array([[sd**2]]))
assert abs(np.log(integral)-scalar)<1e-9
# Repeated measurements share a calibration floor; they are not separate nuisances.
sigma_fractional=.4*np.log(10)*np.sqrt(v[0,0]);repeated=[]
for n in [1,10,100]:
    repeated.append({'measurements':n,'correct_mean_sigma_fractional':float(np.sqrt(.02**2/n+sigma_fractional**2)),'incorrect_independent_calibration_sigma':float(np.sqrt((.02**2+sigma_fractional**2)/n))})
result={'scope':'Normalized linearized Gaussian joint-flux likelihood checks, artificial observations and source only.','dense_comparisons':errors,'scalar_offset_integral_log_error':float(abs(np.log(integral)-scalar)),'repeated_measurement_example':repeated,'observations_fitted':False,'nonlinear_magnitude_marginalization_verified':False}
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'dense_max_error':max(v['absolute_log_likelihood_error'] for v in errors),'scalar_integral_error':result['scalar_offset_integral_log_error'],'repeated_example':repeated}))
