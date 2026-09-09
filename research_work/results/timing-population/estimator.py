"""Integrated event-shape and population likelihood, development only."""
import numpy as np
from scipy.special import log_ndtr, logsumexp, ndtr
from scipy.stats import qmc
from scipy.optimize import minimize


def batched_log_flux(flux, error, shapes):
    y=np.asarray(flux,dtype=float);err=np.asarray(error,dtype=float);s=np.asarray(shapes,dtype=float)
    if y.ndim!=1 or err.shape!=y.shape or s.ndim<1 or s.shape[-1]!=len(y):
        raise ValueError('Invalid observation/shape dimensions')
    if not all(np.all(np.isfinite(v)) for v in (y,err,s)) or np.any(err<=0):
        raise ValueError('Invalid observations')
    if len(y)==0:return np.zeros(s.shape[:-1])
    w=1/err**2
    sw=np.einsum('...i,i->...',s,w)
    sy=np.einsum('...i,i->...',s,w*y)
    ss=np.einsum('...i,...i,i->...',s,s,w)
    aa=4+ss;bb=sw;cc=100+w.sum();det=aa*cc-bb**2
    r0=sy-ss;r1=np.dot(w,y)-sw
    quadratic=np.dot(w,y*y)-2*sy+ss-(cc*r0*r0-2*bb*r0*r1+aa*r1*r1)/det
    logdet=np.log(err**2).sum()+np.log(.25*.01)+np.log(det)
    mean_a=1+(cc*r0-bb*r1)/det;var_a=cc/det
    return -.5*(len(y)*np.log(2*np.pi)+logdet+quadratic)+log_ndtr(mean_a/np.sqrt(var_a))-log_ndtr(2.)


def event_log_likelihood(time,flux,error,log_width,settings,protocol):
    unit=qmc.Sobol(4,scramble=True,seed=settings['seed']).random_base2(settings['sobol_power'])
    priors=protocol['shape_priors']
    ranges=np.array([priors[k] for k in ['peak_days_uniform','rise_fraction_uniform','rise_power_uniform','fall_power_uniform']])
    pars=ranges[:,0]+unit*(ranges[:,1]-ranges[:,0])
    peak,rise,pr,pf=pars.T
    u=np.asarray(time)[None,:]-peak[:,None]
    frac=np.where(u<0,rise[:,None],1-rise[:,None])
    power=np.where(u<0,pr[:,None],pf[:,None])
    base=np.abs(u)/frac
    logs=[]
    for start in range(0,len(log_width),8):
        width=np.exp(log_width[start:start+8])[:,None,None]
        shapes=np.exp(-np.log(2.)*(base[None,:,:]/width)**power[None,:,:])
        logs.extend(logsumexp(batched_log_flux(flux,error,shapes),axis=1)-np.log(len(unit)))
    return np.asarray(logs)


def log_width_weights(log_width,redshift,parameters):
    a,b,log_sigma=parameters;sigma=np.exp(log_sigma)
    mu=a+b*np.log1p(redshift)
    x=np.asarray(log_width)
    weights=np.empty_like(x);weights[0]=(x[1]-x[0])/2;weights[-1]=(x[-1]-x[-2])/2
    weights[1:-1]=(x[2:]-x[:-2])/2
    logs=-.5*((x[None,:]-mu[:,None])/sigma)**2+np.log(weights)[None,:]
    # Normalized finite-grid population model. Boundary mass is audited separately.
    return logs-logsumexp(logs,axis=1)[:,None]


def population_log_likelihood(log_width,event_logs,redshift,parameters):
    return float(logsumexp(event_logs+log_width_weights(log_width,redshift,parameters),axis=1).sum())


def fit_population(log_width,event_logs,redshift,protocol):
    bounds=protocol['bounds']
    limits=[np.log(bounds['mean_width_at_zero_redshift_days']),bounds['b'],np.log(bounds['intrinsic_sigma'])]
    centered=event_logs-event_logs.max(axis=1)[:,None]
    objective=lambda p:-population_log_likelihood(log_width,centered,redshift,p)
    fits=[minimize(objective,[np.log(30),b,np.log(sig)],method='L-BFGS-B',bounds=limits,
                   options={'maxiter':500,'ftol':1e-11,'gtol':1e-6}) for b,sig in [(0,.1),(1,.1),(.5,.3)]]
    valid=[f for f in fits if f.success]
    fit=min(valid or fits,key=lambda f:f.fun)
    a,b,ls=fit.x;sigma=np.exp(ls);mu=a+b*np.log1p(redshift)
    outside=ndtr((log_width[0]-mu)/sigma)+ndtr((mu-log_width[-1])/sigma)
    interior=all(lo+1e-5<value<hi-1e-5 for value,(lo,hi) in zip(fit.x,limits))
    return {'a':float(a),'b':float(b),'sigma':float(sigma),'negative_centered_log_likelihood':float(fit.fun),
            'optimizer_success':bool(fit.success),'interior_solution':bool(interior),
            'maximum_population_mass_outside_width_bounds':float(outside.max()),
            'starts':[{'success':bool(f.success),'objective':float(f.fun),'parameters':f.x.tolist()} for f in fits]}
