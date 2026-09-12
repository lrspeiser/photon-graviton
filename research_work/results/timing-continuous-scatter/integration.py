"""Exact normal averaging of a piecewise-linear likelihood, including sigma=0.
New development component; does not replace the frozen calibration estimator.
"""
import numpy as np
from scipy.special import ndtr

def averaged_likelihood(x, likelihood, mu, sigma):
    x=np.asarray(x,float);L=np.asarray(likelihood,float);mu=np.asarray(mu,float)
    if x.ndim!=1 or len(x)<2 or np.any(np.diff(x)<=0): raise ValueError('Invalid grid')
    if L.shape!=(len(mu),len(x)) or sigma<0 or not np.isfinite(sigma): raise ValueError('Invalid inputs')
    if not all(np.all(np.isfinite(v)) for v in [x,L,mu]) or np.any(L<0): raise ValueError('Nonfinite or negative inputs')
    if sigma==0:
        return np.array([np.interp(m,x,row) if x[0]<=m<=x[-1] else 0. for m,row in zip(mu,L)])
    a=(x[:-1][None,:]-mu[:,None])/sigma;b=(x[1:][None,:]-mu[:,None])/sigma
    P=np.where(a>0,ndtr(-a)-ndtr(-b),ndtr(b)-ndtr(a))
    phi_a=np.exp(-a*a/2)/np.sqrt(2*np.pi);phi_b=np.exp(-b*b/2)/np.sqrt(2*np.pi)
    dx=np.diff(x)[None,:]
    left=((x[1:][None,:]-mu[:,None])*P-sigma*(phi_a-phi_b))/dx
    right=((mu[:,None]-x[:-1][None,:])*P+sigma*(phi_a-phi_b))/dx
    if min(left.min(),right.min()) < -1e-12: raise ArithmeticError('Negative integration mass')
    left=np.maximum(left,0);right=np.maximum(right,0)
    mass=P.sum(axis=1)
    if np.any(mass<=0): raise ArithmeticError('Population outside representable support')
    return (left*L[:,:-1]+right*L[:,1:]).sum(axis=1)/mass

def log_likelihood(x,event_logs,z,a,b,sigma):
    offset=event_logs.max(axis=1)
    value=averaged_likelihood(x,np.exp(event_logs-offset[:,None]),a+b*np.log1p(z),sigma)
    if np.any(value<=0): return -np.inf
    return float(np.sum(np.log(value)+offset))
