"""Normalized Gaussian likelihood with shared linearized zero-point errors."""
import numpy as np

def log_likelihood(observed,measurement_sigma,predicted,band_index,offset_covariance):
    y=np.asarray(observed,dtype=float);e=np.asarray(measurement_sigma,dtype=float);m=np.asarray(predicted,dtype=float);index=np.asarray(band_index)
    v=np.asarray(offset_covariance,dtype=float)
    if y.ndim!=1 or e.shape!=y.shape or m.shape!=y.shape or index.shape!=y.shape or index.dtype.kind not in 'iu':raise ValueError('Invalid measurement arrays')
    if not all(np.isfinite(a).all() for a in [y,e,m,v]) or np.any(e<=0):raise ValueError('Nonfinite values or nonpositive errors')
    if v.ndim!=2 or v.shape[0]!=v.shape[1] or not np.allclose(v,v.T,atol=1e-14,rtol=0):raise ValueError('Invalid offset covariance')
    if np.any(index<0) or np.any(index>=len(v)):raise ValueError('Invalid band index')
    ev,vec=np.linalg.eigh(v)
    if ev.min() < -1e-14:raise ValueError('Offset covariance is not positive semidefinite')
    # Native FLUXCAL decreases with increasing synthetic-magnitude offset.
    jac=np.zeros((len(y),len(v)));jac[np.arange(len(y)),index]=-.4*np.log(10)*m
    factor=jac@(vec*np.sqrt(np.maximum(ev,0)))
    weighted=factor/e[:,None];r=(y-m)/e
    chol=np.linalg.cholesky(np.eye(len(v))+weighted.T@weighted)
    projected=np.linalg.solve(chol,weighted.T@r)
    quadratic=float(r@r-projected@projected)
    logdet=float(2*np.log(e).sum()+2*np.log(np.diag(chol)).sum())
    return -.5*(len(y)*np.log(2*np.pi)+logdet+quadratic)
