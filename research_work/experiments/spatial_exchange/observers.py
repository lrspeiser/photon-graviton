"""Explicit action-clock convention; not a derived material-clock model."""
import numpy as np
from scipy.optimize import root
from probes import evaluate


def clock_rate(q,p,mass,radius,length,arrays):
    alpha,z,beta=arrays;n=alpha.shape[0];line=np.arange(n)*length/n-length/2
    xyz=np.stack(np.meshgrid(line,line,line,indexing='ij'),axis=-1)
    w=np.maximum(1-np.sum((xyz-q)**2,axis=-1)/radius**2,0)**3
    if w.sum()<=0:raise ValueError('Unresolved clock kernel')
    w/=w.sum();E=np.sqrt(mass*mass+z*np.dot(p,p))
    return float(np.sum(w*alpha*mass/E))


def observer(q,velocity,radius,length,arrays):
    velocity=np.asarray(velocity)
    solution=root(lambda p:evaluate(q,p,1.,radius,length,*arrays)['velocity']-velocity,np.zeros(3),tol=1e-11)
    value=evaluate(q,solution.x,1.,radius,length,*arrays)
    residual=float(np.linalg.norm(value['velocity']-velocity))
    if residual>=1e-10:raise RuntimeError('Observer velocity inversion failed')
    rate=clock_rate(q,solution.x,1.,radius,length,arrays)
    if rate<=0:raise ValueError('Nonpositive candidate clock rate')
    return dict(momentum=solution.x,velocity=value['velocity'],rate=rate,energy=value['energy'],residual=residual,solver_success=bool(solution.success))


def photon_energy(q,k,obs,radius,length,arrays):
    H=evaluate(q,k,0.,radius,length,*arrays)['energy']
    return float((H-np.dot(k,obs['velocity']))/obs['rate'])
