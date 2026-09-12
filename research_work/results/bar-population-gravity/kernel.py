"""Known Newtonian Green functions integrated over source ages and symmetries."""
import numpy as np
from scipy.special import roots_legendre

POINTS=np.array([[.1,0,0],[1,0,0],[1,0,1],[0,1,1],[3,0,0]])
SYMMETRIES=np.array([[1,1,1],[-1,-1,1],[1,1,-1],[-1,-1,-1]])

def integrate_kernel(position,steps,T,points=POINTS,nodes=8):
    boundaries=np.r_[0,np.asarray(steps)[(np.asarray(steps)>0)&(np.asarray(steps)<T)],T]
    z,w=roots_legendre(nodes)
    half=np.diff(boundaries)/2
    ages=(boundaries[:-1,None]+half[:,None]*(1+z)).ravel()
    weights=(half[:,None]*w/T).ravel()
    x=np.asarray(position(ages))
    assert x.shape==(len(ages),3)
    potential=np.zeros(len(points));force=np.zeros((len(points),3))
    for symmetry in SYMMETRIES:
        dx=points[:,None,:]-x[None,:,:]*symmetry
        r=np.linalg.norm(dx,axis=2)
        if np.any(r==0):raise ValueError('Point lies exactly on a quadrature trajectory')
        potential-=np.sum(weights[None,:]/r,axis=1)/4
        force-=np.sum(weights[None,:,None]*dx/r[:,:,None]**3,axis=1)/4
    return dict(potential=potential.tolist(),acceleration=force.tolist(),nodes=nodes,
                age_weight_sum=float(sum(weights)))

def compare(a,b):
    dp=np.abs(np.array(a['potential'])-b['potential'])/np.maximum(np.abs(b['potential']),.01)
    da=np.linalg.norm(np.array(a['acceleration'])-b['acceleration'],axis=1)/np.maximum(np.linalg.norm(b['acceleration'],axis=1),.01)
    return dict(potential_scaled_changes=dp.tolist(),acceleration_scaled_changes=da.tolist(),
                passes=bool(max(dp.max(),da.max())<1e-4))
