"""Cache quadrature nodes for the unchanged synthetic deposit geometries."""
import numpy as np
from numpy.polynomial.legendre import leggauss

class CachedCompanion:
    def __init__(self,geometry,amplitude=1000.,soft=.3):
        self.amplitude=amplitude;self.soft=soft
        p=np.arange(128)*2*np.pi/128
        if geometry=='equatorial':
            self.centres=np.c_[1.5*np.cos(p),1.5*np.sin(p),np.zeros(len(p))]
            self.w=np.ones(len(p))/len(p)
        elif geometry=='caps':
            self.centres=np.vstack([np.c_[.5*np.cos(p),.5*np.sin(p),np.full(len(p),z)] for z in [-1.5,1.5]])
            self.w=np.ones(len(self.centres))/len(self.centres)
        elif geometry=='shell':
            mu,ww=leggauss(32)
            self.centres=np.concatenate([np.c_[1.5*np.sqrt(1-u*u)*np.cos(p),1.5*np.sqrt(1-u*u)*np.sin(p),np.full(len(p),1.5*u)] for u in mu])
            self.w=np.repeat(ww/2/len(p),len(p))
        else:raise ValueError(geometry)

    def evaluate(self,xyz):
        xyz=np.atleast_2d(xyz)
        delta=xyz[:,None,:]-self.centres[None,:,:]
        q=np.sum(delta*delta,axis=2)+self.soft**2
        return (-self.amplitude*np.sum(self.w[None,:]/np.sqrt(q),axis=1),
                -self.amplitude*np.sum(self.w[None,:,None]*delta/q[:,:,None]**1.5,axis=1))
