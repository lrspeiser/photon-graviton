"""Separate ordinary components and re-solve the frozen nonlinear field."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear, minimize
from scipy.special import roots_legendre

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SOURCE=HERE.parent/'conservative-field-completion'
spec=importlib.util.spec_from_file_location('frozen',SOURCE/'run.py')
f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
GROUPS=['stellar_disks','gas_disks','central_stars']
PARAMS=[(1.332e9,2.,.3,2.7,'exp'),(8.97e8,2.8,.9,2.7,'exp'),(5.81e7,7.,.085,4.,'sech2'),(2.68e9,1.5,.045,12.,'sech2')]

class DiskPart(f.Expansion):
    def __init__(self,refined,indices):
        r=np.geomspace(1e-5,500,3072 if refined else 1536)
        ells=np.arange(0,257 if refined else 129,2)
        mu,w=roots_legendre(1024 if refined else 512)
        leg,_=f.basis(ells,mu);angular=leg*(w[:,None]*(2*ells+1)[None,:]/2)
        coeff=np.empty((len(r),len(ells)))
        for start in range(0,len(r),32):
            rr=r[start:start+32,None];R=rr*np.sqrt(1-mu**2);z=rr*mu;rho=np.zeros_like(R)
            for i in indices:
                sigma,rd,h,hole,kind=PARAMS[i];surface=sigma*np.exp(-hole/R-R/rd)
                if kind=='exp':rho+=surface*np.exp(-abs(z)/h)/(2*h)
                else:rho+=surface*np.exp(-2*np.logaddexp(z/(2*h),-z/(2*h))+2*np.log(2))/(4*h)
            coeff[start:start+len(rr)]=rho@angular
        inn,out=f.scaled_integrals(r,coeff*r[:,None],coeff*r[:,None],ells)
        fac=4*np.pi*f.G/(2*ells+1)
        super().__init__(r,ells,-fac*(inn+out),fac*((ells+1)*inn-ells*out)/r[:,None])
        self.mass=float(4*np.pi*inn[-1,0]*r[-1])

class Components:
    def __init__(self,refined):
        self.stellar=DiskPart(refined,[0,1]);self.gas=DiskPart(refined,[2,3])
        cache=ROOT/'research_work/data-cache/bar-field'
        self.bar=f.CachedAxisymmetric(cache/'bar-L64.npz');self.nuclei=f.CachedAxisymmetric(cache/'nuclei-L16.npz')
        self.cache={}
    def parts(self,r,mu):
        r,mu=np.broadcast_arrays(np.atleast_1d(r),np.atleast_1d(mu))
        key=hashlib.sha256(r.tobytes()+mu.tobytes()).digest()
        if key not in self.cache:
            central=np.array(self.bar.evaluate(r,mu))+np.array(self.nuclei.evaluate(r,mu))
            q=r*r+.001**2
            bh=np.array([-f.G*4.1e6/np.sqrt(q),f.G*4.1e6*r/q**1.5,np.zeros_like(r)])
            self.cache[key]=np.array([self.stellar.evaluate(r,mu),self.gas.evaluate(r,mu),central,bh])
        return self.cache[key]

class Model:
    def __init__(self,components,scales):self.components=components;self.scales=np.r_[scales,1.]
    def evaluate(self,r,mu):return np.einsum('i,ijk->jk',self.scales,self.components.parts(r,mu))

