"""Exact polar-axis limit for the even-m bar expansion, otherwise original evaluator."""
import numpy as np
def evaluate(bar,xyz):
    x=np.atleast_2d(xyz).astype(float);axis=np.hypot(x[:,0],x[:,1])==0
    potential=np.empty(len(x));acc=np.empty_like(x)
    if np.any(~axis):potential[~axis],acc[~axis]=bar.evaluate(x[~axis])
    if np.any(axis):
        z=x[axis,2];r=abs(z)
        if np.any((r<bar.r[0])|(r>bar.r[-1])):raise ValueError('Axis outside radial cache domain')
        assert np.all(bar.m%2==0),'Polar transverse force needs m=1 for a general field'
        zero=bar.m==0;l=bar.l[zero]
        Y=np.sqrt((2*l+1)/(4*np.pi))[None,:]*np.sign(z)[:,None]**l
        coeff=bar.spline(np.log(r))[:,zero];dr=bar.spline(np.log(r),1)[:,zero]/r[:,None]
        potential[axis]=np.sum(coeff*Y,axis=1)
        acc[axis]=np.c_[np.zeros(len(z)),np.zeros(len(z)),-np.sign(z)*np.sum(dr*Y,axis=1)]
    return potential,acc
