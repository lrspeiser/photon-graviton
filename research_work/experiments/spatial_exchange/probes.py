"""Regularized test-body Hamiltonian shared by massive and massless probes."""
import numpy as np


def coefficients(fields):
    u=.08*fields[0];alpha=np.exp(u);z=np.exp(2*u);A=fields[1:4]
    beta=.04*z*A/np.sqrt(1+.08**2*np.sum(A*A,axis=0))
    return alpha,z,beta


def evaluate(q,p,mass,radius,length,alpha,z,beta):
    n=alpha.shape[0];h=length/n
    if mass<0 or (mass==0 and np.dot(p,p)==0):raise ValueError('Invalid probe mass/momentum')
    if np.any(abs(q)+radius>=length/2):raise ValueError('Probe kernel reaches box edge')
    low=np.maximum(0,np.floor((q-radius+length/2)/h).astype(int))
    high=np.minimum(n,np.ceil((q+radius+length/2)/h).astype(int)+1)
    sl=tuple(slice(int(a),int(b)) for a,b in zip(low,high))
    xyz=np.stack(np.meshgrid(*[np.arange(a,b)*h-length/2 for a,b in zip(low,high)],indexing='ij'),axis=-1)
    delta=xyz-q;t=np.maximum(1-np.sum(delta*delta,axis=-1)/radius**2,0)
    raw=t**3;norm=raw.sum()
    if norm<=0:raise ValueError('Unresolved probe kernel')
    draw=6*t[...,None]**2*delta/radius**2
    w=raw/norm;dw=draw/norm-raw[...,None]*draw.sum(axis=(0,1,2))/norm**2
    al=alpha[sl];zz=z[sl];bb=beta[(slice(None),)+sl]
    E=np.sqrt(mass*mass+zz*np.dot(p,p));density=al*E+np.einsum('jxyz,j->xyz',bb,p)
    bbar=np.sum(w*bb,axis=(1,2,3));cbar=float(np.sum(w*zz))
    velocity=p*np.sum(w*al*zz/E)+bbar
    force=-np.sum(dw*density[...,None],axis=(0,1,2))
    return dict(energy=float(np.sum(w*density)),velocity=velocity,force=force,bbar=bbar,cbar=cbar)


def derivative(state,masses,radius,length,coefficient_arrays):
    result=np.empty_like(state)
    for i,row in enumerate(state):
        value=evaluate(row[:3],row[3:],masses[i],radius,length,*coefficient_arrays)
        result[i,:3]=value['velocity'];result[i,3:]=value['force']
    return result
