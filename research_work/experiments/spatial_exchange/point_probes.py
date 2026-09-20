"""Positive cubic interpolation of the candidate point-particle Hamiltonian."""
import numpy as np


def kernel(q,n,length):
    h=length/n;coordinate=(q+length/2)/h;i=np.floor(coordinate).astype(int);t=coordinate-i
    nodes=[np.arange(j-1,j+3) for j in i]
    if any(j[0]<0 or j[-1]>=n for j in nodes):raise ValueError('Point interpolation reaches box edge')
    w=np.array([(1-t)**3,3*t**3-6*t*t+4,-3*t**3+3*t*t+3*t+1,t**3]).T/6
    dw=np.array([-3*(1-t)**2,9*t*t-12*t,-9*t*t+6*t+3,3*t*t]).T/(6*h)
    weights=np.einsum('i,j,k->ijk',*w)
    derivatives=np.array([np.einsum('i,j,k->ijk',dw[0],w[1],w[2]),np.einsum('i,j,k->ijk',w[0],dw[1],w[2]),np.einsum('i,j,k->ijk',w[0],w[1],dw[2])])
    return np.ix_(*nodes),weights,derivatives


def evaluate(q,p,mass,length,alpha,z,beta):
    if mass<0 or (mass==0 and np.dot(p,p)==0):raise ValueError('Invalid mass/momentum')
    ix,w,dw=kernel(q,alpha.shape[0],length);al=alpha[ix];zz=z[ix];bb=np.array([b[ix] for b in beta])
    E=np.sqrt(mass*mass+zz*np.dot(p,p));density=al*E+np.einsum('jxyz,j->xyz',bb,p)
    bbar=np.sum(w*bb,axis=(1,2,3));cbar=float(np.sum(w*zz))
    return dict(energy=float(np.sum(w*density)),velocity=p*np.sum(w*al*zz/E)+bbar,
                force=-np.sum(dw*density,axis=(1,2,3)),bbar=bbar,cbar=cbar)


def derivative(state,masses,length,arrays):
    result=np.empty_like(state)
    for i,row in enumerate(state):
        value=evaluate(row[:3],row[3:],masses[i],length,*arrays)
        result[i,:3]=value['velocity'];result[i,3:]=value['force']
    return result
