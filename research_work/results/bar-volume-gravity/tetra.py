"""Unsoftened Newtonian gravity of constant-density tetrahedra, G=1."""
import numpy as np

def gravity(vertices,masses,points):
    v=np.asarray(vertices,dtype=float);mass=np.broadcast_to(masses,len(v)).astype(float)
    volume=np.abs(np.einsum('ij,ij->i',v[:,1]-v[:,0],np.cross(v[:,2]-v[:,0],v[:,3]-v[:,0])))/6
    scale=np.max(np.linalg.norm(v-v[:,:1],axis=2),axis=1)
    if np.any(volume<=1e-15*scale**3):raise ValueError('Degenerate cell; no mass may be discarded')
    rho=mass/volume;P=[];A=[]
    for point in np.atleast_2d(points):
        p=np.zeros(len(v));a=np.zeros((len(v),3))
        for opposite in range(4):
            ids=[j for j in range(4) if j!=opposite];face=v[:,ids].copy()
            n=np.cross(face[:,1]-face[:,0],face[:,2]-face[:,0])
            flip=np.einsum('ij,ij->i',n,v[:,opposite]-face[:,0])>0
            face[flip]=face[flip][:,[0,2,1]]
            n=np.cross(face[:,1]-face[:,0],face[:,2]-face[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
            r=face-point;length=np.linalg.norm(r,axis=2)
            d=np.einsum('ij,ij->i',n,r[:,0])
            numerator=np.einsum('ij,ij->i',r[:,0],np.cross(r[:,1],r[:,2]))
            denominator=np.prod(length,axis=1)
            for i,j,k in ((0,1,2),(1,2,0),(2,0,1)):
                denominator+=np.einsum('ij,ij->i',r[:,i],r[:,j])*length[:,k]
            omega=2*np.arctan2(numerator,denominator)
            integral=-d*omega
            for i,j in ((0,1),(1,2),(2,0)):
                edge=face[:,j]-face[:,i];ell=np.linalg.norm(edge,axis=1)
                outward=np.cross(edge/ell[:,None],n)
                h=np.einsum('ij,ij->i',outward,r[:,i])
                den=length[:,i]+length[:,j]-ell
                on_edge=(abs(h)<1e-14*scale)&(abs(d)<1e-14*scale)&(den<=1e-14*scale)
                if np.any((den<=0)&~on_edge):raise FloatingPointError('Unresolved edge distance')
                term=np.zeros(len(v));use=~on_edge
                term[use]=h[use]*np.log1p(2*ell[use]/den[use])
                integral+=term
            p-=.5*rho*d*integral
            a-=rho[:,None]*n*integral[:,None]
        P.append(p.sum());A.append(a.sum(axis=0))
    return np.array(P),np.array(A)
