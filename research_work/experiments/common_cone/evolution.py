"""CC-2 exact Hamiltonian discretization with a positive dissipative ledger."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from model import spherical_weights, rotation


class Evolution:
    def __init__(self,source='rotating',probe='photon',eta=.08,n=32,length=12.,radius=.9,angle=0.):
        self.source,self.probe,self.eta=source,probe,eta
        self.n,self.length,self.radius=n,length,radius
        self.g,self.kappa,self.omega=.08,.5,.2
        self.dx=length/n;self.dv=self.dx**3;self.nf=4*n**3
        c=np.arange(n)*self.dx-length/2
        self.x=np.stack(np.meshgrid(c,c,c,indexing='ij'),axis=-1)
        self.gamma=2*np.maximum(np.max(abs(self.x),axis=-1)-(length/2-1),0)**2
        self.m=np.ones(6 if probe=='none' else 7)
        if probe!='none':self.m[-1]=0 if probe=='photon' else 1e-4
        self.count=len(self.m);self.rot=rotation(angle)

    def unpack(self,y):
        f=y[:self.nf].reshape((4,)+(self.n,)*3)
        pi=y[self.nf:2*self.nf].reshape(f.shape)
        q=y[2*self.nf:2*self.nf+3*self.count].reshape(-1,3)
        p=y[2*self.nf+3*self.count:-1].reshape(-1,3)
        return f,pi,q,p

    def initial(self):
        theta=np.arange(6)*np.pi/3
        q=.7*np.column_stack((np.cos(theta),np.sin(theta),np.zeros(6)))
        direction={'rest':0,'rotating':1,'reverse':-1}[self.source]
        p=.2*direction*np.column_stack((-np.sin(theta),np.cos(theta),np.zeros(6)))
        if self.probe!='none':
            photon=self.probe=='photon'
            q=np.vstack((q,[-1.7 if photon else -1,.6,0]))
            p=np.vstack((p,[1e-4 if photon else 1e-4*.4/np.sqrt(.84),0,0]))
        return np.concatenate((np.zeros(2*self.nf),(q@self.rot.T).ravel(),(p@self.rot.T).ravel(),[0.]))

    def plus(self,f,j):return (np.roll(f,-1,axis=j)-f)/self.dx
    def minus(self,f,j):return (f-np.roll(f,1,axis=j))/self.dx
    def central(self,f,j):return (np.roll(f,-1,axis=j)-np.roll(f,1,axis=j))/(2*self.dx)

    def ingredients(self,y):
        f,pi,q,p=self.unpack(y)
        alpha=np.exp(self.g*f[0]);a=f[1:]
        scale=np.sqrt(1+self.eta**2*np.sum(a*a,axis=0))
        b=self.kappa*self.eta*a/scale;beta=alpha*b
        gp=np.stack([self.plus(f,j+1) for j in range(3)])
        gm=np.stack([self.minus(f,j+1) for j in range(3)])
        gc=(gp+gm)/2
        kinetic=np.sum(pi*pi,axis=0)/2+np.sum(gp*gp+gm*gm,axis=(0,1))/4
        current=np.sum(pi[None]*gc,axis=1)
        w=[];dw=[]
        for position in q:
            wi,dwi=spherical_weights(self.x,position,self.radius)
            w.append(wi);dw.append(dwi)
        w=np.array(w);dw=np.array(dw)
        e=np.sqrt(self.m*self.m+np.sum(p*p,axis=1))
        abar=np.einsum('kxyz,xyz->k',w,alpha)
        bbar=np.einsum('kxyz,jxyz->kj',w,beta)
        velocity=abar[:,None]*p/e[:,None]+bbar
        force=-e[:,None]*np.einsum('kxyzj,xyz->kj',dw,alpha)-np.einsum('kxyzj,ixyz,ki->kj',dw,beta,p)
        density=alpha*kinetic-np.sum(beta*current,axis=0)+self.omega**2*np.sum(f*f,axis=0)/2
        return f,pi,q,p,alpha,a,scale,b,beta,gp,gm,gc,kinetic,current,w,e,abar,bbar,velocity,force,density

    def rhs(self,y,omit_self=False):
        f,pi,q,p,alpha,a,s,b,beta,gp,gm,gc,k,j,w,e,abar,bbar,v,force,density=self.ingredients(y)
        fd=alpha*pi-np.sum(beta[:,None]*gc,axis=0)
        pd=-self.omega**2*f
        for axis in range(3):
            pd+=.5*(self.minus(alpha*gp[axis],axis+1)+self.plus(alpha*gm[axis],axis+1))
            pd-=self.central(beta[axis]*pi,axis+1)
        if not omit_self:
            pd[0]-=self.g*(alpha*k-np.sum(beta*j,axis=0))
            pd[1:]+=alpha*self.kappa*self.eta*(j/s-self.eta**2*a*np.sum(a*j,axis=0)/s**3)
        esource=np.einsum('k,kxyz->xyz',e,w)/self.dv
        psource=np.einsum('ki,kxyz->ixyz',p,w)/self.dv
        pd[0]-=self.g*alpha*(esource+np.sum(b*psource,axis=0))
        pd[1:]-=alpha*self.kappa*self.eta*(psource/s-self.eta**2*a*np.sum(a*psource,axis=0)/s**3)
        pd-=self.gamma*fd
        loss=np.sum(self.gamma*fd*fd)*self.dv
        result=np.concatenate((fd.ravel(),pd.ravel(),v.ravel(),force.ravel(),[loss]))
        if not np.all(np.isfinite(result)):raise FloatingPointError('Nonfinite CC-2 state derivative')
        return result

    def metrics(self,y):
        f,pi,q,p,alpha,a,s,b,beta,gp,gm,gc,k,j,w,e,abar,bbar,v,force,density=self.ingredients(y)
        field=np.sum(density)*self.dv
        matter=np.sum(abar*e+np.sum(bbar*p,axis=1))
        momentum=p.sum(axis=0)-np.sum(j,axis=(1,2,3))*self.dv
        angular=np.cross(q,p).sum(axis=0)-np.cross(self.x,j.transpose(1,2,3,0)).sum(axis=(0,1,2))*self.dv
        angular+=np.cross(f[1:].transpose(1,2,3,0),pi[1:].transpose(1,2,3,0)).sum(axis=(0,1,2))*self.dv
        relative=np.linalg.norm(v-bbar,axis=1)/abar
        cone=max(0,float(relative.max()-1))
        if self.probe=='photon':cone=max(cone,abs(float(relative[-1]-1)))
        curlz=gc[0,2]-gc[1,1]
        plane=self.n//2;mask=np.sum(self.x[:,:,plane,:2]**2,axis=-1)<4
        circulation=np.sum(curlz[:,:,plane][mask])*self.dx**2
        edge=max(np.max(abs(f[:,0])),np.max(abs(f[:,-1])),np.max(abs(f[:,:,0])),np.max(abs(f[:,:,-1])),np.max(abs(f[:,:,:,0])),np.max(abs(f[:,:,:,-1])))
        return dict(total=field+matter,field=field,matter=matter,absorbed=y[-1],ledger=field+matter+y[-1],
                    momentum=momentum,angular=angular,cone_residual=cone,max_characteristic=float(np.max(alpha+np.linalg.norm(beta,axis=0))),
                    amplitude=float(np.max(abs(f))),circulation_z=circulation,edge_amplitude=float(edge),
                    probe_angle=float(np.arctan2(v[-1,1],v[-1,0])) if self.probe!='none' else None)

    def energy(self,y):return self.metrics(y)['total']


def rk4(model,y,dt):
    k1=model.rhs(y);k2=model.rhs(y+dt*k1/2);k3=model.rhs(y+dt*k2/2);k4=model.rhs(y+dt*k3)
    return y+dt*(k1+2*k2+2*k3+k4)/6
