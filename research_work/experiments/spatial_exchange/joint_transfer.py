"""SE-JT1 Hamiltonian: reciprocal scalar/local-transfer/particle dynamics."""
import numpy as np
from local_rotor import curl,derivative,laplacian
from point_probes import kernel


class JointTransfer:
    def __init__(self,n=8,length=8.,masses=(1.,.5,0.),epsilon=2.,direct=1.,g=.08,eta=.08):
        self.n=n;self.length=length;self.h=length/n;self.dv=self.h**3
        self.masses=np.array(masses);self.count=len(masses);self.epsilon=epsilon;self.direct=direct;self.g=g;self.eta=eta
        self.cq=.5;self.kappa=.5;self.omega=.2;self.Omega=1.;self.nf=7*n**3

    def unpack(self,y):
        f=y[:self.nf].reshape(7,self.n,self.n,self.n);p=y[self.nf:2*self.nf].reshape(f.shape)
        particles=y[2*self.nf:].reshape(2,self.count,3)
        return f,p,particles[0],particles[1]

    def evaluate(self,y,flow=True):
        f,mom,x,p=self.unpack(y);phi=f[0];pi=mom[0]
        A=np.moveaxis(f[1:4],0,-1);Q=np.moveaxis(f[4:7],0,-1)
        Pi=np.moveaxis(mom[1:4],0,-1);P=np.moveaxis(mom[4:7],0,-1)
        C=np.exp(2*self.g*phi);a=C*C;alpha=np.sqrt(C);scale=np.sqrt(1+self.eta**2*np.sum(A*A,axis=-1))
        factor=C*self.kappa*self.cq*self.eta
        beta=factor[...,None]*A/scale[...,None]
        B=curl(A,self.h);K=P-self.direct*A-self.epsilon*np.cross(B,Q)
        kin=pi*pi+np.sum(Pi*Pi+K*K,axis=-1);current=np.empty_like(A)
        density=.5*a*kin+.5*(self.omega**2*(phi*phi+np.sum(A*A,axis=-1))+self.Omega**2*np.sum(Q*Q,axis=-1))
        fd0=a*pi;fdA=a[...,None]*Pi;fdQ=a[...,None]*K
        pd0=laplacian(phi,self.h)-self.omega**2*phi
        pdA=laplacian(A,self.h)-self.omega**2*A
        pdQ=self.cq**2*laplacian(Q,self.h)-self.Omega**2*Q
        for j in range(3):
            dp=derivative(phi,j,self.h);dA=derivative(A,j,self.h);dQ=derivative(Q,j,self.h)
            current[...,j]=pi*dp+np.sum(Pi*dA+K*dQ,axis=-1)
            for values,weight in ((phi,1.),(A,1.),(Q,self.cq**2)):
                plus=(np.roll(values,-1,j)-values)/self.h;minus=(values-np.roll(values,1,j))/self.h
                squares=plus*plus+minus*minus
                density+=weight*.25*(squares if values.ndim==3 else np.sum(squares,axis=-1))
            fd0-=beta[...,j]*dp;fdA-=beta[...,j,None]*dA;fdQ-=beta[...,j,None]*dQ
            pd0-=derivative(beta[...,j]*pi,j,self.h)
            pdA-=derivative(beta[...,j,None]*Pi,j,self.h)
            pdQ-=derivative(beta[...,j,None]*K,j,self.h)
        density-=np.sum(beta*current,axis=-1)
        pd0+=-2*self.g*a*kin+2*self.g*np.sum(beta*current,axis=-1)
        def beta_derivative(vector):
            return factor[...,None]*(vector/scale[...,None]-self.eta**2*A*np.sum(A*vector,axis=-1)[...,None]/scale[...,None]**3)
        pdA+=beta_derivative(current)+self.direct*fdQ+self.epsilon*curl(np.cross(Q,fdQ),self.h)
        pdQ-=self.epsilon*np.cross(B,fdQ)
        field=float(np.sum(density)*self.dv);matter=0.;velocity=np.zeros_like(p);force=np.zeros_like(x);cone=[]
        for i,mass in enumerate(self.masses):
            pp=float(p[i]@p[i])
            if mass<0 or (mass==0 and pp==0):raise ValueError('Invalid particle mass/momentum')
            ix,w,dw=kernel(x[i],self.n,self.length);cc=C[ix];al=alpha[ix];bb=beta[ix];aa=A[ix];ss=scale[ix];ff=factor[ix]
            E=np.sqrt(mass*mass+cc*pp);bp=np.sum(bb*p[i],axis=-1);hp=al*E+bp
            matter+=float(np.sum(w*hp));bbar=np.sum(w[...,None]*bb,axis=(0,1,2));cbar=float(np.sum(w*cc))
            velocity[i]=p[i]*np.sum(w*al*cc/E)+bbar
            force[i]=-np.sum(dw*hp[None],axis=(1,2,3))
            pd0[ix]-=w*self.g*(al*(E+cc*pp/E)+2*bp)/self.dv
            response=ff[...,None]*(p[i]/ss[...,None]-self.eta**2*aa*np.sum(aa*p[i],axis=-1)[...,None]/ss[...,None]**3)
            pdA[ix]-=w[...,None]*response/self.dv
            speed=float(np.linalg.norm(velocity[i]-bbar));cone.append(abs(speed-cbar) if mass==0 else max(0.,speed-cbar))
        if not flow:return dict(energy=field+matter,field=field,matter=matter,minimum_density=float(np.min(density)),cone_error=max(cone),max_C=float(np.max(C)))
        fd=np.concatenate((fd0[None],np.moveaxis(fdA,-1,0),np.moveaxis(fdQ,-1,0)))
        pd=np.concatenate((pd0[None],np.moveaxis(pdA,-1,0),np.moveaxis(pdQ,-1,0)))
        result=np.concatenate((fd.ravel(),pd.ravel(),velocity.ravel(),force.ravel()))
        if not np.isfinite(result).all():raise FloatingPointError('Nonfinite joint derivative')
        return result
