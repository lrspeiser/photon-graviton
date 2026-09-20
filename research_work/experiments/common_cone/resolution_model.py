"""SR-2 compact reciprocal sources and positive higher-order field energy."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from spatial_model import SpatialEvolution
from model import spherical_weights


class ResolutionEvolution(SpatialEvolution):
    def __init__(self,order=4,**kwargs):
        if kwargs.get('lam',0)!=0:raise ValueError('SR-2 implements lambda=0 only')
        if order not in (2,4):raise ValueError('Unsupported derivative order')
        super().__init__(**kwargs)
        self.order=order;self.sigma=1/64;self.max_characteristic_seen=0.

    def d4(self,f,axis):
        return (-np.roll(f,-2,axis)+8*np.roll(f,-1,axis)-8*np.roll(f,1,axis)+np.roll(f,2,axis))/(12*self.dx)

    def fourth(self,f,axis):
        return (np.roll(f,-2,axis)-4*np.roll(f,-1,axis)+6*f-4*np.roll(f,1,axis)+np.roll(f,2,axis))/self.dx**4

    def block(self,q):
        if np.any(q-self.radius<=-self.length/2) or np.any(q+self.radius>=self.length/2):
            raise ValueError('Source support reaches domain edge')
        low=np.maximum(0,np.floor((q-self.radius+self.length/2)/self.dx).astype(int))
        high=np.minimum(self.n,np.ceil((q+self.radius+self.length/2)/self.dx).astype(int)+1)
        sl=tuple(slice(int(a),int(b)) for a,b in zip(low,high))
        w,dw=spherical_weights(self.x[sl],q,self.radius)
        return sl,w,dw

    def evaluate(self,y,derivatives=True,omit_self=False):
        f,pi,q,p=self.unpack(y);s=self.s;g=self.g;eta=self.eta
        U=g*f[0];alpha=np.exp(U);z=np.exp(2*s*U)
        a=np.exp((1+3*s)*U);d=np.exp((1-s)*U);c=np.exp((1+s)*U)
        A=f[1:];scale=np.sqrt(1+eta**2*np.sum(A*A,axis=0))
        beta=c*self.kappa*eta*A/scale
        kinetic=np.sum(pi*pi,axis=0)/2
        gradient=np.zeros_like(U);current=np.zeros((3,)+U.shape)
        curlz=np.zeros_like(U)
        if derivatives:fd=a*pi;pd=-self.omega**2*f
        for j in range(3):
            axis=j+1
            if self.order==2:
                gp=self.plus(f,axis);gm=self.minus(f,axis);gc=(gp+gm)/2
                gradient+=np.sum(gp*gp+gm*gm,axis=0)/4
                if derivatives:
                    pd+=.5*(self.minus(d*gp,axis)+self.plus(d*gm,axis))-self.central(beta[j]*pi,axis)
            else:
                gc=self.d4(f,axis);fourth=self.fourth(f,axis)
                gradient+=np.sum(gc*gc,axis=0)/2+self.sigma*self.dx**6*np.sum(fourth*fourth,axis=0)/2
                if derivatives:
                    pd+=self.d4(d*gc,axis)-self.sigma*self.dx**6*self.fourth(d*fourth,axis)-self.d4(beta[j]*pi,axis)
            current[j]=np.sum(pi*gc,axis=0)
            if derivatives:fd-=beta[j]*gc
            if j==0:curlz+=gc[2]
            if j==1:curlz-=gc[1]
        if derivatives and not omit_self:
            pd[0]-=g*((1+3*s)*a*kinetic+(1-s)*d*gradient-(1+s)*np.sum(beta*current,axis=0))
            pd[1:]+=c*self.kappa*eta*(current/scale-eta**2*A*np.sum(A*current,axis=0)/scale**3)
        field=float(np.sum(a*kinetic+d*gradient-np.sum(beta*current,axis=0)+self.omega**2*np.sum(f*f,axis=0)/2)*self.dv)
        matter=0.;velocity=np.zeros_like(p);force=np.zeros_like(q);bbar=np.zeros_like(p);cbar=np.zeros(self.count)
        for i in range(self.count):
            sl,w,dw=self.block(q[i]);vsl=(slice(None),)+sl
            al=alpha[sl];zz=z[sl];bb=beta[vsl];aa=A[vsl];ss=scale[sl]
            pp=np.dot(p[i],p[i]);e=np.sqrt(self.m[i]**2+zz*pp)
            bp=np.einsum('jxyz,j->xyz',bb,p[i]);h=al*e+bp
            matter+=float(np.sum(w*h))
            bbar[i]=np.sum(w*bb,axis=(1,2,3));cbar[i]=np.sum(w*c[sl])
            velocity[i]=p[i]*np.sum(w*al*zz/e)+bbar[i]
            force[i]=-np.sum(dw*h[...,None],axis=(0,1,2))
            if derivatives:
                pd[0][sl]-=g*w*(al*(e+s*zz*pp/e)+(1+s)*bp)/self.dv
                source=c[sl]*self.kappa*eta*(p[i,:,None,None,None]/ss-eta**2*aa*np.einsum('jxyz,j->xyz',aa,p[i])/ss**3)
                pd[(slice(1,None),)+sl]-=source*w/self.dv
        max_characteristic=float(np.max(c+np.linalg.norm(beta,axis=0)))
        self.max_characteristic_seen=max(self.max_characteristic_seen,max_characteristic)
        if derivatives:
            pd-=self.gamma*fd
            loss=float(np.sum(self.gamma*fd*fd)*self.dv)
            rhs=np.concatenate((fd.ravel(),pd.ravel(),velocity.ravel(),force.ravel(),[loss]))
            if not np.all(np.isfinite(rhs)):raise FloatingPointError('Nonfinite SR-2 derivative')
            return rhs
        return dict(f=f,pi=pi,q=q,p=p,field=field,matter=matter,current=current,curlz=curlz,
                    velocity=velocity,bbar=bbar,cbar=cbar,max_characteristic=max_characteristic)

    def rhs(self,y,omit_self=False):return self.evaluate(y,True,omit_self)

    def metrics(self,y):
        t=self.evaluate(y,False);f=t['f'];pi=t['pi'];q=t['q'];p=t['p'];j=t['current']
        relative=np.linalg.norm(t['velocity']-t['bbar'],axis=1)/t['cbar']
        cone=max(0.,float(relative.max()-1))
        if self.probe=='photon':cone=max(cone,abs(float(relative[-1]-1)))
        momentum=p.sum(axis=0)-np.sum(j,axis=(1,2,3))*self.dv
        angular=np.cross(q,p).sum(axis=0)-np.cross(self.x,j.transpose(1,2,3,0)).sum(axis=(0,1,2))*self.dv
        angular+=np.cross(f[1:].transpose(1,2,3,0),pi[1:].transpose(1,2,3,0)).sum(axis=(0,1,2))*self.dv
        edge=max(np.max(abs(np.take(f,i,axis=axis))) for axis in (1,2,3) for i in (0,-1))
        plane=self.n//2;mask=np.sum(self.x[:,:,plane,:2]**2,axis=-1)<4
        # Rotate back before interpreting a signed ray bend.
        vback=t['velocity'][-1]@self.rot
        return dict(total=t['field']+t['matter'],field=t['field'],matter=t['matter'],absorbed=float(y[-1]),ledger=t['field']+t['matter']+y[-1],
                    momentum=momentum,angular=angular,cone_residual=cone,max_characteristic=t['max_characteristic'],
                    edge_amplitude=float(edge),source_extent=float(np.max(abs(q))+self.radius),
                    circulation_z=float(np.sum(t['curlz'][:,:,plane][mask])*self.dx**2),probe_velocity=t['velocity'][-1],
                    probe_angle=float(np.arctan2(vback[1],vback[0])),probe_out_of_plane=float(np.arctan2(vback[2],np.linalg.norm(vback[:2]))))
