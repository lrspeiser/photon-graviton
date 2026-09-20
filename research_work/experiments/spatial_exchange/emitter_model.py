"""SE-2: frozen SE-1 equations extended by constant two-channel emission.

Copied from pinned SE-1; keep that numerical source unchanged.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import sys
import importlib.util
import numpy as np

# Avoid the sibling CC-1 module's generic name colliding with this module.
COMMON=Path(__file__).resolve().parents[1]/'common_cone'
spec=importlib.util.spec_from_file_location('se1_cc1_helpers',COMMON/'model.py')
helpers=importlib.util.module_from_spec(spec);spec.loader.exec_module(helpers)


class ExchangeEvolution:
    def __init__(self,n=32,length=16.,radius=.9,chi=0.,mix=.5,emission=.4,excitation=.001,angle=0.,emitter_angle=0.):
        self.n=n;self.length=length;self.radius=radius;self.chi=chi;self.mix=mix;self.emission=emission;self.excitation=excitation
        self.emitter_direction=np.array([np.cos(emitter_angle),np.sin(emitter_angle)])
        self.g=.08;self.eta=.08;self.kappa=.5;self.omega=.2;self.mu=.75;self.U0=.001;self.Omega=1.
        self.dx=length/n;self.dv=self.dx**3;self.nf=6*n**3;self.count=6;self.m=np.ones(6);self.rot=helpers.rotation(angle)
        c=np.arange(n)*self.dx-length/2;self.x=np.stack(np.meshgrid(c,c,c,indexing='ij'),axis=-1)
        self.gamma=2*np.maximum(np.max(abs(self.x),axis=-1)-(length/2-1),0)**2
        self.maximum_characteristic=0.

    def unpack(self,y):
        f=y[:self.nf].reshape(6,self.n,self.n,self.n);pi=y[self.nf:2*self.nf].reshape(f.shape)
        start=2*self.nf;q=y[start:start+18].reshape(6,3);p=y[start+18:start+36].reshape(6,3)
        Q=y[start+36:start+42];P=y[start+42:start+48]
        return f,pi,q,p,Q,P

    def initial(self):
        theta=np.arange(6)*np.pi/3;q=.7*np.column_stack((np.cos(theta),np.sin(theta),np.zeros(6)))
        p=.2*np.column_stack((-np.sin(theta),np.cos(theta),np.zeros(6)))
        return np.concatenate((np.zeros(2*self.nf),(q@self.rot.T).ravel(),(p@self.rot.T).ravel(),np.full(6,np.sqrt(2*self.excitation)/self.Omega),np.zeros(6),[0.]))

    def plus(self,f,axis):return (np.roll(f,-1,axis)-f)/self.dx
    def minus(self,f,axis):return (f-np.roll(f,1,axis))/self.dx
    def central(self,f,axis):return (np.roll(f,-1,axis)-np.roll(f,1,axis))/(2*self.dx)

    def block(self,q):
        if np.any(q-self.radius<=-self.length/2) or np.any(q+self.radius>=self.length/2):raise ValueError('Source reaches box edge')
        low=np.maximum(0,np.floor((q-self.radius+self.length/2)/self.dx).astype(int))
        high=np.minimum(self.n,np.ceil((q+self.radius+self.length/2)/self.dx).astype(int)+1)
        sl=tuple(slice(int(a),int(b)) for a,b in zip(low,high));w,dw=helpers.spherical_weights(self.x[sl],q,self.radius)
        return sl,w,dw

    def evaluate(self,y,derivatives=True):
        f,pi,q,p,Q,P=self.unpack(y);U=self.g*f[0];alpha=np.exp(U);z=np.exp(2*U);a=np.exp(4*U);A=f[1:4]
        scale=np.sqrt(1+self.eta**2*np.sum(A*A,axis=0));beta=z*self.kappa*self.eta*A/scale
        tangent=np.tanh(U/self.U0);mix=self.mix*tangent;mix_U=self.mix/self.U0*(1-tangent*tangent)
        mass2=self.mu**2*np.exp(2*self.chi*U);residual=f[5]-mix*f[4]
        potential_mix=mass2*residual*residual/2
        kinetic=np.sum(pi*pi,axis=0)/2;gradient=np.zeros_like(U);current=np.zeros((3,)+U.shape)
        # Per-channel propagation terms retain signed shift energy and are positive.
        propagation=a[None]*pi*pi/2
        if derivatives:
            fd=a*pi;pd=np.zeros_like(f);pd[:4]-=self.omega**2*f[:4]
            pd[0]-=self.g*mass2*(self.chi*residual*residual-mix_U*f[4]*residual)
            pd[4]+=mass2*mix*residual;pd[5]-=mass2*residual
        for j in range(3):
            axis=j+1;gp=self.plus(f,axis);gm=self.minus(f,axis);gc=(gp+gm)/2
            gradient+=np.sum(gp*gp+gm*gm,axis=0)/4;current[j]=np.sum(pi*gc,axis=0)
            propagation+=(gp*gp+gm*gm)/4-beta[j]*pi*gc
            if derivatives:
                fd-=beta[j]*gc
                pd+=.5*(self.minus(gp,axis)+self.plus(gm,axis))-self.central(beta[j]*pi,axis)
        if derivatives:
            pd[0]-=self.g*(4*a*kinetic-2*np.sum(beta*current,axis=0))
            pd[1:4]+=z*self.kappa*self.eta*(current/scale-self.eta**2*A*np.sum(A*current,axis=0)/scale**3)
        density=np.sum(propagation,axis=0)+self.omega**2*np.sum(f[:4]**2,axis=0)/2+potential_mix
        field=float(np.sum(density)*self.dv);matter=0.;internal=0.
        velocity=np.zeros_like(p);force=np.zeros_like(q);Qdot=np.zeros(6);Pdot=np.zeros(6);bbar=np.zeros_like(p);cbar=np.zeros(6);masses=np.zeros(6)
        for i in range(6):
            sl,w,dw=self.block(q[i]);vsl=(slice(None),)+sl
            X=self.emitter_direction[0]*f[4][sl]+self.emitter_direction[1]*f[5][sl];xbar=float(np.sum(w*X));difference=Q[i]-self.emission*xbar
            m=self.m[i]+(P[i]**2+self.Omega**2*difference**2)/2;masses[i]=m;internal+=m-self.m[i]
            al=alpha[sl];zz=z[sl];bb=beta[vsl];avec=A[vsl];ss=scale[sl]
            pp=np.dot(p[i],p[i]);e=np.sqrt(m*m+zz*pp);bp=np.einsum('jxyz,j->xyz',bb,p[i]);h=al*e+bp
            response=float(np.sum(w*al*m/e));reaction=response*self.Omega**2*self.emission*difference
            matter+=float(np.sum(w*h));bbar[i]=np.sum(w*bb,axis=(1,2,3));cbar[i]=np.sum(w*zz)
            velocity[i]=p[i]*np.sum(w*al*zz/e)+bbar[i]
            force[i]=-np.sum(dw*h[...,None],axis=(0,1,2))+reaction*np.sum(dw*X[...,None],axis=(0,1,2))
            Qdot[i]=response*P[i];Pdot[i]=-response*self.Omega**2*difference
            if derivatives:
                pd[0][sl]-=self.g*w*(al*(e+zz*pp/e)+2*bp)/self.dv
                source=zz*self.kappa*self.eta*(p[i,:,None,None,None]/ss-self.eta**2*avec*np.einsum('jxyz,j->xyz',avec,p[i])/ss**3)
                pd[(slice(1,4),)+sl]-=source*w/self.dv
                pd[4][sl]+=self.emitter_direction[0]*reaction*w/self.dv
                pd[5][sl]+=self.emitter_direction[1]*reaction*w/self.dv
        characteristic=float(np.max(z+np.linalg.norm(beta,axis=0)));self.maximum_characteristic=max(self.maximum_characteristic,characteristic)
        if derivatives:
            pd-=self.gamma*fd;loss=float(np.sum(self.gamma*fd*fd)*self.dv)
            result=np.concatenate((fd.ravel(),pd.ravel(),velocity.ravel(),force.ravel(),Qdot,Pdot,[loss]))
            if not np.all(np.isfinite(result)):raise FloatingPointError('Nonfinite SE-1 derivative')
            return result
        return dict(f=f,pi=pi,q=q,p=p,current=current,field=field,matter=matter,internal_rest=internal,masses=masses,
                    velocity=velocity,bbar=bbar,cbar=cbar,characteristic=characteristic,
                    radiation_propagation=float(np.sum(propagation[4])*self.dv),companion_propagation=float(np.sum(propagation[5])*self.dv),mixing_potential=float(np.sum(potential_mix)*self.dv))

    def rhs(self,y):return self.evaluate(y,True)

    def metrics(self,y):
        t=self.evaluate(y,False);f=t['f'];pi=t['pi'];q=t['q'];p=t['p'];j=t['current']
        momentum=p.sum(axis=0)-np.sum(j,axis=(1,2,3))*self.dv
        angular=np.cross(q,p).sum(axis=0)-np.cross(self.x,j.transpose(1,2,3,0)).sum(axis=(0,1,2))*self.dv
        angular+=np.cross(f[1:4].transpose(1,2,3,0),pi[1:4].transpose(1,2,3,0)).sum(axis=(0,1,2))*self.dv
        cone=max(0.,float(np.max(np.linalg.norm(t['velocity']-t['bbar'],axis=1)/t['cbar'])-1))
        edge=max(np.max(abs(np.take(f,k,axis=axis))) for axis in (1,2,3) for k in (0,-1))
        result={key:t[key] for key in ('field','matter','internal_rest','radiation_propagation','companion_propagation','mixing_potential')}
        result.update(total=t['field']+t['matter'],ledger=t['field']+t['matter']+y[-1],absorbed=float(y[-1]),momentum=momentum,angular=angular,
                      minimum_mass=float(t['masses'].min()),cone_error=cone,max_characteristic=t['characteristic'],edge_amplitude=float(edge),source_extent=float(np.max(abs(q))+self.radius))
        for index,name in ((4,'radiation'),(5,'companion')):
            weight=f[index]**2;norm=float(weight.sum());result[name+'_amplitude']=float(np.max(abs(f[index])))
            if norm>0:
                center=np.sum(weight[...,None]*self.x,axis=(0,1,2))/norm
                radius=float(np.sqrt(np.sum(weight*np.sum((self.x-center)**2,axis=-1))/norm))
                result[name+'_center']=center;result[name+'_rms_radius']=radius
            else:result[name+'_center']=None;result[name+'_rms_radius']=None
        return result

    def energy(self,y):return self.metrics(y)['total']


def rk4(model,y,dt):
    k1=model.rhs(y);k2=model.rhs(y+dt*k1/2);k3=model.rhs(y+dt*k2/2);k4=model.rhs(y+dt*k3)
    return y+dt*(k1+2*k2+2*k3+k4)/6
