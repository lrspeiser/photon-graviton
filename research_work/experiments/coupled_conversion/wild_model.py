"""Explicitly changed effective laws, with reciprocal Hamiltonian derivatives."""
import numpy as np
from model import Model

class WildModel(Model):
    def __init__(self,power=1,mu2=.16,rest=0.,**kwargs):
        super().__init__(**kwargs)
        self.power=power;self.mu2=mu2;self.rest=rest
        self.offset=mu2**2/(4*self.lam) if mu2<0 else 0.
    def F(self,phi):return self.g*phi**self.power
    def Fp(self,phi):return self.g*self.power*phi**(self.power-1)
    def U(self,phi):return .5*self.mu2*phi**2+self.lam/4*phi**4+self.offset
    def Up(self,phi):return self.mu2*phi+self.lam*phi**3
    def internal_weights(self):return self.I+self.M*self.rest
    def sectors(self):
        a=np.exp(-self.F(self.phi))
        em=.5*np.sum(a*self.D**2)*self.dv
        for d in range(self.dim):
            da=(np.roll(self.A,-1,d)-self.A)/self.dx
            em+=.25*np.sum((a+np.roll(a,-1,d))*da**2)*self.dv
        receiving=np.sum(self.receiving_density())*self.dv
        ph,_,_=self.interpolation(self.phi,self.X)
        material=np.sum(self.internal_weights()*np.exp(-self.b*self.F(ph)))
        return np.array([em,receiving,material,np.sum(self.P**2/(2*self.M[:,None]))])
    def receiving_density(self):
        ans=.5*self.Pi**2+self.U(self.phi)
        for d in range(self.dim):
            ff=(np.roll(self.phi,-1,d)-self.phi)/self.dx
            ans+=.25*self.v**2*(ff**2+np.roll(ff,1,d)**2)
        return ans
    def forces(self):
        a=np.exp(-self.F(self.phi));fp0=self.Fp(self.phi)
        fd=np.zeros(self.shape)
        fpi=self.v**2*self.laplacian(self.phi)-self.Up(self.phi)
        for d in range(self.dim):
            slope=(np.roll(self.A,-1,d)-self.A)/self.dx
            flux=.5*(a+np.roll(a,-1,d))*slope
            fd+=(flux-np.roll(flux,1,d))/self.dx
            fpi+=.25*fp0*a*(slope*slope+np.roll(slope,1,d)**2)
        ph,grad,stencil=self.interpolation(self.phi,self.X)
        omega=np.exp(-self.b*self.F(ph))
        source=self.b*self.Fp(ph)*self.internal_weights()*omega
        for idx,weight in stencil:np.add.at(fpi,idx,source*weight/self.dv)
        return fd,fpi,source[:,None]*grad,omega
    def em_kinetic_step(self,dt):
        a=np.exp(-self.F(self.phi))
        self.A+=dt*a*self.D
        self.Pi+=dt*.5*self.Fp(self.phi)*a*self.D**2
    def rhs(self):
        fd,fpi,fp,omega=self.forces()
        a=np.exp(-self.F(self.phi))
        fpi+=.5*self.Fp(self.phi)*a*self.D**2
        return dict(A=a*self.D,D=fd,phi=self.Pi.copy(),Pi=fpi,
                    X=self.P/self.M[:,None],P=fp,theta=omega)
