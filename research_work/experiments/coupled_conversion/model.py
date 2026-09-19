"""CWC-1 independent lattice implementation of credited Hamiltonian mathematics.
See provenance.md: matched electromagnetic media, scalar fields and symmetric
splitting are prior art. This is not a spin-2 or universal-gravity derivation.
"""
import itertools
import numpy as np

class Model:
    def __init__(self,n=512,dim=1,L=48,g=.2,b=2,m=.4,v=1.,lam=.1,clock_energy=1.):
        self.n=n;self.dim=dim;self.L=L;self.dx=L/n;self.dv=self.dx**dim
        self.g=g;self.b=b;self.m=m;self.v=v;self.lam=lam
        self.shape=(n,)*dim
        self.axis=np.arange(n)*self.dx-L/2
        self.coords=np.meshgrid(*([self.axis]*dim),indexing='ij')
        self.A=np.zeros(self.shape);self.D=self.A.copy()
        self.phi=self.A.copy();self.Pi=self.A.copy()
        if dim==1:self.X=np.array([[-1.],[-1.],[1.],[1.]])
        else:self.X=np.array([[-1.,0.],[1.,0.],[0.,-1.],[0.,1.]])
        self.P=np.zeros_like(self.X);self.M=np.ones(4)
        self.I=np.array([.0005,.001,.0005,.001])*clock_energy
        self.theta=np.zeros(4)
        self.omit_em_reciprocity=False

    def interpolation(self,field,points):
        points=np.atleast_2d(points)
        s=(points+self.L/2)/self.dx
        cell=np.floor(s).astype(int);f=s-cell
        val=np.zeros(len(points));grad=np.zeros((len(points),self.dim));stencil=[]
        for bits in itertools.product([0,1],repeat=self.dim):
            idx=tuple((cell[:,d]+bits[d])%self.n for d in range(self.dim))
            factors=[f[:,d] if bits[d] else 1-f[:,d] for d in range(self.dim)]
            weight=np.prod(factors,axis=0)
            vv=field[idx];val+=weight*vv
            for d in range(self.dim):
                other=np.prod([factors[k] for k in range(self.dim) if k!=d],axis=0) if self.dim>1 else np.ones(len(points))
                grad[:,d]+=vv*(1 if bits[d] else -1)*other/self.dx
            stencil.append((idx,weight))
        return val,grad,stencil

    def laplacian(self,field):
        return sum((np.roll(field,-1,d)-2*field+np.roll(field,1,d))/self.dx**2 for d in range(self.dim))

    def grad(self,field):
        return [(np.roll(field,-1,d)-np.roll(field,1,d))/(2*self.dx) for d in range(self.dim)]

    def sectors(self):
        a=np.exp(-self.g*self.phi)
        kinetic_em=.5*np.sum(a*self.D**2)*self.dv
        magnetic=0.;gradient=0.
        for d in range(self.dim):
            da=(np.roll(self.A,-1,d)-self.A)/self.dx
            dp=(np.roll(self.phi,-1,d)-self.phi)/self.dx
            magnetic+=.25*np.sum((a+np.roll(a,-1,d))*da**2)*self.dv
            gradient+=.5*self.v**2*np.sum(dp*dp)*self.dv
        receiving=(.5*np.sum(self.Pi**2+self.m**2*self.phi**2)+self.lam/4*np.sum(self.phi**4))*self.dv+gradient
        ph,_,_=self.interpolation(self.phi,self.X)
        material=np.sum(self.I*np.exp(-self.b*self.g*ph))
        motion=np.sum(self.P**2/(2*self.M[:,None]))
        return np.array([kinetic_em+magnetic,receiving,material,motion])

    def energy(self):return self.sectors().sum()

    def receiving_density(self):
        ans=.5*(self.Pi**2+self.m**2*self.phi**2)+self.lam/4*self.phi**4
        for d in range(self.dim):
            ff=(np.roll(self.phi,-1,d)-self.phi)/self.dx
            ans+=.25*self.v**2*(ff**2+np.roll(ff,1,d)**2)
        return ans

    def momentum(self):
        return np.array([-np.sum(self.D*dA+self.Pi*dp)*self.dv
                         for dA,dp in zip(self.grad(self.A),self.grad(self.phi))])+self.P.sum(axis=0)

    def forces(self):
        a=np.exp(-self.g*self.phi)
        fd=np.zeros(self.shape)
        fpi=self.v**2*self.laplacian(self.phi)-self.m**2*self.phi-self.lam*self.phi**3
        for d in range(self.dim):
            slope=(np.roll(self.A,-1,d)-self.A)/self.dx
            flux=.5*(a+np.roll(a,-1,d))*slope
            fd+=(flux-np.roll(flux,1,d))/self.dx
            if not self.omit_em_reciprocity:
                fpi+=.25*self.g*a*(slope*slope+np.roll(slope,1,d)**2)
        ph,grad,stencil=self.interpolation(self.phi,self.X)
        omega=np.exp(-self.b*self.g*ph)
        source=self.b*self.g*self.I*omega
        for idx,weight in stencil:np.add.at(fpi,idx,source*weight/self.dv)
        fp=source[:,None]*grad
        return fd,fpi,fp,omega

    def potential_step(self,dt):
        fd,fpi,fp,omega=self.forces()
        self.D+=dt*fd;self.Pi+=dt*fpi;self.P+=dt*fp;self.theta+=dt*omega

    def em_kinetic_step(self,dt):
        a=np.exp(-self.g*self.phi)
        self.A+=dt*a*self.D
        if not self.omit_em_reciprocity:self.Pi+=dt*.5*self.g*a*self.D**2

    def step(self,dt):
        self.potential_step(dt/2)
        self.em_kinetic_step(dt/2)
        self.phi+=dt*self.Pi
        self.X+=dt*self.P/self.M[:,None]
        self.X=(self.X+self.L/2)%self.L-self.L/2
        self.em_kinetic_step(dt/2)
        self.potential_step(dt/2)

    def rhs(self):
        fd,fpi,fp,omega=self.forces()
        a=np.exp(-self.g*self.phi)
        if not self.omit_em_reciprocity:fpi+=.5*self.g*a*self.D**2
        return dict(A=a*self.D,D=fd,phi=self.Pi.copy(),Pi=fpi,
                    X=self.P/self.M[:,None],P=fp,theta=omega)

    def initialize_light(self,energy=1.,wavelength=2.,reverse=False):
        k=2*np.pi/wavelength
        if self.dim==1:
            x=self.coords[0]
            for center in [-14.,-10.]:
                u=x-center
                envelope=np.exp(-u*u/(2*.8**2))
                self.A+=envelope*np.cos(k*u)
                self.D+=envelope*(u/.8**2*np.cos(k*u)+k*np.sin(k*u))
            if reverse:
                # Mirror about receiver/source midpoint: equivalent reversed fixture.
                self.A=np.roll(self.A[::-1],1)
                self.D=-np.roll(self.D[::-1],1)
                self.X=-self.X
        else:
            x,y=self.coords
            for axis in [0,1]:
                longitudinal=[x,y][axis];transverse=[y,x][axis]
                for center in [-6.,6.]:
                    u=longitudinal-center;direction=-np.sign(center)
                    envelope=np.exp(-u*u/(2*.8**2)-transverse**2/(2*1.5**2))
                    self.A+=envelope*np.cos(k*u)
                    self.D+=direction*envelope*(u/.8**2*np.cos(k*u)+k*np.sin(k*u))
        actual=self.sectors()[0]
        scale=np.sqrt(energy/actual) if energy>0 else 0.
        self.A*=scale;self.D*=scale

    def state(self):
        return {k:getattr(self,k).copy() for k in ['A','D','phi','Pi','X','P','theta']}

    def set_state(self,state):
        for k,v in state.items():setattr(self,k,v.copy())
