"""JR-5 reversible ray conversion and axisymmetric shared-potential solver.

Length kpc, density Msun/kpc^3. Source/capture scaffold inherited from R10,
not derived microscopic gravity. No data changed by these functions.
"""
from __future__ import annotations
import numpy as np
from scipy.special import eval_legendre, expit
from scipy.integrate import cumulative_trapezoid
from numba import njit
G=4.30091727003628e-6
C=299792.458
RHOREF=1e7

@njit(cache=True)
def transfer(r, gas, alpha, beta, quadratic=True):
    """Exact transfer with midpoint density in each spherical radial interval.
    gas is (Nr,Nmu) at cell midpoints, r gives cell endpoints, begins at 0.
    """
    nr,nmu=gas.shape
    f=np.empty((nr,nmu))
    prev=np.full(nmu,.5)
    for i in range(nr):
        ds=r[i+1]-r[i]
        for j in range(nmu):
            rho=gas[i,j]
            a=alpha*rho/1e6
            b=beta*rho/1e6*(rho/RHOREF if quadratic else 1.)
            rate=a+b
            val=prev[j]
            if rate>0:
                eq=a/rate
                val=eq+(val-eq)*np.exp(-rate*ds)
            f[i,j]=val
            prev[j]=val
    return f


def gas_density(r,mu,mass,scale,height_ratio=.1,kind='disk'):
    """Analytic density normalized to the stated gas mass in infinite space."""
    R=r[:,None]*np.sqrt(1-mu[None,:]**2)
    z=r[:,None]*mu[None,:]
    h=height_ratio*scale
    if kind=='sphere':
        return np.broadcast_to(mass/(8*np.pi*scale**3)*np.exp(-r[:,None]/scale),R.shape).copy()
    if kind=='disk':
        return mass/(4*np.pi*scale**2*h)*np.exp(-R/scale-z/h)
    if kind=='ring':
        from scipy.special import erf
        w=.35*scale;R0=2*scale
        integ=w*w*np.exp(-R0*R0/(2*w*w))+R0*w*np.sqrt(np.pi/2)*(1+erf(R0/(np.sqrt(2)*w)))
        return mass/(4*np.pi*h*integ)*np.exp(-.5*((R-R0)/w)**2-z/h)
    raise ValueError(kind)


def baseline(r,A,rc,rt,q):
    """Analytic inward force, density and finite total mass of R10 envelope."""
    fac=expit(q*np.log(r/rc))
    g=A/r*fac/np.hypot(1,r/rt)
    slope=rt*rt/(r*r+rt*rt)+q*expit(-q*np.log(r/rc))
    rho=g/(4*np.pi*G*r)*slope
    return g,rho,A*rt/G

class AxisModel:
    def __init__(self,A,rc,rt,q,Re=1.,nr=256,nmu=32,lmax=12):
        self.A,self.rc,self.rt,self.q,self.Re=A,rc,rt,q,Re
        self.r=Re*np.geomspace(1e-5,1e4,nr)
        self.edges=np.r_[0.,self.r]
        self.mid=(self.edges[:-1]+self.edges[1:])/2
        x,w=np.polynomial.legendre.leggauss(nmu)
        self.mu=(x+1)/2;self.w=w/2
        self.ls=np.arange(0,lmax+1,2)
        self.P=np.array([eval_legendre(l,self.mu) for l in self.ls])
        self.P0=np.array([eval_legendre(l,0.) for l in self.ls])
        self.g0,self.rho0,self.M0=baseline(self.r,A,rc,rt,q)
        self.coeff=None
    def solve(self,mass,scale,alpha,beta,quadratic=True,height=.1,kind='disk'):
        self.gas=gas_density(self.mid,self.mu,mass,scale,height,kind)
        self.f=transfer(self.edges,self.gas,float(alpha),float(beta),quadratic)
        self.delta=self.rho0[:,None]*(2*self.f-1)
        self.coeff=(self.delta*self.w)@self.P.T*(2*self.ls+1)[None,:]
        self._potential()
        return self
    def _potential(self):
        r=self.r
        U=[];gr=[]
        for j,l in enumerate(self.ls):
            density=self.coeff[:,j]
            inner=cumulative_trapezoid(density*r**(l+2),r,initial=0.)
            outer=-cumulative_trapezoid((density*r**(1-l))[::-1],r[::-1],initial=0.)[::-1]
            # r_min core correction is negligible at finite evaluated radii;
            # convergence covers this declared approximation.
            pref=4*np.pi*G/(2*l+1)
            U.append(-pref*(inner/r**(l+1)+r**l*outer))
            gr.append(pref*((l+1)*inner/r**(l+2)-l*r**(l-1)*outer))
        self.U=np.array(U).T;self.gr=np.array(gr).T
    def radial(self,r,mu=0.,monopole_only=False):
        r=np.asarray(r,float)
        base=baseline(r,self.A,self.rc,self.rt,self.q)[0]
        count=1 if monopole_only else len(self.ls)
        for j,l in enumerate(self.ls[:count]):
            base=base+np.interp(np.log(r),np.log(self.r),self.gr[:,j])*eval_legendre(l,mu)
        return base
    def gradient(self,R,z):
        R,z=np.broadcast_arrays(np.asarray(R,float),np.asarray(z,float));r=np.sqrt(R*R+z*z)
        mu=z/r;s=R/r
        gr=baseline(r,self.A,self.rc,self.rt,self.q)[0];Um=np.zeros_like(r)
        for j,l in enumerate(self.ls):
            p=eval_legendre(l,mu)
            gr+=np.interp(np.log(r),np.log(self.r),self.gr[:,j])*p
            if l>0:
                derivative=np.polynomial.legendre.Legendre.basis(l).deriv()(mu)
                Um+=np.interp(np.log(r),np.log(self.r),self.U[:,j])*derivative
        return s*gr-mu*s*Um/r,mu*gr+(1-mu*mu)*Um/r
    def bend(self,b,orientation='face',order=192):
        x,w=np.polynomial.legendre.leggauss(order);t=(x+1)*np.pi/4;wt=w*np.pi/4
        ell=b*np.tan(t);jac=b/np.cos(t)**2
        if orientation=='face':force=self.gradient(np.full_like(ell,b),ell)[0]
        elif orientation=='edge':force=self.radial(np.hypot(b,ell),0.)*b/np.hypot(b,ell)
        elif orientation=='edge_vertical':force=self.gradient(ell,np.full_like(ell,b))[1]
        else:raise ValueError(orientation)
        return float(4/C**2*np.sum(wt*force*jac))
    def diagnostic(self,points=(.5,1.,2.,4.,8.)):
        r=np.asarray(points)*self.Re
        gbase=baseline(r,self.A,self.rc,self.rt,self.q)[0]
        m=4*np.pi*np.trapezoid(self.r**2*(self.rho0+self.coeff[:,0]),self.r)
        polar=np.array([np.interp(np.log(x),np.log(self.r),self.f[:,-1]) for x in r])
        equat=np.array([np.interp(np.log(x),np.log(self.r),self.f[:,0]) for x in r])
        return dict(radii_kpc=r.tolist(),equatorial_force_ratio=(self.radial(r,0)/gbase).tolist(),
          polar_force_ratio=(self.radial(r,1)/gbase).tolist(),monopole_force_ratio=(self.radial(r,0,True)/gbase).tolist(),
          polar_companion_fraction=polar.tolist(),equatorial_companion_fraction=equat.tolist(),
          integrated_mass_ratio=float(m/self.M0),minimum_fraction=float(self.f.min()),maximum_fraction=float(self.f.max()))
