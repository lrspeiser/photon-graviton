"""Conservative Newtonian reference fields and proposed companion additions.

Units: kpc, km/s, Msun. Density and potential mathematics are established.
The companion interpretation/coupling is a hypothesis, not derived here.
"""
from pathlib import Path
import numpy as np
from scipy.special import sph_harm_y, gamma, gammainc, gammaincc
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from numpy.polynomial.legendre import leggauss
from published_bar import makeBarDensity

G=4.30091727003628e-6
RO=8.;VO=220.
H=Path(__file__).resolve().parent
RHO_BAR=makeBarDensity()

def bar(xyz):
    with np.errstate(over='ignore',divide='ignore',invalid='ignore'):
        rho=RHO_BAR(np.asarray(xyz))
    if not np.isfinite(rho).all():raise ValueError('Non-finite bar density')
    return rho

def nuclei(xyz):
    """Hunter24 NSC and NSD components; no original 2e9-Msun CMC added."""
    R2=xyz[:,0]**2+xyz[:,1]**2;z=xyz[:,2]
    m=np.sqrt(R2+(z/.73)**2)
    a=.0059
    norm=6.1e7/(4*np.pi*.73*quad(lambda r:r*r*(r/a)**(-.71)*(1+r/a)**(-3.29)*np.exp(-(r/.1)**2),0,2,epsabs=1e-12)[0])
    rho=norm*(m/a)**(-.71)*(1+m/a)**(-3.29)*np.exp(-(m/.1)**2)
    m=np.sqrt(R2+(z/.37)**2)
    for density,cut,p in [(2.00583e12,.00506,.72),(1.53e12,.0246,.79)]:
        rho+=density*np.exp(-(m/cut)**p)
    return rho

class Multipole:
    """Real even-l, even-m cosine harmonics for reflection-symmetric sources.

    phi_lm = -4piG/(2l+1) [r^(-l-1) int_0^r rho_lm s^(l+2)ds
                           +r^l int_r^infty rho_lm s^(1-l)ds].
    A known spherical-harmonic solution of Poisson's equation.
    """
    def __init__(self,r,coeff,l,m):
        self.r=r;self.coeff=coeff;self.l=l;self.m=m
        self.spline=CubicSpline(np.log(r),coeff,axis=0)

    @classmethod
    def build(cls,density,lmax=16,nr=520,ntheta=64,nphi=96,rmin=1e-6,rmax=150):
        r=np.geomspace(rmin,rmax,nr)
        l=np.array([ll for ll in range(0,lmax+1,2) for mm in range(0,ll+1,2)])
        m=np.array([mm for ll in range(0,lmax+1,2) for mm in range(0,ll+1,2)])
        mu,w=leggauss(ntheta);phi=2*np.pi*np.arange(nphi)/nphi
        mu=np.repeat(mu,nphi);theta=np.arccos(mu);phi=np.tile(phi,ntheta)
        unit=np.c_[np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),mu]
        weights=np.repeat(w,nphi)*2*np.pi/nphi
        basis=sph_harm_y(l[:,None],m[:,None],theta[None,:],phi[None,:]).real
        basis*=np.where(m[:,None]>0,np.sqrt(2),1)
        angular=(basis*weights).T
        rho_lm=[]
        for start in range(0,nr,20):
            rr=r[start:start+20]
            vals=density((rr[:,None,None]*unit[None,:,:]).reshape(-1,3)).reshape(len(rr),-1)
            rho_lm.append(vals@angular)
        rho_lm=np.vstack(rho_lm)
        # Recursively accumulate already-scaled integrals. All radius ratios are
        # <=1, avoiding overflow from r**(-l) at high angular order near the centre.
        inside=np.zeros_like(rho_lm);outside=np.zeros_like(rho_lm)
        for i in range(1,nr):
            ratio=(r[i-1]/r[i])**(l+1)
            inside[i]=ratio*inside[i-1]+.5*(r[i]-r[i-1])*(rho_lm[i-1]*r[i-1]*ratio+rho_lm[i]*r[i])
        for i in range(nr-2,-1,-1):
            ratio=(r[i]/r[i+1])**l
            outside[i]=ratio*outside[i+1]+.5*(r[i+1]-r[i])*(rho_lm[i]*r[i]+rho_lm[i+1]*r[i+1]*ratio)
        coeff=-4*np.pi*G/(2*l+1)*(inside+outside)
        obj=cls(r,coeff,l,m)
        obj.integrated_mass=float(np.sqrt(4*np.pi)*inside[-1,0]*r[-1])
        return obj

    def save(self,path):
        np.savez_compressed(path,r=self.r,coeff=self.coeff,l=self.l,m=self.m,mass=getattr(self,'integrated_mass',np.nan))

    @classmethod
    def load(cls,path):
        f=np.load(path);obj=cls(f['r'],f['coeff'],f['l'],f['m']);obj.integrated_mass=float(f['mass']);return obj

    def evaluate(self,xyz):
        xyz=np.atleast_2d(xyz).astype(float);x,y,z=xyz.T
        R=np.hypot(x,y);r=np.sqrt(R*R+z*z)
        if np.any((r<self.r[0])|(r>self.r[-1])|(R==0)):
            raise ValueError('Outside numerical domain or exactly on polar axis')
        theta=np.arctan2(R,z);phi=np.arctan2(y,x)
        Y,der=sph_harm_y(self.l[:,None],self.m[:,None],theta[None,:],phi[None,:],diff_n=1)
        norm=np.where(self.m[:,None]>0,np.sqrt(2),1)
        Y=Y.real*norm;der=der.real*norm[:,:,None]
        coef=self.spline(np.log(r));dr=self.spline(np.log(r),1)/r[:,None]
        potential=np.sum(coef*Y.T,axis=1)
        gr=np.sum(dr*Y.T,axis=1)
        gt=np.sum(coef*der[:,:,0].T,axis=1)
        gp=np.sum(coef*der[:,:,1].T,axis=1)
        aR=-gr*R/r-gt*z/r**2;az=-gr*z/r+gt*R/r**2;ap=-gp/R
        acc=np.c_[aR*x/R-ap*y/R,aR*y/R+ap*x/R,az]
        return potential,acc

def point_potential(xyz,mass=4.1e6,soft=.001):
    xyz=np.atleast_2d(xyz);q=np.sum(xyz*xyz,axis=1)+soft*soft
    return -G*mass/np.sqrt(q),-G*mass*xyz/q[:,None]**1.5

def halo(xyz):
    """Published Hunter24 spherical Einasto halo; comparison only."""
    xyz=np.atleast_2d(xyz);r=np.linalg.norm(xyz,axis=1)
    rho0=2.774e11;a=8.682e-6;p=.1704;x=(r/a)**p
    M=4*np.pi*rho0*a**3/p*gamma(3/p)*gammainc(3/p,x)
    outer=4*np.pi*rho0*a*a/p*gamma(2/p)*gammaincc(2/p,x)
    return -G*(M/r+outer),-G*M[:,None]*xyz/r[:,None]**3

def companion(xyz,geometry='equatorial',amplitude=1000.,soft=.3):
    """Frozen earlier synthetic deposits; amplitude is response, not supplied energy."""
    p=np.arange(128)*2*np.pi/128
    if geometry=='equatorial':
        centres=np.c_[1.5*np.cos(p),1.5*np.sin(p),np.zeros(len(p))];w=np.ones(len(p))/len(p)
    elif geometry=='caps':
        centres=np.vstack([np.c_[.5*np.cos(p),.5*np.sin(p),np.full(len(p),z)] for z in [-1.5,1.5]])
        w=np.ones(len(centres))/len(centres)
    elif geometry=='shell':
        mu,ww=leggauss(32)
        centres=np.concatenate([np.c_[1.5*np.sqrt(1-u*u)*np.cos(p),1.5*np.sqrt(1-u*u)*np.sin(p),np.full(len(p),1.5*u)] for u in mu])
        w=np.repeat(ww/2/len(p),len(p))
    else:raise ValueError(geometry)
    xyz=np.atleast_2d(xyz);pots=[];acc=[]
    for point in xyz:
        delta=point-centres;q=np.sum(delta*delta,axis=1)+soft**2
        pots.append(-amplitude*np.sum(w/np.sqrt(q)))
        acc.append(-amplitude*np.sum(w[:,None]*delta/q[:,None]**1.5,axis=0))
    return np.array(pots),np.array(acc)

def disks(order=24):
    """Two holed stellar disks and two gas disks; Hunter24 normalization.

    This replaces, rather than supplements, the original Sormani example disk.
    No spiral perturbation is included in this baseline.
    """
    from galpy.potential import DiskSCFPotential
    params=[(1.332e9,2.,.3,2.7,'exp'),(8.97e8,2.8,.9,2.7,'exp'),
            (5.81e7,7.,.085,4.,'sech2'),(2.68e9,1.5,.045,12.,'sech2')]
    unit=VO**2/(G*RO)
    def density(R,z):
        R=np.maximum(np.asarray(R),1e-14);z=np.asarray(z)
        value=0.
        for sigma,rd,h,rhole,kind in params:
            sig=sigma/unit*np.exp(-rhole/RO/R-R/(rd/RO))
            if kind=='exp':value+=sig/(2*h/RO)*np.exp(-np.abs(z)/(h/RO))
            else:value+=sig/(4*h/RO)*np.exp(-2*np.logaddexp(z/(2*h/RO),-z/(2*h/RO))+2*np.log(2))
        return value
    return DiskSCFPotential(dens=density,
        Sigma=[dict(type='exp',amp=sigma/unit,h=rd/RO,Rhole=rhole/RO) for sigma,rd,h,rhole,kind in params],
        hz=[dict(type=kind,h=h/RO) for sigma,rd,h,rhole,kind in params],a=2.5,N=order,L=order,ro=RO,vo=VO)

def disk_evaluate(pot,xyz):
    values=[];forces=[]
    for x,y,z in np.atleast_2d(xyz):
        R=np.hypot(x,y)
        values.append(pot(R/RO,z/RO,use_physical=False)*VO**2)
        ar=pot.Rforce(R/RO,z/RO,use_physical=False)*VO**2/RO
        az=pot.zforce(R/RO,z/RO,use_physical=False)*VO**2/RO
        forces.append([ar*x/R,ar*y/R,az])
    return np.array(values),np.array(forces)
