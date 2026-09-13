"""Known spherical Jeans and lensing machinery, conditional empirical force."""
import numpy as np
from scipy.integrate import cumulative_trapezoid,quad
from scipy.stats import ncx2
from scipy.optimize import brentq
G=4.30091727003628e-6;C=299792.458;ARCSEC=206264.80624709636

def aperture_coefficients(a,rap,psf_sigma,A,p,astar,cut,n=4001,order=128):
 r=a*np.geomspace(1e-5,1e5,n);nu=a/(r*(r+a)**3)
 if np.isinf(rap):W=np.ones_like(r)
 elif psf_sigma==0:
  q=np.minimum(1.,(rap/r)**2);W=q/(1+np.sqrt(1-q))
 else:
  nodes,weights=np.polynomial.legendre.leggauss(order);u=(nodes+1)/2
  radius=r[:,None]*np.sqrt(1-u[None,:]**2)
  W=ncx2.cdf((rap/psf_sigma)**2,2,(radius/psf_sigma)**2)@(weights/2)
 gb=G*1e11/(r+a)**2;rt=cut*a;gct=A*astar*(G*1e11/(rt+a)**2/astar)**p
 gc=np.where(r<=rt,A*astar*(gb/astar)**p,gct*(rt/r)**2)
 den=np.trapezoid(r*r*nu*W,r)
 result=[]
 for g in [gb,gc]:
  pressure=-cumulative_trapezoid((nu*g)[::-1],r[::-1],initial=0)[::-1]
  result.append(float(np.trapezoid(r*r*pressure*W,r)/den))
 return result

def mass_from_sigma(sigma,cb,cc,p):
 return 1e11*brentq(lambda y:cb*y+cc*y**p-sigma**2,1e-4,1e3,xtol=1e-12)

def angle(mass,a,Dl,ratio,A,p,astar,cut):
 rt=cut*a;gct=A*astar*(G*mass/(rt+a)**2/astar)**p
 def force(r):
  gb=G*mass/(r+a)**2
  gc=A*astar*(gb/astar)**p if r<=rt else gct*(rt/r)**2
  return gb+gc
 def deflection(b):
  edge=np.arccos(b/rt) if b<rt else None
  val=quad(lambda t:force(b/np.cos(t))*b/np.cos(t),0,np.pi/2,points=[edge] if edge is not None else None,epsabs=1e-6,epsrel=1e-8)[0]
  return 4*val/C**2
 residual=lambda b:ratio*deflection(b)-b/Dl
 lo=a*1e-7;hi=a*1e5
 if residual(lo)*residual(hi)>=0:return None
 b=brentq(residual,lo,hi,xtol=1e-10)
 return b/Dl*ARCSEC
