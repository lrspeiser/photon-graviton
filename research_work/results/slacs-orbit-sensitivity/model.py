"""Constant-anisotropy Jeans aperture kernel; known mathematics."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.stats import ncx2
G=4.30091727003628e-6

def coefficients(a,rap,psf,A,p,astar,cut,betas,n=4001,order=128):
 r=a*np.geomspace(1e-5,1e5,n);nu=a/(r*(r+a)**3)
 nodes,weights=np.polynomial.legendre.leggauss(order);u=(nodes+1)/2;weights=weights/2
 if np.isinf(rap):W=np.ones_like(r);T=np.full_like(r,2/3)
 else:
  if psf<=0:raise ValueError('Finite aperture requires positive seeing here')
  radius=r[:,None]*np.sqrt(1-u[None,:]**2)
  response=ncx2.cdf((rap/psf)**2,2,(radius/psf)**2)
  W=response@weights;T=response@(weights*(1-u*u))
 gb=G*1e11/(r+a)**2;rt=cut*a;gct=A*astar*(G*1e11/(rt+a)**2/astar)**p
 gc=np.where(r<=rt,A*astar*(gb/astar)**p,gct*(rt/r)**2)
 denominator=np.trapezoid(r*r*nu*W,r);out=[]
 for beta in betas:
  assert beta<1
  result=[]
  factor=(r/a)**(2*beta)
  for g in [gb,gc]:
   pressure=-cumulative_trapezoid((nu*g*factor)[::-1],r[::-1],initial=0)[::-1]/factor
   result.append(float(np.trapezoid(r*r*pressure*(W-beta*T),r)/denominator))
  out.append(result)
 return np.array(out)
