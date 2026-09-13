import numpy as np
from scipy.stats import ncx2
from scipy.integrate import cumulative_trapezoid
G=4.30091727003628e-6
class AnnularModel:
 def __init__(self,a,edges,psf,A,p,astar,cut,n=4001,order=128):
  self.a=a;self.p=p;self.r=r=a*np.geomspace(1e-5,1e5,n);self.nu=nu=a/(r*(r+a)**3)
  nodes,weights=np.polynomial.legendre.leggauss(order);u=(nodes+1)/2;weights=weights/2
  radius=r[:,None]*np.sqrt(1-u[None,:]**2);responses=[];press=[]
  for edge in edges:
   response=ncx2.cdf((edge/psf)**2,2,(radius/psf)**2) if edge>0 else np.zeros_like(radius)
   responses.append(response@weights);press.append(response@(weights*(1-u*u)))
  self.W=np.diff(responses,axis=0);self.T=np.diff(press,axis=0)
  assert self.W.min()>-1e-12
  self.W=np.maximum(self.W,0)
  self.den=np.trapezoid(r*r*nu*self.W,r,axis=1);assert np.all(self.den>0)
  gb=G*1e11/(r+a)**2;rt=cut*a;gct=A*astar*(G*1e11/(rt+a)**2/astar)**p
  self.forces=np.array([gb,np.where(r<=rt,A*astar*(gb/astar)**p,gct*(rt/r)**2)])
 def coefficients(self,beta):
  r=self.r;factor=(r/self.a)**(2*beta)
  pressure=-cumulative_trapezoid((self.nu*self.forces*factor)[:,::-1],r[::-1],initial=0,axis=1)[:,::-1]/factor
  return np.array([np.trapezoid(r*r*q*(self.W-beta*self.T),r,axis=1)/self.den for q in pressure])
 def predict(self,mass,beta,extra=True):
  cb,cc=self.coefficients(beta);m=mass/1e11
  return np.sqrt(cb*m+(cc*m**self.p if extra else 0))
