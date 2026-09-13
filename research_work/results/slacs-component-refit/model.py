"""Spherical Abel deprojection of fixed image components, with shared force."""
from pathlib import Path
import importlib.util
import numpy as np
from scipy.special import gamma
from scipy.integrate import cumulative_trapezoid,quad
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
spec=importlib.util.spec_from_file_location('annular',Path(__file__).resolve().parent.parent/'slacs-resolved-fit/model.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
G=base.G;C=299792.458;ARCSEC=206264.80624709636

def deproject(r,components,order=256):
 nodes,weights=np.polynomial.legendre.leggauss(order);density=np.zeros_like(r);total=0.
 for comp in components:
  Re,n,amp,bn=comp['R'],comp['n'],comp['amp'],comp['bn']
  Rmax=Re*(1+100/bn)**n
  upper=np.arccosh(np.maximum(1.,Rmax/r));u=upper[:,None]*(nodes+1)/2;R=r[:,None]*np.cosh(u)
  derivative=amp*np.exp(-bn*((R/Re)**(1/n)-1))*bn/(n*Re)*(R/Re)**(1/n-1)
  density+=np.sum(derivative*weights,axis=1)*upper/(2*np.pi)
  total+=2*np.pi*amp*Re*Re*n*np.exp(bn)*gamma(2*n)/bn**(2*n)
 return density,total

class ComponentModel(base.AnnularModel):
 def __init__(self,a,edges,psf,A,p,astar,cut,components,n=4001,order=128,deproj_order=256):
  super().__init__(a,edges,psf,A,p,astar,cut,n,order)
  r=self.r;nu,total=deproject(r,components,deproj_order);self.nu=nu
  central_slope=np.log(nu[1]/nu[0])/np.log(r[1]/r[0]);initial=4*np.pi*nu[0]*r[0]**3/(3+central_slope)
  enclosed=initial+4*np.pi*cumulative_trapezoid(r*r*nu,r,initial=0)
  self.total_mass_check=float(enclosed[-1]/total);assert abs(self.total_mass_check-1)<2e-4
  # Normalize the tiny finite-quadrature luminosity discrepancy consistently.
  frac=enclosed/enclosed[-1];self.mass_interp=PchipInterpolator(np.log(r),frac,extrapolate=False);self.inner_power=3+central_slope
  self.den=np.trapezoid(r*r*nu*self.W,r,axis=1)
  self.A=A;self.astar=astar;self.rt=cut*a
  gb=G*1e11*frac/r**2;gct=A*astar*(G*1e11*self.mass_fraction(self.rt)/self.rt**2/astar)**p
  self.forces=np.array([gb,np.where(r<=self.rt,A*astar*(gb/astar)**p,gct*(self.rt/r)**2)])
 def mass_fraction(self,r):
  if r<self.r[0]:return float(self.mass_interp(np.log(self.r[0])))*(r/self.r[0])**self.inner_power
  if r>self.r[-1]:return 1.
  return float(self.mass_interp(np.log(r)))
 def angle(self,mass,Dl,ratio,extra):
  rt=self.rt;p=self.p;astar=self.astar;A=self.A if extra else 0
  gct=A*astar*(G*mass*self.mass_fraction(rt)/rt**2/astar)**p
  def g(r):
   gb=G*mass*self.mass_fraction(r)/r**2
   return gb+(A*astar*(gb/astar)**p if r<=rt else gct*(rt/r)**2)
  def residual(b):
   edge=np.arccos(b/rt) if b<rt else None
   alpha=4/C**2*quad(lambda t:g(b/np.cos(t))*b/np.cos(t),0,np.pi/2,points=[edge] if edge is not None else None,epsabs=1e-6,epsrel=1e-7,limit=200)[0]
   return ratio*alpha-b/Dl
  lo=self.a*1e-7;hi=self.a*1e5
  if residual(lo)*residual(hi)>=0:return None
  return brentq(residual,lo,hi,xtol=1e-9)/Dl*ARCSEC
