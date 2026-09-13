from pathlib import Path
import ast,json
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.special import gamma,gammainc
from numpy.polynomial.legendre import leggauss
P=Path(__file__).resolve().parent;BASE=P.parent/'isotropic-galaxy-transfer';G=4.30091727003628e-6
cp=json.loads((BASE/'third-radiation-retention-results.json').read_text())['models']['attenuated']
components={'I':[(460*2.32e7,0,.3),(1700*2.32e7,5.3,.25),(1700*2.32e7,2.6,.8)],'II':[(1600*2.32e7,4.8,.25),(1700*2.32e7,2.,.8)]}
gas=[(53.1,7.,4.),(2180.,1.5,12.)]
tree=ast.parse((BASE/'milky-way-current.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='capture'],type_ignores=[]),str(BASE/'milky-way-current.py'),'exec'))
inputs=json.loads((P.parent/'milky-way-capture/inputs.json').read_text());obs=inputs['eilers']['rows'];Robs=np.array([v['R_kpc'] for v in obs]);Vobs=np.array([v['vc_kms'] for v in obs]);inner=np.arange(len(obs))<20
zobs=inputs['bovy']['rows'];Rz=np.array([v['R_kpc'] for v in zobs]);Kobs=np.array([v['Kz_over_2piG_Msun_pc2'] for v in zobs])

def spherical_binding(r,m):
 dm=np.diff(m);mid=(r[1:]+r[:-1])/2
 inv=dm/mid;outside=np.r_[np.cumsum(inv[::-1])[::-1],0.]
 return G*(np.divide(m,r,out=np.zeros_like(m),where=r>0)+outside)

class Galaxy:
 def __init__(self,variant):
  self.variant=variant;self.r,self.md,self.eta,_=capture(2.6,sum(v[0] for v in components[variant])/.5,8192,192)
  r=self.r;safe=np.maximum(r,1e-30);self.sigma=sum(S*1e6*np.exp(-hole/safe-safe/h) for S,h,hole in gas)
  self.mg=2*np.pi*cumulative_trapezoid(r*self.sigma,r,initial=0)
  # Avoid reciprocal overflow in PCHIP on numerically irrelevant gas tails.
  self.mg[self.mg<1e-100]=0
  self.massD=PchipInterpolator(np.log(r[1:]),self.md[1:]);self.massG=PchipInterpolator(np.log(r[1:]),self.mg[1:])
  self.bindD=PchipInterpolator(np.log(r[1:]),spherical_binding(r,self.md)[1:]);self.bindG=PchipInterpolator(np.log(r[1:]),spherical_binding(r,self.mg)[1:])
  rho=np.gradient(self.md[1:],r[1:])/(4*np.pi*r[1:]**2)
  self.logrho=PchipInterpolator(np.log(r[1:]),np.log(np.maximum(rho,1e-100)))
 def field(self,R,z):
  R,z=np.broadcast_arrays(R,z);r=np.maximum(np.hypot(R,z),1e-5);lr=np.log(r)
  phi=-self.bindD(lr)-self.bindG(lr);fr=G*(self.massD(lr)+self.massG(lr))*R/r**3;fz=G*(self.massD(lr)+self.massG(lr))*z/r**3
  rhoD=np.exp(self.logrho(lr));rhoB=sum(S*1e6*np.exp(-hole/r-r/h) for S,h,hole in gas)/(2*r)
  for M,a,b in components[self.variant]:
   B=np.sqrt(z*z+b*b);d=R*R+(a+B)**2
   phi-=G*M/np.sqrt(d);fr+=G*M*R/d**1.5;fz+=G*M*(a+B)*z/(B*d**1.5)
   rhoB+=b*b*M/(4*np.pi)*(a*R*R+(a+3*B)*(a+B)**2)/(d**2.5*B**3)
  return phi,fr,fz,rhoB,rhoD

def scores(v):
 return {k:dict(n=int(mask.sum()),RMSE_kms=float(np.sqrt(np.mean((v[mask]-Vobs[mask])**2))),bias_kms=float(np.mean(v[mask]-Vobs[mask]))) for k,mask in [('inner',inner),('outer',~inner),('all',np.ones(len(inner),dtype=bool))]}
