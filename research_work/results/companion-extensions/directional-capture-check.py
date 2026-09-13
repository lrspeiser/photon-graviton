"""Independent ray quadrature and oblate homoeoid force checks."""
from pathlib import Path
import ast,json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_legendre
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator

P=Path(__file__).resolve().parent
tree=ast.parse((P/'directional-capture.py').read_text())
fieldtree=ast.Module(body=[v for v in tree.body if isinstance(v,ast.FunctionDef) and v.name=='field'],type_ignores=[])
exec(compile(fieldtree,'directional-capture.py','exec'))
ray=[]
for q in [.25,.5,1.]:
 for R,z in [(0.,0.),(1.,0.),(1.,.5),(0.,2.),(10.,0.)]:
  for direction in [(1.,0.,0.),(0.,0.,1.),(-.6,0.,.8)]:
   nx,ny,nz=direction;A=nx*nx+ny*ny+nz*nz/q**2;b=R*nx+z*nz/q**2;c=1+R*R+z*z/q**2;h=c-b*b/A;y=b/np.sqrt(A*h)
   analytic=(np.arctan2(1.,y)-y/(1+y*y))/(2*np.sqrt(A)*h**1.5)
   numerical=quad(lambda t:(1+(R+t*nx)**2+(t*ny)**2+(z+t*nz)**2/q**2)**-2,0,np.inf,epsabs=1e-12,epsrel=1e-11)[0]
   ray.append(abs(analytic-numerical)/max(analytic,1e-14))
assert max(ray)<1e-8,max(ray)
forces=[]
for refined in [False,True]:
 nr,nt,lmax=(768,48,24) if refined else (384,24,12)
 x=np.geomspace(1e-6,3000,nr);u,uw=leggauss(nt);ls=np.arange(0,lmax+1,2)
 basis=np.array([eval_legendre(l,u)*(2*l+1)*uw/2 for l in ls]);p0=np.array([eval_legendre(l,0) for l in ls])
 for q in [.25,.5,1.]:
  rho=(1+x[:,None]**2*(1+u*u*(q**-2-1)))**-2.5
  lut,mass=field(rho)
  for R in [.1,1.,3.,10.]:
   # Classical oblate homeoid integration; dimensionless v^2/(4 pi G rho0 a^2).
   exact=q*R*R/2*quad(lambda t:(1+R*R/(1+t))**-2.5/((1+t)**2*np.sqrt(q*q+t)),0,np.inf,epsabs=1e-12,epsrel=1e-11)[0]
   approx=float(lut(np.log(R)));forces.append(dict(refined=refined,q=q,R_over_a=R,relative_error=abs(approx/exact-1)))
out=dict(ray_cases=len(ray),max_ray_integral_relative_error=max(ray),ellipsoidal_force_cases=forces,
 scope='Known ray quadrature and smooth oblate density controls; these do not certify arbitrary physical capture or convergence of every production profile.')
(P/'directional-capture-check-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print('ray',max(ray));print('force coarse/fine',[max(r['relative_error'] for r in forces if r['refined']==f) for f in [False,True]])
assert max(r['relative_error'] for r in forces if r['refined'])<.003
