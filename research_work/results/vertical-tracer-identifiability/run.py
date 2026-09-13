"""Conditional 1D equilibrium counterexample for velocity-only identification."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import cumulative_trapezoid
H=Path(__file__).resolve().parent;R=H.parents[2];p=H.parent/'component-vertical-transfer/results.json';d=json.loads(p.read_text(encoding='utf-8'))
source=H.parent/'component-loading-forces/results.json';assert hashlib.sha256(source.read_bytes()).hexdigest()==d['source_sha256']
rows=sorted([v for v in d['rows'] if v['R_kpc']==8 and v['phi_rad']==0 and v['z_kpc']<=2],key=lambda v:v['z_kpc']);z=np.r_[0.,[v['z_kpc'] for v in rows]]
old=np.r_[0.,[v['target_total_vertical_pull'] for v in rows]];new=np.r_[0.,[v['predicted_total_vertical_pull'] for v in rows]]
# Piecewise-linear force surrogate; its integrated potential is exact on each segment.
grid=np.linspace(0,2,2001)
def potential(force):
 values=np.zeros_like(grid)
 for j in range(len(z)-1):
  width=np.clip(grid-z[j],0,z[j+1]-z[j]);slope=(force[j+1]-force[j])/(z[j+1]-z[j]);values+=force[j]*width+.5*slope*width**2
 return values
P0=potential(old);P1=potential(new);out=[]
for sigma in [20.,40.,60.]:
 logn0=-P0/sigma**2;logn1=-P1/sigma**2;n0=np.exp(logn0);n1=np.exp(logn1)
 velocity=np.linspace(-5*sigma,5*sigma,501);g=np.exp(-.5*(velocity/sigma)**2);g/=np.trapezoid(g,velocity)
 # Admissible illustrative position selection functions make observed densities proportional.
 shift=float(np.max(logn0-logn1));selection0=np.full_like(grid,.5);selection1=.5*np.exp(logn0-logn1-shift)
 assert np.all((selection1>0)&(selection1<=.5+1e-14))
 obs0=n0*selection0;obs1=n1*selection1;obs0/=np.trapezoid(obs0,grid);obs1/=np.trapezoid(obs1,grid)
 error=float(np.max(abs(obs0-obs1)));assert error<1e-12
 # Stationary collisionless equation cancels analytically: v*d_z f - K*d_v f=0.
 probe=np.linspace(.001,1.999,99);K=np.interp(probe,z,new);v=17.;F=np.exp(-(v*v/2+np.interp(probe,grid,P1))/sigma**2)
 residual=v*(-K/sigma**2*F)-K*(-v/sigma**2*F);assert np.max(abs(residual))<1e-13
 out.append(dict(sigma_kms=sigma,density_ratio_new_to_old_z1=float(np.exp(np.interp(1.,grid,logn1-logn0))),density_ratio_z2=float(n1[-1]/n0[-1]),selection_new_range=[float(selection1.min()),float(selection1.max())],maximum_normalized_observed_density_difference=error,collisionless_identity_residual=float(np.max(abs(residual)))))
outcome=dict(scope='Exact conditional 1D identifiability example using piecewise-linear force surrogates; not a Galactic equilibrium or observed selection fit',hashes={str(v.relative_to(R)):hashlib.sha256(v.read_bytes()).hexdigest() for v in [Path(__file__),p]},force_knots=dict(z_kpc=z.tolist(),reference=old.tolist(),candidate=new.tolist()),cases=out)
(H/'results.json').write_text(json.dumps(outcome,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(out))
