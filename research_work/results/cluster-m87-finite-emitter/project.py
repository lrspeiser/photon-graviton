"""Resolved spherical M87 deposition projection, checked against ray accounting."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator

OUT=Path(__file__).resolve().parent
source=OUT/'results.json'; prior=json.loads(source.read_text())
weights=np.array(prior['gaussian_weights']); scales=np.array(prior['gaussian_sigma_Mpc'])
L=prior['luminosity_Lsun']; a=.0002488993265191759;k=10.

def volume(r):
    if r==0: return prior['center_volume_power_Lsun_Mpc3']
    total=0.
    for w,sigma in zip(weights,scales):
        # s is distance from a deposition point to an emission point.
        # Gaussian angular average evaluated without overflowing sinh.
        def integrand(x):
            s=r+sigma*x
            t=2*r*s/sigma**2
            angular=-np.expm1(-t)/t if t>1e-12 else 1-t/2
            companion=np.exp(-a*s)*(-np.expm1(-(k-a)*s))
            return np.exp(-x*x/2)*angular*companion
        v=quad(integrand,max(-12.,-r/sigma),12.,epsabs=1e-20,epsrel=2e-9)[0]
        total+=w*sigma*v/(2*np.pi*sigma*sigma)**1.5
    return L*k*a/(k-a)*total

def calculate(n):
    radii=np.geomspace(1e-8,1.,n)
    values=np.array([volume(r) for r in radii])
    assert np.all(values>0) and np.all(np.isfinite(values))
    interp=PchipInterpolator(np.log(radii),np.log(values))
    def q(r):
        return float(np.exp(interp(np.log(max(r,1e-8)))))
    def project(b):
        h=np.sqrt(max(0.,1-b*b))
        return 2*quad(lambda z:q(np.hypot(b,z)),0,h,epsabs=.01,epsrel=2e-7,points=[v for v in [1e-5,1e-4,.001,.01,.1] if v<h],limit=300)[0]
    total=4*np.pi*quad(lambda r:r*r*q(r),0,1,epsabs=.01,epsrel=2e-7,points=[.0001,.001,.01,.1],limit=300)[0]
    bs=np.array([0,.00001,.0001,.001,.003,.01,.03,.1,.3,.5,.9,.99,1.])
    surface=np.array([project(b) for b in bs])
    return radii,values,bs,surface,total,q,project

coarse=calculate(300);fine=calculate(600)
rel=np.max(np.abs(coarse[3][:-1]/fine[3][:-1]-1))
ray=prior['extended_emitter_power_Lsun'];total=fine[4]
assert rel<2e-4
assert abs(total/ray-1)<2e-5
assert fine[3][0]<prior['projected_power_global_upper_bound_Lsun_Mpc2']
assert np.all(np.diff(fine[3])<=0)
# Integrate projected power over the image as another conservation check.
projected_total=2*np.pi*quad(lambda b:b*fine[6](b),0,1,epsabs=1.,epsrel=3e-6,points=[.001,.01,.1],limit=150)[0]
assert abs(projected_total/ray-1)<2e-5
rows=[dict(projected_radius_kpc=float(b*1000),surface_power_Lsun_Mpc2=float(v)) for b,v in zip(fine[2],fine[3])]
for row,b in zip(rows,fine[2]):
    sphere=4*np.pi*quad(lambda r:r*r*fine[5](r),0,b,epsabs=1e-10,epsrel=2e-8,limit=200)[0] if b else 0.
    # Fraction of a spherical shell falling inside a projected circular aperture.
    caps=4*np.pi*quad(lambda r:fine[5](r)*b*b/(1+np.sqrt(max(0.,1-(b/r)**2))),b,1,epsabs=1e-10,epsrel=2e-8,limit=200)[0] if 0<b<1 else 0.
    row['sphere_enclosed_power_Lsun']=sphere
    row['projected_aperture_power_Lsun']=sphere+caps
    assert 0<=sphere<=sphere+caps<=ray*(1+2e-5)
result=dict(scope='Resolved M87-only spherical deposit POWER projection; no accumulation time, gravitational coupling, other emitters or lensing fit.',input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
 central_surface_power_Lsun_Mpc2=float(fine[3][0]),radial_grid_refinement_max_relative=float(rel),
 volume_integrated_power_Lsun=total,projected_integrated_power_Lsun=projected_total,independent_ray_power_Lsun=ray,
 volume_vs_ray_relative=total/ray-1,projection_vs_ray_relative=projected_total/ray-1,profile=rows)
(OUT/'projection.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(result,indent=2))
