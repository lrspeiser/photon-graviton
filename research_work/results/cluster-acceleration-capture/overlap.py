"""Nonadditive weak-capture area for two fixed Plummer gravitational wells."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre
OUT=Path(__file__).resolve().parent
single=3*np.pi**2/4 # integral |g|^2 dV, G=M=a=1

def cross_potential(d):
    # Integration by parts: integral g1.g2 dV = -4*pi*G*integral rho1*Phi2 dV.
    def f(r):
        average=2/(np.sqrt(1+(r+d)**2)+np.sqrt(1+(r-d)**2))
        return 3*r*r/(1+r*r)**2.5*average
    return 4*np.pi*quad(f,0,np.inf,epsabs=1e-10,epsrel=1e-10)[0]

def cross_direct(d,n):
    mu,w=roots_legendre(n)
    def f(r):
        gdot=(r*r-r*d*mu)/((1+r*r)**1.5*(1+r*r+d*d-2*r*d*mu)**1.5)
        return 2*np.pi*r*r*np.dot(w,gdot)
    breaks=sorted(set([0.,max(0.,d-5),d,d+5]))+[np.inf]
    return sum(quad(f,l,h,epsabs=1e-9,epsrel=1e-9,limit=200)[0] for l,h in zip(breaks[:-1],breaks[1:]))

rows=[]
for d in [0.,.5,1.,2.,5.,10.,20.]:
    potential=cross_potential(d);coarse=cross_direct(d,128);fine=cross_direct(d,256)
    err=abs(fine/potential-1)
    assert err<1e-7 and abs(fine/coarse-1)<1e-6
    combined=2*single+2*potential
    assert 2*single<combined<=4*single*(1+1e-12)
    rows.append(dict(separation_over_a=d,cross_field_integral=potential,direct_vs_potential_relative=err,
      angular_refinement_relative=fine/coarse-1,combined_weak_capture_area_coefficient=combined,
      combined_over_sum_of_isolated_areas=combined/(2*single)))
assert abs(rows[0]['combined_over_sum_of_isolated_areas']-2)<1e-12
result=dict(scope='Static identical ordinary-matter Plummer wells, kappa=chi*|g1+g2|^2, optically thin external isotropic capture. Not a merger simulation, full absorption calculation or lensing fit.',
    area_unit='chi*G^2*M^2/a; each well has mass M and scale a',cases=rows,
    midpoint_result='For nonzero separation, g1=-g2 at the midpoint: combined local opacity is zero while the sum of isolated opacities is positive.')
(OUT/'overlap-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(rows,indent=2))
