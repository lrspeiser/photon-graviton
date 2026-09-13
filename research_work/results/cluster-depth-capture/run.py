"""Depth-to-fourth opacity and central illumination curvature in Plummer well."""
from pathlib import Path
import json
from functools import lru_cache
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre,gamma,gammainc
from scipy.optimize import brentq
OUT=Path(__file__).resolve().parent
nodes=lru_cache(None)(roots_legendre)

def column(r,mu,A):
    z=r*mu;B=1+r*r*(1-mu*mu)
    return A*(z/(2*B*(B+z*z))+np.arctan(z/np.sqrt(B))/(2*B**1.5)+np.pi/(4*B**1.5))

def profile(r,A,n=256):
    mu,w=nodes(n)
    return A/(1+r*r)**2*np.dot(w,np.exp(-column(r,mu,A)))/2

threshold=brentq(lambda A:-2+np.pi*A/4+A*A/6,0,10)
rows=[]
for A in [.1,1.,threshold,3.,10.]:
    T=np.pi*A/2;q=1.5
    sigma=np.pi*(T**(1/q)*gamma(1-1/q)*gammainc(1-1/q,T)+np.expm1(-T))
    total=4*np.pi*quad(lambda r:r*r*profile(r,A),0,np.inf,epsabs=1e-8,epsrel=1e-8)[0]
    assert abs(total/sigma-1)<2e-6
    center=A*np.exp(-np.pi*A/4)
    coefficient=-2+np.pi*A/4+A*A/6
    step=1e-4
    numeric=(profile(step,A)/center-1)/step**2
    assert abs(numeric-coefficient)<2e-5
    median=brentq(lambda r:4*np.pi*quad(lambda s:s*s*profile(s,A),0,r,epsabs=1e-8)[0]/sigma-.5,.001,1000)
    rows.append(dict(A=A,effective_area_over_a2=sigma,area_volume_relative_error=total/sigma-1,
       central_deposition_density=center,relative_central_r2_coefficient=coefficient,numeric_curvature=numeric,
       half_deposition_radius_over_a=median))
midpoints=[]
for d in [.5,1.,2.,5.]:
    phi_single=1/np.sqrt(1+(d/2)**2)
    midpoints.append(dict(separation_over_a=d,depth_based_midpoint_opacity_at_unit_coefficient=(2*phi_single)**4,
         combined_over_sum_of_individual_midpoint_opacities=8.,acceleration_squared_midpoint_opacity=0.))
result=dict(scope='Optional isolated-well escape-energy-to-fourth capture, external bath and fixed Plummer potential. No fitted constants, self-consistent storage or observed data.',
    rule='kappa=chi*(-Phi)^4, Phi(infinity)=0; A=chi*(GM)^4/a^3; a*kappa=A/(1+(r/a)^2)^2',
    central_curvature_threshold_A=threshold,cases=rows,midpoint_comparisons=midpoints)
(OUT/'results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(result,indent=2))
