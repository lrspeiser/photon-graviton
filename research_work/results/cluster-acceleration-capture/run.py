"""Candidate kappa proportional to g^2 in a fixed Plummer ordinary-matter well."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre
from scipy.optimize import brentq
OUT=Path(__file__).resolve().parent

def column(r,mu,A):
    z=r*mu;B=1+r*r*(1-mu*mu)
    f2=z/(2*B*(B+z*z))+np.arctan(z/np.sqrt(B))/(2*B**1.5)
    f3=z/(4*B*(B+z*z)**2)+3*f2/(4*B)
    return A*(f2-f3+np.pi/(4*B**1.5)-3*np.pi/(16*B**2.5))

def profile(r,A,n):
    mu,w=roots_legendre(n)
    return A*r*r/(1+r*r)**3*np.dot(w,np.exp(-column(r,mu,A)))/2

def full_column(b,A):
    return A*np.pi*(1+4*b*b)/(8*(1+b*b)**2.5)

rows=[]
for A in [.1,1.,10.,100.]:
    sigma=2*np.pi*quad(lambda b:b*(-np.expm1(-full_column(b,A))),0,np.inf,epsabs=1e-9,epsrel=1e-9)[0]
    totals=[]
    for n in [128,256]:
        total=4*np.pi*quad(lambda r:r*r*profile(r,A,n),0,np.inf,epsabs=1e-8,epsrel=1e-8,limit=250)[0]
        totals.append(total)
    assert abs(totals[-1]/sigma-1)<2e-5
    assert abs(totals[-1]/totals[0]-1)<2e-5
    def cumulative(r):return 4*np.pi*quad(lambda t:t*t*profile(t,A,256),0,r,epsabs=1e-8,epsrel=1e-8)[0]/sigma
    median=brentq(lambda r:cumulative(r)-.5,.01,1000)
    # Independent direct integration of upstream opacity for noncentral locations.
    errors=[]
    for r,mu in [(.5,-.8),(1.,0.),(3.,.9)]:
        b2=r*r*(1-mu*mu);z=r*mu
        num=quad(lambda t:A*(b2+t*t)/(1+b2+t*t)**3,-np.inf,z,epsabs=1e-11)[0]
        errors.append(abs(num-column(r,mu,A)))
    assert max(errors)<1e-8
    rows.append(dict(dimensionless_capture_strength=A,effective_area_over_a2=sigma,effective_radius_over_a=np.sqrt(sigma/np.pi),
        spatial_total_over_a2=totals[-1],area_vs_volume_relative=totals[-1]/sigma-1,angular_refinement_relative=totals[-1]/totals[0]-1,
        maximum_upstream_column_absolute_error=max(errors),half_deposition_radius_over_a=median,
        deposition_fraction_outside_a=1-cumulative(1),deposition_fraction_outside_10a=1-cumulative(10)))
result=dict(scope='Synthetic fixed ordinary-matter Plummer well, hypothetical acceleration-squared opacity, external isotropic companions only. No fit, deposit backreaction or microscopic derivation.',
    definitions=['g(r)=GM r/(r^2+a^2)^(3/2)','kappa(r)=g(r)^2/(ell_star*g_star^2)',
    'A=G^2 M^2/(ell_star*g_star^2*a^3)','P_capture=c*u_companion*a^2*sigma_dimensionless'],cases=rows)
(OUT/'results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(result,indent=2))
