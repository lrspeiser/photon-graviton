"""Measured M87 light-shape sensitivity; spherical/bolometric mapping is assumed."""
from pathlib import Path
import json, hashlib
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre, erfcx, gammaincc
from scipy.optimize import brentq

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
profile_path=OUT/'profile.json'
p=json.loads(profile_path.read_text())
data=np.array(p['rows'])
I=10**data[:,0]; angular_sigma=10**data[:,1]; q=data[:,2]
catalog=ROOT/'companion_causal_test/data/themis.dat'
line=next(x for x in catalog.read_text().splitlines() if x[:23].strip()=='NGC4486')
distance=float(line[49:61]); luminosity=float(line[100:109])
# Preserve each Gaussian's projected luminosity and equal-area radial scale.
weights=I*angular_sigma**2*q; weights/=sum(weights)
sigma=distance*angular_sigma*np.sqrt(q)*np.pi/(180*3600)
a=.0002488993265191759; k=10.; radius=1.

def internal_fraction(r,n):
    mu,w=roots_legendre(n)
    path=-r*mu+np.sqrt(radius**2-r*r+r*r*mu*mu)
    # Stable equivalent to 1-P-C for small conversion rate.
    converted=-np.expm1(-a*path)
    companion=a*np.exp(-a*path)*(-np.expm1(-(k-a)*path))/(k-a)
    return float(w@(converted-companion)/2)

def extended_fraction(n):
    # Integrate each Maxwell radial distribution to 12 sigma, not a fitted cutoff.
    totals=[]
    for s in sigma:
        assert 12*s<radius
        totals.append(quad(lambda x: np.sqrt(2/np.pi)*x*x*np.exp(-x*x/2)*internal_fraction(s*x,n),0,12,epsabs=1e-13,epsrel=1e-9)[0])
    return float(weights@totals)

prefactor=luminosity*k*a/(4*np.pi*(k-a))
analytic=prefactor*np.sum(weights*(erfcx(a*sigma/np.sqrt(2))-erfcx(k*sigma/np.sqrt(2)))/sigma**2)
# Independent radial convolution at the origin, using dimensionless source radius.
numeric=prefactor*sum(w*np.sqrt(2/np.pi)/s**2*quad(lambda x: np.exp(-x*x/2-a*s*x)*(-np.expm1(-(k-a)*s*x)),0,12,epsabs=1e-14)[0] for w,s in zip(weights,sigma))
# K(s)<=k*a/(4*pi*s). Gaussian convolution is erf(r/(sqrt(2)*sigma))/r,
# bounded by sqrt(2/pi)/sigma. Hence every projected ray has finite upper bound.
volume_bound=luminosity*k*a/(4*np.pi)*np.sum(weights*np.sqrt(2/np.pi)/sigma)
coarse=extended_fraction(64);fine=extended_fraction(128);point=internal_fraction(0,128)
assert abs(analytic/numeric-1)<1e-9
assert 0<analytic<=volume_bound
assert abs(coarse-fine)<1e-12
assert 0<fine<1
half=brentq(lambda r: float(weights@(1-np.exp(-r*r/(2*sigma*sigma))))-.5,1e-10,1)
result=dict(scope='M87-only finite emitter sensitivity; circularized spherical deprojection, constant bolometric light-shape mapping, unchanged hypothetical transport. Not a lensing fit or full cluster map.',
 profile_sha256=hashlib.sha256(profile_path.read_bytes()).hexdigest(),catalog_sha256=hashlib.sha256(catalog.read_bytes()).hexdigest(),distance_Mpc=distance,luminosity_Lsun=luminosity,
 gaussian_weights=weights.tolist(),gaussian_sigma_Mpc=sigma.tolist(),circularized_projected_half_light_radius_kpc=half*1000,
 radial_tail_omitted_fraction=float(gammaincc(1.5,72)),center_volume_power_Lsun_Mpc3=float(analytic),center_independent_relative_difference=float(abs(analytic/numeric-1)),
 projected_power_global_upper_bound_Lsun_Mpc2=float(2*radius*volume_bound),
 point_emitter_power_Lsun=luminosity*point,extended_emitter_power_Lsun=luminosity*fine,
 relative_total_power_change=fine/point-1,angular_refinement_absolute_fraction=abs(coarse-fine),
 all_six_objectives='open')
(OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
