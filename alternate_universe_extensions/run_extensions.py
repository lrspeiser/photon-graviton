"""New exploratory toy extensions. FIRAS fit is real derived-data reuse;
dispersion tests and transmission targets are hypothetical diagnostics.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares,brentq
from numpy.polynomial.legendre import leggauss
P=Path(__file__).resolve().parent
a=np.loadtxt(P/'data/firas.txt');c=299792458.;h=6.62607015e-34;kb=1.380649e-23
nu=a[:,0]*100*c;pref=2*h*nu**3/c**2/1e-20
def occ(T,n):
    x=h*nu[:,None]/(kb*T*np.asarray(n)[None,:])
    return np.exp(-x)/(-np.expm1(-x))
ref=pref*occ(2.725,np.array([1.]))[:,0]
obs=ref+a[:,2]/1000;err=a[:,3]/1000;gal=a[:,4]/1000
nodes,weights=leggauss(160)
def spectrum(T,r,ni):
    if np.isinf(r):return pref*occ(T,np.array([1.]))[:,0]
    U=min(r*(1-ni),50.)
    u=(nodes+1)*U/2;w=weights*U/2
    integral=occ(T,1-u/r)@(w*np.exp(-u))
    initial=np.exp(-r*(1-ni))*occ(T,np.array([ni]))[:,0]
    return pref*(integral+initial)
def fit(r,ni=.5):
    def residual(v):return (obs-spectrum(v[0],r,ni)-v[1]*gal)/err
    s=least_squares(residual,[2.725*(1+1/r) if np.isfinite(r) else 2.725,0.],
        bounds=([.1,-100.],[100.,100.]),x_scale=[.01,.1],xtol=1e-11,ftol=1e-11,gtol=1e-8)
    return dict(rate_ratio=None if np.isinf(r) else r,initial_n=ni,T_bath_K=float(s.x[0]),
        galaxy_template_coefficient=float(s.x[1]),chi_squared=float(s.fun@s.fun))
base=fit(np.inf);rows=[fit(r) for r in [10.,30.,100.,300.,1000.,10000.]]
limit=brentq(lambda lr:fit(np.exp(lr))['chi_squared']-base['chi_squared']-3.841458820694124,np.log(10),np.log(1e4))
rmin=float(np.exp(limit));gamma=7.7315e-11
checks=[fit(rmin,ni) for ni in [.1,.5,.9]]
u=.1;ks=np.linspace(0,.2,10001)
gravity=[]
for mass in [0.,.003,.005,.007]:
    om2=mass**2-u*u*ks**2+ks**4
    gravity.append(dict(abs_imaginary_sound_speed=u,m_over_M=mass,
        analytic_min_omega2_over_M2=mass**2-u**4/4,
        grid_min_omega2_over_M2=float(min(om2)),
        minimum_k_over_M=u/np.sqrt(2),
        status='marginal' if np.isclose(mass,u*u/2,rtol=0,atol=1e-14) else ('unstable' if mass<u*u/2 else 'positive in this one-mode toy dispersion')))
test_spectrum=spectrum(2.7353056035,rmin,.5)
nodes,weights=leggauss(320)
integration_error=float(np.max(abs(test_spectrum-spectrum(2.7353056035,rmin,.5)))/max(test_spectrum))
assert integration_error<1e-10
out=dict(status='Exploratory toy extensions, not a completed action or independent data validation',
    thermal_model='dn/dt=gamma; df_k/dt=Gamma(f_eq(k/n,T_bath)-f_k); initial f=f_eq at n_i; n_o=1',
    thermal_parameters='r=Gamma/gamma constant, same rate for all modes in first branch; T_bath and residual Galaxy amplitude fitted',
    blackbody_reference=base,thermal_rate_scan=rows,quadrature_160_vs_320_peak_relative_difference=integration_error,
    conditional_delta_chi_3_841_lower_rate_ratio=rmin,
    memory_sensitivity=checks,
    thermalization_timescale_at_threshold_years=1/(rmin*gamma),
    achromatic_excess_photon_survival_z1=float(np.exp(-rmin*.5)),
    achromatic_excess_photon_log10_survival_z1=-rmin*.5/np.log(10),
    optical_to_microwave_rate_ratio_for_99pct_survival_z1=-np.log(.99)/(rmin*.5),
    transmission_caveat='99% is an illustrative design target, not a measured cosmic transmission. Excess occupation, hence source contrast, damps at Gamma in this model.',
    gravity_dispersion='omega^2=m^2-u^2 k^2+k^4/M^2, positive kinetic normalization assumed; no full gravity constraints derived',
    gravity=gravity,
    limitations=['FIRAS diagonal errors only; temperature is fitted, not predicted.',
       'Constant reservoir temperature assumes an unmodeled energy supply/backreaction balance.',
       'Rate selectivity is an extra hypothesis, not an atomic-protection derivation.',
       'A restoring mass and quartic spatial gradient can stabilize the toy mode, but background compatibility and all gravity modes remain unproved.',
       'No new supernova, angular-distance or CMB anisotropy likelihood.'])
(P/'extension_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
