"""Exact effective capture area and energy requirement of the fitted model."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma,gammainc
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
fit=json.loads((HERE/'results.json').read_text())['models']['attenuated']
curves=[r for r in json.loads((HERE/'predictions.json').read_text()) if r['model']=='attenuated']
meta={}
for line in (ROOT/'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt').read_text().splitlines():
    f=line.split()
    if len(f)!=19:continue
    try:meta[f[0]]=float(f[11])
    except ValueError:pass
C=fit['C_Msun_kpc3'];k0=fit['k0_per_kpc'];s=fit['scale_to_disk']
Msun=1.98847e30;kpc=3.085677581491367e19;c=299792458.
surface=C/k0
fluence=surface*Msun/kpc**2*c*c
mu,w=leggauss(128);rows=[]
def sig(a):
    T=np.pi*k0*a/2
    return np.pi*a*a*(T**(2/3)*gamma(1/3)*gammainc(1/3,T)+np.expm1(-T))
for d in curves:
    a=s*meta[d['galaxy']];T=np.pi*k0*a/2;area=sig(a)
    # Known impact-plane integration; compactified integral via y=1/(1+x^2).
    independent=np.pi*a*a*quad(lambda y:-np.expm1(-T*y**1.5)/y**2,0,1,epsabs=1e-8,epsrel=1e-9)[0]
    assert abs(independent/area-1)<1e-7
    rmax=max(d['R_kpc'])
    def density(r):
        x=r/a;t=x*mu;B2=1+x*x*(1-mu*mu);B=np.sqrt(B2)
        tau=k0*a*(t/(2*B2*(B2+t*t))+(np.arctan(t/B)+np.pi/2)/(2*B**3))
        J=.5*np.sum(np.exp(-np.maximum(tau,0))*w)
        return C*J/(1+x*x)**2
    within=4*np.pi*quad(lambda r:density(r)*r*r,0,rmax,epsabs=1e-2,epsrel=1e-8)[0]
    total=surface*area
    assert 0<within<total
    rows.append(dict(galaxy=d['galaxy'],split=d['split'],capture_scale_kpc=a,chord_center_optical_depth=T,effective_capture_area_kpc2=area,predicted_total_deposit_Msun=total,last_rotation_radius_kpc=rmax,predicted_deposit_inside_last_radius_Msun=within,fraction_mass_outside_measured_rotation=1-within/total,cross_section_relative_check=abs(independent/area-1)))
out=dict(required_mass_exposure_Msun_pc2=surface/1e6,required_integral_c_u_dt_J_m2=fluence,example_years_at_u_1e_minus14_J_m3=fluence/(c*1e-14)/(365.25*86400),example_status='Scaling example only; no assumed bath density or universe age',median_predicted_total_mass_Msun=float(np.median([r['predicted_total_deposit_Msun'] for r in rows])),median_fraction_mass_outside_rotation=float(np.median([r['fraction_mass_outside_measured_rotation'] for r in rows])),max_cross_section_relative_check=max(r['cross_section_relative_check'] for r in rows),rows=rows)
(HERE/'energy-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))
