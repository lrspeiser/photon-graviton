"""Directional consistency of the stipulated universal-clock lapse model."""
from pathlib import Path
import json
import numpy as np
from scipy.special import roots_legendre
from scipy.optimize import brentq

OUT=Path(__file__).resolve().parent
mu,w=roots_legendre(256);w=w/2
rows=[]
for target in [.001,.01,.1,1.]:
    # Constant-gradient ray-integral benchmark: ln(1+z)=-x*mu.
    x=brentq(lambda x:np.sinh(x)/x-1-target,1e-5,10.)
    ell=-x*mu;z=np.expm1(ell)
    mean=float(w@z);meanlog=float(w@ell);variance=float(w@(ell*ell))
    assert abs(mean-target)<1e-11 and abs(meanlog)<1e-14
    assert abs(variance/(x*x/3)-1)<1e-12
    assert abs(float(w[z<0].sum())-.5)<1e-14
    assert np.max(abs(np.exp(ell)*np.exp(ell[::-1])-1))<1e-14
    rows.append(dict(target_mean_redshift=target,gradient_path_product=x,
        mean_log_wavelength_ratio=meanlog,rms_log_wavelength_ratio=variance**.5,
        blueshift_direction_fraction=.5,min_redshift=float(np.expm1(-x)),max_redshift=float(np.expm1(x))))

result=dict(scope='Conditional local directional test for fixed spatial rulers, universal clocks d_tau=dt/n and local photon speed c. Angular benchmark fixes straight-ray integrated gradients; not a self-consistent global ray/field solution or observed sample.',
    provenance='Known lapse kinematics and calculus extended from this project ray-clock identity; no novelty claim for formulas.',
    identities=['ln(S)=-integral ray_direction dot grad(ln n) ds',
    'alpha_local(k_hat)=-k_hat dot grad(ln n); alpha_local(-k_hat)=-alpha_local(k_hat)',
    'isotropic angular mean(alpha_local)=0; angular variance(alpha_local)=|grad(ln n)|^2/3',
    'For prescribed log-shift -x*mu: mean(z)=sinh(x)/x-1 but mean(log(1+z))=0 and half directions blueshift'],
    interpretation='A positive direction-independent local fractional redshift rate cannot arise solely from this lapse/clock rule. Does not prove zero redshift along finite rays in evolving inhomogeneous fields.',
    angular_benchmarks=rows)
(OUT/'directional-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(rows,indent=2))
