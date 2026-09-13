"""Finite effective area for smooth capture profiles without a hard boundary."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma,gammainc

OUT=Path(__file__).resolve().parent

def chord(b,p,k):
    return k*np.sqrt(np.pi)*gamma((p-1)/2)/gamma(p/2)*(1+b*b)**((1-p)/2)

def area(B,p,k):
    return 2*np.pi*quad(lambda b:b*(-np.expm1(-chord(b,p,k))),0,B,epsabs=1e-9,epsrel=1e-9,limit=300)[0]

rows=[]
for p in [2.,3.,4.,6.]:
    for k in [.1,1.,10.]:
        T=chord(0,p,k)
        # Independent line integration checks the analytic optical depth.
        errs=[]
        for b in [0.,1.,10.]:
            integral=2*quad(lambda z:k*(1+b*b+z*z)**(-p/2),0,np.inf,epsabs=1e-11,epsrel=1e-11)[0]
            errs.append(abs(integral/chord(b,p,k)-1))
        assert max(errs)<1e-9
        row=dict(tail_power=p,kappa0_times_scale=k,central_optical_depth=T,area_over_scale2_with_outer_impact_limit={str(B):area(B,p,k) for B in [10.,100.,1000.]},max_chord_relative_error=max(errs))
        if p>3:
            q=(p-1)/2
            analytic=np.pi*(T**(1/q)*gamma(1-1/q)*gammainc(1-1/q,T)+np.expm1(-T))
            numeric=2*np.pi*quad(lambda b:b*(-np.expm1(-chord(b,p,k))),0,np.inf,epsabs=1e-8,epsrel=1e-9,limit=300)[0]
            assert analytic>0 and abs(analytic/numeric-1)<1e-8
            row.update(infinite_area_over_scale2=analytic,effective_radius_over_scale=np.sqrt(analytic/np.pi),independent_area_relative_error=abs(analytic/numeric-1))
        else:row['infinite_area_over_scale2']='diverges'
        rows.append(row)
result=dict(scope='Smooth spherical opacity comparison, straight rays and fixed incoming companion bath. No microphysical gravity-capture relation or observed cluster fit.',
    profile='kappa(r)=kappa0*(1+(r/a)^2)^(-p/2)',cases=rows)
(OUT/'results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(rows,indent=2))
