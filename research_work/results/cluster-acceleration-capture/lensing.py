"""Central shear sign from the calculated hollow deposit profile."""
from pathlib import Path
import json,ast
from functools import lru_cache
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre
OUT=Path(__file__).resolve().parent
path=OUT/'run.py';tree=ast.parse(path.read_text());body=[]
for node in tree.body:
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='rows' for t in node.targets):break
    body.append(node)
ns={'__file__':str(path)};exec(compile(ast.Module(body=body,type_ignores=[]),str(path),'exec'),ns)
ns['roots_legendre']=lru_cache(None)(roots_legendre)
prior=json.loads((OUT/'results.json').read_text())
rows=[]
for case in prior['cases']:
    A=case['dimensionless_capture_strength'];sigma=case['effective_area_over_a2']
    q=lambda r:ns['profile'](r,A,256)
    surf=lambda b:2*quad(lambda z:q(np.hypot(b,z)),0,np.inf,epsabs=1e-10,epsrel=1e-9,limit=250)[0]
    S0=surf(0)
    B=quad(lambda r:q(r)/(r*r),0,np.inf,epsabs=1e-11,epsrel=1e-9,limit=250)[0]
    b=.001
    curvature=(surf(b)-S0)/(b*b)
    assert B>0 and abs(curvature/B-1)<.001
    critical=2*sigma/(np.pi*B)
    samples=[]
    for b in [.1,1.,3.,10.]:
        s=surf(b)
        enclosed3=4*np.pi*quad(lambda r:r*r*q(r),0,b,epsabs=1e-10,epsrel=1e-8)[0]
        caps=4*np.pi*quad(lambda r:q(r)*b*b/(1+np.sqrt(max(0.,1-b*b/r/r))),b,np.inf,epsabs=1e-9,epsrel=1e-8)[0]
        delta=(enclosed3+caps)/(np.pi*b*b)-s
        # Thin-shell aperture geometry and direct image averaging independently agree.
        direct=2*quad(lambda t:t*surf(t),0,b,epsabs=1e-9,epsrel=1e-7)[0]/(b*b)-s
        assert abs(delta-direct)/max(S0,s)<1e-6
        samples.append(dict(radius_over_a=b,normalized_surface_density=s/sigma,normalized_excess_surface_density=delta/sigma,
            deposit_enclosed_mass_fraction=enclosed3/sigma))
    assert samples[0]['normalized_excess_surface_density']<0
    rows.append(dict(A=A,central_surface_power=S0,central_surface_curvature_coefficient=B,
        finite_difference_curvature_relative_error=curvature/B-1,
        deposit_to_baryon_mass_ratio_for_central_total_shear_sign_change=critical,profile=samples))
result=dict(scope='Equal-potential weak-lensing shape diagnostic with stationary density proportional to calculated deposit power. Synthetic fixed Plummer ordinary matter; no observed lens fit or reduced-shear calculation.',
    formula='Sigma_deposit(b)=Sigma0+B*b^2+o(b^2), B=int_0^infinity q(r)/r^2 dr>0; DeltaSigma_deposit=-B*b^2/2+o(b^2)',cases=rows)
(OUT/'lensing-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(rows,indent=2))
