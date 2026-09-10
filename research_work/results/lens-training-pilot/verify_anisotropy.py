"""Independent Jeans pressure and whole-aperture virial verification."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid,quad
H=Path(__file__).resolve().parent
x=np.geomspace(1e-7,1e6,100000);j=1/(x*(1+x)**3);g=1/(1+x)**2
results=[]
for beta in [-.3,0.,.3]:
    pressure=x**(-2*beta)*(-cumulative_trapezoid((j*g*x**(2*beta))[::-1],x[::-1],initial=0)[::-1])
    for ap in [.1,.5,1,3,10]:
        u=np.minimum(1,(ap/x)**2);mu=np.sqrt(1-u);W=u/(1+mu)
        Wb=W-beta*W*W*(mu+2)/3
        K=x**(2*beta)*cumulative_trapezoid(x**(2-2*beta)*Wb,x,initial=0)
        first=np.trapezoid(pressure*x*x*Wb,x);second=np.trapezoid(j*g*K,x)
        err=abs(first/second-1);assert err<1e-5
        results.append(dict(beta=beta,aperture_over_a=ap,kernel_relative_difference=err))
    # For a complete spherical aperture, spherical averaging cancels beta.
    global_sigma=np.trapezoid(pressure*x*x*(1-2*beta/3),x)/np.trapezoid(j*x*x,x)
    assert abs(global_sigma/(1/18)-1)<1e-5
    for mu0 in [0.,.2,.9,.999]:
        W=1-mu0
        expected=W-beta*W*W*(mu0+2)/3
        actual=quad(lambda mu:1-beta*(1-mu*mu),mu0,1)[0]
        assert abs(actual-expected)<1e-12
summary={'maximum_kernel_relative_difference':max(r['kernel_relative_difference'] for r in results),
    'angular_projection_identity_verified':True,'whole_aperture_virial_limit_verified':True,'cases':results,
    'scope':'Numerical moment equations only, not a complete positive-distribution-function or stability proof.'}
for beta in [-.3,.3]:
    fine=json.loads((H/f'predictions-refined-updated-profile-beta{beta:+.1f}.json').read_text())
    coarse=json.loads((H/f'predictions-updated-profile-beta{beta:+.1f}.json').read_text())
    assert len(fine)==len(coarse)==264 and all(r['role']=='training' and r['beta']==beta for r in fine)
    key=lambda r:(r['Name'],r['model'],r['seeing_fwhm_arcsec'],r['cutoff_over_Re'])
    assert [key(r) for r in fine]==[key(r) for r in coarse]
    change=max(abs(r['theta_pred_arcsec']-s['theta_pred_arcsec']) for r,s in zip(fine,coarse));assert change<.001
    summary[f'beta{beta:+.1f}_refinement_max_angle_change_arcsec']=change
(H/'anisotropy-verification.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
