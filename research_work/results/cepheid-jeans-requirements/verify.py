"""Verify Jeans signs, inverse requirements, bounds and source integrity."""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
s=json.loads((HERE/'results.json').read_text())
rows=json.loads((HERE/'requirements.json').read_text())
bins=json.loads((HERE.parent/'cepheid-common-frame/training-bins.json').read_text())
for name,digest in s['input_sha256'].items():
    assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest

# Independent finite differences on positive analytic fields:
# density=exp(-R/rd-z^2/(2h^2)), sR2=S*exp(-2R/rs), q=a*z.
cases=0;max_error=0.
for R in [6.,10.,17.]:
    for rd,rs,a in [(4.,27.3,0.),(3.,15.,-200.),(8.,40.,700.)]:
        h=.25;S=600.;p2=240.**2;eps=1e-4
        def nu(r,z):return np.exp(-r/rd-z*z/(2*h*h))
        def radial(r):return S*np.exp(-2*r/rs)
        dR=(nu(R+eps,0)*radial(R+eps)-nu(R-eps,0)*radial(R-eps))/(2*eps)
        dZ=(nu(R,eps)*a*eps-nu(R,-eps)*a*(-eps))/(2*eps)
        numeric=p2-radial(R)-R/nu(R,0)*(dR+dZ)
        analytic=p2-radial(R)*(1-R/rd-2*R/rs)-R*a
        err=abs(numeric-analytic);max_error=max(max_error,err)
        assert err<.002
        cases+=1

max_inverse_error=0.;bound_checks=0
for row in rows:
    b=bins[row['bin']];R=b['R_mean_kpc'];sr=b['corrected_radial_second_moment'];sp=b['corrected_azimuthal_second_moment']
    derivative=row['required_dln_density_dR_if_Rsigma_27_3_per_kpc']-2/27.3
    reconstructed=sp-sr*(1+R*derivative)
    err=abs(reconstructed-row['model_vc_kms']**2);max_inverse_error=max(max_inverse_error,err)
    assert err<1e-8
    # Any nonpositive logarithmic pressure slope makes Vc^2 >= sp-sr at T=0.
    slopes=-np.geomspace(1e-8,100,201)
    assert np.all(sp-sr*(1+slopes)>=sp-sr)
    bound_checks+=len(slopes)
tilt=json.loads((HERE/'tilt-training.json').read_text())
assert sum(r['n'] for r in tilt)==s['training_rows']==542
assert all(r['above']+r['below']==r['n'] for r in tilt)
out=dict(analytic_finite_difference_cases=cases,max_finite_difference_error_kms2=max_error,
    max_inverse_reconstruction_error_kms2=max_inverse_error,bound_grid_checks=bound_checks,
    input_hashes_unchanged=True,tilt_counts_verified=True,
    qualification='Checks the algebra and implementation, not equilibrium, selection, or the astrophysical completeness of the error model.')
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out,indent=2))
