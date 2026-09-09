"""Exploratory one-parameter curvature fits to already inspected Pantheon+ data.
Keep kappa, Cepheid calibration and photon energy/time factors fixed.
"""
from pathlib import Path
import sys,contextlib,io,json
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.linalg import cho_solve
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'shared_interaction_test'))
with contextlib.redirect_stdout(io.StringIO()):
    import brightness as b
x=np.log1p(b.z)
def scale(a,sign):
    y=a*x
    if sign==1:return np.sinc(y/np.pi)
    if a==0:return np.ones_like(y)
    return np.sinh(y)/y
def chi(a,sign):
    f=scale(a,sign)
    if np.any(f<=0):return np.inf
    r=b.r-5*np.log10(f)
    return float(r@cho_solve(b.vf,r))
out=[]
for name,sign,upper in [('positive',1,.999*np.pi/max(x)),('negative',-1,10.)]:
    sol=minimize_scalar(lambda a:chi(a,sign),bounds=(0,upper),method='bounded',options={'xatol':1e-10})
    opts=[(0.,chi(0.,sign)),(float(sol.x),float(sol.fun))]
    a,ch=min(opts,key=lambda x:x[1])
    grid=np.linspace(0,upper,201)
    # Grid minimum brackets optimized global behavior for this scalar family.
    gridbest=min((chi(aa,sign),aa) for aa in grid)
    rr=b.r-5*np.log10(scale(a,sign))
    bins=[]
    for w,base in zip(b.binweights,b.b):
        bins.append(dict(z_min=base['z_min'],z_max=base['z_max'],residual_mag=float(w@rr)))
    out.append(dict(curvature=name,a_inverse_kappa_Rc=a,
        radius_Mpc=None if a==0 else 1/(b.kap*a),
        chi_squared=ch,delta_chi_vs_flat=ch-b.s['chi_squared'],
        fit_parameters_added=1,best_is_flat_boundary=a==0,
        scan_upper_a=upper,coarse_grid_best_chi=float(gridbest[0]),
        bins_using_fixed_previous_GLS_weights=bins))
result=dict(status='Exploratory fits to the same already-inspected data, not held-out validation',
    kappa_per_Mly=b.cfg['kappa_per_Mly'],calibrated_M=b.M,N=len(b.hi),
    definition='D_L=(1+zHEL)*Rc*sin(R/Rc) for positive curvature; sinh for negative; R=ln(1+zHD)/kappa',
    flat_chi_squared=b.s['chi_squared'],fits=out,
    limitations=['Constant-curvature isotropic geometry only, no rotation or full static gravity solution.',
      'Positive-curvature fit restricted before the first conjugate point; no multiple wraps/images.',
      'One radius fitted to high-z data, so improvements are exploratory and parameter-dependent.',
      'Full released covariance and fixed Cepheid calibration retained; no new source luminosity adjustment.',
      'Negative-curvature upper search bound a=10 is a numerical search range, not a physical prior.'])
(P/'curvature_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
