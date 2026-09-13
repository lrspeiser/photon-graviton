"""Conditional Sachs-source inversion of the fitted optical response."""
from pathlib import Path
import hashlib,json
import numpy as np
import sympy as sp

HERE=Path(__file__).resolve().parent
fitpath=HERE/'regular-area-results.json'
energypath=HERE.parent/'isotropic-galaxy-transfer/energy-results.json'
fit=json.loads(fitpath.read_text());energy=json.loads(energypath.read_text())
predictions=json.loads((HERE/'regular-area-predictions.json').read_text())
f,q=sp.symbols('f q',positive=True)
H=-(1-f)*sp.log(1-f)*sp.sqrt(1+f/(1+q*f))
curvature=sp.lambdify((f,q),-sp.diff(H,f,2)/H,'numpy')
series=sp.series(H,f,0,4).removeO().expand()
assert sp.simplify(series.coeff(f,2))==0
assert sp.simplify(-6*series.coeff(f,3)-(3*q+sp.Rational(13,4)))==0
alpha=fit['alpha_per_mpc'];qfit=fit['q']
c=299792458.;G=6.67430e-11;Mpc=3.085677581491367e22;year=365.25*86400
def evaluate(z):
    z=float(z);frac=z/(1+z)
    r=alpha**2*(3*qfit+13/4 if z==0 else float(curvature(frac,qfit)))
    rSI=r/Mpc**2
    effective=c*c*rSI/(4*np.pi*G*(1+z)**2)
    u=.75*c*c*effective
    return dict(z=z,R_opt_per_Mpc2=r,rho_plus_p_over_c2_kg_m3=effective,
                all_companion_u_J_m3=u)
rows=[dict(CID=r['CID'],role=r['role'],**evaluate(r['zHEL'])) for r in predictions]
representative=[evaluate(z) for z in [0,.01,.1,.3,.6,1,1.5,2,3]]
observer=representative[0]
assert abs(observer['R_opt_per_Mpc2']/fit['observer_optical_focusing_per_Mpc2']-1)<1e-12
duration=energy['required_integral_c_u_dt_J_m2']/(c*observer['all_companion_u_J_m3'])/year
out=dict(status='Conditional source inferred from fitted optics; not measured energy or a global field solution',
         input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [fitpath,energypath,HERE/'regular-area-predictions.json']},
         alpha_per_Mpc=alpha,q=qfit,n_observation_rows=len(rows),
         minimum_sample_R_opt_per_Mpc2=min(r['R_opt_per_Mpc2'] for r in rows),
         maximum_sample_R_opt_per_Mpc2=max(r['R_opt_per_Mpc2'] for r in rows),
         negative_sample_source_count=sum(r['R_opt_per_Mpc2']<0 for r in rows),
         observer=observer,constant_bath_equivalent_duration_years=duration,
         required_capture_fluence_J_m2=energy['required_integral_c_u_dt_J_m2'],
         representative_redshifts=representative,rows=rows)
(HERE/'optical-source-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k not in ['rows','input_sha256']},indent=2))
