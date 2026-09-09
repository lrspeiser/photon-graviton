from pathlib import Path
import os
import json,shutil
import numpy as np
from scipy.integrate import quad
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',Path(__file__).resolve().parents[3]/'research_work/generated'))/'gravity-response'
OUT.mkdir(parents=True,exist_ok=True)
# G_* M_source=1; geometric length units; amplitudes kept modest here.
cases=[]
for bb,bd,m in [(0,1,0),(.1,1,0),(.01,100,0),(1,1,1)]:
    amplitude=2*bb*bd
    r=np.array([.1,1.,10.]);step=r*1e-5
    def phi(x):return -(1+amplitude*np.exp(-m*x))/x
    numerical=(phi(r+step)-phi(r-step))/(2*step)
    analytic=(1+amplitude*(1+m*r)*np.exp(-m*r))/r**2
    derivative_error=float(np.max(abs(numerical/analytic-1)))
    assert derivative_error<1e-8
    lens=[]
    for impact in [.5,2.]:
        # Each gradient is evaluated separately before summing: the scalar
        # terms have opposite signs in the two material-frame potentials.
        def integrand(z):
            rr=np.hypot(impact,z)
            phip=(1+amplitude*(1+m*rr)*np.exp(-m*rr))/rr**2
            psip=(1-amplitude*(1+m*rr)*np.exp(-m*rr))/rr**2
            return (phip+psip)*impact/rr
        value=quad(integrand,-np.inf,np.inf,epsabs=1e-11,epsrel=1e-11)[0]
        assert abs(value/(4/impact)-1)<1e-10
        lens.append(dict(impact=impact,deflection_in_c1_units=value,expected=4/impact))
    local=1+2*bb*bb # massless or short-distance calibration limit
    cases.append(dict(beta_b=bb,beta_d=bd,scalar_mass=m,radii=r.tolist(),dynamics_over_bare_Newton=(analytic*r*r).tolist(),dynamics_over_short_distance_Glab=(analytic*r*r/local).tolist(),lensing_over_Glab=1/local,force_derivative_relative_error=derivative_error,lensing_checks=lens))
pairs=[(.01,100),(.1,10),(1,1)]
rank=[]
for bb,bd in pairs:
    cross=2*bb*bd;ordinary=2*bb*bb;deposit=2*bd*bd
    assert np.isclose(cross*cross,ordinary*deposit,rtol=1e-14)
    rank.append(dict(beta_b=bb,beta_d=bd,cross_excess=cross,ordinary_excess=ordinary,deposit_self_excess=deposit))
result=dict(scope='Linear weak-field canonical conformal scalar comparison; negligible scalar stress backreaction, pressureless sources, static point source; no screening, source production or capture derived',cases=cases,coupling_relation_examples=rank,numerical_tests='Four force-gradient cases, eight light-deflection integrals, three scalar-coupling identities',warning='Large response products from earlier energy budgets cannot be substituted into this linear approximation without checking small fields, scalar energy, stability and calibration')
(OUT/'gravity-response-checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
shutil.copy2(__file__,OUT/Path(__file__).name)
print(json.dumps(dict(force_cases=len(cases),lensing_integrals=8,coupling_identities=3,max_force_derivative_error=max(c['force_derivative_relative_error'] for c in cases)),indent=2))
