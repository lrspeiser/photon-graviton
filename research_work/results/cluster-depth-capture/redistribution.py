"""Positive radial-orbit mixture from at-rest capture; fixed-well diagnostic."""
from pathlib import Path
import ast,json
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.special import roots_legendre
from scipy.optimize import brentq

OUT=Path(__file__).resolve().parent
body=[]
for node in ast.parse((OUT/'run.py').read_text()).body:
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='threshold' for t in node.targets):break
    body.append(node)
ns={'__file__':str(OUT/'run.py')}
exec(compile(ast.Module(body=body,type_ignores=[]),'capture definitions','exec'),ns)

def fall_time(a,theta=np.pi/2,n=64):
    nodes,w=ns['nodes'](n);t=(nodes+1)*theta/2
    # s=a*sin(t) cancels the turning-point square-root singularity exactly.
    h=np.sqrt(1+(a*np.sin(t))**2);H=np.hypot(1,a)
    return float(np.dot(w,np.sqrt(h*H*(h+H)/2))*theta/2)

def fraction_inside(r,a,n=64):
    if r>=a:return 1.
    return fall_time(a,np.arcsin(r/a),n)/fall_time(a,n=n)

prior=json.loads((OUT/'results.json').read_text())['cases']
rows=[]
for item in prior:
    A=item['A']
    if A not in [.1,1.,3.,10.]:continue
    sigma=item['effective_area_over_a2']
    weight=lambda a:4*np.pi*a*a*ns['profile'](a,A)/sigma
    def born(r):return quad(weight,0,r,epsabs=1e-9,epsrel=1e-8)[0]
    def mixed(r,n=64):
        return born(r)+quad(lambda a:weight(a)*fraction_inside(r,a,n),r,np.inf,epsabs=1e-9,epsrel=1e-8)[0]
    total=quad(weight,0,np.inf,epsabs=1e-9)[0]
    assert abs(total-1)<1e-7
    checks=[]
    for r in [.001,.1,1.,3.,10.]:
        frozen=born(r);value=mixed(r);fine=mixed(r,128)
        assert frozen<=value+1e-10 and 0<=value<=1+1e-8
        assert abs(value-fine)<1e-7
        checks.append(dict(radius_over_a=r,birth_mass_fraction=frozen,phase_mixed_mass_fraction=fine,
            quadrature_absolute_difference=abs(value-fine)))
    median=brentq(lambda r:mixed(r)-.5,.001,item['half_deposition_radius_over_a'])
    def cusp_integrand(a):
        H=np.hypot(1,a)
        speed0=np.sqrt(2)*a/np.sqrt(H*(H+1))
        return weight(a)/(4*np.pi*fall_time(a)*speed0)
    cusp=quad(cusp_integrand,0,np.inf,epsabs=1e-9,epsrel=1e-8)[0]
    small=checks[0]['phase_mixed_mass_fraction']
    error=abs(small/(4*np.pi*cusp*.001)-1)
    assert error<1e-4
    rows.append(dict(A=A,normalization_error=total-1,birth_half_mass_radius=item['half_deposition_radius_over_a'],
        phase_mixed_half_mass_radius=median,half_mass_radius_ratio=median/item['half_deposition_radius_over_a'],
        formal_central_r_squared_density_coefficient=cusp,central_mass_asymptote_relative_error=error,samples=checks))

orbits=[]
for a in [.3,1.,3.,10.]:
    T=fall_time(a)
    def rhs(t,y):return [y[1],-y[0]/(1+y[0]**2)**1.5]
    events=[]
    for fraction in [.9,.5,.25,0.]:
        def event(t,y,fraction=fraction):return y[0]-fraction*a
        event.direction=-1
        events.append(event)
    sol=solve_ivp(rhs,[0,1.05*T],[a,0.],events=events,rtol=1e-11,atol=1e-12,method='DOP853')
    assert sol.success and all(len(e)==1 for e in sol.t_events)
    actualT=sol.t_events[-1][0]
    discrepancies=[]
    for fraction,ev in zip([.9,.5,.25],sol.t_events[:3]):
        measured=(actualT-ev[0])/actualT
        discrepancies.append(abs(measured-fraction_inside(fraction*a,a)))
    energy=.5*sol.y[1]**2-1/np.sqrt(1+sol.y[0]**2)
    drift=max(abs(energy/energy[0]-1))
    assert abs(actualT/T-1)<1e-8 and max(discrepancies)<1e-8 and drift<1e-8
    orbits.append(dict(apocenter_over_a=a,fall_time_in_sqrt_a_cubed_over_GM=T,
        ode_fall_time_relative_error=actualT/T-1,max_cdf_absolute_difference=max(discrepancies),relative_energy_drift=float(drift)))

result=dict(scope='Hypothetical at-rest massive-particle capture followed by radial collisionless motion, random orbit axes and fully mixed phases, in a fixed Plummer potential. No capture energy/momentum mechanism, self-gravity, finite-age mixing or data fit.',
    provenance='Known Newtonian orbit and residence-time equations applied to optional depth-capture birth distribution.',
    formulas=['v(r|apocenter)=sqrt(2*(Psi(r)-Psi(apocenter)))','P(r<R|apocenter)=integral_0^min(R,apocenter) dr/v divided by integral_0^apocenter dr/v',
    'M_mixed(<R)=integral birth_weight(apocenter)*P(r<R|apocenter) dapocenter'],
    cases=rows,independent_orbit_checks=orbits)
(OUT/'redistribution-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(result,indent=2))
