"""Eddington inversion of the optional depth capture density as a tracer."""
from pathlib import Path
import ast,json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import roots_legendre,beta
from scipy.optimize import brentq

OUT=Path(__file__).resolve().parent
body=[]
for node in ast.parse((OUT/'run.py').read_text()).body:
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='threshold' for t in node.targets):break
    body.append(node)
ns={'__file__':str(OUT/'run.py')}
exec(compile(ast.Module(body=body,type_ignores=[]),'capture definitions','exec'),ns)

def inversion(A,N=2048,angles=256,quadrature=256):
    psi=np.linspace(0,1,N+1)
    rho=np.r_[0.,[ns['profile'](np.sqrt(max(0.,p**-2-1)),A,angles)/A for p in psi[1:]]]
    C=-2+np.pi*A/4+A*A/6
    spline=CubicSpline(psi,rho,bc_type=((1,0.),(1,-2*np.exp(-np.pi*A/4)*C)))
    nodes,weights=roots_legendre(quadrature)
    t=(nodes+1)/2;w=weights/2
    def f(energy):
        E=np.asarray(energy)
        return 2*np.sqrt(E)*np.sum(spline(E[...,None]*(1-t*t),2)*w,axis=-1)/(np.sqrt(8)*np.pi**2)
    return f,spline

energies=np.linspace(.01,1,1001)
rows=[]
for A in [.1,1.,1.5,1.8332747453568032,3.,10.]:
    f,rho=inversion(A)
    fine,rhofine=inversion(A,4096,512,512)
    values=f(energies);refined=fine(energies)
    scale=max(abs(refined))
    difference=max(abs(values-refined))/scale
    assert difference<2e-4
    # Reintegrate the recovered, possibly signed DF without clipping negative values.
    nodes,w=roots_legendre(256);t=(nodes+1)/2;w=w/2
    checks=[]
    for psi in [.1,.3,.6,.9,1.]:
        rebuilt=8*np.pi*np.sqrt(2)*psi**1.5*np.dot(w,t*t*fine(psi*(1-t*t)))
        target=float(rhofine(psi));error=abs(rebuilt/target-1)
        assert error<2e-4
        checks.append(dict(relative_potential=psi,reconstruction_relative_error=error))
    rows.append(dict(A=A,minimum_sampled_df=float(min(refined)),energy_at_sampled_minimum=float(energies[np.argmin(refined)]),
        maximum_sampled_df=float(max(refined)),df_at_maximum_binding=float(fine(1.)),
        negative_sample_count=int(np.sum(refined<0)),tested_energy_range=[.01,1.],
        refinement_scaled_max_difference=float(difference),density_reconstruction=checks))

# This is a zero of the high-binding endpoint, not a proof of global positivity below it.
endpoint_zero=brentq(lambda A:float(inversion(A)[0](1.)),.1,3.,xtol=1e-7)
endpoint_zero_fine=brentq(lambda A:float(inversion(A,4096,512,512)[0](1.)),.1,3.,xtol=1e-7)
assert abs(endpoint_zero-endpoint_zero_fine)<2e-4
weak,_=inversion(1e-6,4096,512,512)
analytic=12*beta(3,.5)/(np.sqrt(8)*np.pi**2)*energies**2.5
weak_error=float(max(abs(weak(energies)/analytic-1)))
assert weak_error<1e-3
result=dict(scope='Nonrelativistic spherical isotropic collisionless tracer in fixed Plummer potential; density q/A, Psi=1/sqrt(1+r^2). No self-gravity, formation, stability or observations.',
    provenance='Known Eddington inversion applied to the optional depth-to-fourth capture density. No halo prior or expansion model imported.',
    formula='f(E)=1/(sqrt(8)*pi^2) integral_0^E rho_second(Psi)/sqrt(E-Psi) dPsi; rho_prime(0)=0',
    high_binding_endpoint_zero_A=endpoint_zero_fine,endpoint_zero_refinement_absolute_difference=abs(endpoint_zero-endpoint_zero_fine),
    thin_limit_analytic_df_max_relative_error=weak_error,cases=rows)
(OUT/'support-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(result,indent=2))
