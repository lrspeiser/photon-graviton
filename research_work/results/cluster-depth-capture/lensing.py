"""Project the optional depth capture profile; no observed data are fitted."""
from pathlib import Path
import ast
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

OUT = Path(__file__).resolve().parent
path = OUT / 'run.py'
body = []
for node in ast.parse(path.read_text()).body:
    if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'threshold' for t in node.targets):
        break
    body.append(node)
ns = {'__file__': str(path)}
exec(compile(ast.Module(body=body, type_ignores=[]), str(path), 'exec'), ns)

def central_coefficient(A, n=256):
    q0 = A*np.exp(-np.pi*A/4)
    C = -2+np.pi*A/4+A*A/6
    def integrand(r):
        if r < 1e-4:
            return q0*C
        return (ns['profile'](r, A, n)-q0)/r**2
    return quad(integrand, 0, np.inf, epsabs=1e-10, epsrel=1e-9)[0]

transition = brentq(central_coefficient, 1., 10., xtol=1e-9)
transition128 = brentq(lambda A:central_coefficient(A,128), 1., 10., xtol=1e-9)
assert abs(transition-transition128) < 1e-7
# Independent optically thin limit: q=A/(1+r^2)^2 and its analytic projection.
weak_errors = []
for b in [0., .1, 1., 3., 10.]:
    tiny = 1e-6
    measured = 2*quad(lambda z:ns['profile'](np.hypot(b,z),tiny)/tiny,0,np.inf,epsabs=1e-10)[0]
    analytic = np.pi/(2*(1+b*b)**1.5)
    weak_errors.append(abs(measured/analytic-1))
assert max(weak_errors) < 2e-6
cases = json.loads((OUT/'results.json').read_text())['cases']
rows = []
for case in cases:
    A = case['A']
    sigma = case['effective_area_over_a2']
    q = lambda r: ns['profile'](r, A)
    surf = lambda b: 2*quad(lambda z:q(np.hypot(b,z)), 0, np.inf, epsabs=1e-10, epsrel=1e-9)[0]
    S0 = surf(0)
    B = central_coefficient(A)
    B128 = central_coefficient(A, 128)
    step = .001
    numeric = (surf(step)-S0)/step**2
    assert abs(numeric-B)/max(abs(B), S0) < 1e-5
    assert abs(B128-B)/max(abs(B), S0) < 1e-7
    samples = []
    for b in [.1, 1., 3., 10.]:
        S = surf(b)
        M3 = 4*np.pi*quad(lambda r:r*r*q(r),0,b,epsabs=1e-10,epsrel=1e-9)[0]
        caps = 4*np.pi*quad(lambda r:q(r)*b*b/(1+np.sqrt(max(0.,1-b*b/r/r))),b,np.inf,epsabs=1e-10,epsrel=1e-9)[0]
        M2 = M3+caps
        delta = M2/(np.pi*b*b)-S
        independent = 2*quad(lambda t:t*surf(t),0,b,epsabs=1e-10,epsrel=1e-8)[0]/b**2-S
        discrepancy = abs(delta-independent)/max(S0,S)
        assert discrepancy < 1e-7
        samples.append(dict(radius_over_a=b,normalized_surface_density=S/sigma,
            normalized_excess_surface_density=delta/sigma,
            enclosed_3d_fraction=M3/sigma,projected_fraction=M2/sigma,
            deflection_over_spherical_circular_speed_squared_over_c_squared=4*M2/M3,
            aperture_verification_scaled_error=discrepancy))
    rows.append(dict(A=A,central_surface_density=S0,central_surface_curvature=B,
        central_curvature_finite_difference_scaled_error=abs(numeric-B)/max(abs(B),S0),
        angular_refinement_scaled_error=abs(B128-B)/max(abs(B),S0),
        central_total_shear_sign_change_deposit_over_baryon_mass=2*sigma/(np.pi*B) if B>0 else None,
        samples=samples))
result = dict(scope='Synthetic isolated fixed Plummer well; stationary stored density proportional to injection, equal gravitational potentials and negligible stress assumed. No observed lens fit, support or backreaction.',
    provenance='Known spherical projection and weak-lensing mathematics applied to optional kappa=chi*(-Phi)^4 postulate.',
    formula='B=int_0^infinity [q(r)-q(0)]/r^2 dr; Sigma(b)=Sigma(0)+B*b^2+o(b^2); DeltaSigma=-B*b^2/2+o(b^2)',
    projected_central_curvature_transition_A=transition,
    transition_angular_refinement_absolute_error=abs(transition-transition128),
    optically_thin_analytic_projection_max_relative_error=max(weak_errors),cases=rows)
(OUT/'lensing-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(result,indent=2))
