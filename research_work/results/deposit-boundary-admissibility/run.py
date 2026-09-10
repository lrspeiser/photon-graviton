"""Conditional particle support of archived lens sources; no held-out scoring."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
inputs = [RESULTS / s for s in (
    'joint-galaxy-audit/results.json',
    'lens-training-pilot/predictions-refined-updated-profile.json',
    'lens-photometric-audit/training-photometry.json',
    'lensing-data-readiness/conditional-geometry.json')]
pars, predictions, photometry, geometry = [json.loads(f.read_text()) for f in inputs]
p = pars['sparc']['parameters']['p']
A = pars['sparc']['parameters']['A']
astar = pars['sparc']['parameters']['a_star_m_s2'] * 3.085677581491367e19 / 1e6
G = 4.30091727003628e-6
photometry = {r['SDSS']: r for r in photometry}
geometry = {r['Name']: r for r in geometry if r['role'] == 'training'}
selected = [r for r in predictions if r['model'] == 'empirical_companion'
            and r['seeing_fwhm_arcsec'] == 1.5]
assert len(selected) == 99 and all(r['role'] == 'training' for r in selected)


def components(x, c):
    """Dimensionless g/B; density is 4*pi*G*a*rho/B."""
    gc = c / (1+x)**(2*p)
    b = 2/x - 2*p/(1+x)
    rho = gc*b
    logrho_prime = -2*p/(1+x) + (-2/x**2 + 2*p/(1+x)**2)/b
    return 1/(1+x)**2 + gc, rho, logrho_prime


def relative_potential(x, xt, c):
    # Potential above its boundary value in units a*B, with stable near-edge subtraction.
    n = 1-2*p
    extra = c*(1+x)**n * np.expm1(n*np.log1p((xt-x)/(1+x)))/n
    return (xt-x)/((1+x)*(1+xt)) + extra


rows = []
max_mass_error = max_potential_error = max_slope_error = 0.
for row in selected:
    name, xt = row['Name'], row['cutoff_over_Re']*1.8153
    a = geometry[name]['conditional_Dl_Mpc']*1000*photometry[name]['Re(I)']*np.pi/(180*3600)/1.8153
    B = G*row['mass_from_sigma_Msun']/a**2
    c = A*astar*(B/astar)**p/B

    def necessary(x):
        g, rho, slope = components(x, c)
        return -2*relative_potential(x, xt, c)*slope/g - 1

    grid = np.geomspace(xt*1e-6, xt*(1-1e-9), 2000)
    values = np.array([necessary(x) for x in grid])
    transitions = np.where(values[:-1]*values[1:] < 0)[0]
    assert len(transitions) == 1 and values[0] > 0 and values[-1] < 0
    i = transitions[0]
    root = brentq(necessary, grid[i], grid[i+1], xtol=1e-12)
    mt = c*xt**2/(1+xt)**(2*p)
    mr = c*root**2/(1+root)**(2*p)
    integrated = quad(lambda x: components(x,c)[1]*x*x, 0, xt, epsabs=1e-10, epsrel=2e-11)[0]
    max_mass_error = max(max_mass_error, abs(integrated/mt-1))
    for fraction in [.1, .5, .9, .999]:
        x = xt*fraction
        exact = relative_potential(x, xt, c)
        numerical = quad(lambda t: components(t,c)[0], x, xt, epsabs=1e-12)[0]
        max_potential_error = max(max_potential_error, abs(numerical/exact-1))
        h = x*1e-5
        slope_fd = (np.log(components(x+h,c)[1])-np.log(components(x-h,c)[1]))/(2*h)
        max_slope_error = max(max_slope_error, abs(slope_fd/components(x,c)[2]-1))
    rows.append({'Name': name, 'role': 'training', 'cutoff_over_Re': row['cutoff_over_Re'],
                 'c_dimensionless': c, 'edge_density_positive': bool(components(xt,c)[1]>0),
                 'necessary_condition_crossing_over_cutoff': root/xt,
                 'equivalent_source_mass_fraction_outside_crossing': 1-mr/mt,
                 'circular_support_speed_at_Re_km_s': float(np.sqrt(a*B*1.8153*components(1.8153,c)[0])),
                 'circular_support_speed_at_edge_km_s': float(np.sqrt(a*B*xt*components(xt,c)[0]))})

# Independent orbit integration checks a singular circular-shell construction.
# This is an existence construction in a fixed spherical potential, not formation/stability.
orbit_checks = []
for cut in [5., 20., 100.]:
    row = next(r for r in rows if r['cutoff_over_Re'] == cut)
    c, xt = row['c_dimensionless'], cut*1.8153
    for fraction in [.1, .9]:
        x0 = xt*fraction
        speed = np.sqrt(x0*components(x0,c)[0])
        period = 2*np.pi*x0/speed
        def rhs(t, y):
            x = np.hypot(y[0], y[1])
            gc = c/(1+x)**(2*p) if x <= xt else c/(1+xt)**(2*p)*(xt/x)**2
            g = 1/(1+x)**2 + gc
            return [y[2], y[3], -g*y[0]/x, -g*y[1]/x]
        sol = solve_ivp(rhs, (0,10*period), [x0,0,0,speed], method='DOP853',
                        rtol=2e-12, atol=2e-13, t_eval=np.linspace(0,10*period,1001))
        assert sol.success
        drift = float(np.max(np.abs(np.hypot(sol.y[0],sol.y[1])/x0-1)))
        orbit_checks.append({'cutoff_over_Re': cut, 'radius_over_cutoff': fraction,
                             'max_fractional_radius_drift_10_orbits': drift})
        assert drift < 1e-8

assert max_mass_error < 1e-9 and max_potential_error < 1e-9 and max_slope_error < 1e-7
summary = {'classification': 'Conditional collisionless particle-source consistency; not new observational validation',
           'p': p, 'training_configurations': len(rows), 'systems': len({r['Name'] for r in rows}),
           'isotropic_sharp_edge_source_admissible': False,
           'circular_shell_equilibrium_exists_conditionally': True,
           'collective_stability_or_capture_formation_established': False,
           'max_density_integral_relative_error': max_mass_error,
           'max_potential_integral_relative_error': max_potential_error,
           'max_density_slope_relative_error': max_slope_error,
           'orbit_checks': orbit_checks,
           'inputs_sha256': {str(f.relative_to(RESULTS)):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
           'cutoff_summaries': []}
for cut in [5.,20.,100.]:
    subset = [r for r in rows if r['cutoff_over_Re']==cut]
    summary['cutoff_summaries'].append({'cutoff_over_Re':cut, **{
        key: {'min': min(r[key] for r in subset), 'median': float(np.median([r[key] for r in subset])),
              'max': max(r[key] for r in subset)} for key in (
            'necessary_condition_crossing_over_cutoff', 'equivalent_source_mass_fraction_outside_crossing',
            'circular_support_speed_at_Re_km_s', 'circular_support_speed_at_edge_km_s')}})
(HERE/'configurations.json').write_text(json.dumps(rows,indent=2)+'\n')
(HERE/'results.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
