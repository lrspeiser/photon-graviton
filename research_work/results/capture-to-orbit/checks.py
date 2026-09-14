"""Verification checks for the capture-to-orbit test. No data files needed.

    python checks.py [--output PATH]

Exits nonzero if any check fails. Checks cover exact kinematics, the uniform
velocity-ball injection law, the line-spectrum formulas, orbit averaging,
sampling invariance and the Plummer universal-K consequence. They verify the
stated calculations, not the physical interaction or any observation.
"""
import argparse
import json
import os
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import kinematics as kin        # noqa: E402
import supported_profile as sp  # noqa: E402

G = sp.G


def event_balance(rng, n=20000):
    M = 1e6
    u = kin.isotropic(n, rng)*rng.uniform(0, .3, n)[:, None]
    E = rng.uniform(.8, 5., n)
    ev = kin.react(E, kin.isotropic(n, rng), u, M, kin.isotropic(n, rng))
    bad = np.abs(ev['k'] + ev['P'] - ev['X'] - ev['R'])[ev['ok']]
    masses = ev['X'][:, 0]**2 - np.sum(ev['X'][:, 1:]**2, axis=1)
    return dict(max_four_momentum_imbalance_over_M=float(bad.max()/M),
                max_product_mass_error=float(np.max(np.abs(masses[ev['ok']] - 1))),
                passed=bool(bad.max()/M < 1e-12 and np.max(np.abs(masses[ev['ok']] - 1)) < 1e-6))


def drag():
    gray, gray_imb = kin.drag_quadrature(.003, (100., 200.), M=1e7)
    th = kin.threshold(1e7)
    thr, thr_imb = kin.drag_quadrature(3e-4, (th, 2*th), M=1e7)
    pred = kin.drag_coefficient(kin.REFERENCE_CHI)
    return dict(gray_coefficient=gray, gray_expected=4/3, threshold_coefficient=thr, threshold_first_order=pred,
                note='threshold value approaches the first-order result as sqrt(beta)',
                passed=bool(abs(gray - 4/3) < 1e-3 and abs(thr - pred) < 5e-3 and max(gray_imb, thr_imb) < 1e-12))


def window_events(rng, beta, E_lo, E_hi, n, M):
    u = np.zeros((n, 3))
    u[:, 2] = beta
    n_c = kin.isotropic(n, rng)
    ev = kin.react(rng.uniform(E_lo, E_hi, n), n_c, u, M, kin.isotropic(n, rng))
    return ev, ev['weight']*(1 - beta*n_c[:, 2])


def uniform_ball(rng, beta=.01, v_esc=.03, n=2_000_000, M=1e6):
    """Smooth bath spectrum flat on [0.5, 2] E_th: bound products fill the galaxy-frame ball."""
    th = kin.threshold(M)
    g = 1/np.sqrt(1 - beta*beta)
    u_max = v_esc + beta
    x_max = 1/np.sqrt(1 - u_max**2) - 1
    lo, hi = th/(g*(1 + beta)), th*(1 + x_max)/(g*(1 - beta))
    ev, w = window_events(rng, beta, lo, hi, n, M)
    v = ev['X'][:, 1:]/ev['X'][:, :1]
    speed = np.linalg.norm(v, axis=1)
    bound = ev['ok'] & (speed < v_esc)
    evt, wt = window_events(rng, beta, lo, 2*th, n, M)
    rate_bound = (hi - lo)*np.mean(w*bound)
    rate_total = (2*th - lo)*np.mean(wt)
    frac = rate_bound/rate_total
    pred = v_esc**3/3/kin.REFERENCE_RATE
    cube = np.sort((speed[bound]/v_esc)**3)
    wb = w[bound][np.argsort((speed[bound]/v_esc)**3)]
    cdf = np.cumsum(wb)/wb.sum()
    ks = float(np.max(np.abs(cdf - cube)))
    mean_v = np.sum(w[bound, None]*v[bound], axis=0)/w[bound].sum()
    return dict(bound_fraction=float(frac), predicted=float(pred), relative_error=float(frac/pred - 1),
                bound_events=int(bound.sum()), ks_uniform_ball=ks, mean_velocity_over_vesc=(mean_v/v_esc).tolist(),
                receiver_velocity_over_vesc=beta/v_esc,
                passed=bool(abs(frac/pred - 1) < .02 and ks < .01 and np.linalg.norm(mean_v) < .02*v_esc))


def line(rng, n=2_000_000, M=1e6):
    beta, v_esc = .003, .02
    out = {}
    for label, delta in [('inside_band', .001), ('edge_tuned', -beta + 5e-5)]:
        u = np.zeros((n, 3))
        u[:, 2] = beta
        n_c = kin.isotropic(n, rng)
        ev = kin.react(np.full(n, kin.threshold(M)*(1 + delta)), n_c, u, M, kin.isotropic(n, rng))
        w = ev['weight']*(1 - beta*n_c[:, 2])
        v = ev['X'][:, 1:]/ev['X'][:, :1]
        speed = np.linalg.norm(v, axis=1)
        react = ev['ok']
        frac = float(np.sum(w*(react & (speed < v_esc)))/np.sum(w*react))
        rel = np.linalg.norm(v[react] - u[react], axis=1)
        out[label] = dict(delta=delta, reacting_direction_fraction=float(react.mean()), bound_fraction=frac,
                          mean_vz_over_receiver=float(np.sum(w[react]*v[react, 2])/np.sum(w[react])/beta),
                          max_speed_relative_to_receiver=float(rel.max()))
    pred = float(kin.line_bound_fraction(.001, beta, v_esc))
    edge = out['edge_tuned']
    out.update(predicted_inside_band=pred, speed_window_for_edge_tuning_kms=(v_esc - beta)**2/2*kin.C_KMS,
               passed=bool(abs(out['inside_band']['bound_fraction']/pred - 1) < .03 and edge['bound_fraction'] > .999
                           and abs(edge['mean_vz_over_receiver'] - 1) < .02
                           and edge['max_speed_relative_to_receiver'] < 1.02*np.sqrt(2*5e-5)))
    return out


def kepler():
    M, a, e = 1e10, 10., .6
    r = np.geomspace(1e-4, 1e5, 20000)
    pot = sp.SphericalPotential(r, np.full(r.size, M))
    E, L = np.array([-G*M/(2*a)]), np.array([np.sqrt(G*M*a*(1 - e*e))])
    rp, ra = sp.apsides(pot, E, L, np.array([a]))
    reval = np.array([4.5, 6., 8., 10., 12., 15., 15.9])
    prof = sp.population_profile(pot, np.array([a]), E, L, np.array([1.]), reval)
    eta = np.arccos((1 - reval/a)/e)
    err = float(np.max(np.abs(prof['mass'] - (eta - e*np.sin(eta))/np.pi)))
    period = 2*np.pi*np.sqrt(a**3/(G*M))*sp.GYR
    return dict(apsis_errors=[float(rp[0] - 4), float(ra[0] - 16)], time_fraction_max_error=err,
                period_relative_error=float(prof['period_gyr'][0]/period - 1),
                passed=bool(abs(rp[0] - 4) < 1e-4 and abs(ra[0] - 16) < 1e-4 and err < 1e-4
                            and abs(prof['period_gyr'][0]/period - 1) < 1e-4))


def plummer_df():
    """Jeans theorem: injecting f(E)~(-E)^(7/2) at every point reproduces the Plummer mass profile."""
    M, a = 1e10, 2.
    r = np.geomspace(1e-4, 1e5, 20000)
    pot = sp.SphericalPotential(r, M*r**3/(r*r + a*a)**1.5)
    edges = np.geomspace(1e-3, 1e3, 121)
    sites = np.sqrt(edges[1:]*edges[:-1])
    inj = sp.uniform_ball_injections(pot, sites, 4*np.pi*sites**3*np.diff(np.log(edges)), 32, 16)
    w = inj['weight']*(-inj['E'])**3.5
    reval = np.geomspace(.01, 100, 60)
    prof = sp.population_profile(pot, inj['r'], inj['E'], inj['L'], w, reval)
    err = float(np.max(np.abs(prof['mass']/w.sum() - reval**3/(reval**2 + a*a)**1.5)))
    return dict(max_enclosed_fraction_error=err, passed=bool(err < 2e-3))


def sampling_invariance():
    M, a = 1e10, 3.
    r = np.geomspace(1e-4, 1e6, 6000)
    pot = sp.SphericalPotential(r, M*r*r/(r + a)**2)                 # Hernquist
    sites = np.geomspace(.1, 30, 24)
    inj = sp.uniform_ball_injections(pot, sites, sites*np.exp(-sites/a), 12, 6)
    req = np.array([.7, 1.3, 2.2, 4., 6.5, 9., 14., 21.])
    full = sp.population_profile(pot, inj['r'], inj['E'], inj['L'], inj['weight'], req)['mass']
    part = sp.population_profile(pot, inj['r'], inj['E'], inj['L'], inj['weight'], req[::3])['mass']
    more = sp.population_profile(pot, inj['r'], inj['E'], inj['L'], inj['weight'], np.sort(np.r_[req, .5, 3., 30.]))['mass']
    keep = np.isin(np.sort(np.r_[req, .5, 3., 30.]), req)
    # Only floating-point summation order may differ between requests.
    diff = max(np.max(np.abs(full[::3] - part)), np.max(np.abs(full - more[keep])))/inj['weight'].sum()
    return dict(max_fraction_difference_at_common_radii=float(diff), passed=bool(diff < 1e-12))


def plummer_support():
    """Universal-K consequence of P=K rho^(6/5) for the isolated Plummer sphere."""
    M, a = 3e14, 400.
    r = np.geomspace(1e-3, 1e5, 400)
    s = sp.plummer_state(M, a, r)
    K = s['pressure']/s['rho']**1.2
    spread = float(np.max(np.abs(K/sp.plummer_K(M, a) - 1)))
    ratio = sp.plummer_scale_for_K(2*M, sp.plummer_K(M, a))/a
    # Independent Lane-Emden n=5 integration: theta = (1+xi^2/3)^(-1/2) and a = sqrt(3) alpha.
    sol = solve_ivp(lambda xi, y: [y[1], -y[0]**5 - 2*y[1]/xi], (1e-6, 50.), [1 - 1e-12/6, -1e-6/3],
                    rtol=1e-11, atol=1e-13, dense_output=True)
    xi = np.linspace(.1, 50, 200)
    le = float(np.max(np.abs(sol.sol(xi)[0]/(1 + xi*xi/3)**-.5 - 1)))
    rho_c = 3*M/(4*np.pi*a**3)
    alpha = np.sqrt(6*sp.plummer_K(M, a)*rho_c**-.8/(4*np.pi*G))
    return dict(K_constant_relative_spread=spread, scale_ratio_for_doubled_mass=float(ratio),
                lane_emden_max_relative_error=le, lane_emden_scale_over_a=float(np.sqrt(3)*alpha/a),
                statement='At fixed K an isolated Plummer sphere has a = [G (4 pi/3)^(1/5) M^(4/5)/(6K)]^(5/2), so a ~ M^2.',
                passed=bool(spread < 1e-12 and abs(ratio - 4) < 1e-12 and le < 1e-6
                            and abs(np.sqrt(3)*alpha/a - 1) < 1e-12))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    rng = np.random.default_rng(20260913)
    results = dict(event_balance=event_balance(rng), drag=drag(), uniform_ball=uniform_ball(rng), line=line(rng),
                   kepler=kepler(), plummer_distribution_function=plummer_df(),
                   sampling_invariance=sampling_invariance(), plummer_universal_K=plummer_support())
    ok = all(v['passed'] for v in results.values())
    results['all_passed'] = ok
    target = args.output
    if target is None and os.environ.get('PHOTON_GRAVITON_RESULTS'):
        target = Path(os.environ['PHOTON_GRAVITON_RESULTS'])/'capture-to-orbit-checks.json'
    text = json.dumps(results, indent=1) + '\n'
    if target:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding='utf-8', newline='\n')
    print(text)
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
