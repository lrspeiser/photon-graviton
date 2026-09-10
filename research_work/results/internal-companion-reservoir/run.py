"""Synthetic internal-state transport audit; no data or holdouts are used."""
from pathlib import Path
import hashlib
import json
import math

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import erf

HERE = Path(__file__).resolve().parent


def reservoir(P, scale, branch, curvature=0.0):
    P = np.asarray(P)
    if branch == 'quadratic':
        return P * P / (2 * scale), P / scale
    # C2, positive, globally bounded-below completion for the ideal branch.
    # The curvature variant is used only on its invariant P >= 0 domain.
    if curvature:
        assert np.min(P) >= -1e-12
        return scale + P + curvature * P * P / (2 * scale), 1 + curvature * P / scale
    neg = np.minimum(P / scale, 0)
    return (np.where(P < 0, scale * (1 + math.sqrt(math.pi) / 2 * erf(neg)), scale + P),
            np.where(P < 0, np.exp(-neg * neg), 1.0))


def simulate(branch, p, scale, a, T=6.0, curvature=0.0, tight=False):
    P0 = scale if branch == 'quadratic' else 0.0
    seed = float(reservoir(P0, scale, branch, curvature)[0])
    times = np.linspace(0, T, 301)

    def rhs(t, y):
        x, s, P = y
        n = 1 + a * s
        assert n > 0
        return [1 / n, float(reservoir(P, scale, branch, curvature)[1]), p * a / n ** 2]

    sol = solve_ivp(rhs, [0, T], [0, 0, P0], t_eval=times,
                    method='DOP853', rtol=2e-13 if tight else 2e-11,
                    atol=2e-14 if tight else 2e-12, max_step=.1 if tight else .3)
    assert sol.success
    x, s, P = sol.y
    nr = 1 + a * times
    n = 1 + a * s
    ew = p / n
    er, rate = reservoir(P, scale, branch, curvature)
    # Subtract the pre-existing seed explicitly; do not hide drift behind it.
    excess_error = float(np.max(np.abs(ew + (er - seed) - p)) / p)
    assert excess_error < 2e-9
    ref_x = math.log1p(a * T) / a if a else T
    independent_error = None
    if branch == 'quadratic':
        def rate_s(q):
            return math.sqrt(1 + 2 * p / scale * (1 - 1 / (1 + a * q)))
        def time_s(q):
            return quad(lambda z: 1 / rate_s(z), 0, q, epsabs=1e-12, epsrel=1e-12)[0]
        sf = brentq(lambda q: time_s(q) - T, 0, T * math.sqrt(1 + 2 * p / scale) + 1)
        xf = quad(lambda q: 1 / ((1 + a * q) * rate_s(q)), 0, sf,
                  epsabs=1e-12, epsrel=1e-12)[0]
        independent_error = max(abs(sf - s[-1]), abs(xf - x[-1]))
        assert independent_error < 2e-9
    elif curvature == 0:
        exact_P = p * (1 - 1 / nr)
        independent_error = max(float(np.max(abs(s - times))),
                                float(np.max(abs(P - exact_P))), abs(x[-1] - ref_x))
        assert independent_error < 2e-9
    result = dict(branch=branch, p=p, reservoir_scale=scale, curvature=curvature,
                  rolling_rate=a, duration=T, initial_reservoir_energy=seed,
                  seed_over_initial_propagation_energy=seed / p,
                  final_internal_time=float(s[-1]), final_internal_rate=float(rate[-1]),
                  final_propagation_energy_fraction=float(ew[-1] / p),
                  reservoir_gain_over_initial_propagation_energy=float((er[-1] - seed) / p),
                  retained_excess_fraction=float((ew[-1] + er[-1] - seed) / p),
                  final_speed_ratio_to_reference_photon=float(nr[-1] / n[-1]),
                  maximum_sampled_fractional_speed_mismatch=float(np.max(abs(nr / n - 1))),
                  companion_displacement=float(x[-1]), reference_photon_displacement=ref_x,
                  endpoint_lag=float(ref_x - x[-1]), energy_error_over_initial_propagation=excess_error,
                  independent_quadrature_or_exact_error=independent_error)
    if a == 0:
        assert result['maximum_sampled_fractional_speed_mismatch'] < 1e-12
    return result


def main():
    runs = []
    for a in [0.0, .05, .1]:
        for scale in [1., 10., 100., 1000.]:
            runs.append(simulate('quadratic', 1., scale, a))
        for p in [.25, 1., 4.]:
            runs.append(simulate('ideal_linear', p, 1., a))
    for eps in [.01, .1]:
        runs.append(simulate('curved_linear', 1., 1., .05, curvature=eps))
    # Same p/scale gives identical trajectories for the quadratic branch.
    base = simulate('quadratic', 1., 10., .05)
    scaled = simulate('quadratic', 4., 40., .05)
    scaling_error = abs(base['companion_displacement'] - scaled['companion_displacement'])
    assert scaling_error < 2e-10
    refined = simulate('quadratic', 1., 1., .1, tight=True)
    coarse = next(r for r in runs if r['branch'] == 'quadratic' and r['reservoir_scale'] == 1 and r['rolling_rate'] == .1)
    convergence = max(abs(refined[k] - coarse[k]) for k in ['companion_displacement', 'final_internal_time', 'final_speed_ratio_to_reference_photon'])
    assert convergence < 2e-9
    ideal_refined = simulate('ideal_linear', 4., 1., .1, tight=True)
    # Verify the completed ideal reservoir derivative on both sides of P=0.
    derivative_error = 0.0
    for P in [-10., -2., -.1, 0., .1, 2.]:
        h = 1e-5
        fd = (reservoir(P + h, 1., 'ideal_linear')[0] - reservoir(P - h, 1., 'ideal_linear')[0]) / (2 * h)
        derivative_error = max(derivative_error, abs(float(fd - reservoir(P, 1., 'ideal_linear')[1])))
    assert derivative_error < 2e-9
    result = dict(scope='Synthetic conditional transport test; no photon creation, capture, gravity or observed fit.',
                  units='Dimensionless; c0=1. No astronomical normalization or measured tolerance.',
                  hamiltonian='H=p/(1+a*s)+R(P), p>0, 1+a*s>0',
                  initial_state='x=s=0; P=scale for quadratic, otherwise P=0',
                  reference='Unloaded clock s=t; photon speed 1/(1+a*t). Reference photon is not a coupled energetic field.',
                  equations=['dx/dt=1/(1+a*s)', 'dp/dt=0', 'ds/dt=R_prime(P)', 'dP/dt=p*a/(1+a*s)^2'],
                  runs=runs, refinement_maximum=convergence, scaling_control_error=scaling_error,
                  ideal_completion_derivative_error=derivative_error,
                  ideal_completion_lower_bound_over_scale=1-math.sqrt(math.pi)/2,
                  ideal_refined=ideal_refined,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  observational_holdouts_opened=False, full_theory_validated=False)
    (HERE / 'results.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    for r in runs:
        if r['rolling_rate'] == .05:
            print(r['branch'], 'p', r['p'], 'scale', r['reservoir_scale'], 'curvature', r['curvature'],
                  'speed ratio', r['final_speed_ratio_to_reference_photon'], 'retained excess', r['retained_excess_fraction'])


if __name__ == '__main__':
    main()
