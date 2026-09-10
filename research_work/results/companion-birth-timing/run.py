"""Birth-state synchronization versus event timing; reused exposed spectral ages."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'electromagnetic-audit/spectral-aging-predictions.csv'


def trajectory(kappa, birth, a, distance):
    # Canonical ideal reservoir, interpreted only as the previously stipulated
    # packet model. Propagation energy starts at 1; pre-existing reservoir at 1.
    sb = kappa * birth
    nb = 1 + a * sb
    p = nb
    S = np.exp(a * distance)
    delay = nb * np.expm1(a * distance) / a
    def rhs(t, y):
        x, s, P = y
        n = 1 + a * s
        return [1 / n, 1., p * a / n ** 2]
    def arrived(t, y):
        return y[0] - distance
    arrived.terminal = True
    arrived.direction = 1
    sol = solve_ivp(rhs, [birth, birth + 1.1 * delay + 1], [0., sb, 0.],
                    events=arrived, method='DOP853', rtol=2e-12, atol=2e-13,
                    max_step=.2)
    assert sol.success and len(sol.t_events[0]) == 1
    arrival = float(sol.t_events[0][0])
    x, s, P = sol.y_events[0][0]
    energy = p / (1 + a * s)
    assert abs(arrival - birth - delay) < 2e-9
    assert abs(energy - 1 / S) < 2e-11
    assert abs(energy + P - 1) < 2e-11
    return dict(kappa=kappa, birth=birth, a=a, distance=distance,
                arrival=arrival, analytic_arrival=float(birth + delay),
                reference_spectral_factor=float(S), wave_energy=float(energy),
                reservoir_gain=float(P), initial_reservoir_energy=1.,
                expected_event_stretch=float(1 + kappa * (S - 1)))


def main():
    rows = []
    interval_checks = []
    for a in [.01, .05]:
        for D in [1., 6., 15.]:
            for k in [0., .5, 1.]:
                pair = [trajectory(k, b, a, D) for b in [0., 2., 5.]]
                rows.extend(pair)
                for left, right in zip(pair, pair[1:]):
                    measured = (right['arrival'] - left['arrival']) / (right['birth'] - left['birth'])
                    assert abs(measured - left['expected_event_stretch']) < 2e-9
                    interval_checks.append(abs(measured - left['expected_event_stretch']))
    data = pd.read_csv(SOURCE)
    assert len(data) == 35 and data['object'].is_unique
    z = data.z.to_numpy()
    obs = data.observed_aging_rate.to_numpy()
    sigma = data.sigma.to_numpy()
    assert np.all(sigma > 0) and np.all(z >= 0)
    def score(k):
        return float(np.sum(((obs - 1 / (1 + k * z)) / sigma) ** 2))
    fit = minimize_scalar(score, bounds=(0., 4.), method='bounded', options={'xatol':1e-12})
    assert fit.success and .01 < fit.x < 3.99
    target = score(fit.x) + 1
    lo = brentq(lambda k: score(k) - target, 0., fit.x)
    hi = brentq(lambda k: score(k) - target, fit.x, 4.)
    predictions = data[['object', 'z', 'observed_aging_rate', 'sigma']].copy()
    predictions['birth_reset_prediction'] = 1.
    predictions['shared_phase_prediction'] = 1 / (1 + z)
    predictions['fitted_kappa_prediction'] = 1 / (1 + fit.x * z)
    predictions.to_csv(HERE / 'predictions.csv', index=False, lineterminator='\n')
    # A simple global source-aging normalization is a sensitivity, not a
    # justified population model or a full covariance treatment.
    def profiled(k):
        g = 1 / (1 + k * z)
        A = float(np.sum(g * obs / sigma ** 2) / np.sum(g * g / sigma ** 2))
        return float(np.sum(((obs - A * g) / sigma) ** 2)), A
    nuisance = minimize_scalar(lambda k: profiled(k)[0], bounds=(0., 4.), method='bounded')
    assert nuisance.success
    # Equal source/detector ordinary-clock rates q=1/(1+a*t) undo both
    # reference observables in the kappa=1 branch. Exact finite intervals.
    clock_checks = []
    for a, D in [(.01, 6.), (.05, 15.)]:
        pair = [trajectory(1., b, a, D) for b in [2., 5.]]
        source_interval = np.log((1 + a * 5) / (1 + a * 2)) / a
        arrival_interval = np.log((1 + a * pair[1]['arrival']) / (1 + a * pair[0]['arrival'])) / a
        duration = arrival_interval / source_interval
        r = pair[0]
        energy_ratio = r['wave_energy'] * (1 + a * r['arrival']) / (1 + a * r['birth'])
        assert abs(duration - 1) < 2e-10 and abs(energy_ratio - 1) < 2e-10
        clock_checks.append(dict(a=a, distance=D, proper_duration_ratio=duration,
                                 locally_measured_energy_ratio=energy_ratio))
    result = dict(scope='Birth synchronization diagnostic and exposed-data refit, not new validation.',
                  formula_provenance='Conditional solution of stipulated packet dynamics; known integration and least squares.',
                  assumptions=['s_at_birth=kappa*t_birth', 'ds/dt=1', 'dx/dt=1/(1+a*s)',
                               'reference wave energy proportional to 1/(1+a*s)',
                               'fixed source-detector separation', 'ordinary endpoint clocks fixed for main aging comparison'],
                  arrival_law='t_arr=[1+kappa*(exp(a*D)-1)]*t_birth+(exp(a*D)-1)/a',
                  duration_law='S_event=1+kappa*(S_spectrum-1)',
                  aging_law='r_age=1/(1+kappa*z) conditional on intrinsic template aging',
                  trajectories=rows, maximum_interval_check_error=max(interval_checks),
                  spectral_aging=dict(rows=len(data), fixed_birth_reset_chi2=score(0),
                      fixed_shared_phase_chi2=score(1), fitted_kappa=float(fit.x), chi2_fit=score(fit.x),
                      formal_delta_chi2_one_interval=[lo, hi],
                      global_aging_normalization_sensitivity=dict(kappa=float(nuisance.x),
                          normalization=profiled(nuisance.x)[1], chi2=profiled(nuisance.x)[0])),
                  common_clock_cancellation=clock_checks,
                  input_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  new_holdouts_opened=False, photon_conversion_derived=False,
                  astronomical_distance_law_validated=False, strong_case_established=False)
    (HERE / 'results.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps(result['spectral_aging'], indent=2))
    print('Maximum arrival interval error:', max(interval_checks))


if __name__ == '__main__':
    main()
