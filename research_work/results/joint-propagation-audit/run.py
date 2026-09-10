"""Reproducible conditional propagation audit. Run with Python + scipy."""
from pathlib import Path
import csv
import hashlib
import json
import math
from scipy.constants import c, astronomical_unit
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
P = json.loads((HERE / 'protocol.json').read_text())
alpha_path = ROOT / P['alpha_source']
alpha = json.loads(alpha_path.read_text())['alpha_per_mpc']
MPC = 1e6 * astronomical_unit * 648000 / math.pi
YEAR = 31557600
D = P['distance_mpc']
A = alpha * D
S = math.exp(A)

def excess(a):
    # expm1(a)/a - 1 without small-a cancellation.
    return sum(a**n / math.factorial(n+1) for n in range(1, 30))

def delay(distance):
    return distance * MPC / c * excess(alpha * distance)

def analytic(start, width, age=0):
    # Times in units D/c, coordinates in units D.
    return 1 + math.expm1(A) * (age + start) + width * excess(A) + age

def numeric(start, width, age=0):
    y = age + start
    out = solve_ivp(lambda x, z: [1 + A/width*z[0]],
                    (start, start+width), [y], method='DOP853',
                    rtol=2e-13, atol=2e-14)
    assert out.success
    return float(out.y[0,-1]) + 1 - start - width

slabs = []
max_error = 0
for start, width in P['slabs']:
    for age in [0, 0.1, 1]:
        err = abs(numeric(start, width, age) - analytic(start, width, age))
        max_error = max(max_error, err)
        assert err < 1e-11
    dt = 1e-3
    derivative = (numeric(start, width, dt)-numeric(start, width, 0))/dt
    assert abs(derivative-S) < 1e-10
    assert analytic(start, width, .1)-.1 > analytic(start, width, 0)
    slabs.append({'start_fraction': start, 'width_fraction': width,
                  'event_stretch': derivative,
                  'extra_delay_years': (math.expm1(A)*start + width*excess(A))*D*MPC/c/YEAR})
assert excess(0) == 0
for lag in [1, 10, 1000]:
    # For two messengers with the same affine map, common intercept cancels.
    intercept = 123.456
    assert abs(((S*lag+intercept)-intercept)-S*lag) < 1e-10

aging_path = ROOT / 'temporal_candidate_audit/data/spectral_aging.json'
aging = json.loads(aging_path.read_text())
chi = {'no_event_stretch': 0., 'event_stretch_one_plus_z': 0.}
for name, z, rate, sigma in aging:
    chi['no_event_stretch'] += ((rate-1)/sigma)**2
    chi['event_stretch_one_plus_z'] += ((rate-1/(1+z))/sigma)**2
radio_path = ROOT / 'research_work/results/electromagnetic-audit/radio-chromaticity.csv'
radio = list(csv.DictReader(radio_path.open(newline='')))
intrinsic = []
for bound in P['illustrative_absolute_intrinsic_lag_bounds_seconds']:
    allowed = P['gamma_minus_gw_seconds'] + P['quoted_lag_error_seconds'] + S*bound
    intrinsic.append({'hypothetical_intrinsic_bound_seconds': bound,
                      'allowed_extra_flight_seconds': allowed,
                      'thin_interaction_max_source_distance_au': allowed*c/math.expm1(A)/astronomical_unit})
result = {
    'alpha_per_mpc': alpha, 'distance_mpc': D,
    'predicted_transfer_redshift_not_total_observed_z': math.expm1(A),
    'event_stretch': S,
    'uniform_minimum_extra_delay_seconds': delay(D),
    'uniform_minimum_extra_delay_years': delay(D)/YEAR,
    'intrinsic_gamma_emission_offset_needed_years': (P['gamma_minus_gw_seconds']-delay(D))/S/YEAR,
    'distance_sensitivities': [{'distance_mpc': d, 'minimum_extra_delay_years': delay(d)/YEAR}
                               for d in P['distance_sensitivity_mpc']],
    'slab_sensitivities': slabs, 'intrinsic_lag_scenarios': intrinsic,
    'reused_aging_rows': len(aging), 'reused_diagonal_chi_square': chi,
    'reused_radio_predictions': radio,
    'verification': {'max_dimensionless_ode_error': max_error,
                     'finite_interval_stretch': 'pass', 'zero_coupling': 'pass',
                     'positive_age': 'pass', 'shared_messenger_cancellation': 'pass'},
    'input_sha256': {str(path.relative_to(ROOT)).replace('\\','/'): hashlib.sha256(path.read_bytes()).hexdigest()
                     for path in [HERE/'protocol.json', alpha_path, aging_path, radio_path]}
}
(HERE/'results.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
print(json.dumps({key: result[key] for key in ['event_stretch', 'uniform_minimum_extra_delay_years',
      'intrinsic_gamma_emission_offset_needed_years', 'intrinsic_lag_scenarios',
      'reused_diagonal_chi_square', 'verification']}, indent=2))
