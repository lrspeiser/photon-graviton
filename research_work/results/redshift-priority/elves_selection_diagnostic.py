"""Illustrate distance-selection effects; no target outcomes used."""
import hashlib
import json
import math
from pathlib import Path
from scipy.integrate import quad
from scipy.special import ndtr

HERE = Path(__file__).resolve().parent


def phi(x):
    return math.exp(-x*x/2)/math.sqrt(2*math.pi)


if __name__ == '__main__':
    rows = []
    sigma, boundary, nsigma = 0.5, 10.0, 2.0
    cutoff = boundary + nsigma*sigma
    # Declared synthetic grid, not candidate distances or fitted parameters.
    for true_d in (8.0, 10.0, 11.0, 12.0):
        a = (cutoff-true_d)/sigma
        acceptance = float(ndtr(a))
        mean_error = -sigma*phi(a)/acceptance
        nominal_coverage = max(0.0, float(ndtr(min(a, 1.96)))-float(ndtr(-1.96)))/acceptance
        integral_acceptance = quad(phi, -math.inf, a, epsabs=1e-12)[0]
        integral_mean = sigma*quad(lambda x: x*phi(x), -math.inf, a, epsabs=1e-12)[0]/integral_acceptance
        assert abs(acceptance-integral_acceptance) < 1e-10
        assert abs(mean_error-integral_mean) < 1e-9
        rows.append({'synthetic_true_distance_mpc': true_d,
                     'acceptance_probability': acceptance,
                     'mean_measured_minus_true_mpc_given_selected': mean_error,
                     'coverage_of_measured_plus_minus_1p96sigma_given_selected': nominal_coverage})
    screen = json.loads((HERE/'elves-quarantine-screen.json').read_text())
    features = {r['target_name']: r for r in json.loads((HERE/'elves-field-feature-audit.json').read_text())['records']}
    method_counts = {}
    for r in screen['records']:
        if r['decision'] == 'pending_freshness_and_host_audit':
            method = features[r['target_name']]['distance_method']
            method_counts[method] = method_counts.get(method, 0)+1
    out = {'status': 'synthetic selection diagnostic only; no candidate distances/redshifts used',
           'formula_provenance': 'established Gaussian truncation algebra, applied to an explicitly simplified selection example; not a new physical law',
           'fixed_toy_distance_sigma_mpc': sigma, 'lower_bound_boundary_mpc': boundary,
           'lower_bound_sigma_multiple': nsigma, 'measured_distance_cutoff_mpc': cutoff,
           'assumptions': ['unbiased Gaussian distance errors before selection', 'fixed error for every synthetic object',
                           'only the distance lower-bound selection is applied',
                           'fixed true distances; no population prior, spectrum-availability or imaging selection'],
           'pending_method_counts': method_counts,
           'input_hashes': {n: hashlib.sha256((HERE/n).read_bytes()).hexdigest()
                            for n in ('elves-quarantine-screen.json', 'elves-field-feature-audit.json')},
           'records': rows}
    (HERE/'elves-selection-diagnostic.json').write_text(json.dumps(out, indent=2)+'\n', newline='\n')
    print(json.dumps(out, indent=2))
