"""Conditional two-channel comparison; uses exposed observations, no new holdout."""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
source = HERE.parent / 'electromagnetic-audit' / 'spectral-aging-predictions.csv'
with source.open(newline='', encoding='utf-8') as handle:
    rows = list(csv.DictReader(handle))

# DES published central estimate; monotone least-distance fit with physical b bounds.
# No probability distribution is assigned to its quoted systematic estimate.
des_b, stat, systematic = 1.003, .005, .010
fit_b = min(1., max(0., des_b))
comparisons = []
predictions = []
for fraction in [0., .012, .027, .1, .5, 1.]:
    b = 1 - fraction
    scores = {'all': 0., 'low_z': 0., 'high_z': 0.}
    for row in rows:
        z, measured, sigma = map(float, (row['z'], row['observed_aging_rate'], row['sigma']))
        predicted = (1 + z)**(-b)
        residual = (measured - predicted) / sigma
        group = 'low_z' if z < .04 else 'high_z'
        scores['all'] += residual**2
        scores[group] += residual**2
        predictions.append(dict(object=row['object'], z=z, observed_aging_rate=measured,
                                sigma=sigma, stationary_log_shift_fraction=fraction,
                                predicted_aging_rate=predicted, standardized_residual=residual))
    comparisons.append(dict(stationary_log_shift_fraction=fraction, b=b,
                            des_exponent_residual=b-des_b, diagonal_chi_square=scores,
                            duration_at_z1=2**b, bolometric_flux_ratio_at_z1=2**(-1-b)))

result = dict(
    status='Conditional phenomenological decomposition, not a derived interaction or blind test',
    input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    count=len(rows), high_z_count=sum(float(r['z']) > .2 for r in rows),
    calibration_source='https://arxiv.org/html/2406.05050v2',
    comparison_source='https://arxiv.org/html/0804.3595v1#S4.T3',
    calibration=dict(published_b=des_b, statistical_scale=stat, systematic_estimate=systematic,
                     constrained_b=fit_b, stationary_fraction=1-fit_b,
                     method='Project central estimate onto [0,1]; no joint likelihood or Gaussian systematic'),
    sensitivity_bands=[dict(multiplier=k, lower_b=des_b-k*(stat+systematic),
                            maximum_stationary_fraction=max(0., 1-(des_b-k*(stat+systematic))))
                       for k in [1,2,3]],
    comparisons=comparisons)
(HERE / 'channel-budget.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
with (HERE / 'channel-budget-predictions.csv').open('w', newline='', encoding='utf-8') as handle:
    writer=csv.DictWriter(handle, fieldnames=list(predictions[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows(predictions)
print(json.dumps(result, indent=2))
