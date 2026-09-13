"""Report the completed grid integration without rerunning the galaxy models."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
data = json.loads((HERE / 'results.json').read_text(encoding='utf-8'))
for name, expected in data['hashes'].items():
    assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
rows = [r for s in data['systems'] for r in s['predictions']]
assert len(data['systems']) == 6 and len(rows) == 24
assert all(r['numerical_pass'] for r in rows)
lines = [
    '# Outer stellar motions with mass and orbit uncertainty', '',
    'The empirical extra-force model still misses three of six outer measurements '
    'at the level of its central 95% conditional posterior predictive intervals. '
    'The baryonic baseline misses four. Neither count is a calibrated rejection '
    'probability for the physical hypothesis. These are previously exposed training '
    'galaxies, not the reserved independent prediction test.', '',
    '## Observed values and conditional predictions', '',
    'All entries are stellar root-mean-square line-of-sight velocities in km/s, '
    'not circular rotation speeds. Intervals include measurement covariance and '
    'mass/orbit parameter uncertainty, with uniform log mass and uniform orbit '
    'anisotropy as the reference prior.', '',
    '| Galaxy | Observed outer bin | Baryonic 95% interval | Extra-force 95% interval |',
    '|---|---:|---:|---:|'
]
for s in data['systems']:
    selected = [r for r in s['predictions'] if r['prior'] == 'uniform_log_mass_and_beta']
    a, b = [r['grids'][-1]['interval95_kms'] for r in selected]
    lines.append(f"| {s['Name']} | {selected[0]['observed_outer_kms']:.2f} | "
                 f'{a[0]:.2f} to {a[1]:.2f} | {b[0]:.2f} to {b[1]:.2f} |')
max_prior = max(abs(s['predictions'][i]['grids'][-1]['mean_kms'] -
                    s['predictions'][i+1]['grids'][-1]['mean_kms'])
                for s in data['systems'] for i in (0, 2))
max_edge = max(r['grids'][-1]['beta_edge_probability'] for r in rows)
lines += ['',
    f'Using uniform mass instead changes any predictive mean by at most {max_prior:.4f} km/s. '
    'It leaves the inside/outside interval classifications unchanged. '
    'J1112+0826, J1402+6321 and J1630+4520 remain above the extra-force intervals.', '',
    '## Calculation and formula provenance', '',
    'The Jeans projection, conditional multivariate Gaussian, Bayes integration and '
    'trapezoidal quadrature are established mathematics. The extra-force profile is '
    'a frozen empirical prescription from the earlier galaxy fit; it is not derived '
    'from companion production, transport or capture.', '',
    'Let t denote the inner bins and h the outer bin. Only the inner measurements '
    'enter the mass-and-anisotropy posterior. At each parameter pair:', '',
    r'\[\mu_{h|t}=\mu_h+C_{ht}C_{tt}^{-1}(y_t-\mu_t),\qquad',
    r'V_{h|t}=C_{hh}-C_{ht}C_{tt}^{-1}C_{th}.\]', '',
    'The final prediction integrates this Gaussian over the inner-data posterior. '
    'Its variance includes both conditional measurement variance and variation of '
    'the conditional mean across that posterior. The outer observation is used '
    'only to evaluate the prediction, never to weight the parameter grid.', '',
    'Mass bounds are 1e7 to 1e14 solar masses; constant anisotropy bounds are '
    '-2 to 0.45. Both priors are uniform in anisotropy. Light profiles, geometry, '
    'seeing, force coefficients, cutoff and covariance remain fixed.', '',
    '## Numerical evidence and limits', '',
    'All 24 model/prior calculations pass the declared coarse/fine grid checks: '
    '801 by 321 versus 1601 by 641 mass/anisotropy nodes, changes no larger than '
    '0.005 in observed predictive CDF and 0.1 km/s in mean or standard deviation. '
    'This checks parameter integration, not physical accuracy. Saved source hashes '
    'are checked when this report is regenerated.', '',
    f'Maximum posterior probability within 0.05 of an anisotropy boundary is '
    f'{100*max_edge:.2f}%. Finite parameter bounds remain an assumption.', '',
    'The six objects were selected by available released light components and the '
    'release use flag, as documented in the preceding audits. J1538 lacks the '
    'required component profile; J0330 has release use flag zero. Their absence '
    'must not be interpreted as successful predictions.', '',
    'These intervals omit light-profile, seeing, geometry, force-law, covariance '
    'estimation and selection uncertainty. The two mass priors are a limited '
    'sensitivity check, not all possible prior choices. No population coverage '
    'claim or hypothesis p-value follows from six already-exposed systems.', '',
    '## Consequence for the research goal', '',
    'Mass and constant-orbit uncertainty alone does not remove the radial-profile '
    'mismatches under the current fixed assumptions. Further adjustment must have '
    'independent physical or observational support, rather than being selected '
    'solely to match these six outer bins. Lensing must then be recomputed from '
    'the same gravitational source. Companion-source derivation and the reserved '
    'prediction test remain unfinished.', '',
    'Reproduce with `python research_work/results/slacs-outer-predictive/run.py`, '
    'then `python research_work/results/slacs-outer-predictive/summarize.py`.', ''
]
for model in ('baryons', 'empirical_extra'):
    counts = []
    for prior in ('uniform_log_mass_and_beta', 'uniform_mass_and_beta'):
        chosen = [r for r in rows if r['model'] == model and r['prior'] == prior]
        counts.append(sum(not (r['grids'][-1]['interval95_kms'][0] <= r['observed_outer_kms'] <=
                              r['grids'][-1]['interval95_kms'][1]) for r in chosen))
    assert counts == ([4, 4] if model == 'baryons' else [3, 3]), counts
(HERE / 'report.md').write_text('\n'.join(lines), encoding='utf-8', newline='\n')
print('Verified hashes, 24 integration gates, and both-prior interval classifications; wrote report.md')
