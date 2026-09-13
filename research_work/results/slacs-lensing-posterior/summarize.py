from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
initial = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
d = json.loads((HERE/'refined-results.json').read_text(encoding='utf-8'))
assert initial['complete'] and d['complete'] and len(d['rows']) == 6
for result in (initial, d):
    for path, expected in result['hashes'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected, path
initial_failed = sum(not p['numerical_pass'] for r in initial['rows'] for p in r['predictions'])
assert all(p['numerical_pass'] for r in d['rows'] for p in r['predictions'])
lines = ['# Lensing uncertainty from the same stellar-motion mass posterior', '',
         'Only inner stellar-motion measurements determine the mass and constant '
         'orbit-anisotropy posterior. The outer stellar bin and the published lens '
         'angle are excluded from its likelihood. The frozen force and light '
         'profiles then map that mass uncertainty into a lens-angle interval.', '',
         'These six galaxies are exposed training systems. The intervals below '
         'are central 95% conditional credible intervals for the model angle, '
         'not measurement-inclusive prediction intervals or validated coverage.', '',
         '| Galaxy | Published SIE angle | Baryonic angle interval | Extra-force angle interval |',
         '|---|---:|---:|---:|']
outside = {'baryons':0, 'empirical_extra':0}
for r in d['rows']:
    chosen = [p for p in r['predictions'] if p['prior']=='uniform_log_mass_and_beta']
    intervals = [p['resolutions'][-1]['angle_quantiles_arcsec'] for p in chosen]
    for p, vals in zip(chosen, intervals):
        outside[p['model']] += not (vals[0] <= r['catalog_SIE_arcsec'] <= vals[-1])
    a, b = intervals
    lines.append(f"| {r['Name']} | {r['catalog_SIE_arcsec']:.3f} | {a[0]:.3f} to {a[-1]:.3f} | {b[0]:.3f} to {b[-1]:.3f} |")
maxchange = max(p['max_angle_quantile_change_arcsec'] for r in d['rows'] for p in r['predictions'])
priorchange = max(abs(r['predictions'][i]['resolutions'][-1]['angle_quantiles_arcsec'][2]-r['predictions'][i+1]['resolutions'][-1]['angle_quantiles_arcsec'][2]) for r in d['rows'] for i in (0,2))
lines += ['', 'All angles are in arcseconds. Reference prior: uniform log mass and '
          'uniform anisotropy. A catalog SIE angle is a fitted nonspherical lens '
          'summary; it is not an exact measured spherical Einstein radius.', '',
          f"The catalog summary falls outside the displayed model interval in {outside['baryons']}/6 "
          f"baryonic cases and {outside['empirical_extra']}/6 extra-force cases. "
          'This descriptive count is not a hypothesis rejection probability. '
          'It omits lens-summary uncertainty and the structural assumptions below.', '',
          f'Changing to a uniform mass prior changes the median model angle by '
          f'at most {priorchange:.5f} arcseconds.', '',
          '## Method and formula provenance', '',
          'Bayesian marginalization, Jeans projection and the weak-field lens '
          'integral are established mathematics. The additional force remains a '
          'frozen empirical prescription, not a photon-derived deposit source.', '',
          'For mass M and orbit anisotropy beta, the Gaussian inner-bin likelihood '
          'uses the full released covariance. Its posterior is integrated over beta '
          'to obtain p(M | inner motions). The same ComponentModel.angle(M) '
          'used in the earlier fixed-mass lens calculations supplies the angle. '
          'No lensing normalization is fitted. Mass quantiles at 2.5, 16, 50, 84 '
          'and 97.5 percent are mapped through that increasing lens-angle function; '
          'their ordering is checked at both numerical resolutions.', '',
          'The prescribed extra force scales with M to the fixed positive power p. '
          'Both force contributions increase with M at fixed radius. The adopted '
          'single Einstein-root branch gives the increasing mapping used here.', '',
          '## Numerical verification', '',
          f'The initial 801/1601 mass-node and 321/641 anisotropy-node comparison '
          f'failed {initial_failed} of 24 interval checks. Its results remain saved. '
          'The repeated 3201/6401 mass-node and 641/1281 anisotropy-node comparison '
          'passes all 24 unchanged gates: at most 0.005 in log-mass quantiles and '
          '0.005 arcseconds in angle quantiles. '
          f'The largest refined angle change is {maxchange:.6f} arcseconds.', '',
          'These checks establish integration stability at the stated tolerance, '
          'not physical or observational agreement. Both runs retain dependency '
          'hashes, checked by this report generator.', '',
          '## Remaining assumptions and consequence', '',
          'Light components, their spherical deprojection, constant mass-to-light '
          'ratio, seeing, covariance, static redshift-derived geometry, equality '
          'of the two metric potentials, force coefficients and outer cutoff are '
          'fixed. Their uncertainties are absent. The same omissions affect the '
          'inner stellar inference and can correlate its errors with lensing.', '',
          'A genuine joint test needs image/shear measurements with uncertainties '
          'and a consistently modeled stellar distribution. These intervals '
          'cannot substitute for that likelihood. They expose which discrepancies '
          'persist when mass and orbit uncertainty is included, without allowing '
          'a separate lensing mass fit. The companion production, capture, and '
          'transport equations still have to supply the gravitational source.', '',
          'Reproduce: run `run.py`, `refine.py`, then `summarize.py` from this folder '
          'or by repository-relative script paths. No reserved outcomes opened.', '']
(HERE/'report.md').write_text('\n'.join(lines), encoding='utf-8', newline='\n')
print(json.dumps(dict(initial_failed=initial_failed, refined_passed=24, outside=outside,
                      max_angle_change_arcsec=maxchange, max_prior_median_change_arcsec=priorchange)))
