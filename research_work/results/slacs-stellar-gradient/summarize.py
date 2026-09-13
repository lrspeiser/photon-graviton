from pathlib import Path
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
old = json.loads((HERE.parent/'slacs-outer-bin-check/results.json').read_text(encoding='utf-8'))
oldrows = {(r['Name'], r['model']):r for r in old['rows']}
lines = ['# Fixed stellar mass-to-light gradient sensitivity', '',
         'This experiment changes the stellar mass profile while preserving the '
         'observed light profile used to weight stellar velocities. Each scenario '
         'refits total mass and constant orbit anisotropy to inner bins only. '
         'Outer stellar bins and lens angles never enter the likelihood.', '',
         'Stellar population gradients are a known modeling concern, not a new '
         'part of our photon hypothesis. For example, '
         '[Mehrgan et al.](https://arxiv.org/abs/2309.15911) investigate radial '
         'mass-to-light variation in other early-type galaxies. That study does '
         'not constrain the gradients of these six objects or validate the '
         'illustrative factors used here.', '',
         '## Effective postulate and method', '',
         'The positive mass-density shape is proportional to the deprojected '
         'light density times:', '',
         '    w(r) = 1 + eta / (1 + (r / R_pivot)^2).', '',
         'R_pivot is the original published angular effective radius converted '
         'with the unchanged conditional distance. Eta=-0.5, 0 and 1 give central '
         'to asymptotic outer mass/light ratios 0.5, 1 and 2. These are illustrative '
         'sensitivity cases, not measured gradients or allowed confidence bounds. '
         'Every resulting mass profile is normalized to its fitted total mass.', '',
         'The bounded shape above is an analyst-chosen empirical parametrization; '
         'no novelty is claimed. Abel deprojection, Jeans projection, Gaussian '
         'covariance and lensing integration are established mathematics. The '
         'extra force is the same frozen empirical prescription as before.', '',
         'Only the gravitating mass density changes. The tracer light density, '
         'aperture weights and luminosity normalization in the velocity projection '
         'remain fixed. Treating the altered mass profile as altered observed '
         'light would be a different experiment.', '',
         '## Results', '',
         '| Central/outer mass-light ratio | Force | Inner-bin summed chi-square | Lens fractional RMS | Orbit boundary fits |',
         '|---|---|---:|---:|---:|']
allrows = []
max_repro = 0.
for eta, tag in [(-.5,'minus0p5'), (0,'0'), (1,'1')]:
    d = json.loads((HERE/f'eta-{tag}.json').read_text(encoding='utf-8'))
    assert len(d['rows']) == 12
    for path, digest in d['hashes'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    for r in d['rows']:
        r['eta'] = eta
        allrows.append(r)
        if eta == 0:
            previous = oldrows[(r['Name'],r['model'])]
            max_repro = max(max_repro, abs(r['chi2']-previous['chi2']),
                            abs(r['predicted_angle_arcsec']-previous['predicted_angle_arcsec']))
    for model in ['baryons','empirical_extra']:
        rows = [r for r in d['rows'] if r['model']==model]
        rms = np.sqrt(np.mean([(r['angle_ratio']-1)**2 for r in rows]))
        lines.append(f"| {1+eta:g} | {model} | {sum(r['chi2'] for r in rows):.2f} | {100*rms:.2f}% | {sum(r['beta_at_boundary'] for r in rows)} |")
assert max_repro < 1e-6
maxfine = max(r['maximum_relative_vrms_refinement'] for r in allrows)
assert maxfine < 1e-4
lines += ['', 'The extra-force improvement in inner chi-square is 18.69 for '
          'constant mass/light, but only 3.83 for the illustrative centrally '
          'reduced profile. The latter also has two extra-force orbit-boundary '
          'fits. Thus the apparent advantage depends materially on stellar '
          'structure; these results do not independently establish the lower '
          'central mass/light profile or eliminate lensing mismatches.', '',
          'The likelihood uses 34 inner bins in six exposed training galaxies. '
          'Each galaxy has its own fitted mass and anisotropy; the gradient is '
          'fixed at the same value across the six for each scenario. The table '
          'does not select a preferred scenario or supply a joint model score. '
          'The lens fractional RMS compares a spherical prediction with a '
          'published SIE summary and omits its measurement uncertainty.', '',
          '## Individual outer and lens predictions', '',
          'The outer residual below uses the conditional Gaussian measurement '
          'covariance at fitted parameters. It excludes parameter and structural '
          'uncertainty and is not a calibrated significance.', '',
          '| Galaxy | eta | Force | Outer conditional standardized residual | Lens prediction / catalog |',
          '|---|---:|---|---:|---:|']
for r in allrows:
    lines.append(f"| {r['Name']} | {r['eta']:g} | {r['model']} | {r['outer_conditional_standardized_residual']:.3f} | {r['angle_ratio']:.3f} |")
lines += ['', '## Checks and limitations', '',
          f'Eta=0 reproduces the earlier fit and lens calculation within '
          f'{max_repro:.3g} in the checked chi-square/angle values. '
          f'The largest relative velocity change on doubling the radial, angular '
          f'and deprojection quadrature is {maxfine:.3g}, below the fixed 1e-4 gate.', '',
          'This is a best-fit sensitivity calculation, not posterior marginalization '
          'over stellar gradients. Four optimizer starts are retained in every '
          'case; boundary solutions are reported, not removed. Positive mass '
          'density alone does not establish a physical stellar distribution function.', '',
          'A gradient that improves exposed residuals must still be constrained '
          'by spatially resolved stellar populations and a consistent orbit model '
          'before adoption. The exercise supplies neither a companion capture '
          'law nor a funded deposit source. No reserved observations were opened.', '',
          'Run `run.py`, then `summarize.py`. The driver reuses the unchanged '
          'earlier inner-bin algorithm and records its source hash along with '
          'the new mass-profile implementation and all input dependencies.', '']
(HERE/'report.md').write_text('\n'.join(lines), encoding='utf-8', newline='\n')
print(json.dumps(dict(rows=len(allrows), max_baseline_reproduction=max_repro, max_relative_refinement=maxfine)))
