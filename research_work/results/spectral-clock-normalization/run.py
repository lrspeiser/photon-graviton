"""Profile a common spectral-clock normalization without altering measured rows."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from scipy.optimize import minimize_scalar, brentq

HERE = Path(__file__).resolve().parent
source = HERE.parent/'electromagnetic-audit/spectral-aging-predictions.csv'
with source.open(encoding='utf-8', newline='') as f:
    rows = list(csv.DictReader(f))
z, y, sigma = [np.array([float(r[k]) for r in rows]) for k in ['z','observed_aging_rate','sigma']]
assert len(rows)==35 and np.all(sigma > 0)
results = []
for name, mask in [('all_35', np.ones(len(z),bool)), ('nearby_22', z<.1), ('distant_13', z>=.1)]:
    x, values, w = np.log1p(z[mask]), y[mask], sigma[mask]**-2
    def profile(b):
        f = np.exp(-b*x)
        a = np.sum(w*f*values)/np.sum(w*f*f)
        residual = values-a*f
        assert abs(np.sum(w*f*residual)) < 1e-8
        return float(np.sum(w*residual**2)), float(a)
    fit = minimize_scalar(lambda b:profile(b)[0], bounds=(-8,8), method='bounded', options={'xatol':1e-11})
    assert fit.success
    chi, a = profile(fit.x)
    def endpoint(edge):
        if profile(edge)[0]-chi < 1:
            return None
        return float(brentq(lambda b:profile(b)[0]-chi-1, *sorted([edge,float(fit.x)])))
    fixed = {str(b):dict(chi2=profile(b)[0], normalization=profile(b)[1]) for b in [0,1]}
    old = {str(b):float(np.sum(w*(values-np.exp(-b*x))**2)) for b in [0,1]}
    results.append(dict(sample=name, n=int(mask.sum()), b=float(fit.x), normalization=a, chi2=chi,
                        formal_delta_chi2_one_interval=[endpoint(-8),endpoint(8)],
                        profile_search_bounds=[-8,8], fixed_b_free_normalization=fixed,
                        fixed_unit_normalization=old,
                        no_stretch_minus_unit_stretch_chi2=fixed['0']['chi2']-fixed['1']['chi2']))
assert results[1]['n']==22 and results[2]['n']==13
assert abs(results[0]['fixed_unit_normalization']['0']-150.569) < .002
assert abs(results[0]['fixed_unit_normalization']['1']-26.949) < .002
out = dict(scope='Shared calibration sensitivity on exposed published spectral-aging estimates',
           formula='aging_rate = a * (1+z)^(-b)', interpretation='Effective fit, not a microscopic derivation',
           source='https://arxiv.org/abs/0804.3595', results=results,
           assumptions='Published diagonal uncertainties; no full shared template covariance or selection likelihood',
           hashes={source.name:hashlib.sha256(source.read_bytes()).hexdigest(),
                   'run.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(results,indent=2))
