"""Can a uniform additive background preserve deposits without changing frozen fits?"""
from pathlib import Path
import hashlib
import io
import json
import zipfile
import numpy as np
from scipy.optimize import brentq, minimize_scalar

P = Path(__file__).resolve().parent
OLD = P.parent / 'isotropic-galaxy-transfer'
BASE = P.parents[2] / 'temporal_candidate_audit/data'
reference = json.loads((OLD/'third-radiation-retention-results.json').read_text())
saved = json.loads((OLD/'third-radiation-retention-predictions.json').read_text())
inputs = json.loads((P/'threshold-galaxies-results.json').read_text())
xs = {r['galaxy']: r['X'] for r in inputs['rows'] if r['half_width_decades'] == 0}
for name, digest in reference['input_sha256'].items():
    assert hashlib.sha256((BASE/name).read_bytes()).hexdigest() == digest
def eta(x):
    a = x**(1/3)
    return a/(1+a)
def required_background(x, p):
    hi = x
    while eta(hi)/eta(x+hi) < p:
        hi *= 10
    b = brentq(lambda b: eta(b)/eta(x+b)-p, 0, hi, xtol=1e-12)
    assert abs(eta(b)/eta(x+b)-p) < 1e-9
    return b

out = dict(scope='Uniform additive rate proxy B; exact-third response; frozen comparison plus explicitly separated training-only amplitude refits; exposed data',
           input_sha256=reference['input_sha256'],
           model_input_sha256={str(f.relative_to(P.parents[2])).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest()
               for f in [OLD/'third-radiation-retention-results.json',OLD/'third-radiation-retention-predictions.json',P/'threshold-galaxies-results.json']},
           requirements=[], scores=[], rows=[])
scenarios = [(0., 0.)]
cache = []
for p in [.5, .9, .99]:
    req = {name: required_background(x, p) for name, x in xs.items()}
    b = max(req.values())
    binding = max(req, key=req.get)
    assert min(eta(b)/eta(x+b) for x in xs.values()) >= p-1e-9
    out['requirements'].append(dict(retained_fraction=p, common_B=b, binding_galaxy=binding,
        individual_required_B_range=[min(req.values()),max(req.values())],
        background_to_local_range=[b/max(xs.values()),b/min(xs.values())]))
    scenarios.append((p,b))
with zipfile.ZipFile(BASE/'Rotmod_LTG.zip') as archive:
    for row in saved:
        name = row['galaxy']
        arr = np.atleast_2d(np.loadtxt(io.BytesIO(archive.read(name+'_rotmod.dat'))))
        vb = arr[:,3]*abs(arr[:,3])+.5*arr[:,4]*abs(arr[:,4])+.7*arr[:,5]*abs(arr[:,5])
        good = np.isfinite(arr).all(axis=1)&(arr[:,0]>0)&(arr[:,1]>0)&(arr[:,2]>0)&(vb>0)
        arr,vb = arr[good],vb[good]
        assert arr[:,0].tolist()==row['R_kpc'] and arr[:,1].tolist()==row['observed_kms']
        pred = np.array(row['predicted_kms'])
        extra = pred**2-vb
        assert extra.min() > -1e-8
        x = xs[name]
        cache.append(dict(galaxy=name, split=row['split'], vb=vb, extra=extra, observed=arr[:,1], X=x))
        for p,b in scenarios:
            ratio = eta(x+b)/eta(x)
            v = np.sqrt(vb+ratio*extra)
            if b == 0:
                assert max(abs(v-pred)) < 1e-10
            out['rows'].append(dict(galaxy=name, split=row['split'], retained_fraction=p,
                X=x, B=b, extra_gravity_multiplier=ratio,
                RMSE_kms=float(np.sqrt(np.mean((v-arr[:,1])**2))),
                log_mse=float(np.mean(np.log10(v/arr[:,1])**2)),
                max_speed_change_kms=float(max(abs(v-pred))), predicted_kms=v.tolist()))
for p,b in scenarios:
    for split in ['train','validation','test']:
        rows = [r for r in out['rows'] if r['retained_fraction']==p and r['split']==split]
        score = dict(retained_fraction=p,B=b,split=split,galaxies=len(rows),
            RMSE_kms=float(np.sqrt(np.mean([r['RMSE_kms']**2 for r in rows]))),
            log_RMS=float(np.sqrt(np.mean([r['log_mse'] for r in rows]))),
            max_speed_change_kms=max(r['max_speed_change_kms'] for r in rows),
            extra_gravity_multiplier_range=[min(r['extra_gravity_multiplier'] for r in rows),max(r['extra_gravity_multiplier'] for r in rows)])
        if b == 0:
            known = reference['models']['attenuated']
            assert min(abs(score['RMSE_kms']-known[k][split]['RMSE_kms']) for k in ['scores','finer_scores']) < 1e-8
        out['scores'].append(score)
out['training_amplitude_adjustments'] = []
for p,b in scenarios:
    def loss(amplitude, split):
        residuals = []
        for r in cache:
            if r['split'] != split:
                continue
            v = np.sqrt(r['vb']+amplitude*eta(r['X']+b)/eta(r['X'])*r['extra'])
            residuals.append(np.mean((v-r['observed'])**2))
        return float(np.mean(residuals))
    fit = minimize_scalar(lambda a: loss(a,'train'), bounds=(0,4), method='bounded', options={'xatol':1e-10})
    assert fit.success and 1e-7 < fit.x < 4-1e-7
    assert loss(fit.x,'train') <= loss(1,'train')+1e-8
    out['training_amplitude_adjustments'].append(dict(retained_fraction=p,B=b,
        C_multiplier=float(fit.x), criterion='Equal-galaxy training velocity MSE; only C refitted, other reference constants frozen; bound 0..4',
        RMSE_kms={split:float(np.sqrt(loss(fit.x,split))) for split in ['train','validation','test']}))
assert len(saved)==149 and sum(len(r['R_kpc']) for r in saved)==3150
(P/'background-retention-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out['requirements'],indent=2))
print(json.dumps(out['scores'],indent=2))
print(json.dumps(out['training_amplitude_adjustments'],indent=2))
