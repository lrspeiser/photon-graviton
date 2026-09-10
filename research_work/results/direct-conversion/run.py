"""Reused-data direct-conversion branch benchmark; not a fresh validation."""
from pathlib import Path
import csv
import hashlib
import json
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import lambertw

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = ROOT / 'redshift_paper/all_164_groups.csv'
EXPECTED = '8a2044337ecfe108e56c9592d03d053d48169a1ef0c34405437a34f69a2844a0'
C = 299792.458
AU = 149597870700.0
MPC = 1e6 * 648000/math.pi * AU
MLY_MPC = 1e6 * C*1000 * 31557600 / MPC
OLD = json.loads((HERE.parent/'conversion-first/results.json').read_text())
assert hashlib.sha256(DATA.read_bytes()).hexdigest() == EXPECTED
with DATA.open(newline='', encoding='utf-8') as f:
    rows = list(csv.DictReader(f))
assert len(rows) == len({r['group_pgc'] for r in rows}) == 164
D = np.array([float(r['catalog_distance_mpc']) for r in rows])
Y = np.array([float(r['observed_cmb_cz_kms']) for r in rows])
Z = np.array([float(r['observed_cmb_z']) for r in rows])
SPLIT = np.array([r['split'] for r in rows])
REGION = np.array([int(float(r['ra_deg'])//90) + 4*int(float(r['dec_deg'])>=0) for r in rows])
assert np.max(abs(C*Z-Y)) < 1e-9
NEAR = sorted(np.flatnonzero(SPLIT=='test'), key=lambda i: (abs(D[i]-100*MLY_MPC), int(rows[i]['pgc'])))[:5]
NODES, WEIGHTS = np.polynomial.legendre.leggauss(32)
MODELS = ('conversion_path', 'expansion_comoving', 'conversion_luminosity', 'expansion_luminosity', 'linear_control')


def hubble_e(z):
    return np.sqrt(.3*(1+z)**3+.7)


def integral(z):
    z = np.asarray(z, dtype=float)
    samples = .5*z[..., None]*(NODES+1)
    return .5*z*np.sum(WEIGHTS/hubble_e(samples), axis=-1)


def redshift(model, scale, distance):
    x = scale*np.asarray(distance)/C
    if model == 'conversion_path':
        return np.expm1(x)
    if model == 'linear_control':
        return x
    if model == 'conversion_luminosity':
        return np.expm1(2*lambertw(x/2).real)
    z = np.array(x, copy=True)
    for _ in range(10):
        iz = integral(z)
        if model == 'expansion_comoving':
            target, slope = iz, 1/hubble_e(z)
        else:
            target, slope = (1+z)*iz, iz+(1+z)/hubble_e(z)
        z = z - (target-x)/slope
    assert np.all(z >= 0)
    return z


def fit(model, ids):
    loss = lambda scale: float(np.sum((C*redshift(model, scale, D[ids])-Y[ids])**2))
    result = minimize_scalar(loss, bounds=(0.,150.), method='bounded', options={'xatol':1e-9})
    assert result.success
    candidates = [0.,float(result.x),150.]
    answer = min(candidates, key=loss)
    assert 0 < answer < 150
    return answer


def stats(residual):
    return dict(count=len(residual), rms_km_s=float(np.sqrt(np.mean(residual**2))),
                mae_km_s=float(np.mean(abs(residual))), bias_km_s=float(np.mean(residual)))


train = np.flatnonzero(SPLIT=='train')
scales = {m:fit(m,train) for m in MODELS}
assert abs(scales['conversion_path']-C*OLD['fitted']['alpha_per_mpc']) < 1e-5
scales['conversion_path'] = C*OLD['fitted']['alpha_per_mpc']
scales['linear_control'] = OLD['fitted']['linear_c_alpha']
predictions = {m:C*redshift(m,scales[m],D) for m in MODELS}
statistics = {part:{m:stats((predictions[m]-Y)[SPLIT==part]) for m in MODELS} for part in ('train','validation','test')}
assert abs(statistics['test']['conversion_path']['rms_km_s']-OLD['statistics']['test']['conversion_exponential']['rmse_km_s']) < 1e-5

inverse_errors = []
for model in ('expansion_comoving','expansion_luminosity'):
    for i, distance in enumerate(D):
        def luminosity_or_comoving(z):
            iz = quad(lambda x:1/math.sqrt(.3*(1+x)**3+.7),0,z,epsabs=1e-13,epsrel=1e-13)[0]
            return C/scales[model]*iz*((1+z) if model=='expansion_luminosity' else 1)
        reference = brentq(lambda z:luminosity_or_comoving(z)-distance,0.,.2,xtol=1e-14,rtol=1e-14)
        inverse_errors.append(abs(reference-predictions[model][i]/C))
assert max(inverse_errors) < 1e-11
zlum = predictions['conversion_luminosity']/C
lum_reconstructed = C/scales['conversion_luminosity']*np.log1p(zlum)*np.sqrt(1+zlum)
lum_error = float(np.max(abs(lum_reconstructed/D-1)))
assert lum_error < 1e-11

oof = {m:np.empty_like(Y) for m in MODELS}
folds = []
for reg in np.unique(REGION):
    test = np.flatnonzero(REGION==reg)
    other = np.flatnonzero(REGION!=reg)
    fold_scales = {m:fit(m,other) for m in MODELS}
    for m in MODELS:
        oof[m][test] = C*redshift(m,fold_scales[m],D[test])
    folds.append(dict(region=int(reg),count=len(test),training_count=len(other),scales=fold_scales))
oof_stats = {m:stats(oof[m]-Y) for m in MODELS}
rng = np.random.default_rng(2026090923)
differences = {'geometric_convention':[], 'luminosity_convention':[]}
for _ in range(2000):
    ids = np.concatenate([np.flatnonzero(REGION==r) for r in rng.choice(np.unique(REGION),8,replace=True)])
    for label,conv,exp in [('geometric_convention','conversion_path','expansion_comoving'),
                           ('luminosity_convention','conversion_luminosity','expansion_luminosity')]:
        delta=stats((oof[conv]-Y)[ids])['rms_km_s']-stats((oof[exp]-Y)[ids])['rms_km_s']
        differences[label].append(delta)
bootstrap = {k:dict(percentile95_rms_conversion_minus_expansion_km_s=np.percentile(v,[2.5,97.5]).tolist(),
                    draws=len(v), scope='paired whole-region resampling of fixed OOF predictions; exploratory')
             for k,v in differences.items()}

target_d = 100*MLY_MPC
target_a = OLD['fitted']['alpha_per_mpc']*target_d
target_z = math.expm1(target_a)
target = dict(distance_million_ly=100., geometric_path_mpc=target_d,
              conversion_z=target_z, cz_reporting_units_km_s=C*target_z,
              wavelength_factor=1+target_z, wavelength_out_for_500nm=500*(1+target_z),
              photon_energy_fraction=math.exp(-target_a), companion_energy_fraction=-math.expm1(-target_a),
              stationary_model_event_duration_factor=1.,
              benchmark_one_plus_z_duration_factor=1+target_z,
              full_models_at_same_numerical_distance={m:float(redshift(m,scales[m],target_d)) for m in MODELS})
near = [dict(pgc=rows[i]['pgc'],group_pgc=rows[i]['group_pgc'],distance_million_ly=float(D[i]/MLY_MPC),
             distance_mpc=float(D[i]),historical_split=rows[i]['split'],observed_cmb_z=float(Z[i]),
             conversion_prediction=float(predictions['conversion_path'][i]/C),
             expansion_comoving_prediction=float(predictions['expansion_comoving'][i]/C),
             conversion_residual_km_s=float(predictions['conversion_path'][i]-Y[i])) for i in NEAR]

output = dict(status='active direct-conversion branch; conditional reused-data comparisons, not a physical completion',
              source_sha256=EXPECTED,protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),
              alpha_path_per_mpc=OLD['fitted']['alpha_per_mpc'],
              expansion_comparison=dict(omega_m=.3,omega_lambda=.7,radiation=0,active_explanation=False),
              scales=scales,statistics=statistics,coarse_oof_statistics=oof_stats,folds=folds,
              paired_bootstrap=bootstrap,target_100_million_ly=target,nearest_five_old_test_rows=near,
              checks=dict(expansion_inverse_count=len(inverse_errors),max_expansion_absolute_z_error=max(inverse_errors),
                          max_luminosity_conversion_distance_relative_error=lum_error,
                          old_conversion_rate_and_test_rms_reproduced=True,source_unchanged=True))
assert hashlib.sha256(DATA.read_bytes()).hexdigest()==EXPECTED
(HERE/'results.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n',encoding='utf-8')
with (HERE/'all_predictions.csv').open('w',newline='',encoding='utf-8') as f:
    fields=['pgc','group_pgc','historical_split','coarse_region','distance_mpc','observed_cmb_z']
    fields += [m+'_training_fit_z' for m in MODELS]+[m+'_coarse_oof_z' for m in MODELS]
    writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
    for i,r in enumerate(rows):
        row=dict(pgc=r['pgc'],group_pgc=r['group_pgc'],historical_split=r['split'],coarse_region=int(REGION[i]),
                 distance_mpc=float(D[i]),observed_cmb_z=float(Z[i]))
        row.update({m+'_training_fit_z':float(predictions[m][i]/C) for m in MODELS})
        row.update({m+'_coarse_oof_z':float(oof[m][i]/C) for m in MODELS})
        writer.writerow(row)
print(json.dumps(output,indent=2))
