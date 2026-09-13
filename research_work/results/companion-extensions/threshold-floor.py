"""Continuing-source requirements after fading; no fits or physical ages."""
from pathlib import Path
import hashlib
import json
import math
from scipy.integrate import quad
from scipy.optimize import brentq

P = Path(__file__).resolve().parent
source = P / 'threshold-history-results.json'
base = [r for r in json.loads(source.read_text())['rows'] if r['duration'] is None]
L = 6 * math.log(10)
def weight(y):
    return math.sin(math.pi / 3) / (2 * math.pi * (math.cosh(y / 3) + .5))
norm = quad(weight, -L, L)[0]
def integrate(fn):
    return quad(lambda y: fn(math.exp(y)) * weight(y), -L, L, epsabs=1e-11)[0] / norm
def equilibrium(x):
    return integrate(lambda t: x / (x + t))

out = dict(scope='Initially equilibrated stores, abrupt fading to a constant floor; prescribed threshold kinetics, no measured source histories',
           input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(), rows=[])
for row in base:
    x = row['X']
    f0 = equilibrium(x)
    assert abs(f0 - row['equilibrium_fraction']) < 1e-8
    for retained in [.5, .9, .99]:
        y = brentq(lambda y: equilibrium(y) - retained * f0, 0, x, xtol=1e-13)
        fy = equilibrium(y)
        capture = y * (1 - fy)
        release = integrate(lambda t: t * y / (y + t))
        assert abs(capture - release) < 1e-9
        assert abs(fy / f0 - retained) < 1e-8
        # Exact infinite-support identity, independently inverted for comparison.
        beta_exact = retained**3 / (1 + (1 - retained) * x**(1/3))**3
        finite_values = []
        for u in [0, .1, 1, 10, 100, 1000]:
            f = integrate(lambda t: y/(y+t) + (x/(x+t)-y/(y+t))*math.exp(-(y+t)*u))
            finite_values.append(f)
            assert fy - 1e-10 <= f <= f0 + 1e-10
        assert all(b <= a + 1e-10 for a, b in zip(finite_values, finite_values[1:]))
        out['rows'].append(dict(galaxy=row['galaxy'], X=x, retained_fraction=retained,
            floor_X=y, floor_fraction=y/x, infinite_support_floor_fraction=beta_exact,
            initial_occupancy=f0, final_occupancy=fy,
            ongoing_capture_per_capacity_per_dimensionless_time=capture,
            ongoing_power_relative_to_initial=capture/(x*(1-f0))))
out['ranges'] = {}
for p in [.5, .9, .99]:
    rows = [r for r in out['rows'] if r['retained_fraction'] == p]
    out['ranges'][str(p)] = {k: [min(r[k] for r in rows), max(r[k] for r in rows)]
        for k in ['floor_fraction','infinite_support_floor_fraction','ongoing_power_relative_to_initial']}
(P/'threshold-floor-results.json').write_text(json.dumps(out, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
print(json.dumps(out['ranges'], indent=2))
