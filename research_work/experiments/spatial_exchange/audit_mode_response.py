"""Independent matrix-derivative audit of archived SE-1M results."""
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
data = json.loads((ROOT/'mode-response-v2/results.json').read_text())
checks = []
errors = []
expected = set(itertools.product([-.005, -.003, -.001, -.0003, 0., .0003],
                                [-50, 0, 50, 200], [0., .5], [.25, .75, 2.],
                                ['massless', 'massive']))
actual = {(r['U'], r['chi'], r['k0'], r['K'], r['branch']) for r in data['rows']}
checks.append(len(data['rows']) == 288 and actual == expected)
for filename, key in [('mode_response.py', 'source_sha256'), ('mode-response-protocol.md', 'protocol_sha256')]:
    checks.append(hashlib.sha256((ROOT/filename).read_bytes()).hexdigest() == data['summary'][key])
for row in data['rows']:
    u, chi, k0, K = (row[k] for k in ('U', 'chi', 'k0', 'K'))
    k = k0*np.tanh(1000*u)
    dk = 1000*k0*(1-np.tanh(1000*u)**2)
    m = .5625*np.exp(2*chi*u)
    V = m*np.array([[k*k, -k], [-k, 1.]])
    dV = 2*chi*V + m*np.array([[2*k*dk, -dk], [-dk, 0.]])
    values, vectors = np.linalg.eigh(V)
    j = int(row['branch'] == 'massive')
    v = vectors[:, j]
    energy = .5*(K*K+values[j])
    # Mean kinetic energy is E/2; coefficient a_U/a=4.
    response = (2*energy + .25*(v@dV@v))/energy
    error = abs(response-row['source_coefficient'])/max(1., abs(response))
    errors.append(float(error))
    checks.append(error < 1e-10 and abs(energy-row['energy']) < 1e-12)
    honest = row['scaled_error'] < 1e-5 and row['eigenvalue_error'] < 1e-12 and row['eigenvector_residual'] < 1e-12 and row['energy'] > 0 and row['group_speed'] <= row['front_speed']*(1+1e-14)
    checks.append(honest == row['passed'])
checks.append(data['summary']['passed'] == all(r['passed'] for r in data['rows']))
result = dict(passed=bool(all(checks)), checks=len(checks), maximum_independent_response_error=max(errors))
(ROOT/'mode-response-audit.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result))
assert result['passed']
