"""Independent endpoint energy and loop reconstruction for SE-LR2 archives."""
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def reconstruct(state, h, eps):
    a, pi, q, p = state
    n = a.shape[0]
    forward = lambda f, j: (np.take(f, (np.arange(n)+1)%n, axis=j)-f)/h
    center = lambda f, j: (np.take(f, (np.arange(n)+1)%n, axis=j)-np.take(f, (np.arange(n)-1)%n, axis=j))/(2*h)
    da = [center(a, j) for j in range(3)]
    dq = [center(q, j) for j in range(3)]
    b = np.stack([da[1][..., 2]-da[2][..., 1], da[2][..., 0]-da[0][..., 2], da[0][..., 1]-da[1][..., 0]], axis=-1)
    k = p-eps*np.cross(b, q)
    ea = (np.sum(pi*pi+.04*a*a)+sum(np.sum(forward(a, j)**2) for j in range(3)))*h**3/2
    eq = (np.sum(k*k+q*q)+.25*sum(np.sum(forward(q, j)**2) for j in range(3)))*h**3/2
    mom = -np.stack([np.sum(pi*da[j]+p*dq[j], axis=-1) for j in range(3)], axis=-1)
    x = np.stack(np.meshgrid(*(np.arange(n)*h-4 for _ in range(3)), indexing='ij'), axis=-1)
    angular = np.sum(np.cross(x, mom)+np.cross(a, pi)+np.cross(q, p), axis=(0, 1, 2))*h**3
    return np.array([ea+eq, ea, eq, *angular, *np.sum(mom, axis=(0, 1, 2))*h**3])


def loop(a, h, radius, count):
    angle = np.arange(count)*2*np.pi/count
    positions = np.stack((radius*np.cos(angle), radius*np.sin(angle), np.zeros(count)), axis=-1)
    grid = (positions+4)/h
    base = np.floor(grid).astype(int)
    fraction = grid-base
    values = np.zeros((count, 3))
    for i in range(2):
        for j in range(2):
            for k in range(2):
                index = base+np.array([i, j, k])
                w = np.prod(np.where(np.array([i, j, k]), fraction, 1-fraction), axis=-1)
                values += w[:, None]*a[index[:, 0], index[:, 1], index[:, 2]]
    return float(np.sum(-radius*np.sin(angle)*values[:, 0]+radius*np.cos(angle)*values[:, 1])*2*np.pi/count)


def main():
    directory = ROOT/'local-rotor-evolution-v1'
    summary = json.loads((directory/'summary.json').read_text())
    manifest = json.loads((directory/'manifest.json').read_text())
    checks = [hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == digest for f, digest in manifest['hashes'].items()]
    maximum = 0.
    field_angular = {}
    for row in summary['runs']:
        archive = np.load(directory/(row['name']+'.npz'))
        trace = archive['trace']
        field_only = archive['final'].copy()
        field_only[2:] = 0
        field_angular[row['name']] = reconstruct(field_only, 8/row['n'], row['epsilon'])[3:6].tolist()
        checks.append(row == json.loads((directory/(row['name']+'.json')).read_text()))
        for key, t in (('initial', 0), ('final', -1)):
            derived = reconstruct(archive[key], 8/row['n'], row['epsilon'])
            error = float(np.max(abs(derived-trace[t, 1:10])))
            maximum = max(maximum, error)
            checks.append(error < 1e-12)
        edrift = float(np.max(abs(trace[:, 1]-trace[0, 1]))/trace[0, 1])
        jabs = float(np.max(np.linalg.norm(trace[:, 4:7]-trace[0, 4:7], axis=1)))
        jrel = jabs/max(np.linalg.norm(trace[0, 4:7]), 1e-12)
        edge = float(np.max(trace[:, 10])/trace[0, 1])
        checks.extend((abs(edrift-row['energy_drift']) < 1e-12, abs(jabs-row['angular_absolute']) < 1e-12,
                       abs(jrel-row['angular_relative']) < 1e-12, abs(edge-row['edge_fraction']) < 1e-12))
        for measure in row['loops']:
            for name, count in (('coarse', 128), ('fine', 256)):
                checks.append(abs(loop(archive['final'][0], 8/row['n'], measure['radius'], count)-measure[name]) < 1e-12)
        passed = bool(np.isfinite(trace).all() and np.isfinite(archive['final']).all() and edrift < 1e-5
                      and (not row['rotating'] or jrel < .01) and edge < 1e-5
                      and max(abs(x['fine']-x['coarse']) for x in row['loops']) < 1e-5)
        checks.append(passed == row['passed'])
    fields = {r['name']:np.load(directory/(r['name']+'.npz'))['final'] for r in summary['runs']}
    errors = [np.max(abs(fields[k][:2])) for k in ('zero', 'nonrotating')]
    negative = fields['negative'].copy(); negative[:2] *= -1
    errors.append(np.max(abs(negative-fields['rotating'])))
    base = loop(fields['rotating'][0], 8/24, 1.5, 256)
    for name, n in (('time', 24), ('space', 32)):
        fine = loop(fields[name][0], 8/n, 1.5, 256)
        errors.append(abs(base-fine)/max(abs(fine), 1e-10))
    for row, value, limit in zip(summary['comparisons'], errors, (1e-12, 1e-12, 1e-12, .001, .05)):
        checks.append(abs(value-row.get('error', row.get('relative_error'))) < 1e-12)
        checks.append(bool(value < limit) == row['passed'])
    checks.append(summary['passed'] == all(r['passed'] for r in summary['runs']+summary['comparisons']))
    result = dict(checks=len(checks), passed=bool(all(checks)), maximum_endpoint_error=maximum,
                  campaign_passed=summary['passed'], final_A_angular_momentum=field_angular,
                  scope='Endpoint fields and loops independently reconstructed; intermediate ledgers checked from recorded traces.')
    (ROOT/'local-rotor-evolution-audit.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))
    assert result['passed']


if __name__ == '__main__':
    main()
