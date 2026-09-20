"""SE-LR3: local direct-plus-curl momentum shift, not a full gravity model."""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess

import numpy as np

import local_rotor as base
from run_local_rotor import diagnostics, circulation

ROOT = Path(__file__).resolve().parent


def shifted(state, direct):
    result = state.copy()
    result[3] -= direct*state[0]
    return result


def energy(state, h, coupling, direct):
    return base.energy(shifted(state, direct), h, coupling)


def rhs(state, h, coupling, direct):
    flow = base.rhs(shifted(state, direct), h, coupling)
    flow[1] += direct*flow[2]
    return flow


def ledger(state, coords, h, coupling, direct):
    canonical = diagnostics(state, coords, h, coupling)
    kinetic = diagnostics(shifted(state, direct), coords, h, coupling)
    canonical[:3] = kinetic[:3]
    canonical[9] = kinetic[9]
    a, pi = state[:2]
    momentum = -np.stack([np.sum(pi*base.derivative(a, j, h), axis=-1) for j in range(3)], axis=-1)
    angular = np.sum(np.cross(coords, momentum)+np.cross(a, pi), axis=(0, 1, 2))*h**3
    return np.concatenate((canonical, angular))


def main():
    out = ROOT/'local-transfer-v1'
    out.mkdir(exist_ok=False)
    manifest = dict(commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                    hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ('local_rotor.py', 'run_local_rotor.py', 'local_transfer.py')},
                    trace_columns=['time', 'energy', 'A_energy', 'KQ_energy', 'Jx', 'Jy', 'Jz', 'Px', 'Py', 'Pz', 'edge_energy', 'A_Jx', 'A_Jy', 'A_Jz'])
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    rng = np.random.default_rng(20261003)
    controls = []
    for n, eps, direct in itertools.product((8, 12), (0., 2.), (0., 1.)):
        for fixture in range(6):
            state = .03*rng.normal(size=(4, n, n, n, 3))
            direction = rng.normal(size=state.shape)
            direction /= np.linalg.norm(direction)
            h = 8/n
            flow = rhs(state, h, eps, direct)
            expected = float(np.sum(np.stack((-flow[1], flow[0], -flow[3], flow[2]))*direction)*h**3)
            numerical = [(energy(state+d*direction, h, eps, direct)-energy(state-d*direction, h, eps, direct))/(2*d) for d in (1e-6, 5e-7)]
            errors = [abs(value-expected)/max(1., abs(expected)) for value in numerical]
            difference = abs(numerical[0]-numerical[1])/max(1., abs(expected))
            controls.append(dict(n=n, epsilon=eps, direct=direct, fixture=fixture, analytic=expected, numerical=numerical,
                                 errors=errors, step_difference=difference, passed=bool(max(errors+[difference]) < 2e-7)))
    (out/'controls.json').write_text(json.dumps(dict(passed=all(r['passed'] for r in controls), cases=controls), indent=2)+'\n')
    assert all(r['passed'] for r in controls)
    cases = [('curl', 24, .02, 2., 0.), ('direct', 24, .02, 0., 1.),
             ('combined', 24, .02, 2., 1.), ('time', 24, .01, 2., 1.), ('space', 32, .01, 2., 1.)]
    rows = []
    for name, n, dt, eps, direct in cases:
        print('start', name, flush=True)
        h = 8/n
        coords = np.stack(np.meshgrid(*(np.arange(n)*h-4 for _ in range(3)), indexing='ij'), axis=-1)
        f = .03*np.maximum(1-np.sum(coords**2, axis=-1)/1.2**2, 0)**3
        state = np.zeros((4, n, n, n, 3)); state[2, ..., 0] = f; state[3, ..., 1] = f
        initial = state.copy()
        trace = [[0., *ledger(state, coords, h, eps, direct)]]
        for step in range(round(1/dt)):
            k1 = rhs(state, h, eps, direct)
            k2 = rhs(state+.5*dt*k1, h, eps, direct)
            k3 = rhs(state+.5*dt*k2, h, eps, direct)
            k4 = rhs(state+dt*k3, h, eps, direct)
            state += dt*(k1+2*k2+2*k3+k4)/6
            trace.append([(step+1)*dt, *ledger(state, coords, h, eps, direct)])
        trace = np.array(trace)
        loops = [dict(radius=r, coarse=circulation(state[0], h, r, 128), fine=circulation(state[0], h, r, 256)) for r in (1., 1.5, 2.)]
        edrift = float(np.max(abs(trace[:, 1]-trace[0, 1]))/trace[0, 1])
        jdrift = float(np.max(np.linalg.norm(trace[:, 4:7]-trace[0, 4:7], axis=1))/np.linalg.norm(trace[0, 4:7]))
        edge = float(np.max(trace[:, 10])/trace[0, 1])
        row = dict(name=name, n=n, dt=dt, epsilon=eps, direct=direct, energy_drift=edrift, angular_drift=jdrift,
                   edge_fraction=edge, initial_energy=float(trace[0, 1]), A_energy_fraction=float(trace[-1, 2]/trace[0, 1]),
                   A_angular_fraction=float(trace[-1, 13]/trace[0, 6]), A_Jz=float(trace[-1, 13]), loops=loops,
                   passed=bool(np.isfinite(state).all() and np.isfinite(trace).all() and edrift < 1e-5 and jdrift < .01
                               and edge < 1e-5 and max(abs(r['coarse']-r['fine']) for r in loops) < 1e-5))
        np.savez_compressed(out/(name+'.npz'), initial=initial, final=state, trace=trace)
        (out/(name+'.json')).write_text(json.dumps(row, indent=2)+'\n'); rows.append(row)
        print(json.dumps(row), flush=True)
    earlier = np.load(ROOT/'local-rotor-evolution-v1'/'rotating.npz')['final']
    current = np.load(out/'curl.npz')['final']
    error = float(np.max(abs(earlier-current)))
    comparisons = [dict(name='zero-direct-replay', error=error, passed=error < 1e-12)]
    reference = rows[2]
    for row, tolerance in ((rows[3], .001), (rows[4], .05)):
        c = row['loops'][1]['fine']; j = row['A_Jz']
        ce = abs(reference['loops'][1]['fine']-c)/max(abs(c), 1e-10)
        je = abs(reference['A_Jz']-j)/max(abs(j), 1e-10)
        comparisons.append(dict(name=row['name'], circulation_relative_difference=ce, angular_relative_difference=je, passed=ce < tolerance and je < tolerance))
    result = dict(runs=rows, comparisons=comparisons, passed=all(r['passed'] for r in rows+comparisons))
    (out/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(passed=result['passed'], comparisons=comparisons)), flush=True)


if __name__ == '__main__':
    main()
