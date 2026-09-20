"""SE-LR2: finite internal excitation, no material-source identification."""
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.ndimage import map_coordinates

from local_rotor import curl, derivative, rhs

ROOT = Path(__file__).resolve().parent


def diagnostics(state, coords, h, coupling):
    a, pi, q, p = state
    k = p-coupling*np.cross(curl(a, h), q)
    ea = .5*np.sum(pi*pi+.04*a*a, axis=-1)
    eq = .5*np.sum(k*k+q*q, axis=-1)
    momentum = np.empty_like(a)
    for j in range(3):
        ea += .5*np.sum(((np.roll(a, -1, j)-a)/h)**2, axis=-1)
        eq += .125*np.sum(((np.roll(q, -1, j)-q)/h)**2, axis=-1)
        momentum[..., j] = -np.sum(pi*derivative(a, j, h)+p*derivative(q, j, h), axis=-1)
    angular = np.cross(coords, momentum)+np.cross(a, pi)+np.cross(q, p)
    dv = h**3
    edge = np.max(abs(coords), axis=-1) >= 3
    return np.array([np.sum(ea+eq)*dv, np.sum(ea)*dv, np.sum(eq)*dv,
                     *np.sum(angular, axis=(0, 1, 2))*dv,
                     *np.sum(momentum, axis=(0, 1, 2))*dv, np.sum((ea+eq)[edge])*dv])


def circulation(a, h, radius, count):
    theta = np.arange(count)*2*np.pi/count
    pos = np.array([radius*np.cos(theta), radius*np.sin(theta), np.zeros(count)])
    values = np.array([map_coordinates(a[..., j], (pos+4)/h, order=1, mode='wrap', prefilter=False) for j in range(3)])
    tangent = np.array([-radius*np.sin(theta), radius*np.cos(theta), np.zeros(count)])
    return float(np.sum(values*tangent)*2*np.pi/count)


def main():
    out = ROOT/'local-rotor-evolution-v1'
    out.mkdir(exist_ok=False)
    manifest = dict(commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                    hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ('local_rotor.py', 'run_local_rotor.py')},
                    trace_columns=['time', 'energy', 'A_energy', 'KQ_energy', 'Jx', 'Jy', 'Jz', 'Px', 'Py', 'Pz', 'edge_energy'])
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    cases = [('zero', 24, .02, 0., True), ('rotating', 24, .02, 2., True),
             ('negative', 24, .02, -2., True), ('nonrotating', 24, .02, 2., False),
             ('time', 24, .01, 2., True), ('space', 32, .01, 2., True)]
    rows, finals = [], {}
    for name, n, dt, coupling, rotating in cases:
        print('start', name, flush=True)
        h = 8/n
        coords = np.stack(np.meshgrid(*(np.arange(n)*h-4 for _ in range(3)), indexing='ij'), axis=-1)
        f = .03*np.maximum(1-np.sum(coords**2, axis=-1)/1.2**2, 0)**3
        state = np.zeros((4, n, n, n, 3))
        state[2, ..., 0] = f
        state[3, ..., 1 if rotating else 0] = f
        initial = state.copy()
        trace = [[0., *diagnostics(state, coords, h, coupling)]]
        for step in range(round(1/dt)):
            k1 = rhs(state, h, coupling)
            k2 = rhs(state+.5*dt*k1, h, coupling)
            k3 = rhs(state+.5*dt*k2, h, coupling)
            k4 = rhs(state+dt*k3, h, coupling)
            state += dt*(k1+2*k2+2*k3+k4)/6
            trace.append([(step+1)*dt, *diagnostics(state, coords, h, coupling)])
        trace = np.array(trace)
        loops = [dict(radius=r, coarse=circulation(state[0], h, r, 128), fine=circulation(state[0], h, r, 256)) for r in (1., 1.5, 2.)]
        energy_error = float(np.max(abs(trace[:, 1]-trace[0, 1]))/trace[0, 1])
        angular_absolute = float(np.max(np.linalg.norm(trace[:, 4:7]-trace[0, 4:7], axis=1)))
        angular_relative = angular_absolute/max(np.linalg.norm(trace[0, 4:7]), 1e-12)
        edge_fraction = float(np.max(trace[:, 10])/trace[0, 1])
        finite = bool(np.isfinite(trace).all() and np.isfinite(state).all())
        row = dict(name=name, n=n, dt=dt, epsilon=coupling, rotating=rotating,
                   initial_energy=float(trace[0, 1]), final_A_energy=float(trace[-1, 2]),
                   energy_drift=energy_error, angular_absolute=angular_absolute, angular_relative=float(angular_relative),
                   edge_fraction=edge_fraction, loops=loops,
                   passed=bool(finite and energy_error < 1e-5 and (not rotating or angular_relative < .01)
                               and edge_fraction < 1e-5 and max(abs(x['fine']-x['coarse']) for x in loops) < 1e-5))
        np.savez_compressed(out/(name+'.npz'), initial=initial, final=state, trace=trace)
        (out/(name+'.json')).write_text(json.dumps(row, indent=2)+'\n')
        rows.append(row)
        finals[name] = state
        print(json.dumps(row), flush=True)
    comparisons = []
    for name in ('zero', 'nonrotating'):
        error = float(np.max(abs(finals[name][:2])))
        comparisons.append(dict(name=name+'-field', error=error, passed=error < 1e-12))
    mirror = finals['negative'].copy()
    mirror[:2] *= -1
    error = float(np.max(abs(mirror-finals['rotating'])))
    comparisons.append(dict(name='sign-reversal', error=error, passed=error < 1e-12))
    base = rows[1]['loops'][1]['fine']
    for row, limit in ((rows[4], .001), (rows[5], .05)):
        fine = row['loops'][1]['fine']
        error = abs(base-fine)/max(abs(fine), 1e-10)
        comparisons.append(dict(name=row['name'], relative_error=error, passed=error < limit))
    summary = dict(runs=rows, comparisons=comparisons, passed=all(r['passed'] for r in rows+comparisons))
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(dict(passed=summary['passed'], comparisons=comparisons)), flush=True)


if __name__ == '__main__':
    main()
