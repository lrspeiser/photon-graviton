"""SE-1M: fixed-canonical-state response diagnostic; not an observation fit."""
import hashlib
import itertools
import json
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def matrix(u, chi, k0):
    k = k0 * np.tanh(u / .001)
    return .75**2 * np.exp(2 * chi * u) * np.array([[k*k, -k], [-k, 1.]])


def fixture(u, chi, k0, wave_number, branch):
    k = k0 * np.tanh(u / .001)
    ku = k0 / .001 / np.cosh(u / .001)**2
    nu = .75**2 * np.exp(2*chi*u) * (1+k*k)
    eigenvalue = 0. if branch == 'massless' else nu
    vector = np.array([1., k] if branch == 'massless' else [-k, 1.]) / np.sqrt(1+k*k)
    a = np.exp(4*u)
    frequency = np.sqrt(a*(wave_number**2+eigenvalue))
    phase = 2*np.pi*np.arange(64)/64
    field = np.cos(phase)[:, None] * vector
    momentum = frequency/a * np.sin(phase)[:, None] * vector
    gradient = -wave_number * np.sin(phase)[:, None] * vector

    def energy(background):
        return float(np.mean(.5*np.exp(4*background)*np.sum(momentum**2, axis=1)
                             + .5*np.sum(gradient**2, axis=1)
                             + .5*np.einsum('ni,ij,nj->n', field, matrix(background, chi, k0), field)))

    base = energy(u)
    numeric = (energy(u+1e-7)-energy(u-1e-7))/(2e-7*base)
    analytic = 2. if branch == 'massless' else 2+(chi+k*ku/(1+k*k))*nu/(wave_number**2+nu)
    error = abs(numeric-analytic)/max(1., abs(analytic))
    eig_error = float(np.max(abs(np.linalg.eigvalsh(matrix(u, chi, k0))-np.array([0., nu]))))
    residual = float(np.linalg.norm(matrix(u, chi, k0)@vector-eigenvalue*vector))
    front = np.exp(2*u)
    group = front*wave_number/np.sqrt(wave_number**2+eigenvalue)
    passed = error < 1e-5 and eig_error < 1e-12 and residual < 1e-12 and base > 0 and group <= front*(1+1e-14)
    return dict(U=u, chi=chi, k0=k0, K=wave_number, branch=branch, energy=base,
                source_coefficient=analytic, numerical_coefficient=numeric, scaled_error=error,
                eigenvalue_error=eig_error, eigenvector_residual=residual,
                group_speed=group, front_speed=front, passed=bool(passed))


def main():
    directory = ROOT/'mode-response-v1'
    directory.mkdir(exist_ok=False)
    rows = [fixture(*args) for args in itertools.product(
        [-.005, -.003, -.001, -.0003, 0., .0003], [-50, 0, 50, 200], [0., .5],
        [.25, .75, 2.], ['massless', 'massive'])]
    massive = [r for r in rows if r['branch'] == 'massive']
    summary = dict(count=len(rows), passed=all(r['passed'] for r in rows),
                   maximum_scaled_error=max(r['scaled_error'] for r in rows),
                   minimum_response=min(r['source_coefficient'] for r in rows),
                   maximum_response=max(r['source_coefficient'] for r in rows),
                   massive_enhanced=sum(r['source_coefficient'] > 2+1e-10 for r in massive),
                   massive_suppressed=sum(r['source_coefficient'] < 2-1e-10 for r in massive),
                   massive_negative=sum(r['source_coefficient'] < 0 for r in massive),
                   commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   protocol_sha256=hashlib.sha256((ROOT/'mode-response-protocol.md').read_bytes()).hexdigest())
    (directory/'results.json').write_text(json.dumps(dict(summary=summary, rows=rows), indent=2)+'\n')
    print(json.dumps(summary))
    if not summary['passed']:
        raise RuntimeError('SE-1M diagnostic failed; preserve this first run')


if __name__ == '__main__':
    main()
