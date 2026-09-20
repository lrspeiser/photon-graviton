"""SE-LR1: isolated local interaction sector; no matter/observation claim."""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parent


def principal(q, n, speed, coupling):
    t = np.dot(n, q) * np.eye(3) - np.outer(n, q)
    z, eye = np.zeros((3, 3)), np.eye(3)
    matrix = np.block([[z, -coupling*t.T, eye, z],
                       [-coupling*t, z, z, speed**2*eye],
                       [eye, z, z, z], [z, eye, z, z]])
    sym = np.diag([1.]*9 + [speed**2]*3)
    sigma = np.linalg.svd(t, compute_uv=False)
    coefficient = 1 + speed**2 + coupling**2*sigma**2
    upper = (coefficient + np.sqrt(np.maximum(coefficient**2-4*speed**2, 0)))/2
    lower = speed**2/upper
    positive = np.sqrt(np.concatenate((upper, lower)))
    return matrix, sym, np.sort(np.concatenate((-positive, positive)))


def derivative(a, axis, h):
    return (np.roll(a, -1, axis)-np.roll(a, 1, axis))/(2*h)


def curl(a, h):
    return np.stack((derivative(a[..., 2], 1, h)-derivative(a[..., 1], 2, h),
                     derivative(a[..., 0], 2, h)-derivative(a[..., 2], 0, h),
                     derivative(a[..., 1], 0, h)-derivative(a[..., 0], 1, h)), axis=-1)


def laplacian(a, h):
    return sum((np.roll(a, 1, j)+np.roll(a, -1, j)-2*a)/h**2 for j in range(3))


def energy(state, h, coupling, speed=.5, omega=.2, Omega=1.):
    a, pi, q, p = state
    k = p-coupling*np.cross(curl(a, h), q)
    value = np.sum(pi*pi+k*k+omega**2*a*a+Omega**2*q*q)
    for j in range(3):
        value += np.sum(((np.roll(a, -1, j)-a)/h)**2)
        value += speed**2*np.sum(((np.roll(q, -1, j)-q)/h)**2)
    return float(.5*h**3*value)


def rhs(state, h, coupling, speed=.5, omega=.2, Omega=1.):
    a, pi, q, p = state
    b = curl(a, h)
    k = p-coupling*np.cross(b, q)
    return np.stack((pi, laplacian(a, h)-omega**2*a+coupling*curl(np.cross(q, k), h),
                     k, speed**2*laplacian(q, h)-Omega**2*q-coupling*np.cross(b, k)))


def main():
    output = ROOT/'local-rotor-v1'
    output.mkdir(exist_ok=False)
    rng = np.random.default_rng(20261002)
    rows = []
    for amplitude, speed, coupling in itertools.product((.03, .1, .3), (.25, .5, 1.), (0., .4, 2., 10.)):
        for fixture in range(8):
            q, n = rng.normal(size=(2, 3))
            q *= amplitude/np.linalg.norm(q)
            n /= np.linalg.norm(n)
            m, sym, predicted = principal(q, n, speed, coupling)
            numerical = np.linalg.eigvals(m)
            rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
            rotation[:, 0] *= np.linalg.det(rotation)
            rotated = principal(rotation@q, rotation@n, speed, coupling)[0]
            errors = dict(symmetry=float(np.max(abs(sym@m-m.T@sym))),
                          imaginary=float(np.max(abs(numerical.imag))),
                          spectrum=float(np.max(abs(np.sort(numerical.real)-predicted))),
                          rotation=float(np.max(abs(np.sort(np.linalg.eigvals(rotated).real)-predicted))))
            passed = (errors['symmetry'] < 1e-12 and all(errors[k] < 1e-10 for k in ('imaginary', 'spectrum', 'rotation'))
                      and np.min(abs(predicted)) > 0 and np.min(np.diag(sym)) > 0)
            rows.append(dict(amplitude=amplitude, c_Q=speed, epsilon=coupling, fixture=fixture,
                             Q=q.tolist(), direction=n.tolist(), speeds=predicted.tolist(),
                             errors=errors, passed=bool(passed)))
    derivatives = []
    for n, coupling in itertools.product((8, 12), (0., .4, 2., 10.)):
        h = 8/n
        for fixture in range(6):
            state = .03*rng.normal(size=(4, n, n, n, 3))
            direction = rng.normal(size=state.shape)
            direction /= np.linalg.norm(direction)
            flow = rhs(state, h, coupling)
            gradient = np.stack((-flow[1], flow[0], -flow[3], flow[2]))*h**3
            expected = float(np.sum(gradient*direction))
            numerical = [(energy(state+delta*direction, h, coupling)-energy(state-delta*direction, h, coupling))/(2*delta)
                         for delta in (1e-6, 5e-7)]
            scale = max(1., abs(expected))
            errors = [abs(value-expected)/scale for value in numerical]
            difference = abs(numerical[0]-numerical[1])/scale
            derivatives.append(dict(n=n, epsilon=coupling, fixture=fixture, analytic=expected,
                                    numerical=numerical, errors=errors, step_difference=difference,
                                    passed=bool(max(errors+[difference]) < 2e-7)))
    result = dict(seed=20261002, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  principal_cases=rows, derivative_cases=derivatives,
                  passed=all(r['passed'] for r in rows+derivatives))
    (output/'results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(passed=result['passed'], principal_cases=len(rows), derivative_cases=len(derivatives),
                         min_speed=min(min(abs(np.array(r['speeds']))) for r in rows),
                         max_speed=max(max(r['speeds']) for r in rows),
                         max_spectrum_error=max(r['errors']['spectrum'] for r in rows),
                         max_gradient_error=max(max(r['errors']) for r in derivatives))))
    assert result['passed']


if __name__ == '__main__':
    main()
