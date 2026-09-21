#!/usr/bin/env python3
"""CM-1: exact finite projection, not a fitted lifetime or novel formalism.

python phase_junction_network/microscopic/check_companion_memory.py --output-dir NEW_DIR

The unchanged source/conversion/capture/receiver Hamiltonian is extracted from
its original constructor. This verifier never modifies its historical archive.
"""
from __future__ import annotations
import argparse
import ast
import hashlib
import importlib.util
import json
import platform
import subprocess
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import eigh, expm
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BRIDGE = HERE / 'check_companion_bridge.py'
EXPECTED_BLOB = '13df72233738167c8a56e39ea0a835c675ee5143'


def original_hamiltonian():
    raw = BRIDGE.read_bytes()
    blob = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
    if blob != EXPECTED_BLOB:
        raise ValueError('Microscopic source changed: declare a new comparison before accepting it')
    spec = importlib.util.spec_from_file_location('cm1_original_bridge', BRIDGE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = ast.parse(raw.decode())
    fn = next(n for n in source.body if isinstance(n, ast.FunctionDef) and n.name == 'common_hamiltonian')
    prefix = []
    for stmt in fn.body:
        prefix.append(stmt)
        if isinstance(stmt, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'H' for t in stmt.targets):
            break
    else:
        raise ValueError('Original Hamiltonian construction not found')
    fn.name = '_cm1_build_original'
    fn.body = prefix + ast.parse('return H, C, names, u').body
    fn.returns = None
    tree = ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[]))
    namespace = dict(vars(module))
    exec(compile(tree, str(BRIDGE), 'exec'), namespace)
    H, components, names, u = namespace['_cm1_build_original']()
    return H, components, names, u, blob


def evolve(H, initial, times):
    e, v = eigh(H)
    coefficient = v[initial].conj()
    return v @ (coefficient[:, None] * np.exp(-1j * np.outer(e, times)))


def verification(H, components):
    rotating = H - 5 * np.eye(5)
    reference = np.diag([.19, .17, .15, .18], 1) + np.diag([.19, .17, .15, .18], -1)
    pidx, qidx = np.array([3, 4]), np.array([0, 1, 2])
    A, B, V = rotating[np.ix_(pidx, pidx)], rotating[np.ix_(qidx, qidx)], rotating[np.ix_(pidx, qidx)]
    a, b, c = .19, .17, .15
    omega = np.hypot(a, b)
    constant, oscillatory = c*c*a*a/(a*a+b*b), c*c*b*b/(a*a+b*b)

    def kernel(t):
        result = np.zeros((2, 2), complex)
        result[0, 0] = constant + oscillatory * np.cos(omega * t)
        return result

    kernel_error = max(float(np.max(np.abs(kernel(t) - V @ expm(-1j * B * t) @ V.conj().T))) for t in np.linspace(0, 100, 41))
    resolvent_errors = []
    for z in [.07+.11j, .33+.05j, -.21+.17j]:
        full = np.linalg.inv(z * np.eye(5) - rotating)[np.ix_(pidx, pidx)]
        reduced = np.linalg.inv(z * np.eye(2) - A - V @ np.linalg.solve(z * np.eye(3) - B, V.conj().T))
        resolvent_errors.append(float(np.linalg.norm(full-reduced) / np.linalg.norm(full)))
    times = [.5, 2., 20., 40., 100.]
    convolutions = []
    for nodes in [128, 256]:
        gx, gw = leggauss(nodes)
        for initial in [0, 3, 4]:
            q0 = np.eye(5, dtype=complex)[initial, qidx]
            for t in times:
                s, weights = .5*t*(gx+1), .5*t*gw
                path = evolve(rotating, initial, s)
                endpoint = evolve(rotating, initial, np.array([t]))[:, 0]
                integral_q, integral_k = np.zeros(3, complex), np.zeros(2, complex)
                for i in range(nodes):
                    prop = expm(-1j * B * (t-s[i]))
                    integral_q += weights[i] * (prop @ V.conj().T @ path[pidx, i])
                    integral_k += weights[i] * (kernel(t-s[i]) @ path[pidx, i])
                q_reconstructed = expm(-1j * B * t) @ q0 - 1j * integral_q
                drive = V @ expm(-1j * B * t) @ q0
                lhs = (rotating @ endpoint)[pidx]
                rhs = A @ endpoint[pidx] + drive - 1j * integral_k
                convolutions.append(dict(nodes=nodes, initial=initial, time=t,
                    eliminated_state_error=float(np.linalg.norm(q_reconstructed-endpoint[qidx])),
                    reduced_equation_error=float(np.linalg.norm(lhs-rhs)),
                    omitted_initial_source_error=float(np.linalg.norm(lhs-(rhs-drive)))))
    errors = dict(component_hamiltonian_sum_error=float(np.max(np.abs(H-sum(components.values())))),
        rotating_matrix_error=float(np.max(np.abs(rotating-reference))), kernel_error=kernel_error,
        max_projected_resolvent_relative_error=max(resolvent_errors),
        max_eliminated_convolution_error=max(r['eliminated_state_error'] for r in convolutions),
        max_reduced_equation_error=max(r['reduced_equation_error'] for r in convolutions),
        omitted_initial_source_control_error=max(r['omitted_initial_source_error'] for r in convolutions if r['initial']==0))
    passed = (errors['component_hamiltonian_sum_error'] < 1e-12 and errors['rotating_matrix_error'] < 1e-12
              and kernel_error < 1e-12 and max(resolvent_errors) < 1e-10
              and errors['max_eliminated_convolution_error'] < 1e-10
              and errors['max_reduced_equation_error'] < 1e-10
              and errors['omitted_initial_source_control_error'] > 1e-4)
    return dict(errors=errors, passed=bool(passed), convolution_fixtures=convolutions,
        kernel=dict(nonzero_element='K[0,0] for retained indices [3,4]',
                    formula='c^2*(a^2+b^2*cos(sqrt(a^2+b^2)*t))/(a^2+b^2)',
                    omega=omega, period=2*np.pi/omega, constant=constant, cosine_amplitude=oscillatory,
                    minimum=constant-oscillatory, maximum=constant+oscillatory,
                    Kprime_at_zero=0., Ksecond_at_zero=-c*c*b*b,
                    decay_time=None, time_units='dimensionless original Hamiltonian time',
                    interpretation='coherent finite memory, not an exponential irreversible bath'))


def retention(H, components):
    rotating = H-5*np.eye(5)
    e, v = eigh(rotating)
    if np.min(np.diff(e)) < 1e-10:
        raise ValueError('Degenerate spectrum requires a block-diagonal time-average calculation')
    times = np.linspace(0, 2200, 44001)
    projector = np.diag([0, 0, 0, 1, 1])
    results, arrays = [], dict(time=times)
    for initial in [0, 3, 4]:
        psi = evolve(rotating, initial, times)
        probs = np.abs(psi)**2
        bound = probs[3]+probs[4]
        energy = np.real(np.sum(psi.conj()*(H@psi), axis=0))
        component_total = sum(np.real(np.sum(psi.conj()*(matrix@psi), axis=0)) for matrix in components.values())
        old = times <= 220
        peak = int(np.argmax(np.where(old, bound, -1)))
        below = np.flatnonzero((times > times[peak]) & (bound < .1))
        later = np.flatnonzero(times >= 20)
        low = int(later[np.argmin(bound[later])])
        returns = np.flatnonzero(times > 220)
        best_return = int(returns[np.argmax(probs[initial, returns])])
        diag_bound = np.real(np.diag(v.conj().T@projector@v))
        avg = float(np.sum(np.abs(v[initial])**2 * diag_bound))
        result = dict(initial_sector=initial, original_window_peak_bound_probability=float(bound[peak]),
            original_window_peak_time=float(times[peak]),
            first_sampled_bound_below_point_one_after_peak=None if not len(below) else dict(time=float(times[below[0]]), probability=float(bound[below[0]])),
            minimum_bound_after_time20=float(bound[low]), minimum_time=float(times[low]),
            infinite_time_average_bound_probability=avg,
            maximum_return_probability_after_time220=float(probs[initial,best_return]), return_time=float(times[best_return]),
            norm_error=float(np.max(np.abs(probs.sum(0)-1))), energy_drift=float(np.ptp(energy)),
            component_ledger_error=float(np.max(np.abs(component_total-energy))))
        result['numerical_passed'] = bool(max(result[k] for k in ['norm_error','energy_drift','component_ledger_error']) < 1e-10)
        results.append(result)
        arrays['probability_initial_'+str(initial)] = probs
    return dict(rows=results, rotating_eigenvalues=e.tolist(),
        bound_projector_commutator_norm=float(np.linalg.norm(rotating@projector-projector@rotating)),
        passed=all(r['numerical_passed'] for r in results),
        scope='sampled finite unitary trajectories; no position-space escape, thermodynamic limit or physical years'), arrays


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError('Use a fresh directory; historical evidence is not overwritten')
    args.output_dir.mkdir(parents=True)
    H, components, names, u, blob = original_hamiltonian()
    checks = verification(H, components)
    dynamics, arrays = retention(H, components)
    try:
        gitsha = subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        gitsha = None
    output = dict(experiment='CM-1', executed_git_sha=gitsha, python=platform.python_version(),
        numpy=np.__version__, scipy=scipy.__version__, original_source_blob_sha=blob,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [BRIDGE,Path(__file__)]},
        basis=names, hamiltonian_real=H.real.tolist(), hamiltonian_imag=H.imag.tolist(),
        projection=checks, retention=dynamics,
        passed=bool(checks['passed'] and dynamics['passed']),
        novelty='standard exact projection applied to unchanged finite Phase Junction construction; no claim of a new memory formalism or gravitational law')
    (args.output_dir/'results.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(args.output_dir/'trajectories.npz', **arrays)
    summary = {k:v for k,v in output.items() if k not in ['projection','hamiltonian_real','hamiltonian_imag']}
    summary['projection'] = {k:v for k,v in checks.items() if k != 'convolution_fixtures'}
    print('CM1_SUMMARY='+json.dumps(summary,allow_nan=False),flush=True)
    if not output['passed']:
        raise SystemExit('CM-1 implementation verification failed; results preserved')


if __name__ == '__main__':
    main()
