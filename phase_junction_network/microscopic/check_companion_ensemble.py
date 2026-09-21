#!/usr/bin/env python3
"""CM-1E: a conditional many-source limit of the unchanged finite Hamiltonian.

Post-CM-1 analytic exploration, not a preregistered physical discovery.
Run alongside check_companion_memory.py and check_companion_bridge.py.
Standard diagonal-ensemble and independent-copy statistics are not new physics.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, expm
from numpy.polynomial.legendre import leggauss
from check_companion_memory import original_hamiltonian


def calculate():
    H, components, names, _, blob = original_hamiltonian()
    e, v = eigh(H)
    c = v[0].conj()
    if np.min(np.diff(e)) < 1e-10:
        raise ValueError('Degenerate energies require a spectral-block dephasing rule')
    rho = (v * np.abs(c)**2) @ v.conj().T
    P = np.diag([0., 0., 0., 1., 1.])
    p = float(np.trace(rho @ P).real)
    variance = p * (1 - p)
    delta = e[:, None] - e[None, :]
    age_rows = []
    for T in [220., 2200., 22000.]:
        window = np.exp(-.5j * delta * T) * np.sinc(delta * T / (2*np.pi))
        age_state = v @ ((c[:, None] * c[None, :].conj()) * window) @ v.conj().T
        age_rows.append(dict(maximum_age=T, bound_probability=float(np.trace(age_state @ P).real),
            distance_to_diagonal_state=float(np.linalg.norm(age_state-rho)),
            energy=float(np.trace(H @ age_state).real)))
    # Independent quadrature of direct matrix-exponential pure states, not the
    # sinc formula or the same eigenphase evolution.
    gx, gw = leggauss(256)
    T = 220.
    qstate = np.zeros_like(H)
    initial = np.eye(5, dtype=complex)[:, 0]
    for t, w in zip((gx+1)*T/2, gw/2):
        psi = expm(-1j*H*t) @ initial
        qstate += w * np.outer(psi, psi.conj())
    w = np.exp(-.5j*delta*T) * np.sinc(delta*T/(2*np.pi))
    exact_age = v @ ((c[:, None]*c[None, :].conj())*w) @ v.conj().T
    # Explicit two-copy density matrix checks the independent-copy variance.
    rho2 = np.kron(rho, rho)
    N2 = np.kron(P, np.eye(5)) + np.kron(np.eye(5), P)
    mean2 = float(np.trace(rho2 @ N2).real)
    var2 = float(np.trace(rho2 @ N2 @ N2).real) - mean2**2
    checks = dict(trace_error=float(abs(np.trace(rho)-1)),
        hermiticity_error=float(np.linalg.norm(rho-rho.conj().T)),
        minimum_eigenvalue=float(np.min(np.linalg.eigvalsh(rho))),
        stationarity_error=float(np.linalg.norm(H@rho-rho@H)),
        initial_energy_error=float(abs(np.trace(H@rho).real-H[0,0].real)),
        component_ledger_error=float(abs(sum(np.trace(C@rho).real for C in components.values())-np.trace(H@rho).real)),
        direct_exponential_age_quadrature_error=float(np.linalg.norm(qstate-exact_age)),
        two_copy_mean_error=abs(mean2-2*p), two_copy_variance_error=abs(var2-2*variance))
    passed = checks['minimum_eigenvalue'] >= -1e-12 and all(value < 1e-10 for key,value in checks.items() if key != 'minimum_eigenvalue')
    illustrative_target = .01
    return dict(experiment='CM-1E', mode='post-CM-1 analytic exploration; parameters unchanged',
        input_hamiltonian_blob=blob, basis=names,
        state_real=rho.real.tolist(), state_imag=rho.imag.tolist(),
        mean_bound_probability=p, occupation_variance_per_copy=variance,
        initial_and_ensemble_energy=float(np.trace(H@rho).real),
        independent_copy_relative_RMS_coefficient=float(np.sqrt((1-p)/p)),
        independent_copy_examples=[dict(copies=N, relative_occupation_RMS=float(np.sqrt((1-p)/(p*N)))) for N in [10000,1000000]],
        illustrative_relative_RMS_target=illustrative_target,
        minimum_independent_copies_for_target=int(np.ceil((1-p)/(p*illustrative_target**2))),
        large_N_uniform_pairwise_correlation_ceiling_for_target=float(illustrative_target**2*p/(1-p)),
        finite_uniform_age_ensembles=age_rows, checks=checks, passed=bool(passed),
        assumptions=['Identical noninteracting copies with independent age/phase preparation',
            'Bound occupation is the stated projector, not already a derived gravitational source',
            'Diagonal state is an ensemble or time average, not irreversible evolution of one pure state'],
        interpretation='A nonzero stationary mean bound occupation is compatible with reversible microscopic exchange. This supplies a conditional collective route, not a spatial deposit, source-density prediction or new statistical formalism.',
        unproved=['Physical source-age distribution and normalization', 'Spatial interactions and shared-frame correlations',
            'Backreaction, capacity, local transport, field fluctuations and gravitational response', 'Microscopic derivation of a0 or the SM-1 law'])


def main():
    ap=argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args=ap.parse_args()
    if args.output_dir.exists():
        raise FileExistsError('Use a fresh output directory')
    out=calculate()
    out['code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output_dir.mkdir(parents=True)
    (args.output_dir/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2,allow_nan=False))
    if not out['passed']:
        raise SystemExit('CM-1E verification failed; results preserved')


if __name__=='__main__':
    main()
