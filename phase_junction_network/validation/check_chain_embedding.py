#!/usr/bin/env python3
"""Exact finite-link sheet-embedding leakage, plus a separate Gaussian control.

The full H V is computed by enumerating every local move from the embedded
states. This is exact even though the entire 3D Hilbert space is not enumerated.
A failed natural embedding does not rule out every possible dressed embedding.
The Gaussian control is explicitly conditional, not a finite-spin phase proof.
"""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path
import sys
import numpy as np
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'microscopic'))
from finite_em_core import all_plaquettes, gauss_from_flat
from finite_em_dynamics_impl import axial_hamiltonian


def sheet_embedding(L: int, u: float = 0.1, v: float = 0.4) -> dict:
    """V|a,b>: E_x=0, E_y(x,y,z)=a_x, E_z(x,y,z)=b_x.

    Both a and b are spin-1 zero-sum chains. V is an isometry into the full
    periodic Gauss-law Hilbert space. K=1 sets the energy unit, not a fit.
    """
    chains = [s for s in itertools.product((-1, 0, 1), repeat=L) if sum(s) == 0]
    states = []
    for a, b in itertools.product(chains, repeat=2):
        flux = np.zeros((L, L, L, 3), dtype=np.int8)
        flux[..., 1] = np.asarray(a)[:, None, None]
        flux[..., 2] = np.asarray(b)[:, None, None]
        states.append(flux.ravel())
    index = {s.tobytes(): i for i, s in enumerate(states)}
    assert len(index) == len(states)
    assert all(np.max(np.abs(gauss_from_flat(s, L))) == 0 for s in states)
    plaquettes = all_plaquettes(L)
    projected = np.zeros((len(states), len(states)))
    outside_index = {}
    rows, columns, amplitudes = [], [], []
    degrees = []
    for column, state in enumerate(states):
        degree = 0
        for steps in plaquettes:
            indices = np.asarray([i for i, _ in steps], dtype=int)
            shifts = np.asarray([d for _, d in steps], dtype=np.int8)
            for orientation in (-1, 1):
                values = state[indices] + orientation * shifts
                if np.any(values < -1) or np.any(values > 1):
                    continue
                moved = state.copy()
                moved[indices] = values
                key = moved.tobytes()
                degree += 1
                target = index.get(key)
                if target is None:
                    target = outside_index.setdefault(key, len(outside_index))
                    rows.append(target)
                    columns.append(column)
                    amplitudes.append(-1.)
                else:
                    projected[target, column] -= 1.
        degrees.append(degree)
        projected[column, column] += 0.5 * u * float(state @ state.astype(float)) + v * degree
    leakage = coo_matrix((amplitudes, (rows, columns)),
                         shape=(len(outside_index), len(states))).tocsr()
    gram = (leakage.T @ leakage).toarray()
    zero = index[np.zeros(3*L**3, dtype=np.int8).tobytes()]
    vacuum_leakage = float(np.sqrt(gram[zero, zero]))
    # The test acts on all embedded columns, not only on the vacuum witness.
    spectral_leakage = float(np.sqrt(np.max(np.linalg.eigvalsh(gram))))
    off = projected - np.diag(np.diag(projected))
    single = axial_hamiltonian(1, L, u, v).toarray()
    two_chain = np.kron(single, np.eye(len(chains))) + np.kron(np.eye(len(chains)), single)
    mismatch = projected - two_chain
    mismatch -= np.trace(mismatch) / len(states) * np.eye(len(states))
    assert abs(vacuum_leakage**2 / L**3 - 6.) < 1.e-10
    assert spectral_leakage > 1.
    assert np.max(np.abs(off)) < 1.e-12
    return {'L': L, 'spatial_dimension': 3, 'embedded_dimension': len(states),
            'action_support_outside_embedding': len(outside_index),
            'gauss_residual': 0, 'isometry_residual': 0,
            'vacuum_leakage_norm': vacuum_leakage,
            'vacuum_leakage_squared_per_site': vacuum_leakage**2 / L**3,
            'operator_leakage_norm': spectral_leakage,
            'projected_off_diagonal_max': float(np.max(np.abs(off))),
            'two_chain_off_diagonal_max': float(np.max(np.abs(two_chain-np.diag(np.diag(two_chain))))),
            'intertwining_residual_after_energy_shift': float(np.linalg.norm(mismatch, ord=2)),
            'exact_invariant_embedding': False}


def curl(field):
    def d(component, axis):
        return np.roll(field[..., component], -1, axis=axis) - field[..., component]
    return np.stack((d(2, 1)-d(1, 2), d(0, 2)-d(2, 0), d(1, 0)-d(0, 1)), axis=-1)


def curl_adjoint(field):
    def dt(component, axis):
        return np.roll(field[..., component], 1, axis=axis)-field[..., component]
    return np.stack((dt(1, 2)-dt(2, 1), dt(2, 0)-dt(0, 2), dt(0, 1)-dt(1, 0)), axis=-1)


def gaussian_control(L, mode):
    """Full 3D real-space curl-adjoint-curl, no transverse projector inserted."""
    k = 2*np.pi*np.asarray(mode, dtype=float)/L
    symbol = np.exp(1j*k)-1.
    # Gauss says d^dagger A=0. An SVD merely finds source polarizations; the
    # evolution is the full real-space local curl operator on all components.
    _, _, vh = np.linalg.svd(symbol.conj()[None, :], full_matrices=True)
    polarizations = vh.conj().T[:, 1:]
    coordinates = np.indices((L, L, L))
    wave = np.exp(1j*np.einsum('i,ixyz->xyz', k, coordinates))/np.sqrt(L**3)
    q2 = float(np.vdot(symbol, symbol).real)
    residuals = []
    for pol in polarizations.T:
        field = wave[..., None]*pol
        evolved = curl_adjoint(curl(field))
        residuals.append(float(np.linalg.norm(evolved-q2*field)))
    assert max(residuals) < 2.e-12
    return {'L': L, 'mode': list(mode), 'frequency': float(np.sqrt(q2)),
            'continuum_k': float(np.linalg.norm(k)),
            'relative_speed_error': float(np.sqrt(q2)/np.linalg.norm(k)-1),
            'embedding_leakage_residual': max(residuals)}


def run(quick=False):
    sheets = [sheet_embedding(L) for L in ((2, 3) if quick else (2, 3, 4))]
    sizes = (8, 12, 16) if quick else (8, 12, 16, 24, 32, 48)
    gaussian = [gaussian_control(L, mode) for L in sizes for mode in ((1,0,0), (1,1,1), (3,0,0), (2,2,1))]
    axis = [r for r in gaussian if r['mode'] == [1,0,0]]
    exponent = float(np.polyfit(np.log([r['L'] for r in axis]),
                                np.log([-r['relative_speed_error'] for r in axis]), 1)[0])
    assert -2.05 < exponent < -1.95
    angular = []
    for L in sizes:
        a = next(r for r in gaussian if r['L'] == L and r['mode'] == [3,0,0])
        b = next(r for r in gaussian if r['L'] == L and r['mode'] == [2,2,1])
        angular.append({'L':L, 'relative_equal_norm_direction_split':abs(a['frequency']-b['frequency'])/a['continuum_k']})
    assert angular[-1]['relative_equal_norm_direction_split'] < angular[0]['relative_equal_norm_direction_split']
    return {'audit_pass': True, 'finite_spin_exact_reduction_pass': False,
            'mode': 'quick' if quick else 'full', 'finite_link_embedding': sheets,
            'gaussian_3d_control': gaussian, 'gaussian_cutoff_error_power': exponent,
            'gaussian_equal_norm_direction_control': angular,
            'decision': 'Relabel Stage 6B chains as reduced models. The natural sheet '
                        'embedding has nonzero exact 3D leakage. Only the separately '
                        'specified Gaussian parent has the demonstrated exact plane-wave reduction.',
            'not_established': ['a controlled dressed embedding for the interacting spin-1 parent',
                                '3D finite-spin dynamical phase or its infinite-volume speed']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = run(args.quick)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps(report, indent=2, sort_keys=True))
