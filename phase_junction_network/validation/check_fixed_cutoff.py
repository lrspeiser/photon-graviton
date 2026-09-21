#!/usr/bin/env python3
"""Compare Q_N H_full Q_N against rebuilding the completion after truncation.

Only the former measures state-cutoff error for a single fixed Hamiltonian.
Full and quick outputs are separate; quick never claims four-pair convergence.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'microscopic'))
import check_fermionic_multipair_vacuum as model


def maximum(matrix):
    return float(np.max(np.abs(matrix.data))) if matrix.nnz else 0.0


def fit_controls():
    cases = [([1.], [1.]), ([1., 1.], [1., 2.]),
             ([0., 1.], [1., 2.]), ([1., 2.], [-1., 2.]),
             ([1., np.nan], [1., 2.]), ([1., 2.], [1.])]
    rejected = 0
    for xs, ys in cases:
        try:
            model.logarithmic_fit(xs, ys)
        except ValueError:
            rejected += 1
    fit = model.logarithmic_fit([1., 1.25, 1.5], [3., 3./1.25**3, 3./1.5**3])
    assert rejected == len(cases)
    assert abs(fit['power'] + 3.) < 1.e-12
    assert abs(fit['prefactor'] - 3.) < 1.e-12
    return {'invalid_cases_rejected': rejected, 'known_power': fit['power']}


def run(quick=False):
    controls = fit_controls()
    frozen_path = ROOT / 'microscopic/fermionic_multipair_vacuum_results.json'
    original_bytes = frozen_path.read_bytes()
    archived = json.loads(original_bytes)
    reference = archived['full_multipair_result']['photon_gap']
    basis, counts, full, plus, minus, pairs = model.build_system(4)
    assert len(basis) == 6336
    assert counts == {0: 115, 1: 1484, 2: 3138, 3: 1484, 4: 115}
    pb, _, pure, *_ = model.build_system(0)
    bare = model.spectral_summary(pure, model.electric_values(pb), 40)['photon_gap']
    assert abs(bare - model.BARE_PHOTON_GAP) < 1.e-10
    rows = []
    for cutoff in ((1, 2) if quick else (1, 2, 3, 4)):
        selected = np.asarray([i for i, state in enumerate(basis)
                               if state[0].bit_count() <= cutoff], dtype=int)
        subbasis = [basis[i] for i in selected]
        fixed = full[selected][:, selected].tocsr()
        rb, _, rebuilt, *_ = model.build_system(cutoff)
        assert subbasis == rb
        defect = fixed - rebuilt
        off_diagonal = defect.copy()
        off_diagonal.setdiag(0)
        off_diagonal.eliminate_zeros()
        assert maximum(off_diagonal) < 1.e-12
        assert np.min(defect.diagonal().real) > -1.e-12
        result = model.spectral_summary(fixed, model.electric_values(subbasis), 48)
        energies, vectors = result['energies'], result['vectors']
        residual = float(np.max(np.linalg.norm(fixed @ vectors - vectors * energies, axis=0)))
        assert residual < 1.e-8
        gap = result['photon_gap']
        rows.append({'maximum_pairs': cutoff, 'dimension': len(subbasis),
                     'photon_gap': gap, 'photon_residue': result['photon_residue'],
                     'ground_energy': result['ground_energy'],
                     'loop_gap_shift': gap - bare,
                     'loop_shift_relative_error_to_archived_full':
                         (gap-reference)/(reference-bare),
                     'maximum_eigen_residual': residual,
                     'completion_boundary_diagonal_max': maximum(defect),
                     'completion_boundary_off_diagonal_max': maximum(off_diagonal)})
    # Variational ground-energy monotonicity is independent of the pole label.
    assert all(b['ground_energy'] <= a['ground_energy'] + 1.e-10
               for a, b in zip(rows, rows[1:]))
    assert rows[0]['completion_boundary_diagonal_max'] > 0
    convergence = None
    if not quick:
        assert abs(rows[-1]['photon_gap'] - reference) < 1.e-10
        assert abs(rows[1]['loop_shift_relative_error_to_archived_full']) < 0.002
        assert abs(rows[2]['loop_shift_relative_error_to_archived_full']) < 2.e-5
        assert rows[-1]['completion_boundary_diagonal_max'] < 1.e-12
        convergence = True
    assert frozen_path.read_bytes() == original_bytes
    return {'status': 'PASS', 'mode': 'quick' if quick else 'full',
            'method': 'Build the four-pair Hamiltonian once, then Q_N H_full Q_N',
            'fit_controls': controls, 'sector_counts': counts,
            'pure_gauge_gap_recomputed': bare, 'fixed_cutoff_rows': rows,
            'full_four_pair_convergence': convergence,
            'not_tested': ['full_four_pair_convergence'] if quick else [],
            'archived_full_result_sha256': hashlib.sha256(original_bytes).hexdigest(),
            'claim_boundary': 'Finite 2x2 spinless charged-endpoint patch only; '
                              'pole shifts are not infrared speed shifts.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = run(args.quick)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
