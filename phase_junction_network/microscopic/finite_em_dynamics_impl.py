"""Stage 6B: detuned finite electromagnetic dynamics.

This script tests the pure-gauge dynamical continuation of the finite spin-link
Hamiltonian used in Stage 6A.  It keeps three claims separate:

1. a reduced principal-axis chain built from the same finite link amplitudes;
2. a full 2^3 exact-cube anchor for the detuned ground state;
3. the still-open dynamical-matter and interacting-QED completion.
"""
from __future__ import annotations
import argparse
import functools
import itertools
import json
import math
from pathlib import Path
from typing import Any
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import eigsh
from finite_em_core import apply_steps_tuple, enumerate_zero_flux_component
TOL = 1e-09
FULL_U = (0.0, 0.1, 0.2)
FULL_V = (0.2, 0.4, 0.6)
FULL_SIZES = {1: (6, 8, 10), 2: (4, 5, 6)}
QUICK_U = (0.0, 0.2)
QUICK_V = (0.2, 0.6)
QUICK_SIZES = {1: (6, 8, 10), 2: (4, 6, 8)}

def ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [ready(item) for item in value]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    return value

def normalized_raise(m: int, spin: int) -> float:
    if m >= spin:
        return 0.0
    denominator = math.sqrt(spin * (spin + 1.0))
    return math.sqrt(spin * (spin + 1.0) - m * (m + 1.0)) / denominator

def normalized_lower(m: int, spin: int) -> float:
    if m <= -spin:
        return 0.0
    denominator = math.sqrt(spin * (spin + 1.0))
    return math.sqrt(spin * (spin + 1.0) - m * (m - 1.0)) / denominator

@functools.lru_cache(maxsize=None)
def axial_components(spin: int, lattice_size: int):
    """Build one reduced chain in its zero-total-flux sector."""
    values = range(-spin, spin + 1)
    states = [state for state in itertools.product(values, repeat=lattice_size) if sum(state) == 0]
    index = {state: position for position, state in enumerate(states)}
    dimension = len(states)
    electric = np.zeros(dimension, dtype=float)
    flippability = np.zeros(dimension, dtype=float)
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for row, state in enumerate(states):
        electric[row] = 0.5 * sum((value * value for value in state))
        mutable = list(state)
        for site in range(lattice_size):
            neighbor = (site + 1) % lattice_size
            forward = normalized_lower(state[site], spin) * normalized_raise(state[neighbor], spin)
            backward = normalized_raise(state[site], spin) * normalized_lower(state[neighbor], spin)
            flippability[row] += forward + backward
            if forward > 0.0:
                moved = mutable.copy()
                moved[site] -= 1
                moved[neighbor] += 1
                rows.append(row)
                cols.append(index[tuple(moved)])
                data.append(-forward)
            if backward > 0.0:
                moved = mutable.copy()
                moved[site] += 1
                moved[neighbor] -= 1
                rows.append(row)
                cols.append(index[tuple(moved)])
                data.append(-backward)
    hopping = coo_matrix((data, (rows, cols)), shape=(dimension, dimension), dtype=float).tocsr()
    translation = np.asarray([index[(state[-1],) + state[:-1]] for state in states], dtype=np.int64)
    state_array = np.asarray(states, dtype=np.int16)
    return (states, state_array, electric, flippability, hopping, translation)

def axial_hamiltonian(spin: int, lattice_size: int, u_over_t: float, v_over_t: float) -> csr_matrix:
    _, _, electric, flippability, hopping, _ = axial_components(spin, lattice_size)
    return diags(u_over_t * electric + v_over_t * flippability) + hopping

def resolve_translation(energies, vectors, permutation, lattice_size, tolerance=1e-07):
    inverse = np.empty_like(permutation)
    inverse[permutation] = np.arange(len(permutation))
    resolved = []
    start = 0
    while start < len(energies):
        stop = start + 1
        while stop < len(energies) and abs(energies[stop] - energies[start]) < tolerance:
            stop += 1
        block = vectors[:, start:stop]
        phases, rotation = np.linalg.eig(block.conj().T @ block[inverse, :])
        rotated = block @ rotation
        for column, phase in enumerate(phases):
            momentum = int(np.rint(np.angle(phase) % (2.0 * np.pi) * lattice_size / (2.0 * np.pi))) % lattice_size
            vector = rotated[:, column]
            vector /= np.linalg.norm(vector)
            resolved.append({'energy': float(energies[start]), 'momentum_index': momentum, 'translation_modulus': float(abs(phase)), 'vector': vector})
        start = stop
    return sorted(resolved, key=lambda row: row['energy'])

def axial_spectrum(spin: int, lattice_size: int, u_over_t: float, v_over_t: float, eigenpairs=24):
    hamiltonian = axial_hamiltonian(spin, lattice_size, u_over_t, v_over_t)
    _, states, _, _, _, translation = axial_components(spin, lattice_size)
    count = min(eigenpairs, hamiltonian.shape[0] - 2)
    energies, vectors = eigsh(hamiltonian, k=count, which='SA', tol=1e-09, maxiter=30000)
    order = np.argsort(energies)
    resolved = resolve_translation(energies[order], vectors[:, order], translation, lattice_size)
    ground = min(resolved, key=lambda row: row['energy'])
    positive = [row for row in resolved if row['momentum_index'] == 1 and row['energy'] - ground['energy'] > 1e-08]
    negative = [row for row in resolved if row['momentum_index'] == lattice_size - 1 and row['energy'] - ground['energy'] > 1e-08]
    if not positive or not negative:
        raise RuntimeError(f'missing signed momentum branch for S={spin}, L={lattice_size}, u={u_over_t}, v={v_over_t}')
    plus = min(positive, key=lambda row: row['energy'])
    minus = min(negative, key=lambda row: row['energy'])
    plus_gap = plus['energy'] - ground['energy']
    minus_gap = minus['energy'] - ground['energy']
    distinct_plus = sorted({round(row['energy'] - ground['energy'], 10) for row in positive if row['energy'] - ground['energy'] > 1e-08})
    next_gap = distinct_plus[1] if len(distinct_plus) > 1 else None
    lattice_momentum = 2.0 * math.sin(math.pi / lattice_size)
    return {'spin': spin, 'lattice_size': lattice_size, 'hilbert_dimension': hamiltonian.shape[0], 'u_over_t': u_over_t, 'v_over_t': v_over_t, 'lattice_momentum': lattice_momentum, 'positive_momentum_gap': plus_gap, 'negative_momentum_gap': minus_gap, 'signed_momentum_split': abs(plus_gap - minus_gap), 'single_polarization_gap': 0.5 * (plus_gap + minus_gap), 'constructed_copy_count': 2, 'next_same_momentum_gap': next_gap, 'next_branch_over_photon_gap': float(next_gap / (0.5 * (plus_gap + minus_gap))) if next_gap is not None else None, 'ground_vector': ground['vector'], 'photon_vector': plus['vector'], 'states': states, 'resolved': resolved}

def public_spectrum_row(row):
    keys = ('spin', 'lattice_size', 'hilbert_dimension', 'u_over_t', 'v_over_t', 'lattice_momentum', 'positive_momentum_gap', 'negative_momentum_gap', 'signed_momentum_split', 'single_polarization_gap', 'constructed_copy_count', 'next_same_momentum_gap', 'next_branch_over_photon_gap')
    return {key: row[key] for key in keys}

def fit_branch(rows):
    momenta = np.asarray([row['lattice_momentum'] for row in rows], dtype=float)
    gaps = np.asarray([row['single_polarization_gap'] for row in rows], dtype=float)
    power, log_prefactor = np.polyfit(np.log(momenta), np.log(gaps), 1)
    design = np.column_stack((momenta, momenta ** 3))
    speed, cubic = np.linalg.lstsq(design, gaps, rcond=None)[0]
    prediction = design @ np.asarray([speed, cubic])
    ratios = gaps / momenta
    return {'gap_proportional_to_lattice_momentum_power': float(power), 'log_fit_prefactor': float(np.exp(log_prefactor)), 'linear_speed_c_gamma': float(speed), 'cubic_cutoff_coefficient': float(cubic), 'maximum_relative_linear_plus_cubic_fit_residual': float(np.max(np.abs(prediction - gaps) / gaps)), 'gap_over_momentum_coefficient_of_variation': float(np.std(ratios) / np.mean(ratios)), 'minimum_next_branch_over_photon_gap': float(min(values)) if (values := [row['next_branch_over_photon_gap'] for row in rows if row['next_branch_over_photon_gap'] is not None]) else None, 'maximum_signed_momentum_split': float(max((row['signed_momentum_split'] for row in rows)))}

def spectral_residue(spin: int, lattice_size: int, u_over_t: float, v_over_t: float):
    row = axial_spectrum(spin, lattice_size, u_over_t, v_over_t, eigenpairs=40)
    states = row['states'].astype(float)
    phase = np.exp(-2j * np.pi * np.arange(lattice_size) / lattice_size)
    electric_operator = states @ phase
    scalar_operator = states * states @ phase
    ground = row['ground_vector']
    photon_gap = row['single_polarization_gap']
    def weights(operator):
        vector = operator * ground
        norm = float(np.vdot(vector, vector).real)
        entries = []
        for state in row['resolved']:
            if state['momentum_index'] not in (1, lattice_size - 1):
                continue
            residue = float(abs(np.vdot(state['vector'], vector)) ** 2)
            entries.append({'gap': state['energy'] - row['resolved'][0]['energy'], 'momentum_index': state['momentum_index'], 'residue': residue, 'fraction': residue / norm if norm > 0.0 else 0.0})
        return (norm, sorted(entries, key=lambda entry: entry['residue'], reverse=True))
    electric_norm, electric = weights(electric_operator)
    scalar_norm, scalar = weights(scalar_operator)
    electric_low_fraction = sum((entry['fraction'] for entry in electric if abs(entry['gap'] - photon_gap) < 1e-07))
    scalar_low_fraction = sum((entry['fraction'] for entry in scalar if abs(entry['gap'] - photon_gap) < 1e-07))
    scalar_candidates = [entry for entry in scalar if entry['fraction'] > 0.0001]
    scalar_gap = min((entry['gap'] for entry in scalar_candidates))
    return {'spin': spin, 'lattice_size': lattice_size, 'u_over_t': u_over_t, 'v_over_t': v_over_t, 'photon_gap': photon_gap, 'electric_operator_norm': electric_norm, 'electric_residue_in_first_photon_branch': electric_low_fraction, 'scalar_operator_norm': scalar_norm, 'scalar_residue_in_first_photon_branch': scalar_low_fraction, 'first_scalar_residue_gap': scalar_gap, 'scalar_gap_over_photon_gap': scalar_gap / photon_gap, 'top_electric_residues': electric[:5], 'top_scalar_residues': scalar[:5]}

def run_detuned_scan(quick: bool):
    u_values = QUICK_U if quick else FULL_U
    v_values = QUICK_V if quick else FULL_V
    sizes = QUICK_SIZES if quick else FULL_SIZES
    rows = []
    fits = {}
    for spin in (1, 2):
        for u_over_t in u_values:
            for v_over_t in v_values:
                block = [axial_spectrum(spin, lattice_size, u_over_t, v_over_t, eigenpairs=8) for lattice_size in sizes[spin]]
                rows.extend((public_spectrum_row(row) for row in block))
                fits[f'S{spin}:u={u_over_t:.2f}:v={v_over_t:.2f}'] = fit_branch(block)
    controls = {}
    for spin in (1, 2):
        block = [axial_spectrum(spin, lattice_size, 0.0, 1.0, eigenpairs=8) for lattice_size in sizes[spin]]
        controls[f'S{spin}:RK'] = {'rows': [public_spectrum_row(row) for row in block], 'fit': fit_branch(block)}
    residue_sizes = {1: 10, 2: 6 if quick else 8}
    residues = [spectral_residue(spin, residue_sizes[spin], 0.1, 0.4) for spin in (1, 2)]
    return (rows, fits, controls, residues)

def exact_cube_components():
    basis = enumerate_zero_flux_component(L=2, spin=1)
    dimension = len(basis.states)
    rows: list[int] = []
    cols: list[int] = []
    degree = np.zeros(dimension, dtype=float)
    electric = np.zeros(dimension, dtype=float)
    for row, state in enumerate(basis.states):
        electric[row] = 0.5 * sum((value * value for value in state))
        count = 0
        for steps in basis.plaquettes:
            for sign in (-1, 1):
                moved = apply_steps_tuple(state, steps, sign, 1)
                if moved is None:
                    continue
                rows.append(row)
                cols.append(basis.index[moved])
                count += 1
        degree[row] = count
    adjacency = coo_matrix((np.ones(len(rows), dtype=float), (rows, cols)), shape=(dimension, dimension)).tocsr()
    states = np.asarray(basis.states, dtype=np.int8).reshape(dimension, 2, 2, 2, 3)
    phase = np.asarray([1.0, -1.0])[:, None, None]
    operators = {'transverse_y': np.sum(states[..., 1] * phase, axis=(1, 2, 3)), 'transverse_z': np.sum(states[..., 2] * phase, axis=(1, 2, 3)), 'longitudinal_x': np.sum(states[..., 0] * phase, axis=(1, 2, 3)), 'scalar_electric_energy': np.sum(states * states * phase[..., None], axis=(1, 2, 3, 4))}
    return (basis, electric, degree, adjacency, operators)

def exact_cube_scan(quick: bool):
    basis, electric, degree, adjacency, operators = exact_cube_components()
    points = [(0.0, 1.0), (0.1, 0.4)] if quick else [(0.0, 1.0), (0.0, 0.6), (0.1, 0.4), (0.2, 0.2)]
    uniform = np.ones(len(basis.states), dtype=float) / math.sqrt(len(basis.states))
    results = []
    for u_over_t, v_over_t in points:
        hamiltonian = diags(u_over_t * electric + v_over_t * degree) - adjacency
        energies, vectors = eigsh(hamiltonian, k=8, which='SA', tol=1e-09, maxiter=20000, v0=uniform)
        order = np.argsort(energies)
        energies = energies[order]
        vectors = vectors[:, order]
        ground = vectors[:, 0]
        gaps = energies - energies[0]
        first_gap = gaps[1]
        low_mask = np.isclose(gaps, first_gap, rtol=1e-08, atol=1e-08)
        operator_results = {}
        for name, operator in operators.items():
            created = operator * ground
            norm = float(np.vdot(created, created).real)
            residues = np.abs(vectors.conj().T @ created) ** 2
            operator_results[name] = {'norm': norm, 'low_multiplet_residue_fraction': float(np.sum(residues[low_mask]) / norm) if norm > 0.0 else 0.0, 'first_eight_residues': [float(value) for value in residues]}
        results.append({'u_over_t': u_over_t, 'v_over_t': v_over_t, 'gauge_reduced_dimension': len(basis.states), 'lowest_gaps': [float(value) for value in gaps], 'first_gap_multiplicity_in_computed_window': int(np.count_nonzero(low_mask)), 'ground_overlap_with_RK_uniform_state': float(abs(np.vdot(uniform, ground)) ** 2), 'ground_electric_energy_expectation': float(np.sum(np.abs(ground) ** 2 * electric)), 'ground_flippability_expectation': float(np.sum(np.abs(ground) ** 2 * degree)), 'operators': operator_results})
    return results

def full_three_dimensional_symbol():
    tested = []
    maximum_longitudinal = 0.0
    maximum_split = 0.0
    for lattice_size in (12, 16, 24, 32, 48, 64):
        for mode in ((1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 0)):
            momentum = 2.0 * np.sin(np.pi * np.asarray(mode, dtype=float) / lattice_size)
            norm = float(np.linalg.norm(momentum))
            unit = momentum / norm
            projector = np.eye(3) - np.outer(unit, unit)
            eigenvalues = np.linalg.eigvalsh(norm * norm * projector)
            positive = eigenvalues[eigenvalues > 1e-12]
            longitudinal = eigenvalues[0]
            maximum_longitudinal = max(maximum_longitudinal, abs(float(longitudinal)))
            maximum_split = max(maximum_split, float(abs(positive[1] - positive[0]) / (norm * norm)))
            tested.append({'lattice_size': lattice_size, 'mode': list(mode), 'lattice_momentum_norm': norm, 'kernel_eigenvalues': [float(value) for value in eigenvalues], 'transverse_rank': int(len(positive))})
    anisotropy = []
    for lattice_size in (12, 16, 24, 32, 48, 64, 96, 128):
        axis = np.linalg.norm(2.0 * np.sin(np.pi * np.asarray((3, 0, 0)) / lattice_size))
        mixed = np.linalg.norm(2.0 * np.sin(np.pi * np.asarray((2, 2, 1)) / lattice_size))
        relative = abs(axis - mixed) / (0.5 * (axis + mixed))
        anisotropy.append({'lattice_size': lattice_size, 'relative_anisotropy': float(relative)})
    slope, intercept = np.polyfit(np.log([row['lattice_size'] for row in anisotropy]), np.log([row['relative_anisotropy'] for row in anisotropy]), 1)
    return {'method': 'Gauss-reduced cubic Maxwell symbol inferred from the exact axial transfer blocks', 'tested_modes': tested, 'maximum_longitudinal_eigenvalue_abs': maximum_longitudinal, 'maximum_relative_transverse_split': maximum_split, 'equal_continuum_norm_anisotropy': anisotropy, 'anisotropy_proportional_to_L_power': float(slope), 'anisotropy_prefactor': float(np.exp(intercept))}

def build_report(quick: bool):
    rows, fits, controls, residues = run_detuned_scan(quick)
    cube = exact_cube_scan(quick)
    symbol = full_three_dimensional_symbol()
    fit_values = list(fits.values())
    powers = [row['gap_proportional_to_lattice_momentum_power'] for row in fit_values]
    residuals = [row['maximum_relative_linear_plus_cubic_fit_residual'] for row in fit_values]
    separations = [row['minimum_next_branch_over_photon_gap'] for row in fit_values if row['minimum_next_branch_over_photon_gap'] is not None]
    signed_splits = [row['maximum_signed_momentum_split'] for row in fit_values]
    rk_powers = [entry['fit']['gap_proportional_to_lattice_momentum_power'] for entry in controls.values()]
    acceptance = {'finite_detuned_rectangle_not_single_point': len(fits) >= (4 if quick else 18), 'all_detuned_gap_powers_between_0p94_and_1p16': min(powers) >= 0.94 and max(powers) <= 1.16, 'all_linear_plus_cubic_fit_residuals_below_0p002': max(residuals) < 0.002, 'signed_momentum_pair_degenerate': max(signed_splits) < 1e-08, 'two_identical_copies_constructed_not_3d_degeneracy_measurement': all((row['constructed_copy_count'] == 2 for row in rows)), 'electric_operator_residue_above_0p98': min((row['electric_residue_in_first_photon_branch'] for row in residues)) > 0.98, 'scalar_residue_in_photon_branch_below_1e_minus_8': max((row['scalar_residue_in_first_photon_branch'] for row in residues)) < 1e-08, 'scalar_branch_separated_by_factor_gt_2': min((row['scalar_gap_over_photon_gap'] for row in residues)) > 2.0, 'rk_point_is_z2_not_z1_control': min(rk_powers) > 1.75, 'exact_cube_longitudinal_operator_null': max((point['operators']['longitudinal_x']['norm'] for point in cube)) < 1e-12, 'exact_cube_scalar_low_multiplet_residue_below_1e_minus_8': max((point['operators']['scalar_electric_energy']['low_multiplet_residue_fraction'] for point in cube)) < 1e-08, 'exact_cube_detuned_ground_continuously_overlaps_rk_state': min((point['ground_overlap_with_RK_uniform_state'] for point in cube)) > 0.6, 'full_symbol_has_two_transverse_and_zero_longitudinal': symbol['maximum_longitudinal_eigenvalue_abs'] < 1e-12 and symbol['maximum_relative_transverse_split'] < 1e-12 and all((row['transverse_rank'] == 2 for row in symbol['tested_modes'])), 'cubic_anisotropy_vanishes_as_L_minus_2': abs(symbol['anisotropy_proportional_to_L_power'] + 2.0) < 0.1}
    stage_pass = all(acceptance.values())
    return {'module': 'phase_junction_network/microscopic/check_finite_em_dynamics.py', 'status': 'Stage 6B reduced-chain PASS; finite-spin 3D phase, interacting matter and QED remain open' if stage_pass else 'Stage 6B FAIL', 'run_mode': 'quick' if quick else 'full', 'model': {'parent_hamiltonian': 'H_A/t = (u/t)/2 sum E^2 - sum_p(W_p+W_p^dagger) + (v/t) sum_p D_p', 'axial_transverse_block': 'H_perp/t = (u/t)/2 sum_z E_z^2 - sum_z(W_z+W_z^dagger) + (v/t) sum_z D_z, W_z=U_{z+1}U_z^dagger', 'gauss_sector': 'zero total chain flux; not a proof of 3D Gauss reduction', 'polarizations': 'two identical reduced-chain copies; 3D embedding not established', 'rk_control': 'u/t=0, v/t=1', 'detuned_rectangle': {'u_over_t': list(QUICK_U if quick else FULL_U), 'v_over_t': list(QUICK_V if quick else FULL_V)}}, 'detuned_spectrum_rows': rows, 'detuned_fits': fits, 'rk_z2_controls': controls, 'spectral_residue_controls': residues, 'exact_periodic_cube_anchor': cube, 'full_three_dimensional_symbol': symbol, 'summary': {'minimum_detuned_power': min(powers), 'maximum_detuned_power': max(powers), 'maximum_linear_plus_cubic_fit_residual': max(residuals), 'minimum_higher_branch_separation_where_resolved': min(separations) if separations else None, 'maximum_signed_momentum_split': max(signed_splits), 'minimum_electric_photon_residue': min((row['electric_residue_in_first_photon_branch'] for row in residues)), 'maximum_scalar_residue_in_photon_branch': max((row['scalar_residue_in_photon_branch'] for row in residues)), 'minimum_scalar_gap_over_photon_gap': min((row['scalar_gap_over_photon_gap'] for row in residues)), 'rk_control_powers': rk_powers, 'minimum_exact_cube_rk_overlap': min((point['ground_overlap_with_RK_uniform_state'] for point in cube))}, 'acceptance': acceptance, 'stage_pass': stage_pass, 'physical_3d_phase_established': False, 'claim_boundary': {'established': ['an exactly diagonalized reduced chain using Stage-6A finite-link amplitudes', 'a finite two-dimensional bare-coupling rectangle with z approximately 1 for spin 1 and spin 2', 'two degenerate transverse copies with positive electric-field spectral residue', 'negligible scalar residue in the photon branch and a separated scalar control branch', 'the RK point is a z approximately 2 negative control rather than the relativistic phase', 'an exact 2^3 detuned full-cube anchor with continuous RK overlap and no longitudinal operator', 'a separately assumed Gaussian cubic symbol has rank two and L^-2 anisotropy recovery'], 'not_established': ['the complete infinite-volume three-dimensional many-body spectrum away from RK', 'dynamical issue-4 charged matter inside the detuned phase', 'finite interacting Ward identities, vacuum polarization, or charge renormalization', 'a regulator-independent QED fixed point or precision Lorentz bounds', 'a common photon, companion, matter, and graviton cone', 'the physical normalization of c_gamma or its relation to alpha and G']}, 'next_gate': 'Embed the |q|=1 issue-4 endpoint in the detuned finite photon interval and test current continuity, a finite Ward identity, photon dressing, and the issue-7 reversible conversion vertex using actual photon eigenmodes.'}

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('phase_junction_network/microscopic/finite_em_dynamics_revalidated.json'))
    parser.add_argument('--quick', action='store_true')
    args = parser.parse_args()
    report = build_report(bool(args.quick))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(ready(report), indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({'status': report['status'], 'stage_pass': report['stage_pass'], 'summary': report['summary'], 'failed_acceptance': [key for key, value in report['acceptance'].items() if not value]}, indent=2, sort_keys=True))
    return 0 if report['stage_pass'] else 1
if __name__ == '__main__':
    raise SystemExit(main())

# Validation scope: chains are reduced models, not established exact 3D sectors.
