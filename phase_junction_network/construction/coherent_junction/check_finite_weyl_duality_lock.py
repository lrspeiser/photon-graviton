#!/usr/bin/env python3
"""Stage 5J: finite Weyl self-duality lock for the doubled bosonic candidate.

Stage 5I showed that the common chiral pole is physically viable only when it
acts as the symplectic structure of a doubled positive bosonic system. This
script implements the corresponding self-duality constraints with finite Weyl
operators on a local one-dimensional propagation block.

The exact lock has the expected finite physical-band dimension. The projected
positive magnetic Hamiltonian is then tested under local-dimension refinement
and spatial-size scaling. The outcome is a controlled rejection of the naive
*energetic* lock as the continuum regulator: its first gap does not retain the
required z=1 scaling.
"""
from __future__ import annotations

import argparse
import builtins
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, eigsh

HERE = Path(__file__).resolve().parent
import sys
sys.path.insert(0, str(HERE))
from check_common_chiral_generator import build_report as chiral_report

FULL_CASES = ((3, 3), (5, 3), (7, 3), (11, 3), (3, 5))
QUICK_CASES = ((3, 3), (5, 3), (3, 5))


def ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [ready(item) for item in value]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return value


def inverse_mod(value: int, prime: int) -> int:
    return builtins.pow(int(value) % prime, -1, prime)


def reduced_self_duality_form(prime: int, lattice_size: int) -> np.ndarray:
    inverse_two = inverse_mod(2, prime)
    derivative = np.zeros((lattice_size, lattice_size), dtype=int)
    for site in range(lattice_size):
        derivative[site, (site + 1) % lattice_size] = inverse_two
        derivative[site, (site - 1) % lattice_size] = (-inverse_two) % prime
    polarization_rotation = np.asarray([[0, -1], [1, 0]], dtype=int) % prime
    curl = np.kron(polarization_rotation, derivative) % prime
    duality_rotation = np.asarray([[0, 1], [-1, 0]], dtype=int) % prime
    full = (inverse_two * np.kron(duality_rotation, curl)) % prime

    columns = []
    for field in range(4):
        for site in range(1, lattice_size):
            item = np.zeros(4 * lattice_size, dtype=int)
            item[field * lattice_size + site] = 1
            columns.append(item)
    embedding = np.column_stack(columns)
    return (embedding.T @ full @ embedding) % prime


def symplectic_basis(form: np.ndarray, prime: int):
    dimension = form.shape[0]
    vectors = [np.eye(dimension, dtype=int)[:, index] for index in range(dimension)]
    output = []
    while vectors:
        first = vectors.pop(0) % prime
        partner_index = None
        for index, candidate in enumerate(vectors):
            if int(first @ form @ candidate) % prime:
                partner_index = index
                break
        if partner_index is None:
            raise ValueError("The reduced self-duality form is degenerate.")
        second = vectors.pop(partner_index) % prime
        pairing = int(first @ form @ second) % prime
        second = (second * inverse_mod(pairing, prime)) % prime
        remaining = []
        for vector in vectors:
            alpha = int(first @ form @ vector) % prime
            beta = int(second @ form @ vector) % prime
            remaining.append((vector - alpha * second + beta * first) % prime)
        output.extend((first, second))
        vectors = remaining
    transformation = np.column_stack(output) % prime
    canonical = (transformation.T @ form @ transformation) % prime
    return transformation, canonical


def modular_rank(matrix: np.ndarray, prime: int) -> int:
    work = np.asarray(matrix, dtype=int).copy() % prime
    row = 0
    for column in range(work.shape[1]):
        pivot = next(
            (index for index in range(row, work.shape[0]) if work[index, column]),
            None,
        )
        if pivot is None:
            continue
        work[[row, pivot]] = work[[pivot, row]]
        work[row] = (work[row] * inverse_mod(work[row, column], prime)) % prime
        for index in range(work.shape[0]):
            if index != row and work[index, column]:
                work[index] = (
                    work[index] - work[index, column] * work[row]
                ) % prime
        row += 1
        if row == work.shape[0]:
            break
    return row


def pair_ground_band(prime: int):
    dimension = prime * prime
    states = np.arange(dimension, dtype=np.int64)
    first = states % prime
    second = (states // prime) % prime
    omega = np.exp(2j * np.pi / prime)
    form = np.asarray([[0, 1], [-1, 0]], dtype=int) % prime
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    coordinates = np.column_stack((first, second))
    for constraint in range(2):
        phase_vector = (-form[constraint]) % prime
        for state in range(dimension):
            phase = omega ** (int(phase_vector @ coordinates[state]) % prime)
            target = coordinates[state].copy()
            target[constraint] = (target[constraint] + 1) % prime
            target_index = int(target[0] + prime * target[1])
            hamiltonian[target_index, state] -= phase
            hamiltonian[state, target_index] -= np.conj(phase)
        hamiltonian += 2.0 * np.eye(dimension)
    hamiltonian = 0.5 * (hamiltonian + hamiltonian.conj().T)
    eigenvalues, eigenvectors = eigh(hamiltonian)
    mask = np.isclose(eigenvalues, eigenvalues[0], atol=1.0e-10)
    return {
        "eigenvalues": eigenvalues,
        "ground_basis": eigenvectors[:, mask],
        "ground_dimension": int(mask.sum()),
        "lock_gap": float(eigenvalues[int(mask.sum())] - eigenvalues[0]),
    }


def projected_phase_operator(
    prime: int,
    ground_basis: np.ndarray,
    first_coefficient: int,
    second_coefficient: int,
):
    states = np.arange(prime * prime, dtype=np.int64)
    first = states % prime
    second = (states // prime) % prime
    phases = np.exp(
        2j
        * np.pi
        / prime
        * ((first_coefficient * first + second_coefficient * second) % prime)
    )
    return ground_basis.conj().T @ (phases[:, None] * ground_basis)


def magnetic_terms(prime: int, lattice_size: int):
    form = reduced_self_duality_form(prime, lattice_size)
    transformation, canonical = symplectic_basis(form, prime)
    pair = pair_ground_band(prime)
    pair_count = form.shape[0] // 2
    identity = np.eye(prime, dtype=complex)

    expected = np.kron(
        np.eye(pair_count, dtype=int),
        np.asarray([[0, 1], [-1, 0]], dtype=int),
    ) % prime
    if np.max((canonical - expected) % prime):
        raise AssertionError("Symplectic decomposition did not reach canonical form.")

    terms = []
    for field in range(4):
        def coordinate(site: int):
            if site == 0:
                return None
            return field * (lattice_size - 1) + site - 1

        for site in range(lattice_size):
            phase = np.zeros(form.shape[0], dtype=int)
            left = coordinate(site)
            right = coordinate((site + 1) % lattice_size)
            if left is not None:
                phase[left] -= 1
            if right is not None:
                phase[right] += 1
            transformed = (transformation.T @ phase) % prime
            operators = []
            for pair_index in range(pair_count):
                first = int(transformed[2 * pair_index])
                second = int(transformed[2 * pair_index + 1])
                if first == 0 and second == 0:
                    operators.append(None)
                else:
                    operators.append(
                        projected_phase_operator(
                            prime,
                            pair["ground_basis"],
                            first,
                            second,
                        )
                    )
            terms.append(operators)
    return form, pair, terms


def apply_axis(tensor: np.ndarray, operator: np.ndarray, axis: int):
    output = np.tensordot(operator, tensor, axes=(1, axis))
    return np.moveaxis(output, 0, axis)


def projected_hamiltonian(prime: int, lattice_size: int):
    form, pair, terms = magnetic_terms(prime, lattice_size)
    pair_count = form.shape[0] // 2
    dimension = prime**pair_count

    def matvec(vector):
        tensor = vector.reshape((prime,) * pair_count)
        output = np.zeros_like(tensor, dtype=complex)
        for operators in terms:
            forward = tensor
            backward = tensor
            for axis, operator in enumerate(operators):
                if operator is None:
                    continue
                forward = apply_axis(forward, operator, axis)
                backward = apply_axis(backward, operator.conj().T, axis)
            output += tensor - 0.5 * (forward + backward)
        return output.reshape(-1)

    return (
        LinearOperator(
            (dimension, dimension),
            matvec=matvec,
            rmatvec=matvec,
            dtype=complex,
        ),
        form,
        pair,
    )


def spectrum(prime: int, lattice_size: int, static_coefficient: float, residue: float):
    hamiltonian, form, pair = projected_hamiltonian(prime, lattice_size)
    count = min(10, hamiltonian.shape[0] - 2)
    eigenvalues = np.sort(
        eigsh(
            hamiltonian,
            k=count,
            which="SA",
            tol=1.0e-9,
            maxiter=60000,
            return_eigenvectors=False,
        )
    )
    gaps = eigenvalues - eigenvalues[0]
    first_gap = float(gaps[1])
    first_multiplicity = int(np.count_nonzero(np.isclose(gaps, first_gap, atol=1.0e-7)))
    canonical_scale = prime**3 / (2.0 * np.pi) ** 2
    scaled_gap = first_gap * canonical_scale
    physical_gap = static_coefficient * scaled_gap
    central_momentum = abs(np.sin(2.0 * np.pi / lattice_size))
    expected = static_coefficient / residue * central_momentum
    term_count = 4 * lattice_size
    mixing_bound = (
        static_coefficient
        * canonical_scale
        * 2.0
        * term_count
        / pair["lock_gap"]
    )
    return {
        "prime": prime,
        "lattice_size": lattice_size,
        "constraint_dimension": int(form.shape[0]),
        "constraint_rank_mod_p": modular_rank(form, prime),
        "symplectic_pair_count": int(form.shape[0] // 2),
        "projected_hilbert_dimension": int(hamiltonian.shape[0]),
        "expected_projected_dimension": int(prime ** (form.shape[0] // 2)),
        "pair_lock_ground_dimension": pair["ground_dimension"],
        "pair_lock_gap": pair["lock_gap"],
        "lowest_projected_energies": [float(value) for value in eigenvalues],
        "first_dimensionless_gap": first_gap,
        "first_gap_multiplicity": first_multiplicity,
        "canonical_p_cubed_scaled_gap": scaled_gap,
        "physical_gap": physical_gap,
        "central_lattice_momentum": central_momentum,
        "duality_target_frequency": expected,
        "relative_frequency_error": physical_gap / expected - 1.0,
        "rigorous_lock_mixing_bound": mixing_bound,
        "positive_projected_hamiltonian": bool(eigenvalues[0] >= -1.0e-9),
    }


def build_report(quick: bool):
    parent = chiral_report()
    static = parent["photon"]["static_odd_coefficient"]
    residue = parent["photon"]["projected_residue"]
    cases = QUICK_CASES if quick else FULL_CASES
    rows = [spectrum(prime, size, static, residue) for prime, size in cases]

    p3_l3 = next(row for row in rows if row["prime"] == 3 and row["lattice_size"] == 3)
    p3_l5 = next(row for row in rows if row["prime"] == 3 and row["lattice_size"] == 5)
    spatial_power = math.log(
        p3_l5["first_dimensionless_gap"] / p3_l3["first_dimensionless_gap"]
    ) / math.log(
        (2.0 * np.sin(np.pi / 5.0)) / (2.0 * np.sin(np.pi / 3.0))
    )

    l3_rows = [row for row in rows if row["lattice_size"] == 3]
    small_prime_doublets = [
        row for row in l3_rows if row["prime"] in (3, 5, 7)
    ]
    high_prime = next(
        (row for row in l3_rows if row["prime"] == 11),
        None,
    )
    checks = {
        "self_duality_form_is_nondegenerate": all(
            row["constraint_rank_mod_p"] == row["constraint_dimension"]
            for row in rows
        ),
        "finite_lock_has_expected_ground_band": all(
            row["pair_lock_ground_dimension"] == row["prime"]
            and row["projected_hilbert_dimension"]
            == row["expected_projected_dimension"]
            for row in rows
        ),
        "projected_hamiltonians_are_positive": all(
            row["positive_projected_hamiltonian"] for row in rows
        ),
        "small_prime_first_gap_is_twofold": all(
            row["first_gap_multiplicity"] == 2
            for row in small_prime_doublets
        ),
        "lock_mixing_is_perturbatively_controlled": max(
            row["rigorous_lock_mixing_bound"] for row in rows
        ) < 0.01,
        "spatial_dispersion_is_not_linear": abs(spatial_power - 1.0) > 1.0,
        "local_dimension_refinement_has_no_stable_plateau": (
            high_prime is None
            or abs(high_prime["relative_frequency_error"]) > 0.5
        ),
    }
    return {
        "module": (
            "phase_junction_network/construction/coherent_junction/"
            "check_finite_weyl_duality_lock.py"
        ),
        "status": (
            "Stage 5J NO-GO/PARTIAL: the finite Weyl self-duality lock has the "
            "right physical-band dimension and positivity, but its projected "
            "magnetic dynamics does not retain z=1 spatial scaling"
        ),
        "run_mode": "quick" if quick else "full",
        "construction": {
            "fields": "two duality copies times two transverse polarizations",
            "spatial_derivative": "local periodic central difference",
            "constraints": "p-(1/2) epsilon_duality tensor curl q",
            "global_modes": "one coordinate per field/polarization fixed",
            "finite_lock": "sum_i [2-W_i-W_i^dagger]",
            "positive_dynamics": "projected sum_links [1-cos(Delta q)]",
            "new_continuous_coefficients": 0,
        },
        "ring_input": {
            "static_odd_coefficient": static,
            "projected_residue": residue,
            "duality_speed": static / residue,
        },
        "rows": rows,
        "p3_spatial_gap_power_from_L3_L5": spatial_power,
        "checks": checks,
        "execution_pass": all(checks.values()),
        "decision": (
            "Energetically penalizing the noncommuting self-duality constraints "
            "does produce a finite positive physical band, but it changes the "
            "infrared derivative order. The p=3 spatial gap scales roughly as "
            "k^3.82 between L=3 and L=5, not k. The correct first-order "
            "symplectic dynamics must therefore be implemented as a local "
            "unitary/canonical update rather than as a static constraint-lock "
            "Hamiltonian."
        ),
        "claim_boundary": {
            "established": [
                "an exact finite Weyl realization of the local self-duality constraints",
                "the expected finite physical-band dimension",
                "a positive projected Hamiltonian with a twofold first gap at small primes",
                "perturbative separation from the nonphysical lock band",
                "failure of the naive energetic lock to preserve z=1 spatial scaling",
            ],
            "not_established": [
                "a finite unitary realization of the first-order symplectic evolution",
                "a spatial continuum limit for the doubled finite model",
                "nonlinear duality-symmetric gravity",
                "matter coupling in the revised variables",
            ],
        },
        "next_gate": (
            "Implement the same local curl equations as an exact finite Clifford "
            "symplectic automaton. A leapfrog composition of two local shears is "
            "unitary on odd-prime qudits and has omega proportional to |k|; it "
            "must preserve Gauss/frame constraints and reproduce two modes."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build_report(args.quick)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(ready(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "execution_pass": report["execution_pass"],
                "spatial_power": report["p3_spatial_gap_power_from_L3_L5"],
                "rows": report["rows"],
                "failed_checks": [
                    key
                    for key, value in report["checks"].items()
                    if not value
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if report["execution_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
