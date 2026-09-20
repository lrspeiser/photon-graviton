#!/usr/bin/env python3
"""Constrained local reduction of the linear Phase Junction frame/connection system.

The calculation starts from the local extended variables
  (h_ij, p^ij, C^a_bc, P_a^bc)
and performs:
  1. second-class elimination of the auxiliary connection,
  2. first-class symplectic reduction by the three vector constraints and
     one scalar constraint,
  3. extraction of the two physical frequencies without constructing a
     transverse-traceless projector.

All nonzero modes on the requested odd periodic lattices are checked.
"""
from __future__ import annotations

import argparse
import itertools
import json
from dataclasses import asdict, dataclass

import numpy as np
from scipy.linalg import null_space


def symmetric_basis() -> np.ndarray:
    basis: list[np.ndarray] = []
    for i in range(3):
        element = np.zeros((3, 3), dtype=float)
        element[i, i] = 1.0
        basis.append(element)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        element = np.zeros((3, 3), dtype=float)
        element[i, j] = element[j, i] = 1.0 / np.sqrt(2.0)
        basis.append(element)
    return np.asarray(basis)


def connection_basis() -> np.ndarray:
    basis: list[np.ndarray] = []
    for a in range(3):
        for b in range(3):
            for c in range(b, 3):
                element = np.zeros((3, 3, 3), dtype=float)
                element[a, b, c] = 1.0
                if b != c:
                    element[a, c, b] = 1.0
                basis.append(element)
    return np.asarray(basis)


H_BASIS = symmetric_basis()
C_BASIS = connection_basis()


def h_vector(matrix: np.ndarray) -> np.ndarray:
    return np.asarray([np.sum(element * matrix) for element in H_BASIS])


def h_matrix(vector: np.ndarray) -> np.ndarray:
    return np.tensordot(vector, H_BASIS, axes=1)


def c_vector(tensor: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            np.sum(element * tensor) / np.sum(element * element)
            for element in C_BASIS
        ]
    )


def quadratic_matrix(basis: np.ndarray, value_function) -> np.ndarray:
    size = len(basis)
    matrix = np.zeros((size, size), dtype=float)
    diagonal = [value_function(element) for element in basis]
    for i in range(size):
        matrix[i, i] = diagonal[i]
        for j in range(i):
            cross = 0.5 * (
                value_function(basis[i] + basis[j])
                - diagonal[i]
                - diagonal[j]
            )
            matrix[i, j] = matrix[j, i] = cross
    return matrix


def connection_quadratic(connection: np.ndarray) -> float:
    value = 0.0
    for i, a, b in itertools.product(range(3), repeat=3):
        value += (
            connection[a, b, i] * connection[b, a, i]
            - connection[a, i, i] * connection[b, a, b]
        )
    return float(value)


CONNECTION_METRIC = quadratic_matrix(C_BASIS, connection_quadratic)


def gamma_tensor(h: np.ndarray, k: np.ndarray) -> np.ndarray:
    connection = np.zeros((3, 3, 3), dtype=float)
    for a, b, c in itertools.product(range(3), repeat=3):
        connection[a, b, c] = 0.5 * (
            k[b] * h[a, c]
            + k[c] * h[a, b]
            - k[a] * h[b, c]
        )
    return connection


def gamma_map(k: np.ndarray) -> np.ndarray:
    return np.column_stack(
        [c_vector(gamma_tensor(element, k)) for element in H_BASIS]
    )


def gauge_map(k: np.ndarray) -> np.ndarray:
    columns: list[np.ndarray] = []
    for direction in range(3):
        unit = np.eye(3)[direction]
        columns.append(
            h_vector(np.outer(k, unit) + np.outer(unit, k))
        )
    return np.column_stack(columns)


def scalar_constraint_vector(k: np.ndarray) -> np.ndarray:
    k_squared = float(k @ k)
    return h_vector(np.outer(k, k) - k_squared * np.eye(3))


def fierz_pauli_value(h: np.ndarray, k: np.ndarray) -> float:
    k_squared = float(k @ k)
    trace = float(np.trace(h))
    k_dot_h = k @ h
    return float(
        k_squared * np.sum(h * h)
        - 2.0 * np.dot(k_dot_h, k_dot_h)
        + 2.0 * np.dot(k_dot_h, k) * trace
        - k_squared * trace * trace
    )


def fierz_pauli_matrix(k: np.ndarray) -> np.ndarray:
    return quadratic_matrix(
        H_BASIS,
        lambda h: fierz_pauli_value(h, k),
    )


TRACE_VECTOR = np.asarray([np.trace(element) for element in H_BASIS])
KINETIC_METRIC = np.eye(6) - 0.5 * np.outer(TRACE_VECTOR, TRACE_VECTOR)


def reduced_poisson() -> np.ndarray:
    return np.block(
        [
            [np.zeros((6, 6)), np.eye(6)],
            [-np.eye(6), np.zeros((6, 6))],
        ]
    )


REDUCED_POISSON = reduced_poisson()
REDUCED_SYMPLECTIC = -REDUCED_POISSON


@dataclass(frozen=True)
class ModeResult:
    second_class_rank: int
    first_class_rank: int
    constraint_surface_dimension: int
    physical_phase_dimension: int
    minimum_physical_hessian_eigenvalue: float
    frequencies: list[float]
    relative_frequency_error: float
    maximum_growth_rate: float
    first_second_bracket_residual: float
    first_first_bracket_residual: float
    connection_pullback_residual: float
    gauge_constraint_residual: float
    symplectic_gauge_residual: float
    hamiltonian_descent_residual: float


def extended_system(k: np.ndarray):
    connection_map = gamma_map(k)
    gauge = gauge_map(k)
    scalar = scalar_constraint_vector(k)

    poisson = np.zeros((48, 48), dtype=float)
    poisson[:6, 6:12] = np.eye(6)
    poisson[6:12, :6] = -np.eye(6)
    poisson[12:30, 30:48] = np.eye(18)
    poisson[30:48, 12:30] = -np.eye(18)

    hamiltonian = np.zeros((48, 48), dtype=float)
    hamiltonian[6:12, 6:12] = KINETIC_METRIC
    hamiltonian[12:30, 12:30] = 4.0 * CONNECTION_METRIC
    hamiltonian[:6, 12:30] = (
        -4.0 * connection_map.T @ CONNECTION_METRIC
    )
    hamiltonian[12:30, :6] = (
        -4.0 * CONNECTION_METRIC @ connection_map
    )

    # Second-class constraints:
    #   primary: P_C = 0
    #   secondary: C - Gamma[h] = 0
    second_class = np.zeros((36, 48), dtype=float)
    second_class[:18, 30:48] = np.eye(18)
    second_class[18:, 12:30] = np.eye(18)
    second_class[18:, :6] = -connection_map

    # Extended first-class vector generators transform both h and C.
    first_class = np.zeros((4, 48), dtype=float)
    first_class[:3, 6:12] = gauge.T
    first_class[:3, 30:48] = (connection_map @ gauge).T
    first_class[3, :6] = scalar

    # Embed the second-class surface into the 12-dimensional frame phase space.
    embedding = np.zeros((48, 12), dtype=float)
    embedding[:12, :12] = np.eye(12)
    embedding[12:30, :6] = connection_map

    reduced_hamiltonian = embedding.T @ hamiltonian @ embedding
    reduced_constraints = np.zeros((4, 12), dtype=float)
    reduced_constraints[:3, 6:12] = gauge.T
    reduced_constraints[3, :6] = scalar

    return (
        poisson,
        hamiltonian,
        second_class,
        first_class,
        embedding,
        reduced_hamiltonian,
        reduced_constraints,
    )


def analyze_mode(k: np.ndarray, tolerance: float = 1.0e-9) -> ModeResult:
    k = np.asarray(k, dtype=float)
    k_norm = float(np.linalg.norm(k))
    if k_norm <= tolerance:
        raise ValueError("The propagating-mode reduction excludes k=0.")

    (
        extended_poisson,
        _,
        second_class,
        first_class,
        _,
        reduced_hamiltonian,
        reduced_constraints,
    ) = extended_system(k)

    second_bracket = (
        second_class @ extended_poisson @ second_class.T
    )
    first_second = (
        first_class @ extended_poisson @ second_class.T
    )
    first_first = first_class @ extended_poisson @ first_class.T

    expected_reduced = np.block(
        [
            [fierz_pauli_matrix(k), np.zeros((6, 6))],
            [np.zeros((6, 6)), KINETIC_METRIC],
        ]
    )

    constraint_basis = null_space(reduced_constraints)
    gauge_directions = (
        REDUCED_POISSON @ reduced_constraints.T
    )
    pulled_symplectic = (
        constraint_basis.T
        @ REDUCED_SYMPLECTIC
        @ constraint_basis
    )
    left, singular_values, _ = np.linalg.svd(pulled_symplectic)
    physical_dimension = int(np.count_nonzero(singular_values > tolerance))
    physical_basis = constraint_basis @ left[:, :physical_dimension]

    physical_symplectic = (
        physical_basis.T
        @ REDUCED_SYMPLECTIC
        @ physical_basis
    )
    physical_hamiltonian = (
        physical_basis.T
        @ reduced_hamiltonian
        @ physical_basis
    )
    dynamical_matrix = np.linalg.solve(
        physical_symplectic,
        physical_hamiltonian,
    )
    dynamical_eigenvalues = np.linalg.eigvals(dynamical_matrix)
    frequencies = sorted(
        float(value.imag)
        for value in dynamical_eigenvalues
        if value.imag > tolerance
    )

    return ModeResult(
        second_class_rank=int(
            np.linalg.matrix_rank(second_bracket, tol=tolerance)
        ),
        first_class_rank=int(
            np.linalg.matrix_rank(first_class, tol=tolerance)
        ),
        constraint_surface_dimension=int(constraint_basis.shape[1]),
        physical_phase_dimension=physical_dimension,
        minimum_physical_hessian_eigenvalue=float(
            np.min(np.linalg.eigvalsh(physical_hamiltonian))
        ),
        frequencies=frequencies,
        relative_frequency_error=float(
            max(abs(value - k_norm) for value in frequencies) / k_norm
        ),
        maximum_growth_rate=float(
            np.max(np.abs(dynamical_eigenvalues.real))
        ),
        first_second_bracket_residual=float(
            np.max(np.abs(first_second))
        ),
        first_first_bracket_residual=float(
            np.max(np.abs(first_first))
        ),
        connection_pullback_residual=float(
            np.max(np.abs(reduced_hamiltonian - expected_reduced))
        ),
        gauge_constraint_residual=float(
            np.max(
                np.abs(reduced_constraints @ gauge_directions)
            )
        ),
        symplectic_gauge_residual=float(
            np.max(
                np.abs(
                    constraint_basis.T
                    @ REDUCED_SYMPLECTIC
                    @ gauge_directions
                )
            )
        ),
        hamiltonian_descent_residual=float(
            np.max(
                np.abs(
                    constraint_basis.T
                    @ reduced_hamiltonian
                    @ gauge_directions
                )
            )
        ),
    )


def centered_modes(lattice_size: int):
    half = lattice_size // 2
    for mode in itertools.product(
        range(-half, half + 1),
        repeat=3,
    ):
        if mode != (0, 0, 0):
            yield mode


def lattice_momentum(mode, lattice_size: int) -> np.ndarray:
    return 2.0 * np.sin(
        np.pi * np.asarray(mode, dtype=float) / lattice_size
    )


def analyze_lattice(lattice_size: int):
    maxima = {
        "relative_frequency_error": 0.0,
        "maximum_growth_rate": 0.0,
        "first_second_bracket_residual": 0.0,
        "first_first_bracket_residual": 0.0,
        "connection_pullback_residual": 0.0,
        "gauge_constraint_residual": 0.0,
        "symplectic_gauge_residual": 0.0,
        "hamiltonian_descent_residual": 0.0,
    }
    minimum_hessian = float("inf")
    failures = 0
    count = 0
    representative = None

    for mode in centered_modes(lattice_size):
        result = analyze_mode(lattice_momentum(mode, lattice_size))
        count += 1
        if representative is None:
            representative = {
                "mode": list(mode),
                **asdict(result),
            }

        passed = (
            result.second_class_rank == 36
            and result.first_class_rank == 4
            and result.constraint_surface_dimension == 8
            and result.physical_phase_dimension == 4
            and len(result.frequencies) == 2
            and result.minimum_physical_hessian_eigenvalue > 0.0
        )
        if not passed:
            failures += 1

        minimum_hessian = min(
            minimum_hessian,
            result.minimum_physical_hessian_eigenvalue,
        )
        for key in maxima:
            maxima[key] = max(maxima[key], getattr(result, key))

    return {
        "lattice_size": lattice_size,
        "nonzero_modes": count,
        "failures": failures,
        "minimum_physical_hessian_eigenvalue": minimum_hessian,
        "maximum_residuals": maxima,
        "representative": representative,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sizes",
        default="3,5,7",
        help="comma-separated odd periodic lattice sizes",
    )
    parser.add_argument("--output")
    args = parser.parse_args()

    sizes = [int(value) for value in args.sizes.split(",")]
    if any(size < 3 or size % 2 == 0 for size in sizes):
        raise ValueError("All lattice sizes must be odd and at least 3.")

    lattices = [analyze_lattice(size) for size in sizes]
    report = {
        "status": (
            "linear constrained local-to-physical reduction verified "
            "without a TT projector"
        ),
        "phase_space_count": {
            "extended_phase_dimension": 48,
            "second_class_constraints": 36,
            "first_class_constraints": 4,
            "physical_phase_dimension": 48 - 36 - 2 * 4,
            "physical_configuration_modes": 2,
        },
        "method": (
            "Eliminate P_C=0 and C-Gamma[h]=0 as a full-rank second-class "
            "set; then quotient the four first-class frame constraints by "
            "the null directions of the pulled-back symplectic form."
        ),
        "uses_transverse_traceless_projector": False,
        "zero_mode_policy": (
            "k=0 is excluded from the propagating-mode test and retained "
            "as the previously identified periodic global frame sector."
        ),
        "lattices": lattices,
        "decision": (
            "The local linear frame/connection system reduces to a positive "
            "four-dimensional physical phase space with two equal frequencies "
            "|k_hat|. The next unresolved step is a finite local implementation "
            "of this reduction and its nonlinear constraint algebra."
        ),
    }

    assert all(entry["failures"] == 0 for entry in lattices)
    assert max(
        entry["maximum_residuals"]["relative_frequency_error"]
        for entry in lattices
    ) < 1.0e-12
    assert max(
        entry["maximum_residuals"]["maximum_growth_rate"]
        for entry in lattices
    ) < 1.0e-12

    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")


if __name__ == "__main__":
    main()
