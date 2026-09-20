#!/usr/bin/env python3
"""Real-space constrained reduction of the linear Phase Junction frame sector.

The reduction is performed on full periodic real-space matrices, with no
Fourier block diagonalization and no transverse-traceless projector. Fourier
symbols are used only after the reduction to compare the resulting spectrum
with the analytic lattice expectation.
"""
from __future__ import annotations

import argparse
import itertools
import json

import numpy as np
from scipy.linalg import null_space

from check_constrained_local_reduction import (
    CONNECTION_METRIC,
    H_BASIS,
    KINETIC_METRIC,
)


def periodic_central_derivative(lattice_size: int, axis: int) -> np.ndarray:
    sites = lattice_size**3
    derivative = np.zeros((sites, sites), dtype=float)

    def index(x: int, y: int, z: int) -> int:
        return (x * lattice_size + y) * lattice_size + z

    for x, y, z in itertools.product(range(lattice_size), repeat=3):
        coordinate = [x, y, z]
        plus = coordinate.copy()
        minus = coordinate.copy()
        plus[axis] = (plus[axis] + 1) % lattice_size
        minus[axis] = (minus[axis] - 1) % lattice_size
        row = index(x, y, z)
        derivative[row, index(*plus)] = 0.5
        derivative[row, index(*minus)] = -0.5
    return derivative


def real_space_maps(lattice_size: int):
    sites = lattice_size**3
    derivatives = [
        periodic_central_derivative(lattice_size, axis)
        for axis in range(3)
    ]
    root_two = np.sqrt(2.0)

    gauge = np.zeros((6 * sites, 3 * sites), dtype=float)
    gauge[0 * sites : 1 * sites, 0 * sites : 1 * sites] = (
        2.0 * derivatives[0]
    )
    gauge[1 * sites : 2 * sites, 1 * sites : 2 * sites] = (
        2.0 * derivatives[1]
    )
    gauge[2 * sites : 3 * sites, 2 * sites : 3 * sites] = (
        2.0 * derivatives[2]
    )
    gauge[3 * sites : 4 * sites, 1 * sites : 2 * sites] = (
        root_two * derivatives[0]
    )
    gauge[3 * sites : 4 * sites, 0 * sites : 1 * sites] = (
        root_two * derivatives[1]
    )
    gauge[4 * sites : 5 * sites, 2 * sites : 3 * sites] = (
        root_two * derivatives[0]
    )
    gauge[4 * sites : 5 * sites, 0 * sites : 1 * sites] = (
        root_two * derivatives[2]
    )
    gauge[5 * sites : 6 * sites, 2 * sites : 3 * sites] = (
        root_two * derivatives[1]
    )
    gauge[5 * sites : 6 * sites, 1 * sites : 2 * sites] = (
        root_two * derivatives[2]
    )

    second = [derivative @ derivative for derivative in derivatives]
    scalar = np.zeros((sites, 6 * sites), dtype=float)
    scalar[:, 0 * sites : 1 * sites] = -(second[1] + second[2])
    scalar[:, 1 * sites : 2 * sites] = -(second[0] + second[2])
    scalar[:, 2 * sites : 3 * sites] = -(second[0] + second[1])
    scalar[:, 3 * sites : 4 * sites] = (
        root_two * derivatives[0] @ derivatives[1]
    )
    scalar[:, 4 * sites : 5 * sites] = (
        root_two * derivatives[0] @ derivatives[2]
    )
    scalar[:, 5 * sites : 6 * sites] = (
        root_two * derivatives[1] @ derivatives[2]
    )

    connection = np.zeros((18 * sites, 6 * sites), dtype=float)
    component = 0
    for a in range(3):
        for b in range(3):
            for c in range(b, 3):
                for frame_component, basis in enumerate(H_BASIS):
                    block = 0.5 * (
                        basis[a, c] * derivatives[b]
                        + basis[a, b] * derivatives[c]
                        - basis[b, c] * derivatives[a]
                    )
                    connection[
                        component * sites : (component + 1) * sites,
                        frame_component * sites : (frame_component + 1) * sites,
                    ] = block
                component += 1

    return gauge, scalar, connection


def independent_row_basis(matrix: np.ndarray, tolerance: float = 1.0e-9):
    _, singular_values, right = np.linalg.svd(matrix, full_matrices=False)
    rank = int(np.count_nonzero(singular_values > tolerance))
    return right[:rank, :], rank


def expected_central_frequencies(lattice_size: int) -> np.ndarray:
    half = lattice_size // 2
    frequencies: list[float] = []
    for mode in itertools.product(range(-half, half + 1), repeat=3):
        if mode == (0, 0, 0):
            continue
        momentum = np.sin(
            2.0 * np.pi * np.asarray(mode, dtype=float) / lattice_size
        )
        frequencies.extend([float(np.linalg.norm(momentum))] * 2)
    return np.sort(np.asarray(frequencies))


def inertia(values: np.ndarray, tolerance: float = 1.0e-10):
    return {
        "negative": int(np.count_nonzero(values < -tolerance)),
        "zero": int(np.count_nonzero(np.abs(values) <= tolerance)),
        "positive": int(np.count_nonzero(values > tolerance)),
    }


def analyze_lattice(lattice_size: int, tolerance: float = 1.0e-9):
    sites = lattice_size**3
    gauge, scalar, connection = real_space_maps(lattice_size)

    connection_metric = np.kron(CONNECTION_METRIC, np.eye(sites))
    frame_potential = -4.0 * connection.T @ connection_metric @ connection
    frame_kinetic = np.kron(KINETIC_METRIC, np.eye(sites))
    zero = np.zeros_like(frame_potential)
    hamiltonian = np.block(
        [
            [frame_potential, zero],
            [zero, frame_kinetic],
        ]
    )

    phase_dimension = 12 * sites
    poisson = np.block(
        [
            [np.zeros((6 * sites, 6 * sites)), np.eye(6 * sites)],
            [-np.eye(6 * sites), np.zeros((6 * sites, 6 * sites))],
        ]
    )
    symplectic = -poisson

    constraints = np.zeros((4 * sites, phase_dimension), dtype=float)
    constraints[: 3 * sites, 6 * sites :] = gauge.T
    constraints[3 * sites :, : 6 * sites] = scalar

    # Remove the six constant frame coordinates and six constant momenta.
    # This isolates the propagating sector while keeping the reduction itself
    # in real space.
    site_mean_zero = null_space(np.ones((1, sites)))
    frame_mean_zero = np.kron(np.eye(6), site_mean_zero)
    phase_mean_zero = np.block(
        [
            [frame_mean_zero, np.zeros_like(frame_mean_zero)],
            [np.zeros_like(frame_mean_zero), frame_mean_zero],
        ]
    )

    reduced_hamiltonian = phase_mean_zero.T @ hamiltonian @ phase_mean_zero
    reduced_symplectic = phase_mean_zero.T @ symplectic @ phase_mean_zero
    reduced_constraints = constraints @ phase_mean_zero

    row_basis, constraint_rank = independent_row_basis(
        reduced_constraints,
        tolerance,
    )
    constraint_surface = null_space(row_basis)
    pulled_symplectic = (
        constraint_surface.T
        @ reduced_symplectic
        @ constraint_surface
    )
    left, singular_values, _ = np.linalg.svd(pulled_symplectic)
    physical_dimension = int(np.count_nonzero(singular_values > tolerance))
    physical_basis = constraint_surface @ left[:, :physical_dimension]

    physical_symplectic = (
        physical_basis.T
        @ reduced_symplectic
        @ physical_basis
    )
    physical_hamiltonian = (
        physical_basis.T
        @ reduced_hamiltonian
        @ physical_basis
    )
    dynamics = np.linalg.solve(physical_symplectic, physical_hamiltonian)
    dynamical_eigenvalues = np.linalg.eigvals(dynamics)
    positive_frequencies = np.sort(
        dynamical_eigenvalues.imag[
            dynamical_eigenvalues.imag > tolerance
        ]
    )
    expected = expected_central_frequencies(lattice_size)

    return {
        "lattice_size": lattice_size,
        "sites": sites,
        "extended_connection_second_class_constraints": 36 * sites,
        "mean_zero_phase_dimension": int(phase_mean_zero.shape[1]),
        "first_class_constraint_rank": constraint_rank,
        "expected_first_class_constraint_rank": 4 * (sites - 1),
        "constraint_surface_dimension": int(constraint_surface.shape[1]),
        "physical_phase_dimension": physical_dimension,
        "expected_physical_phase_dimension": 4 * (sites - 1),
        "positive_frequency_count": int(len(positive_frequencies)),
        "expected_positive_frequency_count": 2 * (sites - 1),
        "minimum_physical_hessian_eigenvalue": float(
            np.min(np.linalg.eigvalsh(physical_hamiltonian))
        ),
        "maximum_relative_frequency_error": float(
            np.max(np.abs(positive_frequencies - expected) / expected)
        ),
        "maximum_growth_rate": float(
            np.max(np.abs(dynamical_eigenvalues.real))
        ),
        "scalar_vector_closure_residual": float(
            np.max(np.abs(scalar @ gauge))
        ),
        "minimum_nonzero_symplectic_singular_value": float(
            singular_values[physical_dimension - 1]
        ),
        "maximum_gauge_null_singular_value": float(
            singular_values[physical_dimension]
        ),
        "minimum_frequency": float(positive_frequencies[0]),
        "maximum_frequency": float(positive_frequencies[-1]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sizes",
        default="3,5",
        help="comma-separated odd periodic lattice sizes",
    )
    parser.add_argument("--output")
    args = parser.parse_args()

    sizes = [int(value) for value in args.sizes.split(",")]
    if any(size < 3 or size % 2 == 0 for size in sizes):
        raise ValueError("All lattice sizes must be odd and at least 3.")

    lattices = [analyze_lattice(size) for size in sizes]
    global_kinetic_eigenvalues = np.linalg.eigvalsh(KINETIC_METRIC)
    report = {
        "status": (
            "full real-space propagating-sector reduction verified without "
            "Fourier block reduction or a TT projector"
        ),
        "uses_fourier_during_reduction": False,
        "uses_transverse_traceless_projector": False,
        "derivative_regulator": (
            "periodic central difference with symbol i*sin(2*pi*n/L)"
        ),
        "lattices": lattices,
        "periodic_global_sector": {
            "phase_dimension": 12,
            "hamiltonian_inertia": {
                "negative": 1,
                "zero": 6,
                "positive": 5,
            },
            "kinetic_eigenvalues": [
                float(value) for value in global_kinetic_eigenvalues
            ],
            "interpretation": (
                "The six constant frame coordinates have zero potential and "
                "the homogeneous trace momentum has negative kinetic sign. "
                "The propagating test therefore fixes/removes the global frame "
                "sector. A fundamental model must gauge, constrain, lift, or "
                "otherwise explain this global conformal direction."
            ),
        },
        "decision": (
            "The local real-space constraints reduce to exactly two positive "
            "propagating modes for every nonzero lattice momentum. The next "
            "gate is a finite local implementation of the second-class "
            "connection constraints and a principled treatment of the global "
            "volume/conformal mode."
        ),
    }

    assert all(
        entry["first_class_constraint_rank"]
        == entry["expected_first_class_constraint_rank"]
        for entry in lattices
    )
    assert all(
        entry["physical_phase_dimension"]
        == entry["expected_physical_phase_dimension"]
        for entry in lattices
    )
    assert all(
        entry["positive_frequency_count"]
        == entry["expected_positive_frequency_count"]
        for entry in lattices
    )
    assert max(
        entry["maximum_relative_frequency_error"] for entry in lattices
    ) < 1.0e-12
    assert max(entry["maximum_growth_rate"] for entry in lattices) < 1.0e-12
    assert all(
        entry["minimum_physical_hessian_eigenvalue"] > 0.0
        for entry in lattices
    )
    assert inertia(global_kinetic_eigenvalues) == {
        "negative": 1,
        "zero": 0,
        "positive": 5,
    }

    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")


if __name__ == "__main__":
    main()
