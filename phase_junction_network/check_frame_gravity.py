#!/usr/bin/env python3
"""Zero-data checks for the linear Phase-Junction Frame gravity prototype.

The script verifies, for a symmetric spatial junction field h_ij and its
conjugate pi^ij:
  1. Three momentum constraints plus one curvature constraint are first-class
     at nonzero momentum and leave exactly two configuration modes.
  2. The Fierz-Pauli spatial stiffness is positive and degenerate on the two
     transverse-traceless modes.
  3. Omitting the scalar curvature constraint leaves one negative-energy
     scalar mode.
  4. The static lattice Green function approaches A/r + B.

This is an internal consistency test, not an observational fit.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass

import numpy as np
from numpy.linalg import eigvalsh, matrix_rank
from scipy.linalg import null_space


@dataclass(frozen=True)
class SpectrumResult:
    samples: int
    failures: int
    max_constraint_gauge_residual: float
    max_momentum_scalar_gauge_residual: float
    max_relative_tt_eigenvalue_error: float
    representative_k_squared: float
    representative_tt_eigenvalues: list[float]
    representative_unconstrained_scalar_eigenvalues: list[float]


@dataclass(frozen=True)
class GreenFunctionResult:
    lattice_size: int
    fit_r_min: int
    fit_r_max: int
    coefficient_A: float
    offset_B: float
    normalized_rms_residual: float


def symmetric_basis() -> np.ndarray:
    """Return an orthonormal basis for real symmetric 3x3 matrices."""
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


BASIS = symmetric_basis()


def vectorize(matrix: np.ndarray) -> np.ndarray:
    return np.asarray([np.sum(element * matrix) for element in BASIS])


def matricize(vector: np.ndarray) -> np.ndarray:
    return np.tensordot(vector, BASIS, axes=1)


def gauge_matrix(k: np.ndarray) -> np.ndarray:
    """Map a spatial relabeling xi_i to delta h_ij = k_i xi_j + k_j xi_i."""
    columns: list[np.ndarray] = []
    for index in range(3):
        xi = np.zeros(3, dtype=float)
        xi[index] = 1.0
        columns.append(vectorize(np.outer(k, xi) + np.outer(xi, k)))
    return np.column_stack(columns)


def curvature_constraint(k: np.ndarray) -> np.ndarray:
    """C = (k_i k_j - k^2 delta_ij) h_ij."""
    k_squared = float(k @ k)
    operator = np.outer(k, k) - k_squared * np.eye(3)
    return vectorize(operator)[None, :]


def momentum_constraints(k: np.ndarray) -> np.ndarray:
    """G_i = k_j pi_ij."""
    result = np.zeros((3, 6), dtype=float)
    for column, element in enumerate(BASIS):
        result[:, column] = element @ k
    return result


def trace_constraint() -> np.ndarray:
    return np.asarray([[np.trace(element) for element in BASIS]])


def spatial_stiffness(k: np.ndarray) -> np.ndarray:
    """Quadratic Fierz-Pauli spatial stiffness matrix.

    V[h] = k^2 h_ij h_ij - 2(k_i h_ij)^2
           + 2(k_i h_ij) k_j h - k^2 h^2.
    """
    k_squared = float(k @ k)

    def value(vector: np.ndarray) -> float:
        h = matricize(vector)
        trace = float(np.trace(h))
        k_dot_h = k @ h
        return float(
            k_squared * np.sum(h * h)
            - 2.0 * np.dot(k_dot_h, k_dot_h)
            + 2.0 * np.dot(k_dot_h, k) * trace
            - k_squared * trace * trace
        )

    matrix = np.zeros((6, 6), dtype=float)
    for i in range(6):
        ei = np.zeros(6, dtype=float)
        ei[i] = 1.0
        matrix[i, i] = value(ei)
        for j in range(i):
            ej = np.zeros(6, dtype=float)
            ej[j] = 1.0
            cross = 0.5 * (value(ei + ej) - value(ei) - value(ej))
            matrix[i, j] = matrix[j, i] = cross
    return matrix


def run_spectrum_checks(samples: int, seed: int) -> SpectrumResult:
    rng = np.random.default_rng(seed)
    failures = 0
    max_cg = 0.0
    max_md = 0.0
    max_relative_error = 0.0
    representative: tuple[float, list[float], list[float]] | None = None

    for _ in range(samples):
        k = rng.normal(size=3)
        if np.linalg.norm(k) < 1.0e-12:
            continue

        gauge = gauge_matrix(k)
        curvature = curvature_constraint(k)
        momentum = momentum_constraints(k)
        stiffness = spatial_stiffness(k)
        k_squared = float(k @ k)

        # TT tensors are transverse and traceless.
        tt_basis = null_space(np.vstack((momentum, trace_constraint())))
        tt_eigenvalues = eigvalsh(tt_basis.T @ stiffness @ tt_basis)

        # If the scalar constraint is omitted, the transverse sector has
        # one negative scalar eigenvalue and two positive tensor eigenvalues.
        transverse_basis = null_space(momentum)
        unconstrained_eigenvalues = eigvalsh(
            transverse_basis.T @ stiffness @ transverse_basis
        )

        max_cg = max(max_cg, float(np.max(np.abs(curvature @ gauge))))
        max_md = max(
            max_md,
            float(np.max(np.abs(momentum @ curvature.ravel()))),
        )
        max_relative_error = max(
            max_relative_error,
            float(np.max(np.abs(tt_eigenvalues - k_squared)) / k_squared),
        )

        coordinate_dof = (6 - matrix_rank(curvature)) - matrix_rank(gauge)
        momentum_dof = (6 - matrix_rank(momentum)) - 1
        passed = (
            coordinate_dof == 2
            and momentum_dof == 2
            and tt_basis.shape[1] == 2
            and np.all(tt_eigenvalues > 0.0)
            and unconstrained_eigenvalues[0] < 0.0
        )
        if not passed:
            failures += 1

        if representative is None:
            representative = (
                k_squared,
                tt_eigenvalues.tolist(),
                unconstrained_eigenvalues.tolist(),
            )

    if representative is None:
        raise RuntimeError("No nonzero wavevectors were generated.")

    return SpectrumResult(
        samples=samples,
        failures=failures,
        max_constraint_gauge_residual=max_cg,
        max_momentum_scalar_gauge_residual=max_md,
        max_relative_tt_eigenvalue_error=max_relative_error,
        representative_k_squared=representative[0],
        representative_tt_eigenvalues=representative[1],
        representative_unconstrained_scalar_eigenvalues=representative[2],
    )


def run_green_function_check(
    lattice_size: int,
    fit_r_min: int,
    fit_r_max: int,
) -> GreenFunctionResult:
    if lattice_size < 16:
        raise ValueError("lattice_size must be at least 16")
    if not (1 < fit_r_min < fit_r_max < lattice_size // 2):
        raise ValueError("Fit range must lie inside the periodic half-box.")

    n = lattice_size
    center = n // 2
    source = np.zeros((n, n, n), dtype=float)
    source[center, center, center] = 1.0
    source -= source.mean()  # Remove the periodic zero mode.

    wave = 2.0 * np.pi * np.fft.fftfreq(n)
    kx, ky, kz = np.meshgrid(wave, wave, wave, indexing="ij")
    lattice_laplacian = 4.0 * (
        np.sin(kx / 2.0) ** 2
        + np.sin(ky / 2.0) ** 2
        + np.sin(kz / 2.0) ** 2
    )

    source_k = np.fft.fftn(source)
    potential_k = np.zeros_like(source_k, dtype=complex)
    nonzero = lattice_laplacian > 1.0e-14
    potential_k[nonzero] = source_k[nonzero] / lattice_laplacian[nonzero]
    potential = np.fft.ifftn(potential_k).real

    coordinate = np.arange(n) - center
    x, y, z = np.meshgrid(coordinate, coordinate, coordinate, indexing="ij")
    radius = np.sqrt(x * x + y * y + z * z)

    radii: list[int] = []
    shell_values: list[float] = []
    for integer_radius in range(2, n // 4):
        shell = (radius >= integer_radius - 0.5) & (radius < integer_radius + 0.5)
        if np.any(shell):
            radii.append(integer_radius)
            shell_values.append(float(np.mean(potential[shell])))

    radii_array = np.asarray(radii, dtype=float)
    values_array = np.asarray(shell_values, dtype=float)
    selected = (radii_array >= fit_r_min) & (radii_array <= fit_r_max)
    design = np.column_stack(
        (1.0 / radii_array[selected], np.ones(np.count_nonzero(selected)))
    )
    coefficient_a, offset_b = np.linalg.lstsq(
        design, values_array[selected], rcond=None
    )[0]
    prediction = coefficient_a / radii_array[selected] + offset_b
    value_range = float(np.ptp(values_array[selected]))
    normalized_rms = float(
        np.sqrt(np.mean((values_array[selected] - prediction) ** 2)) / value_range
    )

    return GreenFunctionResult(
        lattice_size=n,
        fit_r_min=fit_r_min,
        fit_r_max=fit_r_max,
        coefficient_A=float(coefficient_a),
        offset_B=float(offset_b),
        normalized_rms_residual=normalized_rms,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--lattice-size", type=int, default=96)
    parser.add_argument("--fit-r-min", type=int, default=4)
    parser.add_argument("--fit-r-max", type=int, default=20)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    spectrum = run_spectrum_checks(args.samples, args.seed)
    green = run_green_function_check(
        args.lattice_size,
        args.fit_r_min,
        args.fit_r_max,
    )
    report = {
        "spectrum": asdict(spectrum),
        "static_green_function": asdict(green),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if spectrum.failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
