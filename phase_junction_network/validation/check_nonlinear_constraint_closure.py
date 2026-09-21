#!/usr/bin/env python3
"""Stage 5C: nonlinear local constraint closure and finite-regulator obstruction.

Use the curvature coefficient derived from the coherent frame-area vertex and
the inherited frame kinetic residue. The continuum/smooth spectral Hamiltonian
is

  H[N] = int N [ a/sqrt(g) (pi^ij pi_ij - pi^2/2)
                 - b sqrt(g) R ],

with a=U_frame/2 and b=kappa/4. The code verifies the exact nonlinear D-D,
D-H and H-H brackets, evolves nontrivial constrained data without manual
projection, and counts the nonlinear principal modes.

The same calculation is repeated with a local central-difference regulator.
That regulator converges quadratically but does not close exactly at finite
spacing because the ordinary product rule is broken. This failure is retained
as a physical gate rather than hidden by projection.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

X = 0.13554178509861228
U_FRAME = 0.006789913753821769
A_KINETIC = U_FRAME / 2.0


def ring_parameters(x: float = X):
    d = 1.0 + 2.0 * x * x
    alpha = (d + np.sqrt(d * d - 4.0 * x * x)) / 2.0
    r = x / alpha
    return d, alpha, r, r**4


_, ALPHA, _, RHO = ring_parameters()
KAPPA = X * X * RHO / (2.0 * np.sqrt(ALPHA) * (1.0 - RHO) ** 2)
B_CURVATURE = KAPPA / 4.0


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


def symmetric_basis() -> np.ndarray:
    basis = []
    for index in range(3):
        matrix = np.zeros((3, 3), dtype=float)
        matrix[index, index] = 1.0
        basis.append(matrix)
    for first, second in ((0, 1), (0, 2), (1, 2)):
        matrix = np.zeros((3, 3), dtype=float)
        matrix[first, second] = matrix[second, first] = 1.0 / np.sqrt(2.0)
        basis.append(matrix)
    return np.asarray(basis)


H_BASIS = symmetric_basis()
POISSON_12 = np.block(
    [
        [np.zeros((6, 6)), np.eye(6)],
        [-np.eye(6), np.zeros((6, 6))],
    ]
)


def matrix_coordinates(matrix: np.ndarray) -> np.ndarray:
    return np.asarray([np.sum(element * matrix) for element in H_BASIS])


def derivative(field: np.ndarray, axis: int, method: str) -> np.ndarray:
    length = field.shape[axis]
    if method == "spectral":
        frequency = np.fft.fftfreq(length, 1.0 / length)
        shape = [1] * field.ndim
        shape[axis] = length
        return np.fft.ifft(
            1j * frequency.reshape(shape) * np.fft.fft(field, axis=axis),
            axis=axis,
        ).real
    if method != "central":
        raise ValueError(method)
    spacing = 2.0 * np.pi / length
    return (
        np.roll(field, -1, axis=axis) - np.roll(field, 1, axis=axis)
    ) / (2.0 * spacing)


def geometry(metric: np.ndarray, method: str):
    inverse = np.linalg.inv(metric)
    metric_derivative = np.stack(
        [derivative(metric, axis, method) for axis in range(3)], axis=-3
    )
    connection = np.zeros(metric.shape[:-2] + (3, 3, 3), dtype=float)
    for upper in range(3):
        for first in range(3):
            for second in range(3):
                connection[..., upper, first, second] = 0.5 * sum(
                    inverse[..., upper, lower]
                    * (
                        metric_derivative[..., first, second, lower]
                        + metric_derivative[..., second, first, lower]
                        - metric_derivative[..., lower, first, second]
                    )
                    for lower in range(3)
                )
    connection_derivative = np.stack(
        [derivative(connection, axis, method) for axis in range(3)], axis=-4
    )
    ricci = np.zeros_like(metric)
    for first in range(3):
        for second in range(3):
            ricci[..., first, second] = sum(
                connection_derivative[..., k, k, first, second]
                - connection_derivative[..., second, k, first, k]
                + sum(
                    connection[..., k, k, lower]
                    * connection[..., lower, first, second]
                    - connection[..., k, second, lower]
                    * connection[..., lower, first, k]
                    for lower in range(3)
                )
                for k in range(3)
            )
    scalar = np.einsum("...ij,...ij->...", inverse, ricci)
    return inverse, connection, ricci, scalar


def covariant_hessian(lapse: np.ndarray, connection: np.ndarray, method: str):
    gradient = np.stack(
        [derivative(lapse, axis, method) for axis in range(3)], axis=-1
    )
    hessian = np.stack(
        [derivative(gradient, axis, method) for axis in range(3)], axis=-2
    )
    for first in range(3):
        for second in range(3):
            hessian[..., first, second] -= sum(
                connection[..., upper, first, second] * gradient[..., upper]
                for upper in range(3)
            )
    return gradient, hessian


def scalar_gradients(
    metric: np.ndarray,
    momentum: np.ndarray,
    lapse: np.ndarray,
    method: str,
):
    inverse, connection, ricci, scalar = geometry(metric, method)
    determinant = np.linalg.det(metric)
    root = np.sqrt(determinant)
    trace = np.einsum("...ij,...ij->...", metric, momentum)
    norm = np.einsum(
        "...ij,...ik,...jl,...kl->...", momentum, metric, metric, momentum
    )
    lowered = np.einsum("...ik,...kl,...lj->...ij", metric, momentum, metric)
    d_momentum = (
        2.0
        * A_KINETIC
        * lapse[..., None, None]
        / root[..., None, None]
        * (lowered - 0.5 * trace[..., None, None] * metric)
    )

    momentum_metric_momentum = np.einsum(
        "...ik,...kl,...lj->...ij", momentum, metric, momentum
    )
    kinetic_scalar = norm - 0.5 * trace * trace
    d_metric_kinetic = (
        A_KINETIC
        * lapse[..., None, None]
        / root[..., None, None]
        * (
            2.0 * momentum_metric_momentum
            - trace[..., None, None] * momentum
            - 0.5 * kinetic_scalar[..., None, None] * inverse
        )
    )

    _, hessian = covariant_hessian(lapse, connection, method)
    laplacian = np.einsum("...ij,...ij->...", inverse, hessian)
    ricci_up = np.einsum("...ik,...jl,...kl->...ij", inverse, inverse, ricci)
    einstein_up = ricci_up - 0.5 * scalar[..., None, None] * inverse
    hessian_up = np.einsum("...ik,...jl,...kl->...ij", inverse, inverse, hessian)
    d_metric_curvature = B_CURVATURE * root[..., None, None] * (
        lapse[..., None, None] * einstein_up
        + inverse * laplacian[..., None, None]
        - hessian_up
    )
    return d_metric_kinetic + d_metric_curvature, d_momentum


def scalar_value(
    metric: np.ndarray,
    momentum: np.ndarray,
    lapse: np.ndarray,
    method: str,
) -> float:
    _, _, _, scalar = geometry(metric, method)
    root = np.sqrt(np.linalg.det(metric))
    trace = np.einsum("...ij,...ij->...", metric, momentum)
    norm = np.einsum(
        "...ij,...ik,...jl,...kl->...", momentum, metric, metric, momentum
    )
    density = (
        A_KINETIC / root * (norm - 0.5 * trace * trace)
        - B_CURVATURE * root * scalar
    )
    return float(np.mean(lapse * density))


def lie_metric(vector: np.ndarray, metric: np.ndarray, method: str):
    output = sum(
        vector[..., axis, None, None] * derivative(metric, axis, method)
        for axis in range(3)
    )
    vector_derivative = np.stack(
        [derivative(vector, axis, method) for axis in range(3)], axis=-2
    )
    output += np.einsum("...ik,...kj->...ij", vector_derivative, metric)
    output += np.einsum("...jk,...ik->...ij", vector_derivative, metric)
    return output


def lie_momentum(vector: np.ndarray, momentum: np.ndarray, method: str):
    output = sum(
        vector[..., axis, None, None] * derivative(momentum, axis, method)
        for axis in range(3)
    )
    vector_derivative = np.stack(
        [derivative(vector, axis, method) for axis in range(3)], axis=-2
    )
    output -= np.einsum("...ki,...kj->...ij", vector_derivative, momentum)
    output -= np.einsum("...kj,...ik->...ij", vector_derivative, momentum)
    output += np.trace(
        vector_derivative, axis1=-2, axis2=-1
    )[..., None, None] * momentum
    return output


def vector_gradients(
    metric: np.ndarray, momentum: np.ndarray, vector: np.ndarray, method: str
):
    return -lie_momentum(vector, momentum, method), lie_metric(vector, metric, method)


def vector_value(
    metric: np.ndarray, momentum: np.ndarray, vector: np.ndarray, method: str
) -> float:
    return float(
        np.mean(
            np.sum(momentum * lie_metric(vector, metric, method), axis=(-2, -1))
        )
    )


def poisson_bracket(first, second) -> float:
    first_metric, first_momentum = first
    second_metric, second_momentum = second
    return float(
        np.mean(
            np.sum(
                first_metric * second_momentum
                - first_momentum * second_metric,
                axis=(-2, -1),
            )
        )
    )


def vector_commutator(first: np.ndarray, second: np.ndarray, method: str):
    return sum(
        first[..., axis, None] * derivative(second, axis, method)
        - second[..., axis, None] * derivative(first, axis, method)
        for axis in range(3)
    )


def sample_fields(length: int):
    x, y, z = 2.0 * np.pi * np.indices((length, length, length)) / length
    perturbation = np.zeros((length, length, length, 3, 3), dtype=float)
    perturbation[..., 0, 0] = 0.03 * np.sin(x + y)
    perturbation[..., 1, 1] = 0.024 * np.cos(y + z)
    perturbation[..., 2, 2] = 0.018 * np.sin(z + x)
    perturbation[..., 0, 1] = perturbation[..., 1, 0] = 0.009 * np.cos(x - z)
    perturbation[..., 1, 2] = perturbation[..., 2, 1] = 0.006 * np.sin(x + y)
    eigenvalues, eigenvectors = np.linalg.eigh(perturbation)
    metric = np.einsum(
        "...ia,...a,...ja->...ij", eigenvectors, np.exp(eigenvalues), eigenvectors
    )
    momentum = np.zeros_like(metric)
    momentum[..., 0, 0] = 0.05 * np.cos(x) + 0.02 * np.sin(y)
    momentum[..., 1, 1] = -0.03 * np.sin(y + z)
    momentum[..., 2, 2] = 0.02 * np.cos(z + x)
    momentum[..., 0, 1] = momentum[..., 1, 0] = 0.015 * np.sin(x) * np.cos(y)
    momentum[..., 1, 2] = momentum[..., 2, 1] = 0.01 * np.cos(y) * np.sin(z)
    first_lapse = 1.0 + 0.2 * np.sin(x) + 0.1 * np.cos(y) + 0.05 * np.sin(z)
    second_lapse = 1.0 + 0.15 * np.cos(x) + 0.12 * np.sin(2.0 * y) + 0.04 * np.cos(z)
    first_vector = np.stack(
        (
            0.2 * np.sin(y) + 0.1 * np.cos(z),
            0.15 * np.sin(z) + 0.08 * np.cos(x),
            0.12 * np.sin(x) + 0.07 * np.cos(y),
        ),
        axis=-1,
    )
    second_vector = np.stack(
        (
            0.18 * np.cos(z),
            0.11 * np.cos(x) + 0.06 * np.sin(y),
            0.14 * np.cos(y) + 0.05 * np.sin(z),
        ),
        axis=-1,
    )
    return metric, momentum, first_lapse, second_lapse, first_vector, second_vector


def bracket_row(length: int, method: str) -> dict[str, Any]:
    metric, momentum, first_lapse, second_lapse, first_vector, second_vector = sample_fields(length)
    first_d = vector_gradients(metric, momentum, first_vector, method)
    second_d = vector_gradients(metric, momentum, second_vector, method)
    dd = poisson_bracket(first_d, second_d)
    commutator = vector_commutator(first_vector, second_vector, method)
    dd_target = vector_value(metric, momentum, commutator, method)

    first_h = scalar_gradients(metric, momentum, first_lapse, method)
    second_h = scalar_gradients(metric, momentum, second_lapse, method)
    dh = poisson_bracket(first_d, first_h)
    lie_lapse = sum(
        first_vector[..., axis] * derivative(first_lapse, axis, method)
        for axis in range(3)
    )
    dh_target = scalar_value(metric, momentum, lie_lapse, method)

    hh = poisson_bracket(first_h, second_h)
    inverse = np.linalg.inv(metric)
    first_gradient = np.stack(
        [derivative(first_lapse, axis, method) for axis in range(3)], axis=-1
    )
    second_gradient = np.stack(
        [derivative(second_lapse, axis, method) for axis in range(3)], axis=-1
    )
    shift = np.einsum(
        "...ij,...j->...i",
        inverse,
        first_lapse[..., None] * second_gradient
        - second_lapse[..., None] * first_gradient,
    )
    hh_target = (
        A_KINETIC
        * B_CURVATURE
        * vector_value(metric, momentum, shift, method)
    )

    def relative(measured: float, target: float) -> float:
        return abs(measured - target) / max(abs(target), 1.0e-30)

    return {
        "L": length,
        "derivative": method,
        "DD_measured": dd,
        "DD_target": dd_target,
        "DD_relative_error": relative(dd, dd_target),
        "DH_measured": dh,
        "DH_target": dh_target,
        "DH_relative_error": relative(dh, dh_target),
        "HH_measured": hh,
        "HH_target": hh_target,
        "HH_relative_error": relative(hh, hh_target),
    }


def cubic_kinetic_bootstrap(seed: int = 5742) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    design, target = [], []
    for _ in range(160):
        momentum, generator = rng.normal(size=(2, 3, 3))
        momentum = 0.5 * (momentum + momentum.T)
        metric_variation = generator + generator.T
        momentum_variation = (
            -generator.T @ momentum
            - momentum @ generator
            + np.trace(generator) * momentum
        )
        quadratic = np.trace(momentum @ momentum) - 0.5 * np.trace(momentum) ** 2
        quadratic_variation = np.sum(
            (2.0 * momentum - np.trace(momentum) * np.eye(3))
            * momentum_variation
        )
        design.append(
            [
                np.trace(metric_variation) * np.trace(momentum @ momentum),
                np.trace(metric_variation) * np.trace(momentum) ** 2,
                np.trace(metric_variation @ momentum @ momentum),
                np.trace(momentum) * np.trace(metric_variation @ momentum),
            ]
        )
        target.append(np.trace(generator) * quadratic - quadratic_variation)
    design = np.asarray(design)
    target = np.asarray(target)
    coefficients, _, rank, _ = np.linalg.lstsq(design[:80], target[:80], rcond=None)
    residual = float(np.max(np.abs(design[80:] @ coefficients - target[80:])))
    return {
        "basis": [
            "tr(h) tr(pi^2)",
            "tr(h) (tr pi)^2",
            "tr(h pi^2)",
            "tr(pi) tr(h pi)",
        ],
        "coefficients": coefficients.tolist(),
        "rank": int(rank),
        "heldout_covariance_residual": residual,
        "derived_T3": (
            "2 tr(h pi^2)-tr(pi)tr(h pi)-1/2 tr(h)tr(pi^2)"
            "+1/4 tr(h)(tr pi)^2"
        ),
    }


def principal_mode_count(samples: int, seed: int) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(samples):
        raw = rng.normal(size=(3, 3))
        metric = raw @ raw.T + np.eye(3)
        metric /= np.linalg.det(metric) ** (1.0 / 3.0)
        inverse = np.linalg.inv(metric)
        covector = rng.normal(size=3)
        raised = inverse @ covector
        norm_squared = covector @ raised

        constraints = np.zeros((4, 12), dtype=float)
        for spatial in range(3):
            for component, basis in enumerate(H_BASIS):
                constraints[spatial, 6 + component] = -2.0 * np.einsum(
                    "j,k,jk", metric[spatial, :], covector, basis
                )
        scalar_symbol = np.outer(raised, raised) - norm_squared * inverse
        constraints[3, :6] = matrix_coordinates(scalar_symbol)
        bracket = constraints @ POISSON_12 @ constraints.T

        gauge_fixed = np.zeros((6, 12), dtype=float)
        gauge_fixed[:4] = constraints
        gauge_fixed[4, :6] = matrix_coordinates(inverse)
        gauge_fixed[5, 6:] = matrix_coordinates(metric)
        fixed_bracket = gauge_fixed @ POISSON_12 @ gauge_fixed.T

        first_class = int(np.linalg.matrix_rank(constraints, tol=1.0e-10))
        fixed_rank = int(np.linalg.matrix_rank(gauge_fixed, tol=1.0e-10))
        fixed_second = int(np.linalg.matrix_rank(fixed_bracket, tol=1.0e-10))
        fixed_first = fixed_rank - fixed_second
        rows.append(
            {
                "metric_condition_number": float(np.linalg.cond(metric)),
                "constraint_rank": first_class,
                "first_class_bracket_residual": float(np.max(np.abs(bracket))),
                "physical_configurations": int((12 - 2 * first_class) / 2),
                "gauge_fixed_rank": fixed_rank,
                "gauge_fixed_second_class": fixed_second,
                "gauge_fixed_first_class": fixed_first,
                "gauge_fixed_physical_configurations": int(
                    (12 - 2 * fixed_first - fixed_second) / 2
                ),
            }
        )
    return {
        "samples": samples,
        "maximum_condition_number": max(row["metric_condition_number"] for row in rows),
        "maximum_first_class_bracket_residual": max(
            row["first_class_bracket_residual"] for row in rows
        ),
        "constraint_ranks": sorted(set(row["constraint_rank"] for row in rows)),
        "physical_configuration_counts": sorted(
            set(row["physical_configurations"] for row in rows)
        ),
        "local_volume_gauge_classifications": sorted(
            set(
                (
                    row["gauge_fixed_second_class"],
                    row["gauge_fixed_first_class"],
                    row["gauge_fixed_physical_configurations"],
                )
                for row in rows
            )
        ),
        "decision": (
            "The nonlinear principal symbol retains four first-class constraints "
            "and two physical configurations. Local determinant/trace conditions "
            "are gauge fixings: four second-class plus two first-class constraints, "
            "still leaving two configurations."
        ),
    }


def constraint_densities(metric: np.ndarray, momentum: np.ndarray, method: str):
    _, connection, _, scalar = geometry(metric, method)
    root = np.sqrt(np.linalg.det(metric))
    trace = np.einsum("...ij,...ij->...", metric, momentum)
    norm = np.einsum(
        "...ij,...ik,...jl,...kl->...", momentum, metric, metric, momentum
    )
    scalar_constraint = (
        A_KINETIC / root * (norm - 0.5 * trace * trace)
        - B_CURVATURE * root * scalar
    )
    divergence = np.zeros(metric.shape[:-2] + (3,), dtype=float)
    for upper in range(3):
        divergence[..., upper] = sum(
            derivative(momentum[..., upper, lower], lower, method)
            for lower in range(3)
        ) + sum(
            connection[..., upper, first, second] * momentum[..., second, first]
            for first in range(3)
            for second in range(3)
        )
    vector_constraint = -2.0 * np.einsum(
        "...ij,...j->...i", metric, divergence
    )
    return scalar_constraint, vector_constraint


def evolution_rhs(
    metric: np.ndarray, momentum: np.ndarray, lapse: np.ndarray, method: str
):
    d_metric_hamiltonian, d_momentum_hamiltonian = scalar_gradients(
        metric, momentum, lapse, method
    )
    return d_momentum_hamiltonian, -d_metric_hamiltonian


def rk4_step(
    metric: np.ndarray,
    momentum: np.ndarray,
    lapse: np.ndarray,
    step: float,
    method: str,
):
    first_metric, first_momentum = evolution_rhs(metric, momentum, lapse, method)
    second_metric, second_momentum = evolution_rhs(
        metric + 0.5 * step * first_metric,
        momentum + 0.5 * step * first_momentum,
        lapse,
        method,
    )
    third_metric, third_momentum = evolution_rhs(
        metric + 0.5 * step * second_metric,
        momentum + 0.5 * step * second_momentum,
        lapse,
        method,
    )
    fourth_metric, fourth_momentum = evolution_rhs(
        metric + step * third_metric,
        momentum + step * third_momentum,
        lapse,
        method,
    )
    new_metric = metric + step * (
        first_metric + 2.0 * second_metric + 2.0 * third_metric + fourth_metric
    ) / 6.0
    new_momentum = momentum + step * (
        first_momentum
        + 2.0 * second_momentum
        + 2.0 * third_momentum
        + fourth_momentum
    ) / 6.0
    new_metric = 0.5 * (new_metric + np.swapaxes(new_metric, -1, -2))
    new_momentum = 0.5 * (new_momentum + np.swapaxes(new_momentum, -1, -2))
    return new_metric, new_momentum


def evolution_run(
    length: int, method: str, steps: int, time_step: float
) -> dict[str, Any]:
    x, y, z = 2.0 * np.pi * np.indices((length, length, length)) / length
    metric = np.broadcast_to(
        np.eye(3), (length, length, length, 3, 3)
    ).copy()
    momentum = np.zeros_like(metric)
    lapse = (
        1.0 + 0.1 * np.cos(x) + 0.07 * np.sin(y) + 0.04 * np.cos(z)
    ) / np.sqrt(A_KINETIC * B_CURVATURE)

    maximum_scalar = 0.0
    maximum_vector = 0.0
    for index in range(steps + 1):
        scalar_constraint, vector_constraint = constraint_densities(
            metric, momentum, method
        )
        maximum_scalar = max(maximum_scalar, float(np.max(np.abs(scalar_constraint))))
        maximum_vector = max(maximum_vector, float(np.max(np.abs(vector_constraint))))
        if index < steps:
            metric, momentum = rk4_step(
                metric, momentum, lapse, time_step, method
            )
    eigenvalues = np.linalg.eigvalsh(metric)
    return {
        "L": length,
        "derivative": method,
        "steps": steps,
        "time_step": time_step,
        "maximum_scalar_constraint": maximum_scalar,
        "maximum_vector_constraint": maximum_vector,
        "maximum_metric_change": float(np.max(np.abs(metric - np.eye(3)))),
        "maximum_momentum": float(np.max(np.abs(momentum))),
        "minimum_metric_eigenvalue": float(np.min(eigenvalues)),
        "manual_projection_steps": 0,
    }


def build_report(quick: bool) -> dict[str, Any]:
    central_sizes = (8, 12, 16) if quick else (8, 12, 16, 24, 32, 48)
    central_rows = [bracket_row(length, "central") for length in central_sizes]
    spectral = bracket_row(16, "spectral")
    powers = {
        name: float(
            np.polyfit(
                np.log(np.asarray(central_sizes, dtype=float)),
                np.log([row[f"{name}_relative_error"] for row in central_rows]),
                1,
            )[0]
        )
        for name in ("DD", "DH", "HH")
    }
    cubic = cubic_kinetic_bootstrap()
    modes = principal_mode_count(30 if quick else 200, 1505)
    spectral_evolution = evolution_run(
        12, "spectral", 30 if quick else 60, 1.0e-3
    )
    local_evolution_sizes = (8, 12) if quick else (8, 12, 16, 24)
    local_evolution = [
        evolution_run(length, "central", 25 if quick else 50, 1.0e-3)
        for length in local_evolution_sizes
    ]

    checks = {
        "exact_DD_closure": spectral["DD_relative_error"] < 1.0e-12,
        "exact_DH_closure": spectral["DH_relative_error"] < 1.0e-12,
        "exact_HH_closure": spectral["HH_relative_error"] < 1.0e-12,
        "local_regulator_errors_fall_quadratically": all(
            -2.2 < powers[name] < -1.8 for name in powers
        ),
        "strict_finite_local_closure_rejected": max(
            central_rows[-1][f"{name}_relative_error"] for name in ("DD", "DH", "HH")
        ) > 1.0e-8,
        "unprojected_spectral_evolution_preserves_constraints": (
            spectral_evolution["maximum_scalar_constraint"] < 1.0e-12
            and spectral_evolution["maximum_vector_constraint"] < 1.0e-12
            and spectral_evolution["maximum_metric_change"] > 1.0e-5
        ),
        "local_evolution_drift_decreases": all(
            right["maximum_vector_constraint"] < left["maximum_vector_constraint"]
            for left, right in zip(local_evolution, local_evolution[1:])
        ),
        "cubic_kinetic_coefficients_unique": (
            cubic["rank"] == 4
            and cubic["heldout_covariance_residual"] < 1.0e-10
            and np.max(
                np.abs(
                    np.asarray(cubic["coefficients"])
                    - np.asarray((-0.5, 0.25, 2.0, -1.0))
                )
            ) < 1.0e-10
        ),
        "nonlinear_principal_symbol_has_two_modes": (
            modes["constraint_ranks"] == [4]
            and modes["physical_configuration_counts"] == [2]
            and modes["maximum_first_class_bracket_residual"] < 1.0e-10
        ),
        "local_volume_pair_is_gauge_fixing_not_extra_mode_removal": (
            modes["local_volume_gauge_classifications"] == [(4, 2, 2)]
        ),
    }

    report = {
        "module": (
            "phase_junction_network/validation/"
            "check_nonlinear_constraint_closure.py"
        ),
        "status": (
            "Stage 5C PASS at continuum nonlinear-closure and finite-regulator "
            "obstruction scope"
            if all(checks.values())
            else "Stage 5C FAIL"
        ),
        "run_mode": "quick" if quick else "full",
        "derived_coefficients": {
            "frame_kinetic_U": U_FRAME,
            "a_equals_U_over_two": A_KINETIC,
            "coherent_ring_kappa": KAPPA,
            "b_equals_kappa_over_four": B_CURVATURE,
            "scalar_density": (
                "H=a/sqrt(g)(pi^ij pi_ij-pi^2/2)-b sqrt(g) R"
            ),
            "HH_structure_function": (
                "{H[N],H[M]}=a b D[g^-1(N dM-M dN)]"
            ),
            "connection_energy_accounting": (
                "the torsion-free connection is eliminated once into sqrt(g)R; "
                "no second connection-energy term is added"
            ),
        },
        "smooth_spectral_constraint_algebra": spectral,
        "local_central_difference_constraint_algebra": central_rows,
        "local_error_powers": powers,
        "cubic_kinetic_bootstrap": cubic,
        "nonlinear_principal_mode_count": modes,
        "unprojected_spectral_evolution": spectral_evolution,
        "unprojected_local_evolution": local_evolution,
        "checks": checks,
        "stage_pass": all(checks.values()),
        "full_finite_local_nonlinear_closure": False,
        "continuum_nonlinear_constraint_closure": True,
        "decision": (
            "The coherent area vertex plus the uniquely covariant DeWitt kinetic "
            "completion closes the nonlinear continuum constraint algebra and "
            "preserves two modes without projection. A fixed central-difference "
            "lattice does not close exactly; it approaches the algebra with O(a^2) "
            "errors. Do not relabel continuum closure as an exact finite-lattice result."
        ),
        "claim_boundary": {
            "established": [
                "exact nonlinear continuum D-D, D-H and H-H brackets with derived coefficients",
                "frame kinetic and torsion-free connection/curvature self-energy in one scalar constraint",
                "unprojected nontrivial constrained evolution in the smooth spectral regulator",
                "two nonlinear principal configuration modes and no regenerated scalar",
                "the local determinant/trace pair acts as gauge fixing and preserves the two-mode count",
                "quadratic convergence, but not exact closure, for the local central-difference regulator",
            ],
            "not_established": [
                "an exact finite local nonlinear diffeomorphism algebra",
                "a unique microscopic derivation of the temporal kinetic residue",
                "nonperturbative strong-field global existence or black-hole solutions",
                "local vacuum-energy sequestering",
                "a quantum nonlinear constraint algebra or anomaly freedom",
                "Newton's constant or empirical validity",
            ],
        },
        "next_gate": (
            "Either construct a finite local derivative/move algebra with an exact "
            "discrete Leibniz rule, or accept diffeomorphism symmetry as emergent and "
            "demonstrate regulator universality. In parallel, test strong-field and "
            "local vacuum sources with the closed continuum constraints."
        ),
    }
    if not report["stage_pass"]:
        raise AssertionError(json.dumps(ready(checks), indent=2, sort_keys=True))
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "phase_junction_network/validation/results/"
            "nonlinear_constraint_closure_full.json"
        ),
    )
    args = parser.parse_args()
    report = build_report(bool(args.quick))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(ready(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "stage_pass": report["stage_pass"],
                "spectral": report["smooth_spectral_constraint_algebra"],
                "modes": report["nonlinear_principal_mode_count"],
                "failed_checks": [
                    key for key, value in report["checks"].items() if not value
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
