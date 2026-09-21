#!/usr/bin/env python3
"""Stage 5A: exact nonlinear homogeneous fixed-volume reduction.

For a positive-definite homogeneous three-metric g and symmetric momentum pi,
test the canonical pair

    C_V = log det(g) = 0,
    C_P = g_ij pi^ij = 0.

The pair removes the exact conformal momentum direction from the nonlinear
DeWitt kinetic form and makes a pure volume energy constant on the reduced
homogeneous shape manifold.

This is not local nonlinear constraint closure or a cosmological-constant
solution. It is an exact nonlinear extension of the previous linear global
fixed-volume control.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

DIM = 3
DEWITT_LAMBDA = 0.5


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


def symmetric_basis() -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for index in range(DIM):
        matrix = np.zeros((DIM, DIM), dtype=float)
        matrix[index, index] = 1.0
        basis.append(matrix)
    for first, second in ((0, 1), (0, 2), (1, 2)):
        matrix = np.zeros((DIM, DIM), dtype=float)
        matrix[first, second] = matrix[second, first] = 1.0 / math.sqrt(2.0)
        basis.append(matrix)
    return basis


BASIS = symmetric_basis()


def random_metric(rng: np.random.Generator, spread: float = 0.65) -> np.ndarray:
    matrix = rng.normal(size=(DIM, DIM))
    rotation, _ = np.linalg.qr(matrix)
    logarithms = rng.normal(scale=spread, size=DIM)
    logarithms -= np.mean(logarithms)
    metric = rotation @ np.diag(np.exp(logarithms)) @ rotation.T
    metric /= np.linalg.det(metric) ** (1.0 / DIM)
    return 0.5 * (metric + metric.T)


def random_symmetric(rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(DIM, DIM))
    return 0.5 * (matrix + matrix.T)


def trace_momentum(metric: np.ndarray, momentum: np.ndarray) -> float:
    return float(np.einsum("ij,ij->", metric, momentum))


def project_tracefree(metric: np.ndarray, momentum: np.ndarray) -> np.ndarray:
    return momentum - trace_momentum(metric, momentum) / DIM * np.linalg.inv(metric)


def kinetic(metric: np.ndarray, momentum: np.ndarray, lam: float = DEWITT_LAMBDA) -> float:
    determinant = float(np.linalg.det(metric))
    norm = float(np.einsum("ij,ik,jl,kl->", momentum, metric, metric, momentum))
    trace = trace_momentum(metric, momentum)
    return (norm - lam * trace * trace) / math.sqrt(determinant)


def kinetic_matrix(metric: np.ndarray, lam: float = DEWITT_LAMBDA) -> np.ndarray:
    determinant = float(np.linalg.det(metric))
    traces = np.asarray([np.einsum("ij,ij->", metric, item) for item in BASIS])
    output = np.zeros((len(BASIS), len(BASIS)), dtype=float)
    for left, first in enumerate(BASIS):
        for right, second in enumerate(BASIS):
            norm = np.einsum("ij,ik,jl,kl->", first, metric, metric, second)
            output[left, right] = (
                norm - lam * traces[left] * traces[right]
            ) / math.sqrt(determinant)
    return 0.5 * (output + output.T)


def nullspace(vector: np.ndarray) -> np.ndarray:
    _, _, right = np.linalg.svd(vector.reshape(1, -1), full_matrices=True)
    return right[1:].T


def constraint_bracket(metric: np.ndarray) -> float:
    # d log(det g)/dg_ij=(g^-1)_ij and d(g:pi)/dpi^ij=g_ij.
    return float(np.einsum("ij,ij->", np.linalg.inv(metric), metric))


def volume_preservation(metric: np.ndarray, momentum: np.ndarray) -> tuple[float, float]:
    determinant = float(np.linalg.det(metric))
    trace = trace_momentum(metric, momentum)
    expected = 2.0 * (1.0 - DIM * DEWITT_LAMBDA) * trace / math.sqrt(determinant)
    derivative = 2.0 / math.sqrt(determinant) * (
        metric @ momentum @ metric - DEWITT_LAMBDA * trace * metric
    )
    direct = float(np.einsum("ij,ij->", np.linalg.inv(metric), derivative))
    return expected, direct


def tracefree_path_data(
    metric: np.ndarray,
    rng: np.random.Generator,
    epsilon: float = 2.0e-4,
) -> dict[str, float]:
    values, vectors = np.linalg.eigh(metric)
    square_root = vectors @ np.diag(np.sqrt(values)) @ vectors.T
    inverse_square_root = vectors @ np.diag(1.0 / np.sqrt(values)) @ vectors.T
    raw = random_symmetric(rng)
    generator = inverse_square_root @ raw @ inverse_square_root
    generator -= np.trace(generator) / DIM * np.eye(DIM)
    generator = 0.5 * (generator + generator.T)

    def exponential(scale: float) -> np.ndarray:
        eigenvalues, eigenvectors = np.linalg.eigh(generator)
        return eigenvectors @ np.diag(np.exp(scale * eigenvalues)) @ eigenvectors.T

    path = [
        square_root @ exponential(scale) @ square_root
        for scale in (-epsilon, 0.0, epsilon)
    ]
    determinants = [float(np.linalg.det(item)) for item in path]
    volumes = [math.sqrt(item) for item in determinants]
    first = (volumes[2] - volumes[0]) / (2.0 * epsilon)
    second = (volumes[2] - 2.0 * volumes[1] + volumes[0]) / epsilon**2

    conformal = [math.exp(2.0 * scale) * metric for scale in (-epsilon, 0.0, epsilon)]
    conformal_volumes = [math.sqrt(float(np.linalg.det(item))) for item in conformal]
    conformal_first = (conformal_volumes[2] - conformal_volumes[0]) / (2.0 * epsilon)
    return {
        "maximum_determinant_error": max(abs(item - 1.0) for item in determinants),
        "shape_first_derivative": first,
        "shape_second_derivative": second,
        "conformal_first_derivative": conformal_first,
    }


def one_sample(
    metric: np.ndarray,
    momentum: np.ndarray,
    rng: np.random.Generator,
) -> dict[str, Any]:
    full_matrix = kinetic_matrix(metric)
    full_eigenvalues = np.linalg.eigvalsh(full_matrix)
    trace_vector = np.asarray([np.einsum("ij,ij->", metric, item) for item in BASIS])
    reduction = nullspace(trace_vector)
    reduced_matrix = reduction.T @ full_matrix @ reduction
    reduced_eigenvalues = np.linalg.eigvalsh(0.5 * (reduced_matrix + reduced_matrix.T))

    projected = project_tracefree(metric, momentum)
    expected, direct = volume_preservation(metric, momentum)
    path = tracefree_path_data(metric, rng)
    return {
        "determinant": float(np.linalg.det(metric)),
        "condition_number": float(np.linalg.cond(metric)),
        "constraint_bracket": constraint_bracket(metric),
        "full_negative_count": int(np.count_nonzero(full_eigenvalues < -1.0e-10)),
        "full_positive_count": int(np.count_nonzero(full_eigenvalues > 1.0e-10)),
        "reduced_negative_count": int(np.count_nonzero(reduced_eigenvalues < -1.0e-10)),
        "reduced_positive_count": int(np.count_nonzero(reduced_eigenvalues > 1.0e-10)),
        "minimum_reduced_eigenvalue": float(np.min(reduced_eigenvalues)),
        "projected_trace_residual": trace_momentum(metric, projected),
        "projected_kinetic_energy": kinetic(metric, projected),
        "conformal_kinetic_energy": kinetic(metric, np.linalg.inv(metric)),
        "volume_preservation_residual": direct - expected,
        **path,
    }


def build_report(samples: int, seed: int) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    rows = [
        one_sample(random_metric(rng), random_symmetric(rng), rng)
        for _ in range(samples)
    ]
    full_inertia = sorted(
        set((item["full_negative_count"], item["full_positive_count"]) for item in rows)
    )
    reduced_inertia = sorted(
        set((item["reduced_negative_count"], item["reduced_positive_count"]) for item in rows)
    )
    summary = {
        "sample_count": samples,
        "maximum_metric_condition_number": max(item["condition_number"] for item in rows),
        "maximum_determinant_error": max(abs(item["determinant"] - 1.0) for item in rows),
        "maximum_constraint_bracket_error_from_three": max(abs(item["constraint_bracket"] - 3.0) for item in rows),
        "full_kinetic_inertia_counts": full_inertia,
        "reduced_kinetic_inertia_counts": reduced_inertia,
        "minimum_reduced_kinetic_eigenvalue": min(item["minimum_reduced_eigenvalue"] for item in rows),
        "maximum_projected_trace_residual": max(abs(item["projected_trace_residual"]) for item in rows),
        "minimum_projected_kinetic_energy": min(item["projected_kinetic_energy"] for item in rows),
        "maximum_conformal_kinetic_energy": max(item["conformal_kinetic_energy"] for item in rows),
        "maximum_volume_preservation_formula_residual": max(abs(item["volume_preservation_residual"]) for item in rows),
        "maximum_shape_path_determinant_error": max(item["maximum_determinant_error"] for item in rows),
        "maximum_volume_potential_shape_first_derivative": max(abs(item["shape_first_derivative"]) for item in rows),
        "maximum_volume_potential_shape_second_derivative": max(abs(item["shape_second_derivative"]) for item in rows),
        "minimum_volume_potential_conformal_first_derivative": min(abs(item["conformal_first_derivative"]) for item in rows),
    }
    checks = {
        "constraint_pair_has_exact_rank_two": summary["maximum_constraint_bracket_error_from_three"] < 1.0e-12,
        "unreduced_dewitt_metric_has_one_negative_direction": full_inertia == [(1, 5)],
        "reduced_shape_metric_has_five_positive_directions": reduced_inertia == [(0, 5)],
        "reduced_kinetic_energy_positive": summary["minimum_reduced_kinetic_eigenvalue"] > 0.0 and summary["minimum_projected_kinetic_energy"] >= -1.0e-12,
        "conformal_direction_negative": summary["maximum_conformal_kinetic_energy"] < 0.0,
        "trace_projection_exact": summary["maximum_projected_trace_residual"] < 1.0e-12,
        "hamiltonian_preservation_generates_trace_constraint": summary["maximum_volume_preservation_formula_residual"] < 1.0e-12,
        "tracefree_exponential_path_preserves_volume": summary["maximum_shape_path_determinant_error"] < 1.0e-11,
        "pure_volume_term_has_no_shape_force": summary["maximum_volume_potential_shape_first_derivative"] < 1.0e-8 and summary["maximum_volume_potential_shape_second_derivative"] < 1.0e-6,
        "pure_volume_term_has_conformal_force_before_reduction": summary["minimum_volume_potential_conformal_first_derivative"] > 2.9,
    }
    return {
        "module": "phase_junction_network/microscopic/check_nonlinear_volume_reduction.py",
        "status": "Stage 5A PASS: exact nonlinear homogeneous fixed-volume reduction removes the conformal canonical pair" if all(checks.values()) else "Stage 5A FAIL",
        "canonical_constraints": {
            "volume": "C_V=log det(g)=0",
            "trace_momentum": "C_P=g_ij pi^ij=0",
            "poisson_bracket": "{C_V,C_P}=3",
            "classification": "second class",
        },
        "kinetic_form": {
            "formula": "T=(pi^ij g_ik g_jl pi^kl - 1/2 (g_ij pi^ij)^2)/sqrt(det g)",
            "unreduced_inertia": "one negative conformal plus five positive shape directions",
            "reduced_inertia": "five positive shape directions",
        },
        "summary": summary,
        "representative_rows": rows[:3],
        "checks": checks,
        "stage_pass": all(checks.values()),
        "claim_boundary": {
            "established": [
                "the fixed-volume and trace-momentum pair is exactly second class for finite nonlinear homogeneous metrics",
                "the pair removes the one negative DeWitt conformal momentum direction",
                "the five homogeneous tracefree shape directions have positive kinetic energy",
                "Hamiltonian preservation of fixed volume generates the trace-momentum constraint",
                "a pure volume energy has no force along the reduced homogeneous shape manifold",
            ],
            "not_established": [
                "local inhomogeneous nonlinear constraint closure",
                "suppression of local vacuum-energy curvature",
                "a microscopic reason that selects fixed volume",
                "gravitational self-energy sourcing and strong-field solutions",
                "absence of nonlinear scalar modes away from the homogeneous sector",
            ],
        },
        "next_gate": (
            "Promote the exact determinant/trace pair to the local lattice frame and test the scalar/vector constraint algebra through quadratic order. "
            "The local calculation must include frame self-energy and must not impose an Einstein-Hilbert vertex by hand."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=500)
    parser.add_argument("--seed", type=int, default=57)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "phase_junction_network/microscopic/nonlinear_volume_reduction_results.json"
        ),
    )
    args = parser.parse_args()
    report = build_report(args.samples, args.seed)
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
                "summary": report["summary"],
                "failed_checks": [
                    key for key, value in report["checks"].items() if not value
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if report["stage_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
