#!/usr/bin/env python3
"""Stage 5B: derive a coherent frame-area curvature vertex.

The current coherent four-stage junction fixes the full oriented return
resolvent and its absolute gain. This check couples the low entrance port to a
Hermitian frame-area matrix B, keeps the physical and reference rings, and
eliminates the heavy states exactly. The linear-curvature vertex is therefore
obtained from the declared finite parent rather than inserted as an
Einstein-Hilbert interaction.

A separate coframe calculation verifies that the finite Clifford area trace is
the three-dimensional Palatini density sqrt(g) R. The check distinguishes the
exact smooth/spectral identity from the O(a^2) error of a central-difference
regulator.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from solve_coherent_junction import X, adj, return_block, ring_dirac, ring_parameters

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "microscopic"))
import chiral_matter_model as matter


I2 = np.eye(2, dtype=complex)
SPIN = np.asarray([np.kron(sigma, I2) / 2 for sigma in (matter.SX, matter.SY, matter.SZ)])
GAMMA = np.asarray(matter.SPATIAL_GAMMAS)
EPS = np.zeros((3, 3, 3), dtype=float)
EPS[0, 1, 2] = EPS[1, 2, 0] = EPS[2, 0, 1] = 1.0
EPS[0, 2, 1] = EPS[2, 1, 0] = EPS[1, 0, 2] = -1.0


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


def norm(matrix: np.ndarray) -> float:
    return float(np.linalg.norm(matrix, 2))


def bridge_with_area(holonomy: np.ndarray, area: np.ndarray):
    """Physical/reference coherent bridge with a frame-area entrance vertex."""
    n = len(holonomy)
    if holonomy.shape != (n, n) or area.shape != (n, n):
        raise ValueError("matching square holonomy and area matrices required")
    if norm(adj(holonomy) @ holonomy - np.eye(n)) > 1.0e-10:
        raise ValueError("unitary holonomy required")
    if norm(area - adj(area)) > 1.0e-10:
        raise ValueError("Hermitian area required")

    physical = ring_dirac(holonomy)[0]
    reference = ring_dirac(np.eye(n))[0]
    heavy_dimension = len(physical)
    heavy = np.block(
        [
            [physical, np.zeros_like(physical)],
            [np.zeros_like(reference), reference],
        ]
    )
    coupling = np.zeros((2 * n, 2 * heavy_dimension), dtype=complex)
    coupling[:n, :n] = X / np.sqrt(2.0) * area
    coupling[:n, heavy_dimension : heavy_dimension + n] = X / np.sqrt(2.0) * area
    coupling[n:, 4 * n : 5 * n] = X / np.sqrt(2.0) * np.eye(n)
    coupling[n:, heavy_dimension + 4 * n : heavy_dimension + 5 * n] = -X / np.sqrt(2.0) * np.eye(n)
    full = np.block(
        [
            [np.zeros((2 * n, 2 * n), dtype=complex), coupling],
            [adj(coupling), heavy],
        ]
    )
    effective = -coupling @ np.linalg.solve(heavy, adj(coupling))
    return full, effective, coupling, heavy


def projected_area_operator(holonomy: np.ndarray, area: np.ndarray) -> np.ndarray:
    """Project the two low ports onto the +sigma_y transport quadrature."""
    _, effective, _, _ = bridge_with_area(holonomy, area)
    n = len(holonomy)
    cross = effective[:n, n:]
    return 0.5j * (cross - adj(cross))


def coherent_vertex_audit(samples: int, seed: int) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    _, alpha, _, rho = ring_parameters()
    kappa = X * X * rho / (2.0 * np.sqrt(alpha) * (1.0 - rho) ** 2)
    errors = {
        "schur_block": 0.0,
        "diagonal_low_blocks": 0.0,
        "flat_holonomy": 0.0,
        "unitary_covariance": 0.0,
    }
    for _ in range(samples):
        raw = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        generator = 0.5 * (raw + adj(raw))
        generator /= max(norm(generator), 1.0e-30)
        holonomy = expm(0.37j * generator)
        raw_area = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        area = 0.5 * (raw_area + adj(raw_area))
        _, effective, _, _ = bridge_with_area(holonomy, area)
        difference = return_block(holonomy) - return_block(np.eye(4))
        predicted = -0.5 * X * X * area @ adj(difference)
        errors["schur_block"] = max(
            errors["schur_block"], norm(effective[:4, 4:] - predicted)
        )
        errors["diagonal_low_blocks"] = max(
            errors["diagonal_low_blocks"],
            norm(effective[:4, :4]),
            norm(effective[4:, 4:]),
        )
        errors["flat_holonomy"] = max(
            errors["flat_holonomy"], norm(projected_area_operator(np.eye(4), area))
        )

        raw_unitary = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        unitary, _ = np.linalg.qr(raw_unitary)
        transformed = projected_area_operator(
            unitary @ holonomy @ adj(unitary), unitary @ area @ adj(unitary)
        )
        errors["unitary_covariance"] = max(
            errors["unitary_covariance"],
            norm(transformed - unitary @ projected_area_operator(holonomy, area) @ adj(unitary)),
        )

    raw = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    curvature = 0.5 * (raw + adj(raw))
    curvature /= max(norm(curvature), 1.0e-30)
    raw_area = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    area = 0.5 * (raw_area + adj(raw_area))
    target = -0.5 * kappa * (area @ curvature + curvature @ area)
    epsilons = np.asarray((0.002, 0.004, 0.008, 0.016), dtype=float)
    derivative_errors = []
    for epsilon in epsilons:
        plus = projected_area_operator(expm(1j * epsilon * curvature), area)
        minus = projected_area_operator(expm(-1j * epsilon * curvature), area)
        derivative = (plus - minus) / (2.0 * epsilon)
        derivative_errors.append(norm(derivative - target))
    taylor_power = float(
        np.polyfit(np.log(epsilons), np.log(derivative_errors), 1)[0]
    )

    return {
        "x": X,
        "alpha": float(alpha),
        "rho": float(rho),
        "absolute_curvature_gain_kappa": float(kappa),
        "exact_zero_energy_schur_rule": (
            "H_LR=-(x^2/2) B [G_00(W)-G_00(I)]^dagger"
        ),
        "transport_quadrature": "H_area=(i/2)(H_LR-H_LR^dagger)",
        "linear_vertex": (
            "d H_area/d epsilon|0=-(kappa/2){B,F}, W=exp(i epsilon F)"
        ),
        "maximum_errors": errors,
        "central_derivative_epsilons": epsilons.tolist(),
        "central_derivative_errors": derivative_errors,
        "central_derivative_error_power": taylor_power,
    }


def derivative(field: np.ndarray, axis: int, spectral: bool) -> np.ndarray:
    length = field.shape[axis]
    if spectral:
        frequencies = np.fft.fftfreq(length, 1.0 / length)
        shape = [1] * field.ndim
        shape[axis] = length
        return np.fft.ifft(
            1j * frequencies.reshape(shape) * np.fft.fft(field, axis=axis),
            axis=axis,
        ).real
    spacing = 2.0 * np.pi / length
    return (
        np.roll(field, -1, axis=axis) - np.roll(field, 1, axis=axis)
    ) / (2.0 * spacing)


def metric_geometry(coframe: np.ndarray, spectral: bool):
    metric = np.einsum("...ia,...ja->...ij", coframe, coframe)
    inverse = np.linalg.inv(metric)
    metric_derivative = np.stack(
        [derivative(metric, axis, spectral) for axis in range(3)], axis=-3
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
        [derivative(connection, axis, spectral) for axis in range(3)], axis=-4
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
    return metric, inverse, connection, scalar


def spin_curvature(coframe: np.ndarray, spectral: bool):
    metric, _, christoffel, scalar = metric_geometry(coframe, spectral)
    inverse_frame = np.linalg.inv(coframe)
    frame_derivative = np.stack(
        [derivative(coframe, axis, spectral) for axis in range(3)], axis=-3
    )
    omega = np.zeros(coframe.shape[:-2] + (3, 3, 3), dtype=float)
    for direction in range(3):
        covariant = frame_derivative[..., direction, :, :].copy()
        for spatial in range(3):
            for internal in range(3):
                covariant[..., spatial, internal] -= sum(
                    christoffel[..., upper, direction, spatial]
                    * coframe[..., upper, internal]
                    for upper in range(3)
                )
        omega[..., direction, :, :] = -np.einsum(
            "...ja,...cj->...ac", covariant, inverse_frame
        )
    vector_connection = 0.5 * np.einsum("cab,...iab->...ic", EPS, omega)
    connection_derivative = np.stack(
        [derivative(vector_connection, axis, spectral) for axis in range(3)],
        axis=-3,
    )
    curvature = np.zeros(coframe.shape[:-2] + (3, 3, 3), dtype=float)
    for first in range(3):
        for second in range(3):
            curvature[..., first, second, :] = (
                connection_derivative[..., first, second, :]
                - connection_derivative[..., second, first, :]
                - np.einsum(
                    "abc,...b,...c->...a",
                    EPS,
                    vector_connection[..., first, :],
                    vector_connection[..., second, :],
                )
            )
    palatini = np.einsum("ijk,...ia,...jka->...", EPS, coframe, curvature)
    return metric, scalar, omega, curvature, palatini


def matrix_area_density(coframe: np.ndarray, curvature: np.ndarray) -> np.ndarray:
    inverse_frame = np.linalg.inv(coframe)
    volume = np.linalg.det(coframe)
    gamma = np.einsum("...ai,abc->...ibc", inverse_frame, GAMMA)
    total = np.zeros(coframe.shape[:-2], dtype=float)
    for first in range(3):
        for second in range(first + 1, 3):
            area = volume[..., None, None] * (
                gamma[..., first, :, :] @ gamma[..., second, :, :]
                - gamma[..., second, :, :] @ gamma[..., first, :, :]
            ) / (2.0j)
            curvature_matrix = np.einsum(
                "...a,aij->...ij", curvature[..., first, second, :], SPIN
            )
            total += np.trace(
                area @ curvature_matrix, axis1=-2, axis2=-1
            ).real
    return total


def sample_coframe(length: int, amplitude: float = 0.04) -> np.ndarray:
    x, y, z = 2.0 * np.pi * np.indices((length, length, length)) / length
    perturbation = np.zeros((length, length, length, 3, 3), dtype=float)
    perturbation[..., 0, 0] = amplitude * np.sin(x + y)
    perturbation[..., 1, 1] = 0.8 * amplitude * np.cos(y + z)
    perturbation[..., 2, 2] = 0.6 * amplitude * np.sin(z + x)
    perturbation[..., 0, 1] = perturbation[..., 1, 0] = (
        0.3 * amplitude * np.cos(x - z)
    )
    perturbation[..., 1, 2] = perturbation[..., 2, 1] = (
        0.2 * amplitude * np.sin(x + y)
    )
    eigenvalues, eigenvectors = np.linalg.eigh(perturbation)
    return np.einsum(
        "...ia,...a,...ja->...ij",
        eigenvectors,
        np.exp(0.5 * eigenvalues),
        eigenvectors,
    )


def palatini_audit(quick: bool) -> dict[str, Any]:
    sizes = (12, 16, 24) if quick else (12, 16, 24, 32, 48)
    rows = []
    for length in sizes:
        coframe = sample_coframe(length)
        metric, scalar, omega, curvature, palatini = spin_curvature(coframe, False)
        target = np.sqrt(np.linalg.det(metric)) * scalar
        matrix_density = matrix_area_density(coframe, curvature)
        rows.append(
            {
                "L": length,
                "relative_palatini_error": float(
                    np.linalg.norm(palatini - target) / np.linalg.norm(target)
                ),
                "relative_matrix_area_error": float(
                    np.linalg.norm(matrix_density - palatini)
                    / np.linalg.norm(palatini)
                ),
                "spin_connection_antisymmetry_error": float(
                    np.linalg.norm(omega + np.swapaxes(omega, -1, -2))
                    / np.linalg.norm(omega)
                ),
            }
        )
    local_power = float(
        np.polyfit(
            np.log(np.asarray(sizes, dtype=float)),
            np.log([row["relative_palatini_error"] for row in rows]),
            1,
        )[0]
    )

    coframe = sample_coframe(16)
    metric, scalar, omega, curvature, palatini = spin_curvature(coframe, True)
    target = np.sqrt(np.linalg.det(metric)) * scalar
    matrix_density = matrix_area_density(coframe, curvature)
    spectral = {
        "relative_palatini_error": float(
            np.linalg.norm(palatini - target) / np.linalg.norm(target)
        ),
        "relative_matrix_area_error": float(
            np.linalg.norm(matrix_density - palatini) / np.linalg.norm(palatini)
        ),
        "spin_connection_antisymmetry_error": float(
            np.linalg.norm(omega + np.swapaxes(omega, -1, -2))
            / np.linalg.norm(omega)
        ),
    }
    return {
        "identity": (
            "sum_{i<j} tr(B^{ij} F_{ij})="
            "epsilon^{ijk} e_i^a F_{jk}^a=sqrt(det g) R"
        ),
        "central_difference_rows": rows,
        "central_difference_error_power": local_power,
        "smooth_spectral_control": spectral,
    }


def build_report(quick: bool) -> dict[str, Any]:
    coherent = coherent_vertex_audit(6 if quick else 20, 5926)
    palatini = palatini_audit(quick)
    kappa = coherent["absolute_curvature_gain_kappa"]
    curvature_coefficient = kappa / 4.0
    checks = {
        "exact_heavy_elimination": coherent["maximum_errors"]["schur_block"] < 1.0e-12,
        "flat_reference_cancellation": coherent["maximum_errors"]["flat_holonomy"] < 1.0e-12,
        "local_frame_covariance": coherent["maximum_errors"]["unitary_covariance"] < 1.0e-12,
        "linear_area_curvature_vertex": 1.9 < coherent["central_derivative_error_power"] < 2.1,
        "spectral_palatini_identity": palatini["smooth_spectral_control"]["relative_palatini_error"] < 1.0e-10,
        "finite_clifford_area_identity": palatini["smooth_spectral_control"]["relative_matrix_area_error"] < 1.0e-12,
        "central_regulator_converges_quadratically": -2.2 < palatini["central_difference_error_power"] < -1.8,
    }
    report = {
        "module": (
            "phase_junction_network/construction/coherent_junction/"
            "check_frame_area_vertex.py"
        ),
        "status": (
            "Stage 5B PASS: the coherent junction derives a first-order "
            "frame-area curvature vertex with fixed absolute amplitude"
            if all(checks.values())
            else "Stage 5B FAIL"
        ),
        "run_mode": "quick" if quick else "full",
        "coherent_junction": coherent,
        "palatini_geometry": palatini,
        "derived_scalar_curvature_coefficient": {
            "convention": (
                "energy is tr(H_area)/4 in the four-component Clifford block"
            ),
            "b_equals_kappa_over_four": curvature_coefficient,
            "candidate_density": "H_curvature=-b sqrt(det g) R",
        },
        "checks": checks,
        "stage_pass": all(checks.values()),
        "claim_boundary": {
            "established": [
                "an exact finite heavy-state elimination with a frame-area entrance port",
                "a flat-reference cancellation inside the microscopic Hamiltonian",
                "a fixed anticommutator area-curvature vertex at first order",
                "the exact three-dimensional Palatini identity for the same finite Clifford area",
                "an absolute curvature coefficient kappa/4 including the coherent ring residue",
            ],
            "not_established": [
                "uniqueness of the area-port architecture",
                "a complete finite quantum nonlinear gravity Hamiltonian",
                "local finite-lattice diffeomorphism closure",
                "the temporal/kinetic term from the same coherent parent",
                "strong-field solutions, Newton's constant, or empirical validity",
            ],
        },
        "next_gate": (
            "Use the derived b=kappa/4 with the inherited kinetic residue to "
            "test the exact nonlinear scalar/vector constraint algebra, "
            "unprojected evolution, and nonlinear mode count."
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
            "phase_junction_network/construction/coherent_junction/"
            "frame_area_vertex_results.json"
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
                "kappa": report["coherent_junction"][
                    "absolute_curvature_gain_kappa"
                ],
                "b": report["derived_scalar_curvature_coefficient"][
                    "b_equals_kappa_over_four"
                ],
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
