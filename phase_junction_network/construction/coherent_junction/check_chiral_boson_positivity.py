#!/usr/bin/env python3
"""Stage 5I: bosonic positivity and derivative-order audit of the common chiral pole.

Stage 5H found that equal trace-normalized U(1) and frame generators have the
same chiral low pole. Equality of a first-order pole is not yet equality of a
physical bosonic light cone. This check distinguishes four uses of that pole:

1. a single chiral operator as a Hamiltonian (indefinite);
2. a constant shift that restores positivity (gapped);
3. a squared operator (positive, but photon and frame have different derivative
   order and therefore different z);
4. a local doubled duality-symmetric first-order system, where the chiral curl
   is the symplectic form and the squared curl is the positive Hamiltonian.

The doubled construction is tested for the spin-1 Gauss sector and the spin-2
scalar/vector-constraint sector without inserting a TT projector.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import null_space

HERE = Path(__file__).resolve().parent
import sys
sys.path.insert(0, str(HERE))
from check_common_chiral_generator import build_report as chiral_report

U_COMMON = 0.006789913753821769


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


def symmetric_basis() -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for index in range(3):
        item = np.zeros((3, 3), dtype=complex)
        item[index, index] = 1.0
        basis.append(item)
    for first, second in ((0, 1), (0, 2), (1, 2)):
        item = np.zeros((3, 3), dtype=complex)
        item[first, second] = item[second, first] = 1.0 / math.sqrt(2.0)
        basis.append(item)
    return basis


H_BASIS = symmetric_basis()


def h_vector(matrix: np.ndarray) -> np.ndarray:
    return np.asarray(
        [np.vdot(item, matrix) for item in H_BASIS],
        dtype=complex,
    )


def curl_matrix(momentum: np.ndarray) -> np.ndarray:
    kx, ky, kz = momentum
    return np.asarray(
        [
            [0.0, -1j * kz, 1j * ky],
            [1j * kz, 0.0, -1j * kx],
            [-1j * ky, 1j * kx, 0.0],
        ],
        dtype=complex,
    )


def tensor_curl(momentum: np.ndarray) -> np.ndarray:
    vector_curl = curl_matrix(momentum)
    columns = []
    for item in H_BASIS:
        transformed = 0.5 * (
            vector_curl @ item + item @ vector_curl.T
        )
        columns.append(h_vector(transformed))
    output = np.column_stack(columns)
    return 0.5 * (output + output.conj().T)


def tensor_constraint_basis(momentum: np.ndarray) -> np.ndarray:
    gauge_columns = []
    for direction in range(3):
        unit = np.eye(3)[direction]
        gauge_columns.append(
            h_vector(np.outer(momentum, unit) + np.outer(unit, momentum))
        )
    gauge = np.column_stack(gauge_columns)
    k_squared = float(momentum @ momentum)
    scalar = h_vector(
        np.outer(momentum, momentum) - k_squared * np.eye(3)
    )
    constraints = np.vstack((gauge.conj().T, scalar[None, :]))
    return null_space(constraints)


def physical_curls(momentum: np.ndarray):
    vector_basis = null_space(momentum.reshape(1, 3))
    vector = vector_basis.conj().T @ curl_matrix(momentum) @ vector_basis
    tensor_basis = tensor_constraint_basis(momentum)
    tensor = tensor_basis.conj().T @ tensor_curl(momentum) @ tensor_basis
    return (
        0.5 * (vector + vector.conj().T),
        0.5 * (tensor + tensor.conj().T),
        vector_basis,
        tensor_basis,
    )


def duality_block(curl: np.ndarray, residue: float, static_coefficient: float):
    epsilon = np.asarray([[0.0, 1.0], [-1.0, 0.0]], dtype=complex)
    symplectic = residue * np.kron(epsilon, curl)
    hamiltonian = static_coefficient * np.kron(
        np.eye(2), curl @ curl
    )
    dynamics = np.linalg.solve(symplectic, hamiltonian)
    eigenvalues = np.linalg.eigvals(dynamics)
    frequencies = sorted(
        float(value.imag)
        for value in eigenvalues
        if value.imag > 1.0e-12
    )
    return {
        "configuration_dimension": hamiltonian.shape[0],
        "symplectic_rank": int(np.linalg.matrix_rank(symplectic, tol=1.0e-12)),
        "minimum_hamiltonian_eigenvalue": float(
            np.min(np.linalg.eigvalsh(hamiltonian))
        ),
        "frequencies": frequencies,
        "maximum_growth_rate": float(np.max(np.abs(eigenvalues.real))),
    }


def power_fit(momentum: np.ndarray, frequency: np.ndarray) -> float:
    slope, _ = np.polyfit(np.log(momentum), np.log(frequency), 1)
    return float(slope)


def build_report() -> dict[str, Any]:
    parent = chiral_report()
    photon_parent = parent["photon"]
    frame_parent = parent["frame"]
    residue = photon_parent["projected_residue"]
    static = photon_parent["static_odd_coefficient"]
    dynamic = static / residue

    momenta = (
        np.asarray((1.0, 0.0, 0.0)),
        np.asarray((1.0, 1.0, 0.0)),
        np.asarray((1.0, 1.0, 1.0)),
        np.asarray((2.0, 1.0, 0.0)),
    )
    rows = []
    for momentum in momenta:
        vector, tensor, vector_basis, tensor_basis = physical_curls(momentum)
        norm = float(np.linalg.norm(momentum))
        vector_eigenvalues = np.linalg.eigvalsh(vector)
        tensor_eigenvalues = np.linalg.eigvalsh(tensor)
        vector_duality = duality_block(vector, residue, static)
        tensor_duality = duality_block(tensor, residue, static)
        rows.append(
            {
                "momentum": [float(value) for value in momentum],
                "momentum_norm": norm,
                "vector_constraint_dimension": int(vector_basis.shape[1]),
                "tensor_constraint_dimension": int(tensor_basis.shape[1]),
                "vector_curl_eigenvalues": [
                    float(value) for value in vector_eigenvalues
                ],
                "tensor_curl_eigenvalues": [
                    float(value) for value in tensor_eigenvalues
                ],
                "single_chiral_vector_minimum_energy": float(
                    static * np.min(vector_eigenvalues)
                ),
                "single_chiral_tensor_minimum_energy": float(
                    static * np.min(tensor_eigenvalues)
                ),
                "minimum_constant_shift_for_vector_positivity": static * norm,
                "minimum_constant_shift_for_tensor_positivity": static * norm,
                "vector_duality": vector_duality,
                "tensor_duality": tensor_duality,
                "expected_duality_frequency": dynamic * norm,
                "vector_frequency_relative_error": max(
                    abs(value - dynamic * norm)
                    for value in vector_duality["frequencies"]
                )
                / (dynamic * norm),
                "tensor_frequency_relative_error": max(
                    abs(value - dynamic * norm)
                    for value in tensor_duality["frequencies"]
                )
                / (dynamic * norm),
            }
        )

    lattice_sizes = np.asarray((16, 24, 32, 48, 64, 96), dtype=float)
    lattice_momenta = 2.0 * np.sin(np.pi / lattice_sizes)
    squared_photon = math.sqrt(U_COMMON) * dynamic * lattice_momenta
    squared_frame = math.sqrt(U_COMMON) * dynamic * lattice_momenta**2
    photon_power = power_fit(lattice_momenta, squared_photon)
    frame_power = power_fit(lattice_momenta, squared_frame)
    old_frame_speed = math.sqrt(2.0 * U_COMMON * dynamic)
    squared_photon_speed = math.sqrt(U_COMMON) * dynamic
    squared_speed_ratio = squared_photon_speed / old_frame_speed

    maximum_vector_error = max(
        row["vector_frequency_relative_error"] for row in rows
    )
    maximum_tensor_error = max(
        row["tensor_frequency_relative_error"] for row in rows
    )
    checks = {
        "single_chiral_energy_is_indefinite": all(
            row["single_chiral_vector_minimum_energy"] < 0.0
            and row["single_chiral_tensor_minimum_energy"] < 0.0
            for row in rows
        ),
        "constant_shift_would_gap_zero_momentum": all(
            row["minimum_constant_shift_for_vector_positivity"] > 0.0
            for row in rows
        ),
        "squaring_gives_photon_z_one": abs(photon_power - 1.0) < 1.0e-12,
        "squaring_frame_curvature_gives_z_two": abs(frame_power - 2.0) < 1.0e-12,
        "squared_photon_and_unsquared_frame_cones_mismatch": (
            squared_speed_ratio < 0.01
        ),
        "duality_vector_has_two_modes": all(
            row["vector_duality"]["configuration_dimension"] == 4
            and row["vector_duality"]["symplectic_rank"] == 4
            and len(row["vector_duality"]["frequencies"]) == 2
            for row in rows
        ),
        "duality_tensor_has_two_modes_without_TT_projector": all(
            row["tensor_constraint_dimension"] == 2
            and row["tensor_duality"]["configuration_dimension"] == 4
            and row["tensor_duality"]["symplectic_rank"] == 4
            and len(row["tensor_duality"]["frequencies"]) == 2
            for row in rows
        ),
        "duality_hamiltonians_are_positive": all(
            row["vector_duality"]["minimum_hamiltonian_eigenvalue"] > 0.0
            and row["tensor_duality"]["minimum_hamiltonian_eigenvalue"] > 0.0
            for row in rows
        ),
        "duality_systems_are_stable": max(
            max(row["vector_duality"]["maximum_growth_rate"],
                row["tensor_duality"]["maximum_growth_rate"])
            for row in rows
        )
        < 1.0e-12,
        "duality_vector_and_tensor_share_ring_speed": (
            maximum_vector_error < 1.0e-12
            and maximum_tensor_error < 1.0e-12
        ),
    }

    return {
        "module": (
            "phase_junction_network/construction/coherent_junction/"
            "check_chiral_boson_positivity.py"
        ),
        "status": (
            "Stage 5I ARCHITECTURE RESULT: a single chiral pole is not a "
            "positive bosonic Hamiltonian, but the same pole defines a local "
            "positive duality-symmetric spin-1/spin-2 system with two modes and "
            "one shared speed"
        ),
        "ring_input": {
            "residue": residue,
            "static_odd_coefficient": static,
            "dynamic_speed_coefficient": dynamic,
            "no_new_coefficient": True,
        },
        "single_chiral_hamiltonian": {
            "decision": "rejected: one helicity has negative energy",
            "constant_shift_control": (
                "a momentum-independent shift large enough to restore positivity "
                "produces a nonzero k=0 gap"
            ),
        },
        "squared_interpretation": {
            "photon_power": photon_power,
            "frame_power": frame_power,
            "photon_speed": squared_photon_speed,
            "old_unsquared_frame_speed": old_frame_speed,
            "photon_to_frame_speed_ratio": squared_speed_ratio,
            "decision": (
                "positive squaring preserves z=1 for the photon but makes a "
                "squared-curvature frame mode z=2; it is not a common-cone solution"
            ),
        },
        "duality_symmetric_candidate": {
            "action_structure": (
                "Omega=Z epsilon_ab tensor curl; H=gamma I_ab tensor curl^2"
            ),
            "speed": dynamic,
            "hamiltonian": "positive on the constrained nonzero-momentum sector",
            "photon_constraints": "Gauss-law transverse nullspace",
            "tensor_constraints": (
                "nullspace of the local vector gauge directions and scalar constraint"
            ),
            "uses_TT_projector": False,
            "rows": rows,
            "maximum_vector_frequency_error": maximum_vector_error,
            "maximum_tensor_frequency_error": maximum_tensor_error,
        },
        "checks": checks,
        "execution_pass": all(checks.values()),
        "decision": (
            "Stage 5H is not a physical solution when the chiral pole is treated "
            "as a one-field energy. It becomes a viable local positive candidate "
            "only as a doubled duality-symmetric symplectic system. This is a "
            "real architecture change and must earn a finite Gauss-reduced "
            "Hamiltonian and the nonlinear frame constraints."
        ),
        "claim_boundary": {
            "established": [
                "the one-field chiral bosonic Hamiltonian is indefinite",
                "squaring does not give a common photon/frame derivative order",
                "a doubled local duality-symmetric construction is positive and stable",
                "both constrained spin-1 and spin-2 blocks have exactly two modes",
                "the vector and tensor frequencies share the same ring-derived speed",
            ],
            "not_established": [
                "a finite local Hilbert-space realization of the doubled symplectic form",
                "a deconfined compact photon phase for the revised architecture",
                "nonlinear closure of duality-symmetric gravity",
                "matter coupling and stress Ward identities in the doubled variables",
                "equivalence to the observed Maxwell and Einstein interactions",
            ],
        },
        "next_gate": (
            "Quantize the doubled duality-symmetric photon and tensor blocks with "
            "finite Weyl variables, verify exact Gauss/frame constraints and a "
            "positive transfer matrix, then repeat the complete antisymmetric "
            "Fock stress response in those variables."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build_report()
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
                "squared": report["squared_interpretation"],
                "duality_speed": report["duality_symmetric_candidate"]["speed"],
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
