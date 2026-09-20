#!/usr/bin/env python3
"""Finite Weyl implementation of the connection second-class constraint pair."""
from __future__ import annotations

import argparse
import json

import numpy as np


def clock_operators(prime: int):
    omega = np.exp(2j * np.pi / prime)
    coordinate = np.diag([omega**value for value in range(prime)])
    shift = np.zeros((prime, prime), dtype=complex)
    for value in range(prime):
        shift[(value + 1) % prime, value] = 1.0
    return shift, coordinate


def relative_constraint_oscillator(prime: int):
    shift, coordinate = clock_operators(prime)
    identity = np.eye(prime, dtype=complex)
    momentum_cost = 2.0 * identity - shift - shift.conj().T
    coordinate_cost = 2.0 * identity - coordinate - coordinate.conj().T
    scale = prime / (4.0 * np.pi)
    hamiltonian = scale * (momentum_cost + coordinate_cost)
    energies, states = np.linalg.eigh(hamiltonian)
    ground = states[:, 0]
    momentum_expectation = float(
        np.real(ground.conj() @ momentum_cost @ ground)
    )
    coordinate_expectation = float(
        np.real(ground.conj() @ coordinate_cost @ ground)
    )
    return {
        "prime": prime,
        "scale": scale,
        "ground_energy": float(energies[0]),
        "gap": float(energies[1] - energies[0]),
        "gap_error_from_one": float(1.0 - (energies[1] - energies[0])),
        "momentum_constraint_cost": momentum_expectation,
        "coordinate_constraint_cost": coordinate_expectation,
        "total_constraint_cost": momentum_expectation + coordinate_expectation,
        "coordinate_chord_rms": float(np.sqrt(coordinate_expectation)),
        "momentum_chord_rms": float(np.sqrt(momentum_expectation)),
    }


def embed(operator: np.ndarray, index: int, count: int, prime: int):
    result = np.asarray([[1.0]], dtype=complex)
    identity = np.eye(prime, dtype=complex)
    for location in range(count):
        result = np.kron(
            result,
            operator if location == index else identity,
        )
    return result


def coupled_lock_example(prime: int, map_matrix: np.ndarray):
    map_matrix = np.asarray(map_matrix, dtype=int) % prime
    connection_count, frame_count = map_matrix.shape
    variable_count = frame_count + connection_count
    dimension = prime**variable_count
    shift, coordinate = clock_operators(prime)
    shifts = [
        embed(shift, index, variable_count, prime)
        for index in range(variable_count)
    ]
    coordinates = [
        embed(coordinate, index, variable_count, prime)
        for index in range(variable_count)
    ]
    identity = np.eye(dimension, dtype=complex)
    scale = prime / (4.0 * np.pi)
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)

    for connection_index in range(connection_count):
        relative = coordinates[frame_count + connection_index].copy()
        for frame_index in range(frame_count):
            exponent = int(-map_matrix[connection_index, frame_index]) % prime
            relative = relative @ np.linalg.matrix_power(
                coordinates[frame_index],
                exponent,
            )
        connection_shift = shifts[frame_count + connection_index]
        hamiltonian += scale * (
            2.0 * identity
            - connection_shift
            - connection_shift.conj().T
            + 2.0 * identity
            - relative
            - relative.conj().T
        )

    dressed_shifts = []
    for frame_index in range(frame_count):
        dressed = shifts[frame_index].copy()
        for connection_index in range(connection_count):
            exponent = int(map_matrix[connection_index, frame_index]) % prime
            dressed = dressed @ np.linalg.matrix_power(
                shifts[frame_count + connection_index],
                exponent,
            )
        dressed_shifts.append(dressed)

    energies = np.linalg.eigvalsh(hamiltonian)
    ground_energy = energies[0]
    ground_degeneracy = int(
        np.count_nonzero(
            np.isclose(
                energies,
                ground_energy,
                rtol=1.0e-10,
                atol=1.0e-10,
            )
        )
    )
    gap = float(energies[ground_degeneracy] - ground_energy)

    hamiltonian_coordinate_commutators = [
        float(
            np.max(
                np.abs(
                    hamiltonian @ coordinates[index]
                    - coordinates[index] @ hamiltonian
                )
            )
        )
        for index in range(frame_count)
    ]
    hamiltonian_dressed_commutators = [
        float(
            np.max(
                np.abs(
                    hamiltonian @ dressed
                    - dressed @ hamiltonian
                )
            )
        )
        for dressed in dressed_shifts
    ]

    omega = np.exp(2j * np.pi / prime)
    weyl_residuals = []
    for coordinate_index in range(frame_count):
        for shift_index in range(frame_count):
            phase = omega if coordinate_index == shift_index else 1.0
            weyl_residuals.append(
                float(
                    np.max(
                        np.abs(
                            coordinates[coordinate_index]
                            @ dressed_shifts[shift_index]
                            - phase
                            * dressed_shifts[shift_index]
                            @ coordinates[coordinate_index]
                        )
                    )
                )
            )

    dressed_shift_commutators = []
    for left in range(frame_count):
        for right in range(left):
            dressed_shift_commutators.append(
                float(
                    np.max(
                        np.abs(
                            dressed_shifts[left] @ dressed_shifts[right]
                            - dressed_shifts[right] @ dressed_shifts[left]
                        )
                    )
                )
            )

    return {
        "prime": prime,
        "frame_qudits": frame_count,
        "connection_qudits": connection_count,
        "hilbert_dimension": dimension,
        "map_matrix_mod_p": map_matrix.tolist(),
        "ground_energy": float(ground_energy),
        "ground_band_dimension": ground_degeneracy,
        "expected_ground_band_dimension": prime**frame_count,
        "gap_above_ground_band": gap,
        "maximum_H_Zh_commutator": max(
            hamiltonian_coordinate_commutators,
            default=0.0,
        ),
        "maximum_H_dressed_Xh_commutator": max(
            hamiltonian_dressed_commutators,
            default=0.0,
        ),
        "maximum_dressed_frame_weyl_residual": max(
            weyl_residuals,
            default=0.0,
        ),
        "maximum_dressed_shift_commutator": max(
            dressed_shift_commutators,
            default=0.0,
        ),
    }


def power_fit(primes, values):
    slope, intercept = np.polyfit(
        np.log(np.asarray(primes, dtype=float)),
        np.log(np.asarray(values, dtype=float)),
        1,
    )
    return {
        "power": float(slope),
        "prefactor": float(np.exp(intercept)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--primes",
        default="3,5,7,11,17,23,31,47,67,97",
    )
    parser.add_argument("--output")
    args = parser.parse_args()

    primes = [int(value) for value in args.primes.split(",")]
    rows = [relative_constraint_oscillator(prime) for prime in primes]
    asymptotic = [row for row in rows if row["prime"] >= 11]

    constraint_fit = power_fit(
        [row["prime"] for row in asymptotic],
        [row["total_constraint_cost"] for row in asymptotic],
    )
    rms_fit = power_fit(
        [row["prime"] for row in asymptotic],
        [row["coordinate_chord_rms"] for row in asymptotic],
    )
    gap_error_fit = power_fit(
        [row["prime"] for row in asymptotic],
        [row["gap_error_from_one"] for row in asymptotic],
    )

    coupled_examples = [
        coupled_lock_example(
            3,
            np.asarray([[1, 1], [1, 2]], dtype=int),
        ),
        coupled_lock_example(
            5,
            np.asarray([[1], [2]], dtype=int),
        ),
    ]

    report = {
        "status": (
            "finite Weyl relative oscillator implements the connection "
            "second-class pair with an exact dressed frame algebra"
        ),
        "single_relative_qudit": {
            "hamiltonian": (
                "H_sc = p/(4*pi) * [(2-X-X^dagger) + "
                "(2-Z-Z^dagger)]"
            ),
            "rows": rows,
            "asymptotic_total_constraint_cost_fit": constraint_fit,
            "asymptotic_constraint_rms_fit": rms_fit,
            "asymptotic_gap_error_fit": gap_error_fit,
        },
        "coupled_frame_connection_examples": coupled_examples,
        "general_construction": {
            "relative_coordinate": (
                "R_a = Z_Ca * product_i Z_hi^(-A_ai)"
            ),
            "connection_momentum": "M_a = X_Ca",
            "physical_coordinate": "Zbar_i = Z_hi",
            "dressed_physical_shift": (
                "Xbar_i = X_hi * product_a X_Ca^(A_ai)"
            ),
            "factorization": (
                "The invertible change of variables r=C-Ah converts "
                "H_sc into independent relative-qudit oscillators tensor "
                "an untouched frame Hilbert space."
            ),
        },
        "decision": (
            "Use the gapped relative oscillator to implement the finite "
            "second-class connection pair. Build all physical frame moves "
            "from Zbar and dressed Xbar so they commute exactly with the "
            "constraint lock. The next step is expressing the finite "
            "Fierz-Pauli/frame Hamiltonian in this dressed algebra and "
            "testing its full lattice spectrum."
        ),
        "limitations": [
            (
                "The constraints are minimum-uncertainty locks at finite p, "
                "not simultaneous sharp eigenvalue equations."
            ),
            (
                "The chord-cost expectation vanishes as approximately 1/p "
                "and its RMS as approximately 1/sqrt(p)."
            ),
            (
                "The construction is linear in the map A; nonlinear frame "
                "dependence and nonlinear constraint closure remain open."
            ),
            (
                "The finite dressed frame Hamiltonian and the shared "
                "electromagnetic/gravitational stiffness calculation have "
                "not yet been completed."
            ),
        ],
    }

    assert all(row["gap"] > 0.75 for row in rows)
    assert -1.1 < constraint_fit["power"] < -0.8
    assert -0.6 < rms_fit["power"] < -0.4
    assert -1.1 < gap_error_fit["power"] < -0.8
    assert all(
        row["ground_band_dimension"]
        == row["expected_ground_band_dimension"]
        for row in coupled_examples
    )
    assert max(
        row["maximum_H_Zh_commutator"] for row in coupled_examples
    ) < 1.0e-12
    assert max(
        row["maximum_H_dressed_Xh_commutator"]
        for row in coupled_examples
    ) < 1.0e-12
    assert max(
        row["maximum_dressed_frame_weyl_residual"]
        for row in coupled_examples
    ) < 1.0e-12

    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")


if __name__ == "__main__":
    main()
