#!/usr/bin/env python3
"""Finite transfer-matrix benchmark for the reduced two-TT-mode sector.

This verifies that finite local oscillator truncation is not itself an obstacle
once the gravitational constraints have been solved. It does not prove that
the reduced oscillators arise from a local microscopic frame/connection model.
"""
from __future__ import annotations
import argparse
import json
import numpy as np


def truncated_oscillator(dimension, omega, boundary_penalty):
    annihilation = np.zeros((dimension, dimension), dtype=complex)
    for n in range(1, dimension):
        annihilation[n - 1, n] = np.sqrt(n)
    creation = annihilation.conj().T
    q = (annihilation + creation) / np.sqrt(2.0 * omega)
    p = -1j * np.sqrt(omega / 2.0) * (annihilation - creation)
    top = np.zeros((dimension, dimension), dtype=complex)
    top[-1, -1] = 1.0
    hamiltonian = (
        0.5 * (p @ p + omega * omega * q @ q)
        + boundary_penalty * omega * top
    )
    return q, p, hamiltonian, top


def tensor_product_hamiltonian(single):
    dimension = single.shape[0]
    identity = np.eye(dimension, dtype=complex)
    return np.kron(single, identity) + np.kron(identity, single)


def boundary_projector(top):
    dimension = top.shape[0]
    identity = np.eye(dimension, dtype=complex)
    return np.kron(top, identity) + np.kron(identity, top) - np.kron(top, top)


def reduced_spectrum(dimension, omega, boundary_penalty, time_step):
    q, _, single, top = truncated_oscillator(dimension, omega, boundary_penalty)
    hamiltonian = tensor_product_hamiltonian(single)
    energies, states = np.linalg.eigh(hamiltonian)
    energies = energies.real
    gaps = energies - energies[0]
    transfer_eigenvalues = np.exp(-time_step * energies)
    transfer_eigenvalues = transfer_eigenvalues[np.argsort(transfer_eigenvalues)[::-1]]
    transfer_gap = -np.log(transfer_eigenvalues[1] / transfer_eigenvalues[0]) / time_step
    boundary = boundary_projector(top)
    boundary_occupations = [
        float(np.real(states[:, index].conj() @ boundary @ states[:, index]))
        for index in range(min(4, states.shape[1]))
    ]
    return {
        "ground_energy": float(energies[0]),
        "first_six_gaps": [float(x) for x in gaps[:6]],
        "first_gap": float(gaps[1]),
        "first_gap_degeneracy": int(
            np.count_nonzero(np.isclose(gaps, gaps[1], rtol=1e-11, atol=1e-12))
        ),
        "transfer_matrix_gap": float(transfer_gap),
        "low_state_boundary_occupations": boundary_occupations,
        "maximum_low_boundary_occupation": max(boundary_occupations),
        "relative_gap_error": float(abs(gaps[1] - omega) / omega),
        "relative_transfer_gap_error": float(abs(transfer_gap - omega) / omega),
        "single_q": q,
        "single_hamiltonian": single,
    }


def euclidean_correlator(q, hamiltonian, omega, times):
    energies, states = np.linalg.eigh(hamiltonian)
    ground = states[:, 0]
    amplitudes = states.conj().T @ q @ ground
    gaps = energies - energies[0]
    values = []
    exact = []
    for time in times:
        value = np.sum(np.abs(amplitudes) ** 2 * np.exp(-gaps * time))
        values.append(float(np.real(value)))
        exact.append(float(np.exp(-omega * time) / (2.0 * omega)))
    reflection = np.empty((len(times), len(times)), dtype=float)
    for i, left in enumerate(times):
        for j, right in enumerate(times):
            reflection[i, j] = float(
                np.real(
                    np.sum(
                        np.abs(amplitudes) ** 2
                        * np.exp(-gaps * (left + right))
                    )
                )
            )
    reflection_eigenvalues = np.linalg.eigvalsh(reflection)
    return {
        "times": times,
        "values": values,
        "exact_values": exact,
        "maximum_relative_error": float(
            max(abs(value - target) / target for value, target in zip(values, exact))
        ),
        "reflection_kernel_minimum_eigenvalue": float(np.min(reflection_eigenvalues)),
        "reflection_kernel_eigenvalues": [float(x) for x in reflection_eigenvalues],
    }


def power_fit(x, y):
    slope, intercept = np.polyfit(
        np.log(np.asarray(x, dtype=float)),
        np.log(np.asarray(y, dtype=float)),
        1,
    )
    return float(slope), float(np.exp(intercept))


def compact_spectrum_row(row):
    return {
        "L": row["L"],
        "dimension": row["dimension"],
        "omega": row["omega"],
        "first_gap": row["first_gap"],
        "transfer_matrix_gap": row["transfer_matrix_gap"],
        "first_gap_degeneracy": row["first_gap_degeneracy"],
        "relative_gap_error": row["relative_gap_error"],
        "relative_transfer_gap_error": row["relative_transfer_gap_error"],
        "maximum_low_boundary_occupation": row["maximum_low_boundary_occupation"],
    }


def compact_correlator_row(row):
    return {
        "L": row["L"],
        "dimension": row["dimension"],
        "times": row["times"],
        "values": row["values"],
        "exact_values": row["exact_values"],
        "maximum_relative_error": row["maximum_relative_error"],
        "reflection_kernel_minimum_eigenvalue": row[
            "reflection_kernel_minimum_eigenvalue"
        ],
        "reflection_kernel_eigenvalues": row["reflection_kernel_eigenvalues"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", default="12,16,24,32,48,64")
    parser.add_argument("--dimensions", default="3,4,5,6,8")
    parser.add_argument("--boundary-penalty", type=float, default=10.0)
    parser.add_argument("--time-step", type=float, default=0.25)
    parser.add_argument("--output")
    args = parser.parse_args()
    sizes = [int(x) for x in args.sizes.split(",")]
    dimensions = [int(x) for x in args.dimensions.split(",")]
    times = [0.25, 0.5, 0.75, 1.0]
    rows = []
    correlator_rows = []
    maximum_gap_error = 0.0
    maximum_transfer_error = 0.0
    maximum_boundary = 0.0
    minimum_reflection = float("inf")
    for lattice_size in sizes:
        omega = float(2.0 * np.sin(np.pi / lattice_size))
        for dimension in dimensions:
            result = reduced_spectrum(
                dimension, omega, args.boundary_penalty, args.time_step
            )
            q = result.pop("single_q")
            single = result.pop("single_hamiltonian")
            correlator = euclidean_correlator(q, single, omega, times)
            rows.append({"L": lattice_size, "dimension": dimension, "omega": omega, **result})
            correlator_rows.append({"L": lattice_size, "dimension": dimension, **correlator})
            maximum_gap_error = max(maximum_gap_error, result["relative_gap_error"])
            maximum_transfer_error = max(
                maximum_transfer_error, result["relative_transfer_gap_error"]
            )
            maximum_boundary = max(
                maximum_boundary, result["maximum_low_boundary_occupation"]
            )
            minimum_reflection = min(
                minimum_reflection, correlator["reflection_kernel_minimum_eigenvalue"]
            )
    reference_dimension = max(dimensions)
    scaling_rows = [row for row in rows if row["dimension"] == reference_dimension]
    slope, prefactor = power_fit(
        [row["L"] for row in scaling_rows],
        [row["first_gap"] for row in scaling_rows],
    )
    report = {
        "status": (
            "reduced physical TT transfer matrix and Euclidean correlator verified; "
            "local microscopic embedding remains open"
        ),
        "parameters": {
            "sizes": sizes,
            "dimensions": dimensions,
            "boundary_penalty_in_units_of_omega": args.boundary_penalty,
            "time_step": args.time_step,
        },
        "spectrum_scaling_rows": [compact_spectrum_row(row) for row in scaling_rows],
        "dimension_rows_at_reference_L": [
            compact_spectrum_row(row)
            for row in rows
            if row["L"] == sizes[len(sizes) // 2]
        ],
        "correlator_examples": [
            compact_correlator_row(row)
            for row in correlator_rows
            if row["L"] == sizes[len(sizes) // 2]
            and row["dimension"] in (min(dimensions), max(dimensions))
        ],
        "summary": {
            "maximum_relative_energy_gap_error": maximum_gap_error,
            "maximum_relative_transfer_gap_error": maximum_transfer_error,
            "maximum_low_state_boundary_occupation": maximum_boundary,
            "minimum_reflection_kernel_eigenvalue": minimum_reflection,
            "gap_scaling_power": slope,
            "gap_scaling_prefactor": prefactor,
            "target_gap_scaling_power": -1.0,
        },
        "interpretation": (
            "After solving the constraints, two finite truncated oscillators "
            "support a positive transfer matrix, a twofold first excitation, "
            "the expected Euclidean tensor correlator, and a 1/L gap. This is "
            "a benchmark for the physical sector, not proof that a local finite "
            "frame/connection model reduces to it."
        ),
    }
    assert maximum_gap_error < 1.0e-11
    assert maximum_transfer_error < 1.0e-11
    assert maximum_boundary < 1.0e-11
    assert minimum_reflection > -1.0e-11
    assert abs(slope + 1.0) < 0.02
    assert all(row["first_gap_degeneracy"] == 2 for row in rows)
    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")


if __name__ == "__main__":
    main()
