"""Charged endpoint and neutral companion sectors for Stage 3F."""
from __future__ import annotations
from typing import Any
import numpy as np
from all_sector_move_common import (
    X_ROOT, COEFFICIENTS, SMALL_X, canonical_effective,
    completed_hamiltonian, factorized_completion_residual,
    partial_operator, power_fit, two_state_inventory,
)


def matter_model(x: float, coefficient: float, phase: float = 0.0):
    spin = 2
    states = [
        (position, flux, ancilla)
        for position in (0, 1)
        for flux in range(-spin, spin + 1)
        for ancilla in (-1, 0, 1)
    ]
    lookup = {state: index for index, state in enumerate(states)}
    dimension = len(states)
    left_gauss = np.zeros(dimension)
    right_gauss = np.zeros(dimension)
    ancilla_energy = np.zeros(dimension)
    left_number = np.zeros(dimension)

    for index, (position, flux, ancilla) in enumerate(states):
        n_left = 1 - position
        n_right = position
        left_gauss[index] = flux + n_left + ancilla - 1
        right_gauss[index] = -flux + n_right - ancilla
        ancilla_energy[index] = ancilla * ancilla
        left_number[index] = n_left

    moves = (
        partial_operator(states, lookup, lambda state: (1, state[1], 1)
            if state[0] == 0 and state[2] == 0 else None),
        partial_operator(states, lookup, lambda state: (1, state[1], 0)
            if state[0] == 0 and state[2] == -1 else None),
        partial_operator(states, lookup, lambda state: (state[0], state[1] + 1, 0)
            if state[2] == 1 and state[1] < spin else None),
        partial_operator(states, lookup, lambda state: (state[0], state[1] + 1, -1)
            if state[2] == 0 and state[1] < spin else None),
    )
    left_generator = np.diag(left_gauss).astype(complex)
    right_generator = np.diag(right_gauss).astype(complex)
    base = 0.5 * (
        left_generator @ left_generator + right_generator @ right_generator
    ) + np.diag(ancilla_energy)
    hamiltonian = completed_hamiltonian(
        base, moves, x, coefficient, (phase, phase, 0.0, 0.0)
    )
    physical = np.flatnonzero(
        np.isclose(left_gauss, 0.0) & np.isclose(right_gauss, 0.0)
    )
    low = [index for index in physical if states[index][2] == 0]
    physical_hamiltonian = hamiltonian[np.ix_(physical, physical)]
    local_low = [list(physical).index(index) for index in low]
    low_basis = np.eye(len(physical), dtype=complex)[:, local_low]
    effective, energies, eigenvectors = canonical_effective(
        physical_hamiltonian, low_basis
    )

    left_number_operator = np.diag(left_number).astype(complex)
    matter_move = moves[0] + moves[1]
    current = 1j * x * (matter_move - matter_move.conj().T)
    continuity = float(np.max(np.abs(
        1j * (hamiltonian @ left_number_operator - left_number_operator @ hamiltonian)
        + current
    )))

    physical_number = left_number_operator[np.ix_(physical, physical)]
    physical_current = current[np.ix_(physical, physical)]
    number_eigenbasis = eigenvectors.conj().T @ physical_number @ eigenvectors
    current_eigenbasis = eigenvectors.conj().T @ physical_current @ eigenvectors
    energy_difference = energies[:, None] - energies[None, :]
    spectral_ward = float(np.max(np.abs(
        1j * energy_difference * number_eigenbasis + current_eigenbasis
    )))

    derivative = -1j * x * (
        matter_move[np.ix_(physical, physical)]
        - matter_move[np.ix_(physical, physical)].conj().T
    )
    second_derivative = x * (
        matter_move[np.ix_(physical, physical)]
        + matter_move[np.ix_(physical, physical)].conj().T
    )
    ground = eigenvectors[:, 0]
    diamagnetic = float(np.real(ground.conj() @ second_derivative @ ground))
    paramagnetic = 0.0
    for index in range(1, len(energies)):
        amplitude = eigenvectors[:, index].conj() @ derivative @ ground
        paramagnetic += 2.0 * abs(amplitude) ** 2 / (energies[0] - energies[index])

    return {
        "states": states,
        "physical_dimension": len(physical),
        "effective": effective,
        "energies": energies,
        "low_band_hopping": float(abs(effective[0, 1])),
        "low_band_common_diagonal": float(np.real(effective[0, 0])),
        "operator_inventory": two_state_inventory(effective),
        "maximum_factorized_completion_residual": max(
            factorized_completion_residual(move, x) for move in moves
        ),
        "next_band_gap_over_hopping": float(
            (energies[2] - energies[1]) / abs(effective[0, 1])
        ),
        "maximum_left_gauss_commutator": float(np.max(np.abs(
            hamiltonian @ left_generator - left_generator @ hamiltonian
        ))),
        "maximum_right_gauss_commutator": float(np.max(np.abs(
            hamiltonian @ right_generator - right_generator @ hamiltonian
        ))),
        "continuity_residual": continuity,
        "spectral_ward_residual": spectral_ward,
        "diamagnetic_curvature": diamagnetic,
        "paramagnetic_curvature": float(paramagnetic),
        "total_pure_gauge_curvature": float(diamagnetic + paramagnetic),
    }


def matter_sector() -> dict[str, Any]:
    controls = []
    for coefficient in COEFFICIENTS:
        row = matter_model(0.005, coefficient)
        controls.append({
            "coefficient": coefficient,
            "common_diagonal_over_x_squared": row["low_band_common_diagonal"] / 0.005**2,
            "hopping_over_two_x_squared": row["low_band_hopping"] / (2.0 * 0.005**2),
        })
    diagonal = np.asarray([
        abs(matter_model(value, 1.0)["low_band_common_diagonal"])
        for value in SMALL_X
    ])
    hopping = np.asarray([
        matter_model(value, 1.0)["low_band_hopping"] for value in SMALL_X
    ])
    root = matter_model(X_ROOT, 1.0)
    spectra = [matter_model(X_ROOT, 1.0, phase)["energies"]
               for phase in (0.0, 0.2, 0.7, 1.3)]
    phase_spread = max(
        float(np.max(np.abs(spectrum - spectra[0]))) for spectrum in spectra
    )
    return {
        "full_hilbert_dimension": len(root["states"]),
        "physical_gauss_sector_dimension": root["physical_dimension"],
        "coefficient_controls": controls,
        "small_x_common_diagonal_fit": power_fit(SMALL_X, diagonal),
        "small_x_charged_hopping_fit": power_fit(SMALL_X, hopping),
        "root": {key: root[key] for key in (
            "low_band_hopping", "low_band_common_diagonal",
            "next_band_gap_over_hopping", "maximum_left_gauss_commutator",
            "maximum_right_gauss_commutator", "continuity_residual",
            "spectral_ward_residual", "diamagnetic_curvature",
            "paramagnetic_curvature", "total_pure_gauge_curvature",
            "operator_inventory", "maximum_factorized_completion_residual",
        )},
        "phase_twist_maximum_spectral_spread": phase_spread,
        "interpretation": (
            "A charged endpoint moves through two exact Gauss-preserving "
            "elementary steps. The same unit completion cancels its individual "
            "self-energies while leaving second-order gauge-covariant hopping."
        ),
    }


def companion_model(x: float, coefficient: float):
    base = np.diag((0.0, 1.0, 0.0)).astype(complex)
    first = np.zeros((3, 3), dtype=complex)
    second = np.zeros((3, 3), dtype=complex)
    first[1, 0] = 1.0
    second[2, 1] = 1.0
    hamiltonian = completed_hamiltonian(base, (first, second), x, coefficient)
    low_basis = np.eye(3, dtype=complex)[:, (0, 2)]
    effective, energies, eigenvectors = canonical_effective(hamiltonian, low_basis)
    times = np.linspace(0.0, 100.0, 1001)
    forward = []
    reverse = []
    for time in times:
        evolution = eigenvectors @ np.diag(np.exp(-1j * energies * time)) @ eigenvectors.conj().T
        forward.append(abs(evolution[2, 0]) ** 2)
        reverse.append(abs(evolution[0, 2]) ** 2)
    return {
        "effective": effective,
        "energies": energies,
        "conversion_hopping": float(abs(effective[0, 1])),
        "common_diagonal": float(np.real(effective[0, 0])),
        "operator_inventory": two_state_inventory(effective),
        "maximum_factorized_completion_residual": max(
            factorized_completion_residual(first, x),
            factorized_completion_residual(second, x),
        ),
        "next_band_gap_over_hopping": float(
            (energies[2] - energies[1]) / abs(effective[0, 1])
        ),
        "maximum_conversion_probability": float(max(forward)),
        "forward_reverse_probability_residual": float(
            np.max(np.abs(np.asarray(forward) - np.asarray(reverse)))
        ),
    }


def companion_sector() -> dict[str, Any]:
    controls = []
    for coefficient in COEFFICIENTS:
        row = companion_model(0.005, coefficient)
        controls.append({
            "coefficient": coefficient,
            "common_diagonal_over_x_squared": row["common_diagonal"] / 0.005**2,
            "conversion_over_x_squared": row["conversion_hopping"] / 0.005**2,
        })
    diagonal = np.asarray([
        abs(companion_model(value, 1.0)["common_diagonal"]) for value in SMALL_X
    ])
    hopping = np.asarray([
        companion_model(value, 1.0)["conversion_hopping"] for value in SMALL_X
    ])
    root = companion_model(X_ROOT, 1.0)
    return {
        "hilbert_dimension": 3,
        "coefficient_controls": controls,
        "small_x_common_diagonal_fit": power_fit(SMALL_X, diagonal),
        "small_x_conversion_fit": power_fit(SMALL_X, hopping),
        "root": {key: root[key] for key in (
            "conversion_hopping", "common_diagonal",
            "next_band_gap_over_hopping", "maximum_conversion_probability",
            "forward_reverse_probability_residual", "operator_inventory",
            "maximum_factorized_completion_residual",
        )},
        "interpretation": (
            "The neutral companion conversion uses the same two-step completed "
            "move; coefficient one cancels endpoint self-energies and preserves "
            "a reversible second-order conversion channel."
        ),
    }
