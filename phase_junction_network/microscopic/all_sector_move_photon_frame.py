"""Photon and dressed-frame sectors for Stage 3F."""
from __future__ import annotations
from typing import Any
import numpy as np
from check_seagull_completion import matched_coefficients_with_coefficient, photon_result
from all_sector_move_common import (
    X_ROOT, COEFFICIENTS, SMALL_X, canonical_effective,
    completed_hamiltonian, factorized_completion_residual, power_fit,
)


def photon_sector() -> dict[str, Any]:
    root = photon_result(X_ROOT)
    controls = []
    for coefficient in COEFFICIENTS:
        row = matched_coefficients_with_coefficient(0.005, coefficient)
        controls.append({
            "coefficient": coefficient,
            "delta_U_over_x_squared": row["delta_U_over_Delta"] / 0.005**2,
            "ring_over_20_x_fourth": row["K_plaquette_over_Delta"] / (20.0 * 0.005**4),
        })
    photon_move = np.zeros((3, 3), dtype=complex)
    photon_move[1, 0] = photon_move[2, 1] = 1.0
    return {
        "factorized_completion_residual": factorized_completion_residual(photon_move, X_ROOT),
        "root": {
            "x_t_over_Delta": X_ROOT,
            "K_A_over_Delta": root["K_plaquette_over_Delta"],
            "U_over_Delta": root["u_physical_over_Delta"],
            "v_over_K": root["v_over_K"],
            "double_flux_over_K": root["double_flux_over_K"],
            "photon_gap_power": root["photon_fit"]["gap_proportional_to_lattice_momentum_power"],
            "photon_fit_residual": root["photon_fit"]["maximum_relative_linear_plus_cubic_fit_residual"],
            "photon_speed_in_Delta_units": root["photon_speed_in_Delta_units"],
        },
        "coefficient_controls": controls,
        "interpretation": (
            "The finite photon link fixes the unit completion coefficient; "
            "the ring survives and the O(x^2) electric self-energy cancels."
        ),
    }


def clock_operators(prime: int):
    omega = np.exp(2j * np.pi / prime)
    coordinate = np.diag([omega**value for value in range(prime)])
    shift = np.zeros((prime, prime), dtype=complex)
    for value in range(prime):
        shift[(value + 1) % prime, value] = 1.0
    return shift, coordinate


def dressed_frame_model(prime: int, x: float, coefficient: float):
    shift, coordinate = clock_operators(prime)
    identity = np.eye(prime, dtype=complex)
    x_frame = np.kron(shift, identity)
    x_connection = np.kron(identity, shift)
    z_frame = np.kron(coordinate, identity)
    z_connection = np.kron(identity, coordinate)
    relative = z_connection @ z_frame.conj().T
    lock_scale = prime / (4.0 * np.pi)
    lock = lock_scale * (
        4.0 * np.eye(prime * prime)
        - x_connection - x_connection.conj().T
        - relative - relative.conj().T
    )
    dressed_shift = x_frame @ x_connection

    mediator_raise = np.asarray([[0.0, 0.0], [1.0, 0.0]], dtype=complex)
    mediator_lower = mediator_raise.conj().T
    mediator_identity = np.eye(2, dtype=complex)
    move_one = np.kron(np.eye(prime * prime), mediator_raise)
    move_two = np.kron(dressed_shift, mediator_lower)
    base = np.kron(lock, mediator_identity) + np.kron(
        np.eye(prime * prime), np.diag((0.0, 1.0))
    )
    hamiltonian = completed_hamiltonian(base, (move_one, move_two), x, coefficient)

    lock_energies, lock_states = np.linalg.eigh(lock)
    lock_mask = np.isclose(lock_energies, lock_energies[0], atol=1.0e-10)
    ground_band = lock_states[:, lock_mask]
    low_basis = np.kron(ground_band, np.asarray([[1.0], [0.0]], dtype=complex))
    effective, energies, _ = canonical_effective(hamiltonian, low_basis)
    projected_shift = ground_band.conj().T @ dressed_shift @ ground_band
    target = -(projected_shift + projected_shift.conj().T)
    target -= np.trace(target) / prime * np.eye(prime)
    centered = effective - np.trace(effective) / prime * np.eye(prime)
    hopping = float(np.real(np.vdot(target, centered) / np.vdot(target, target)))
    nearest_only_residual = float(
        np.linalg.norm(centered - hopping * target)
        / max(np.linalg.norm(centered), 1.0e-30)
    )

    inventory_operators = [np.eye(prime, dtype=complex)]
    inventory_names = ["identity"]
    for harmonic in range(1, prime // 2 + 1):
        inventory_operators.append(-(
            np.linalg.matrix_power(projected_shift, harmonic)
            + np.linalg.matrix_power(projected_shift.conj().T, harmonic)
        ))
        inventory_names.append(f"cosine_harmonic_{harmonic}")
    design = np.column_stack([operator.reshape(-1) for operator in inventory_operators])
    inventory_coefficients, *_ = np.linalg.lstsq(
        design, effective.reshape(-1), rcond=None
    )
    inventory_prediction = sum(
        value * operator
        for value, operator in zip(inventory_coefficients, inventory_operators)
    )
    inventory_residual = float(
        np.linalg.norm(effective - inventory_prediction)
        / max(np.linalg.norm(effective), 1.0e-30)
    )
    operator_inventory = {
        name: float(np.real(value))
        for name, value in zip(inventory_names, inventory_coefficients)
    }
    lock_on_total = np.kron(lock, mediator_identity)
    projector = ground_band @ ground_band.conj().T
    return {
        "prime": prime,
        "hilbert_dimension": hamiltonian.shape[0],
        "ground_band_dimension": int(lock_mask.sum()),
        "expected_ground_band_dimension": prime,
        "lock_gap": float(lock_energies[lock_mask.sum()] - lock_energies[0]),
        "maximum_lock_commutator": float(
            np.max(np.abs(hamiltonian @ lock_on_total - lock_on_total @ hamiltonian))
        ),
        "dressed_shift_lock_leakage": float(
            np.linalg.norm((np.eye(prime * prime) - projector) @ dressed_shift @ ground_band, 2)
        ),
        "dressed_weyl_residual": float(
            np.max(np.abs(
                z_frame @ dressed_shift
                - np.exp(2j * np.pi / prime) * dressed_shift @ z_frame
            ))
        ),
        "effective_hopping": hopping,
        "nearest_only_fit_residual": nearest_only_residual,
        "complete_circulant_fit_residual": inventory_residual,
        "operator_inventory": operator_inventory,
        "maximum_factorized_completion_residual": max(
            factorized_completion_residual(move_one, x),
            factorized_completion_residual(move_two, x),
        ),
        "common_energy": float(np.trace(effective).real / prime),
        "low_band_separation_over_hopping": float(
            (energies[prime] - energies[prime - 1]) / hopping
        ),
    }


def frame_sector() -> dict[str, Any]:
    rows = [dressed_frame_model(prime, X_ROOT, 1.0) for prime in (3, 5)]
    scaling = {}
    for prime in (3, 5):
        hoppings = np.asarray([
            dressed_frame_model(prime, value, 1.0)["effective_hopping"]
            for value in SMALL_X
        ])
        scaling[str(prime)] = power_fit(SMALL_X, hoppings)
    coefficient_controls = [{
        "coefficient": coefficient,
        "effective_hopping": dressed_frame_model(3, 0.01, coefficient)["effective_hopping"],
        "common_energy": dressed_frame_model(3, 0.01, coefficient)["common_energy"],
    } for coefficient in COEFFICIENTS]
    return {
        "root_rows": rows,
        "small_x_hopping_fits": scaling,
        "coefficient_controls": coefficient_controls,
        "interpretation": (
            "The completed moves commute exactly with the connection lock. "
            "Because the dressed Weyl shift is unitary, the seagull is central "
            "in this frame representation: coefficient one is compatible but "
            "the frame sector alone does not determine it."
        ),
    }
