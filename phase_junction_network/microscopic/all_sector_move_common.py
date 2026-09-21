"""Shared utilities for Stage 3F all-sector completed-move checks."""
from __future__ import annotations
from typing import Any
import numpy as np

X_ROOT = 0.13554178509861228
COEFFICIENTS = (0.8, 0.9, 1.0, 1.1, 1.2)
SMALL_X = np.asarray((0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05))


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


def power_fit(xs: np.ndarray, ys: np.ndarray) -> dict[str, float]:
    power, log_prefactor = np.polyfit(np.log(xs), np.log(np.abs(ys)), 1)
    return {"power": float(power), "prefactor": float(np.exp(log_prefactor))}


def canonical_effective(hamiltonian: np.ndarray, low_basis: np.ndarray):
    count = low_basis.shape[1]
    energies, eigenvectors = np.linalg.eigh(hamiltonian)
    overlap = low_basis.conj().T @ eigenvectors[:, :count]
    left, _, right = np.linalg.svd(overlap)
    polar = left @ right
    effective = polar @ np.diag(energies[:count]) @ polar.conj().T
    return effective, energies, eigenvectors


def partial_operator(states, lookup, transition):
    operator = np.zeros((len(states), len(states)), dtype=complex)
    for state in states:
        target = transition(state)
        if target is not None and target in lookup:
            operator[lookup[target], lookup[state]] = 1.0
    return operator


def completed_hamiltonian(base, moves, x, coefficient=1.0, phases=None):
    output = np.asarray(base, dtype=complex).copy()
    phases = phases or [0.0] * len(moves)
    for move, phase in zip(moves, phases):
        output += -x * (
            np.exp(1j * phase) * move
            + np.exp(-1j * phase) * move.conj().T
        )
        output += coefficient * x * x * (
            move.conj().T @ move + move @ move.conj().T
        )
    return output


def factorized_completion_residual(move: np.ndarray, x: float) -> float:
    """Verify the universal completed move is a symmetric positive square."""
    identity = np.eye(move.shape[0], dtype=complex)
    q = identity - 2.0 * x * move
    factorized = 0.25 * (
        q.conj().T @ q + q @ q.conj().T - 2.0 * identity
    )
    completed = -x * (move + move.conj().T) + x * x * (
        move.conj().T @ move + move @ move.conj().T
    )
    return float(np.max(np.abs(factorized - completed)))


def two_state_inventory(effective: np.ndarray) -> dict[str, float]:
    """Pauli-basis inventory for a two-state low band."""
    identity = 0.5 * np.trace(effective)
    sigma_x = 0.5 * (effective[0, 1] + effective[1, 0])
    sigma_y = (effective[1, 0] - effective[0, 1]) / (2.0j)
    sigma_z = 0.5 * (effective[0, 0] - effective[1, 1])
    return {
        "identity": float(np.real(identity)),
        "sigma_x": float(np.real(sigma_x)),
        "sigma_y": float(np.real(sigma_y)),
        "sigma_z": float(np.real(sigma_z)),
        "maximum_imaginary_coefficient": float(
            max(
                abs(np.imag(identity)),
                abs(np.imag(sigma_x)),
                abs(np.imag(sigma_y)),
                abs(np.imag(sigma_z)),
            )
        ),
    }
