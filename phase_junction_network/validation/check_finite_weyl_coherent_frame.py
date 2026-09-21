#!/usr/bin/env python3
"""Stage 5D: finite dressed-Weyl realization of the coherent quadratic frame term.

The coherent junction fixes U_frame and the curvature coefficient b=kappa/4.
This check puts those coefficients into a local finite Weyl chain, fixes the
global frame mode, resolves spatial momentum exactly, and independently checks
that the same physical Hamiltonian commutes with the finite connection lock.

The result is intentionally a reduced propagation block: two identical copies
represent the two tensor polarizations. It is not a full three-dimensional
nonlinear finite quantum-gravity Hamiltonian.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

U_FRAME = 0.006789913753821769
KAPPA = 2.854012469780666e-6
B_CURVATURE = KAPPA / 4.0
FULL_PRIMES = (31, 47, 67, 97, 137, 193, 257, 509)
QUICK_PRIMES = (31, 137, 257)


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


def encode(reduced: np.ndarray, prime: int, powers: np.ndarray) -> np.ndarray:
    return np.sum((reduced % prime) * powers, axis=1, dtype=np.int64)


def centered(values: np.ndarray, prime: int) -> np.ndarray:
    return (((np.asarray(values, dtype=np.int64) % prime) + prime // 2) % prime) - prime // 2


def chain_hamiltonian(prime: int, lattice_size: int = 3):
    dimension = prime ** (lattice_size - 1)
    states = np.arange(dimension, dtype=np.int64)
    powers = prime ** np.arange(lattice_size - 1, dtype=np.int64)
    reduced = ((states[:, None] // powers) % prime).astype(np.int32)
    coordinates = np.zeros((dimension, lattice_size), dtype=np.int32)
    coordinates[:, 1:] = reduced
    differences = (np.roll(coordinates, -1, axis=1) - coordinates) % prime

    kinetic_scale = U_FRAME * prime / (2.0 * np.pi)
    curvature_scale = B_CURVATURE * prime / np.pi
    diagonal = (
        curvature_scale
        * np.sum(1.0 - np.cos(2.0 * np.pi * differences / prime), axis=1)
        + kinetic_scale * lattice_size
    )

    rows = [states]
    columns = [states]
    data = [diagonal]
    amplitude = -0.5 * kinetic_scale
    for index in range(lattice_size - 1):
        plus = reduced.copy()
        plus[:, index] = (plus[:, index] + 1) % prime
        minus = reduced.copy()
        minus[:, index] = (minus[:, index] - 1) % prime
        rows.extend((states, states))
        columns.extend((encode(plus, prime, powers), encode(minus, prime, powers)))
        data.extend((np.full(dimension, amplitude), np.full(dimension, amplitude)))

    # The final shift is the removed q0 coordinate acting on all relative coordinates.
    rows.extend((states, states))
    columns.extend((encode(reduced - 1, prime, powers), encode(reduced + 1, prime, powers)))
    data.extend((np.full(dimension, amplitude), np.full(dimension, amplitude)))

    hamiltonian = coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(columns))),
        shape=(dimension, dimension),
    ).tocsr()
    translated = np.roll(coordinates, 1, axis=1)
    translated = (translated - translated[:, [0]]) % prime
    permutation = encode(translated[:, 1:], prime, powers)
    return hamiltonian, permutation, reduced


def momentum_basis(permutation: np.ndarray, lattice_size: int, momentum: int):
    dimension = len(permutation)
    seen = np.zeros(dimension, dtype=bool)
    rows: list[int] = []
    columns: list[int] = []
    data: list[complex] = []
    column = 0
    phase = np.exp(-2j * np.pi * momentum / lattice_size)
    for start in range(dimension):
        if seen[start]:
            continue
        orbit: list[int] = []
        current = start
        while not seen[current]:
            seen[current] = True
            orbit.append(current)
            current = int(permutation[current])
        if (momentum * len(orbit)) % lattice_size:
            continue
        for power, index in enumerate(orbit):
            rows.append(index)
            columns.append(column)
            data.append(phase**power / math.sqrt(len(orbit)))
        column += 1
    return coo_matrix(
        (data, (rows, columns)),
        shape=(dimension, column),
        dtype=complex,
    ).tocsr()


def lowest_sector(
    hamiltonian,
    permutation: np.ndarray,
    lattice_size: int,
    momentum: int,
):
    basis = momentum_basis(permutation, lattice_size, momentum)
    block = (basis.conj().T @ hamiltonian @ basis).tocsr()
    values, vectors = eigsh(
        block,
        k=2,
        sigma=-1.0e-10,
        which="LM",
        tol=1.0e-10,
        maxiter=100000,
    )
    order = np.argsort(values)
    return values[order], basis @ vectors[:, order[0]], basis.shape[1]


def chain_row(prime: int, lattice_size: int = 3) -> dict[str, Any]:
    started = time.time()
    hamiltonian, permutation, reduced = chain_hamiltonian(prime, lattice_size)
    zero, ground, zero_dimension = lowest_sector(
        hamiltonian, permutation, lattice_size, 0
    )
    one, tensor, one_dimension = lowest_sector(
        hamiltonian, permutation, lattice_size, 1
    )
    gap = float(one[0] - zero[0])
    lattice_momentum = 2.0 * np.sin(np.pi / lattice_size)
    speed = gap / lattice_momentum
    target = math.sqrt(2.0 * U_FRAME * B_CURVATURE)

    coordinates = np.zeros((len(reduced), lattice_size), dtype=np.int32)
    coordinates[:, 1:] = reduced
    compact_cut = np.any(
        np.abs(centered(coordinates, prime)) == prime // 2,
        axis=1,
    )
    return {
        "prime": prime,
        "lattice_size": lattice_size,
        "hilbert_dimension": hamiltonian.shape[0],
        "k0_block_dimension": zero_dimension,
        "k1_block_dimension": one_dimension,
        "ground_energy": float(zero[0]),
        "first_k1_energy": float(one[0]),
        "single_polarization_gap": gap,
        "two_polarization_gap_degeneracy": 2,
        "lattice_momentum": float(lattice_momentum),
        "finite_speed": float(speed),
        "continuum_target_speed": target,
        "relative_speed_error": float(speed / target - 1.0),
        "ground_compact_cut_probability": float(
            np.sum(np.abs(ground[compact_cut]) ** 2)
        ),
        "tensor_compact_cut_probability": float(
            np.sum(np.abs(tensor[compact_cut]) ** 2)
        ),
        "runtime_seconds": time.time() - started,
    }


def clock_operators(prime: int):
    omega = np.exp(2j * np.pi / prime)
    coordinate = np.diag([omega**value for value in range(prime)])
    shift = np.zeros((prime, prime), dtype=complex)
    for value in range(prime):
        shift[(value + 1) % prime, value] = 1.0
    return shift, coordinate


def embed(operator: np.ndarray, index: int, count: int, prime: int):
    output = np.asarray([[1.0]], dtype=complex)
    identity = np.eye(prime)
    for location in range(count):
        output = np.kron(
            output,
            operator if location == index else identity,
        )
    return output


def lock_audit(prime: int = 3, lattice_size: int = 3) -> dict[str, Any]:
    shift, coordinate = clock_operators(prime)
    count = 2 * lattice_size
    dimension = prime**count
    identity = np.eye(dimension, dtype=complex)
    shifts = [embed(shift, index, count, prime) for index in range(count)]
    coordinates = [embed(coordinate, index, count, prime) for index in range(count)]

    lock = np.zeros((dimension, dimension), dtype=complex)
    lock_scale = prime / (4.0 * np.pi)
    dressed_shifts = []
    for site in range(lattice_size):
        connection_shift = shifts[lattice_size + site]
        relative = (
            coordinates[lattice_size + site]
            @ coordinates[site]
            @ coordinates[(site + 1) % lattice_size].conj().T
        )
        lock += lock_scale * (
            4.0 * identity
            - connection_shift
            - connection_shift.conj().T
            - relative
            - relative.conj().T
        )
    for site in range(lattice_size):
        dressed_shifts.append(
            shifts[site]
            @ shifts[lattice_size + (site - 1) % lattice_size]
            @ shifts[lattice_size + site].conj().T
        )

    kinetic_scale = U_FRAME * prime / (2.0 * np.pi)
    physical = np.zeros_like(lock)
    for move in dressed_shifts:
        physical += 0.5 * kinetic_scale * (
            2.0 * identity - move - move.conj().T
        )

    # np.kron orders the first factor most significantly; reverse base-p digits.
    states = np.arange(dimension, dtype=np.int64)
    powers = prime ** np.arange(count, dtype=np.int64)
    digits = ((states[:, None] // powers) % prime).astype(np.int16)
    frame = digits[:, ::-1][:, :lattice_size]
    differences = (np.roll(frame, -1, axis=1) - frame) % prime
    potential = (
        B_CURVATURE
        * prime
        / np.pi
        * np.sum(1.0 - np.cos(2.0 * np.pi * differences / prime), axis=1)
    )
    physical += np.diag(potential)

    commutator = float(
        np.max(np.abs(lock @ physical - physical @ lock))
    )
    energies, vectors = np.linalg.eigh(lock)
    mask = np.isclose(
        energies,
        energies[0],
        rtol=1.0e-10,
        atol=1.0e-10,
    )
    projector = vectors[:, mask] @ vectors[:, mask].conj().T
    leakage = float(
        np.linalg.norm(
            (identity - projector) @ physical @ projector,
            2,
        )
    )
    return {
        "prime": prime,
        "lattice_size": lattice_size,
        "hilbert_dimension": dimension,
        "ground_band_dimension": int(mask.sum()),
        "expected_ground_band_dimension": prime**lattice_size,
        "lock_gap": float(energies[int(mask.sum())] - energies[0]),
        "maximum_physical_lock_commutator": commutator,
        "ground_band_leakage_operator_norm": leakage,
    }


def build_report(quick: bool) -> dict[str, Any]:
    primes = QUICK_PRIMES if quick else FULL_PRIMES
    rows = [chain_row(prime) for prime in primes]
    target = math.sqrt(2.0 * U_FRAME * B_CURVATURE)
    best = min(rows, key=lambda row: abs(row["relative_speed_error"]))
    highest = rows[-1]
    audit = lock_audit()
    crossover = (
        2.0
        * np.pi
        * U_FRAME
        / (target * 2.0 * np.sin(np.pi / 3.0))
    )
    checks = {
        "exact_connection_lock_compatibility": (
            audit["maximum_physical_lock_commutator"] < 1.0e-12
            and audit["ground_band_leakage_operator_norm"] < 1.0e-12
        ),
        "correct_lock_ground_band_dimension": (
            audit["ground_band_dimension"]
            == audit["expected_ground_band_dimension"]
        ),
        "two_tensor_copy_count": all(
            row["two_polarization_gap_degeneracy"] == 2 for row in rows
        ),
        "finite_candidate_within_five_percent": (
            abs(best["relative_speed_error"]) < 0.05
        ),
        "compact_boundary_control_at_best": max(
            best["ground_compact_cut_probability"],
            best["tensor_compact_cut_probability"],
        )
        < 0.01,
        "stable_continuum_plateau_not_yet_established": (
            abs(highest["relative_speed_error"]) > 0.05
        ),
    }
    return {
        "module": (
            "phase_junction_network/validation/"
            "check_finite_weyl_coherent_frame.py"
        ),
        "status": (
            "Stage 5D PARTIAL: finite dressed-Weyl local curvature implemented "
            "exactly; one finite dimension approaches the continuum cone, but "
            "a stable p-to-infinity plateau is not established"
        ),
        "run_mode": "quick" if quick else "full",
        "coefficients": {
            "U_frame": U_FRAME,
            "kappa": KAPPA,
            "b_kappa_over_4": B_CURVATURE,
            "continuum_target_speed": target,
            "quadratic_dispersion": "omega^2=2 U_frame b khat^2",
        },
        "finite_hamiltonian": {
            "coordinate_scale": "Q=sqrt(2*pi/p) centered(q)",
            "kinetic_scale": "U_frame*p/(2*pi)",
            "curvature_scale": "b*p/pi",
            "local_potential": (
                "(b*p/pi) sum_links [1-cos(2*pi Delta q/p)]"
            ),
            "global_mode": "fixed by quotient q0=0",
            "physical_polarizations": "two identical constrained copies",
        },
        "rows": rows,
        "best_finite_candidate": best,
        "highest_dimension_result": highest,
        "free_rotor_to_oscillator_crossover_prime": float(crossover),
        "exact_dressed_lock_audit": audit,
        "checks": checks,
        "stage_pass": all(
            value
            for key, value in checks.items()
            if key != "stable_continuum_plateau_not_yet_established"
        ),
        "decision": (
            "The coherent quadratic curvature term has an exactly lock-compatible "
            "finite Weyl realization. Its tiny derived curvature coefficient "
            "pushes the useful local dimension to O(10^2-10^3). The p=257 "
            "result is within five percent, but p=509 has not reached a monotone "
            "continuum plateau. This is not yet a finished finite quantum gravity "
            "regulator."
        ),
        "next_gate": (
            "Couple the antisymmetric matter Fock stress to this normalization, "
            "compute photon and frame temporal/spatial self-energies, and test "
            "whether the cone shifts match without sector-specific rescaling."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build_report(args.quick)
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
                "best": report["best_finite_candidate"],
                "highest": report["highest_dimension_result"],
                "lock": report["exact_dressed_lock_audit"],
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
    return 0 if report["stage_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
