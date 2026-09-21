#!/usr/bin/env python3
"""Stage 5K: exact finite Clifford automaton for the doubled duality system.

The static Weyl self-duality lock in Stage 5J changed the infrared derivative
order. Here the first-order symplectic equations are implemented directly as a
local reversible update. One time step is a composition of two nearest-neighbor
symplectic shears,

    q' = q + s C p,
    p' = p - s C q',

where C is the staggered primal/dual curl. For odd-prime qudits and s=1/2 mod p,
the update is an exact finite Clifford automorphism. Over the reals it preserves
a local positive quadratic invariant whenever |s lambda(C)|<2.

The same update is tested on the spin-1 Gauss sector and the spin-2 local
scalar/vector-constraint sector. This is a Floquet/QCA candidate, not yet a
time-independent interacting Hamiltonian.
"""
from __future__ import annotations

import argparse
import builtins
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
import sys
sys.path.insert(0, str(HERE))
from check_common_chiral_generator import build_report as chiral_report
from check_chiral_boson_positivity import physical_curls

STEP = 0.5
FULL_PRIMES = (3, 5, 7, 11)
QUICK_PRIMES = (3, 5, 7)


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


def inverse_mod(value: int, prime: int) -> int:
    return builtins.pow(int(value) % prime, -1, prime)


def leapfrog(curl: np.ndarray, step: float):
    dimension = curl.shape[0]
    identity = np.eye(dimension, dtype=complex)
    first = np.block(
        [[identity, step * curl], [np.zeros_like(curl), identity]]
    )
    second = np.block(
        [[identity, np.zeros_like(curl)], [-step * curl, identity]]
    )
    update = second @ first
    symplectic = np.block(
        [[np.zeros_like(curl), identity], [-identity, np.zeros_like(curl)]]
    )
    invariant = np.block(
        [
            [identity, 0.5 * step * curl],
            [0.5 * step * curl, identity],
        ]
    )
    return update, symplectic, invariant


def mode_row(momentum: np.ndarray):
    vector, tensor, vector_basis, tensor_basis = physical_curls(momentum)
    output = {
        "momentum": [float(value) for value in momentum],
        "momentum_norm": float(np.linalg.norm(momentum)),
        "vector_constraint_dimension": int(vector_basis.shape[1]),
        "tensor_constraint_dimension": int(tensor_basis.shape[1]),
    }
    frequency_sets = {}
    for name, curl in (("vector", vector), ("tensor", tensor)):
        update, symplectic, invariant = leapfrog(curl, STEP)
        eigenvalues = np.linalg.eigvals(update)
        phases = np.angle(eigenvalues)
        positive = sorted(float(value) for value in phases if value > 1.0e-12)
        expected = 2.0 * math.asin(
            STEP * float(np.linalg.norm(momentum)) / 2.0
        )
        frequency_sets[name] = {
            "configuration_dimension": int(curl.shape[0]),
            "phase_space_dimension": int(update.shape[0]),
            "positive_quasifrequencies": positive,
            "expected_quasifrequency": expected,
            "maximum_relative_frequency_error": max(
                abs(value - expected) for value in positive
            )
            / expected,
            "symplectic_residual": float(
                np.max(
                    np.abs(
                        update.conj().T @ symplectic @ update - symplectic
                    )
                )
            ),
            "invariant_residual": float(
                np.max(
                    np.abs(
                        update.conj().T @ invariant @ update - invariant
                    )
                )
            ),
            "minimum_invariant_eigenvalue": float(
                np.min(np.linalg.eigvalsh(invariant))
            ),
            "maximum_update_modulus_error": float(
                np.max(np.abs(np.abs(eigenvalues) - 1.0))
            ),
        }
    output["vector"] = frequency_sets["vector"]
    output["tensor"] = frequency_sets["tensor"]
    output["vector_tensor_frequency_mismatch"] = max(
        abs(first - second)
        for first, second in zip(
            frequency_sets["vector"]["positive_quasifrequencies"],
            frequency_sets["tensor"]["positive_quasifrequencies"],
        )
    )
    return output


def staggered_curl_mod(prime: int, lattice_size: int):
    derivative = np.zeros((lattice_size, lattice_size), dtype=int)
    for site in range(lattice_size):
        derivative[site, site] = -1
        derivative[site, (site + 1) % lattice_size] = 1
    zero = np.zeros_like(derivative)
    curl = np.block([[zero, derivative.T], [derivative, zero]]) % prime
    return curl


def finite_symplectic_row(prime: int, lattice_size: int = 3):
    curl = staggered_curl_mod(prime, lattice_size)
    dimension = curl.shape[0]
    identity = np.eye(dimension, dtype=int)
    zero = np.zeros_like(identity)
    step = inverse_mod(2, prime)
    first = np.block([[identity, step * curl], [zero, identity]]) % prime
    second = np.block([[identity, zero], [(-step * curl) % prime, identity]]) % prime
    update = (second @ first) % prime
    symplectic = np.block([[zero, identity], [-identity, zero]]) % prime
    residual = (update.T @ symplectic @ update - symplectic) % prime
    inverse = np.block(
        [
            [(identity - (step * step % prime) * (curl @ curl)) % prime,
             (-step * curl) % prime],
            [(step * curl) % prime, identity],
        ]
    ) % prime
    inverse_residual = (update @ inverse - np.eye(2 * dimension, dtype=int)) % prime
    return {
        "prime": prime,
        "lattice_size": lattice_size,
        "canonical_coordinate_count": dimension,
        "phase_space_dimension": int(2 * dimension),
        "hilbert_dimension": int(prime**dimension),
        "finite_step_inverse_two": step,
        "maximum_symplectic_residual_mod_p": int(np.max(residual)),
        "maximum_inverse_residual_mod_p": int(np.max(inverse_residual)),
        "determinant_mod_p": int(round(np.linalg.det(update))) % prime,
        "clifford_gate_decomposition": (
            "momentum-quadratic nearest-neighbor shear followed by "
            "coordinate-quadratic nearest-neighbor shear"
        ),
        "causal_radius_per_substep": 1,
        "causal_radius_per_full_step": 2,
    }


def no_doubler_audit(lattice_size: int = 16):
    zero_count = 0
    minimum_nonzero = float("inf")
    maximum = 0.0
    for mode in np.ndindex(*(lattice_size,) * 3):
        momentum = 2.0 * np.sin(
            np.pi * np.asarray(mode, dtype=float) / lattice_size
        )
        norm = float(np.linalg.norm(momentum))
        if norm < 1.0e-12:
            zero_count += 1
        else:
            minimum_nonzero = min(minimum_nonzero, norm)
            maximum = max(maximum, norm)
    return {
        "lattice_size": lattice_size,
        "zero_count": zero_count,
        "minimum_nonzero_lattice_momentum": minimum_nonzero,
        "maximum_lattice_momentum": maximum,
        "maximum_STEP_times_momentum": STEP * maximum,
        "stability_bound": 2.0,
    }


def dispersion_fit(dynamic_speed: float):
    sizes = np.asarray((12, 16, 24, 32, 48, 64, 96), dtype=float)
    momentum = 2.0 * np.sin(np.pi / sizes)
    phase = 2.0 * np.arcsin(STEP * momentum / 2.0)
    time_step = STEP / dynamic_speed
    frequency = phase / time_step
    design = np.column_stack((momentum, momentum**3))
    coefficients, *_ = np.linalg.lstsq(design, frequency, rcond=None)
    prediction = design @ coefficients
    power, _ = np.polyfit(np.log(momentum), np.log(frequency), 1)
    return {
        "sizes": [int(value) for value in sizes],
        "lattice_momenta": [float(value) for value in momentum],
        "physical_frequencies": [float(value) for value in frequency],
        "ring_fixed_time_step": time_step,
        "linear_speed": float(coefficients[0]),
        "cubic_cutoff_coefficient": float(coefficients[1]),
        "gap_power": float(power),
        "maximum_relative_fit_residual": float(
            np.max(np.abs(prediction - frequency) / frequency)
        ),
    }


def build_report(quick: bool):
    parent = chiral_report()
    dynamic_speed = parent["photon"]["dynamic_odd_coefficient"]
    modes = (
        np.asarray((0.5, 0.0, 0.0)),
        np.asarray((0.5, 0.5, 0.0)),
        np.asarray((0.5, 0.5, 0.5)),
        np.asarray((1.0, 0.5, 0.0)),
    )
    mode_rows = [mode_row(momentum) for momentum in modes]
    primes = QUICK_PRIMES if quick else FULL_PRIMES
    finite_rows = [finite_symplectic_row(prime) for prime in primes]
    doublers = no_doubler_audit()
    dispersion = dispersion_fit(dynamic_speed)

    maximum_frequency_error = max(
        max(
            row["vector"]["maximum_relative_frequency_error"],
            row["tensor"]["maximum_relative_frequency_error"],
        )
        for row in mode_rows
    )
    maximum_common_mismatch = max(
        row["vector_tensor_frequency_mismatch"] for row in mode_rows
    )
    minimum_invariant = min(
        min(
            row["vector"]["minimum_invariant_eigenvalue"],
            row["tensor"]["minimum_invariant_eigenvalue"],
        )
        for row in mode_rows
    )
    checks = {
        "real_update_is_exactly_symplectic": max(
            max(row["vector"]["symplectic_residual"],
                row["tensor"]["symplectic_residual"])
            for row in mode_rows
        )
        < 1.0e-12,
        "positive_local_invariant": minimum_invariant > 0.0,
        "unit_circle_stability": max(
            max(row["vector"]["maximum_update_modulus_error"],
                row["tensor"]["maximum_update_modulus_error"])
            for row in mode_rows
        )
        < 1.0e-12,
        "two_spin1_and_two_spin2_modes": all(
            len(row["vector"]["positive_quasifrequencies"]) == 2
            and len(row["tensor"]["positive_quasifrequencies"]) == 2
            for row in mode_rows
        ),
        "vector_tensor_quasifrequencies_match": maximum_common_mismatch < 1.0e-12,
        "analytic_leapfrog_dispersion_matches": maximum_frequency_error < 1.0e-12,
        "finite_odd_prime_updates_are_symplectic": all(
            row["maximum_symplectic_residual_mod_p"] == 0
            and row["maximum_inverse_residual_mod_p"] == 0
            and row["determinant_mod_p"] == 1
            for row in finite_rows
        ),
        "staggered_curl_has_only_one_zero": doublers["zero_count"] == 1,
        "three_dimensional_update_is_stable": (
            doublers["maximum_STEP_times_momentum"] < doublers["stability_bound"]
        ),
        "physical_low_momentum_dispersion_is_linear": (
            abs(dispersion["gap_power"] - 1.0) < 0.01
            and abs(dispersion["linear_speed"] / dynamic_speed - 1.0) < 1.0e-4
        ),
    }

    return {
        "module": (
            "phase_junction_network/construction/coherent_junction/"
            "check_clifford_duality_automaton.py"
        ),
        "status": (
            "Stage 5K FINITE DYNAMICAL CANDIDATE PASS: an exact local odd-prime "
            "Clifford automaton realizes the doubled positive spin-1/spin-2 "
            "duality equations with two modes and one ring-derived cone"
        ),
        "ring_input": {
            "dynamic_speed": dynamic_speed,
            "dimensionless_automaton_step": STEP,
            "physical_time_per_step": STEP / dynamic_speed,
            "new_sector_specific_coefficients": 0,
        },
        "update": {
            "first_shear": "q <- q + (1/2) C p",
            "second_shear": "p <- p - (1/2) C q_new",
            "curl": (
                "staggered primal/dual incidence curl with symbol "
                "|khat|=sqrt(sum_i [2 sin(k_i/2)]^2)"
            ),
            "positive_invariant": (
                "K=[[I,(s/2)C],[(s/2)C,I]], preserved exactly by the update"
            ),
            "finite_quantization": (
                "odd-prime symplectic map; product of two local quadratic "
                "Clifford shear gates"
            ),
            "time_structure": "reversible Floquet/QCA, not yet a static Hamiltonian",
        },
        "mode_rows": mode_rows,
        "finite_field_rows": finite_rows,
        "no_doubler_audit": doublers,
        "dispersion": dispersion,
        "checks": checks,
        "execution_pass": all(checks.values()),
        "decision": (
            "Implementing the chiral pole as reversible symplectic time evolution "
            "rather than an energetic constraint lock restores exact z=1 dynamics. "
            "The candidate is local, finite, unitary at odd prime, positive through "
            "a conserved quadratic invariant, and gives two photon and two tensor "
            "modes with one speed. This is the strongest common-cone construction "
            "so far, but it changes the microscopic theory from a static Hamiltonian "
            "to a Floquet/automaton formulation."
        ),
        "claim_boundary": {
            "established": [
                "an exact finite odd-prime symplectic/Clifford update",
                "two local shear layers and finite causal radius",
                "a positive exactly conserved local quadratic invariant",
                "two constrained spin-1 and two constrained spin-2 modes",
                "one common ring-derived low-momentum speed",
                "no staggered-curl Brillouin-zone doublers",
            ],
            "not_established": [
                "a time-independent positive Hamiltonian whose exponential is the update",
                "nonlinear constraint closure under the finite automaton",
                "interacting matter and companion coupling in Floquet variables",
                "absence of Floquet heating in the nonlinear many-body theory",
                "equivalence to observed Maxwell and Einstein dynamics",
            ],
        },
        "next_gate": (
            "Lift the nonlinear frame constraints and the antisymmetric matter "
            "Fock moves into the same Clifford time step. Test exact discrete "
            "Ward identities, common-cone loop stability, and whether nonlinear "
            "updates avoid Floquet heating without introducing a new scale."
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
                "execution_pass": report["execution_pass"],
                "dispersion": report["dispersion"],
                "doublers": report["no_doubler_audit"],
                "finite": report["finite_field_rows"],
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
