#!/usr/bin/env python3
"""Test auxiliary-field stability alternatives for the first-order frame branch."""
from __future__ import annotations
import argparse
import json
import numpy as np

from check_first_order_frame import J, fp_matrix, gamma_map, tt_basis


def inertia(values, tolerance=1.0e-10):
    values = np.asarray(values)
    return {
        "negative": int(np.count_nonzero(values < -tolerance)),
        "zero": int(np.count_nonzero(np.abs(values) <= tolerance)),
        "positive": int(np.count_nonzero(values > tolerance)),
    }


def positive_auxiliary_no_go(samples, seed):
    """For A>0, -B^T A^-1 B cannot be a positive stiffness."""
    rng = np.random.default_rng(seed)
    maximum_positive_eigenvalue = 0.0
    failures = 0
    for _ in range(samples):
        n_aux = int(rng.integers(3, 9))
        n_frame = int(rng.integers(2, 6))
        raw = rng.normal(size=(n_aux, n_aux))
        auxiliary_metric = raw.T @ raw + np.eye(n_aux)
        coupling = rng.normal(size=(n_aux, n_frame))
        effective = -coupling.T @ np.linalg.solve(auxiliary_metric, coupling)
        largest = float(np.max(np.linalg.eigvalsh(effective)))
        maximum_positive_eigenvalue = max(maximum_positive_eigenvalue, largest)
        if largest > 1.0e-10:
            failures += 1
    return {
        "samples": samples,
        "failures": failures,
        "maximum_largest_eigenvalue": maximum_positive_eigenvalue,
        "theorem": (
            "For H=1/2 C^T A C + C^T B h with A positive definite and no "
            "bare h^2 term, eliminating C gives "
            "H_eff=-1/2 h^T B^T A^-1 B h, which is negative semidefinite."
        ),
    }


def first_order_inertia(samples, seed):
    rng = np.random.default_rng(seed)
    connection_inertia = inertia(np.linalg.eigvalsh(J))
    block_inertias = []
    fp_inertias = []
    for _ in range(samples):
        k = rng.normal(size=3)
        mapping = gamma_map(k)
        block = np.block(
            [
                [np.zeros((6, 6)), -4.0 * mapping.T @ J],
                [-4.0 * J @ mapping, 4.0 * J],
            ]
        )
        block_inertias.append(inertia(np.linalg.eigvalsh(block)))
        fp_inertias.append(inertia(np.linalg.eigvalsh(fp_matrix(k))))
    unique_blocks = sorted(
        {tuple(result[key] for key in ("negative", "zero", "positive"))
         for result in block_inertias}
    )
    unique_fp = sorted(
        {tuple(result[key] for key in ("negative", "zero", "positive"))
         for result in fp_inertias}
    )
    return {
        "samples": samples,
        "connection_quadratic_inertia": connection_inertia,
        "unique_full_first_order_block_inertias": [
            {"negative": x[0], "zero": x[1], "positive": x[2]}
            for x in unique_blocks
        ],
        "unique_unconstrained_fierz_pauli_inertias": [
            {"negative": x[0], "zero": x[1], "positive": x[2]}
            for x in unique_fp
        ],
        "interpretation": (
            "The independent connection is an auxiliary constrained variable, "
            "not an ordinary positive-energy oscillator. The three block zero "
            "modes are spatial gauge directions; the unconstrained Fierz-Pauli "
            "kernel also has one scalar negative direction removed by the "
            "scalar constraint."
        ),
    }


def stable_completion(samples, seed):
    """Show what a bounded real Hamiltonian must add."""
    rng = np.random.default_rng(seed)
    minimum_block_eigenvalue = float("inf")
    maximum_schur_error = 0.0
    failures = 0
    for _ in range(samples):
        n_aux = int(rng.integers(3, 8))
        n_frame = 2
        raw = rng.normal(size=(n_aux, n_aux))
        auxiliary_metric = raw.T @ raw + np.eye(n_aux)
        coupling = rng.normal(size=(n_aux, n_frame))
        scale = float(rng.uniform(0.1, 3.0))
        desired = scale * np.eye(n_frame)
        bare_frame = desired + coupling.T @ np.linalg.solve(
            auxiliary_metric, coupling
        )
        block = np.block(
            [
                [bare_frame, coupling.T],
                [coupling, auxiliary_metric],
            ]
        )
        schur = bare_frame - coupling.T @ np.linalg.solve(
            auxiliary_metric, coupling
        )
        minimum = float(np.min(np.linalg.eigvalsh(block)))
        error = float(np.linalg.norm(schur - desired, ord=2))
        minimum_block_eigenvalue = min(minimum_block_eigenvalue, minimum)
        maximum_schur_error = max(maximum_schur_error, error)
        if minimum < -1.0e-10 or error > 1.0e-10:
            failures += 1
    return {
        "samples": samples,
        "failures": failures,
        "minimum_full_block_eigenvalue": minimum_block_eigenvalue,
        "maximum_schur_error": maximum_schur_error,
        "interpretation": (
            "A bounded real Hamiltonian can retain an auxiliary connection only "
            "if it also contains a sufficiently positive bare frame stiffness. "
            "The connection then localizes or renormalizes the dynamics rather "
            "than generating the entire positive stiffness from a zero h-h block."
        ),
    }


def euclidean_imaginary_route(samples, seed):
    """Check the sign produced by an imaginary Hubbard-Stratonovich coupling."""
    rng = np.random.default_rng(seed)
    minimum_effective_eigenvalue = float("inf")
    failures = 0
    for _ in range(samples):
        n_aux = int(rng.integers(3, 9))
        n_frame = int(rng.integers(2, 6))
        raw = rng.normal(size=(n_aux, n_aux))
        auxiliary_metric = raw.T @ raw + np.eye(n_aux)
        coupling = rng.normal(size=(n_aux, n_frame))
        effective = coupling.T @ np.linalg.solve(auxiliary_metric, coupling)
        minimum = float(np.min(np.linalg.eigvalsh(effective)))
        minimum_effective_eigenvalue = min(minimum_effective_eigenvalue, minimum)
        if minimum < -1.0e-10:
            failures += 1
    return {
        "samples": samples,
        "failures": failures,
        "minimum_effective_eigenvalue": minimum_effective_eigenvalue,
        "identity": (
            "Integrating exp[-1/2 C^T A C - i C^T B h] over real C with A>0 "
            "produces exp[-1/2 h^T B^T A^-1 B h]."
        ),
        "limitation": (
            "The imaginary mixed term gives a complex pre-integration weight. "
            "Reflection positivity, a local transfer matrix, and the required "
            "finite gauge constraints must still be established."
        ),
    }


def physical_tt_check(samples, seed):
    rng = np.random.default_rng(seed)
    minimum_tt = float("inf")
    maximum_relative_error = 0.0
    for _ in range(samples):
        k = rng.normal(size=3)
        k_squared = float(k @ k)
        basis = tt_basis(k)
        physical = basis.T @ fp_matrix(k) @ basis
        values = np.linalg.eigvalsh(physical)
        minimum_tt = min(minimum_tt, float(np.min(values)))
        maximum_relative_error = max(
            maximum_relative_error,
            float(np.max(np.abs(values - k_squared)) / k_squared),
        )
    return {
        "samples": samples,
        "minimum_physical_tt_eigenvalue": minimum_tt,
        "maximum_relative_tt_eigenvalue_error": maximum_relative_error,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=500)
    parser.add_argument("--seed", type=int, default=44)
    parser.add_argument("--output")
    args = parser.parse_args()

    report = {
        "status": (
            "naive positive auxiliary Hamiltonian rejected; constrained-action, "
            "complex Euclidean auxiliary, and positive bare-stiffness routes separated"
        ),
        "positive_auxiliary_no_go": positive_auxiliary_no_go(
            args.samples, args.seed
        ),
        "actual_first_order_inertia": first_order_inertia(
            args.samples, args.seed + 1
        ),
        "physical_tt_sector": physical_tt_check(
            args.samples, args.seed + 2
        ),
        "bounded_real_hamiltonian_completion": stable_completion(
            args.samples, args.seed + 3
        ),
        "euclidean_imaginary_auxiliary_route": euclidean_imaginary_route(
            args.samples, args.seed + 4
        ),
        "decision": (
            "Do not quantize the 18 connection components as ordinary real "
            "positive-energy oscillators with a zero bare h-h block. Advance "
            "three separately labeled candidates: (A) a constrained first-order "
            "phase-space/Palatini transfer matrix, (B) a complex Euclidean "
            "auxiliary representation with a reflection-positivity test, and "
            "(C) a bounded Hamiltonian with explicit bare frame stiffness."
        ),
    }

    assert report["positive_auxiliary_no_go"]["failures"] == 0
    assert report["bounded_real_hamiltonian_completion"]["failures"] == 0
    assert report["euclidean_imaginary_auxiliary_route"]["failures"] == 0
    assert report["physical_tt_sector"]["maximum_relative_tt_eigenvalue_error"] < 1.0e-12

    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")


if __name__ == "__main__":
    main()
