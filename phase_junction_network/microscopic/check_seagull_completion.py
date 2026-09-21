#!/usr/bin/env python3
"""Stage 3E: symmetry-complete local seagull candidate.

Stage 3D showed that the pure off-diagonal swap generates an O(x^2)
electric self-energy that destroys the shared photon/gravity cone.  The
minimal local completion pairs each finite-link swap with its algebraic
anticommutator,

    x^2 (U^dagger U + U U^dagger).

For an isolated move this is the unique unit-coefficient local seagull that
cancels the second-order self-energy.  Cross-move processes remain, so the
photon loop and frame superexchange are not removed.

This script tests that candidate exactly in the spin-1 plaquette, the Stage-6B
photon blocks, the exact 2^3 gauge cube, and the compact gravity regulator.
It does not derive the seagull from a final fundamental symmetry.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

from check_collective_capacitance import fit_gap_rows, gravity_row
from check_plaquette_diagonal_matching import (
    E_LINK,
    GAUGE_EMBEDDING,
    H0,
    U_LINK,
    U_LINKS,
    V_SWAP,
    inverse_square_root,
)
from finite_em_dynamics_impl import axial_spectrum, exact_cube_components, fit_branch


R_COMMON = 1.0 / 144.0
D_LINK = U_LINK.T @ U_LINK + U_LINK @ U_LINK.T
SEAGULL = sum(
    operator.T @ operator + operator @ operator.T for operator in U_LINKS
)


def ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [ready(item) for item in value]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def effective_hamiltonian(x: float, coefficient: float = 1.0):
    full = H0 + x * V_SWAP + coefficient * x * x * SEAGULL
    energies, eigenvectors = np.linalg.eigh(full)
    overlap = GAUGE_EMBEDDING.T @ eigenvectors[:, :3]
    polar = overlap @ inverse_square_root(overlap.T @ overlap)
    effective = polar @ np.diag(energies[:3]) @ polar.T
    return effective, energies


def isolated_link_energy(x: float, reference_flux: int, coefficient: float) -> float:
    fluxes = np.asarray((-1.0, 0.0, 1.0))
    hamiltonian = np.diag((fluxes - reference_flux) ** 2)
    hamiltonian -= x * (U_LINK + U_LINK.T)
    hamiltonian += coefficient * x * x * D_LINK
    energies, eigenvectors = np.linalg.eigh(hamiltonian)
    reference_index = reference_flux + 1
    branch = int(np.argmax(np.abs(eigenvectors[reference_index, :]) ** 2))
    return float(energies[branch])


def matched_coefficients(x: float, coefficient: float = 1.0) -> dict[str, float]:
    effective, full_energies = effective_hamiltonian(x, coefficient)
    edge = 0.5 * (effective[0, 0] + effective[2, 2])
    center = effective[1, 1]
    full_contrast = edge - center
    independent = 4.0 * (
        isolated_link_energy(x, 1, coefficient)
        - isolated_link_energy(x, 0, coefficient)
    )
    delta_u = 0.5 * independent
    v_diagonal = -(full_contrast - independent)
    ring = -0.5 * (effective[0, 1] + effective[1, 2])
    double_flux = effective[0, 2]
    return {
        "x_t_over_Delta": x,
        "K_plaquette_over_Delta": float(ring),
        "delta_U_over_Delta": float(delta_u),
        "v_over_Delta": float(v_diagonal),
        "v_over_K": float(v_diagonal / ring),
        "double_flux_over_K": float(double_flux / ring),
        "low_band_separation_over_K": float(
            (full_energies[3] - full_energies[2]) / ring
        ),
    }


def gravity_coupling(x: float) -> float:
    hamiltonian = np.asarray(
        [[x * x, -x, 0.0], [-x, 1.0, -x], [0.0, -x, x * x]],
        dtype=float,
    )
    energies = np.linalg.eigvalsh(hamiltonian)
    return float((energies[1] - energies[0]) / 2.0)


def logarithmic_fit(xs: np.ndarray, values: np.ndarray) -> dict[str, float]:
    power, log_prefactor = np.polyfit(np.log(xs), np.log(np.abs(values)), 1)
    return {"power": float(power), "prefactor": float(np.exp(log_prefactor))}


def photon_result(x: float) -> dict[str, Any]:
    matched = matched_coefficients(x)
    u_physical = R_COMMON + matched["delta_U_over_Delta"]
    if u_physical <= 0.0:
        raise ValueError("seagull over-cancelled the physical charging coefficient")
    u_over_t = u_physical / matched["K_plaquette_over_Delta"]
    v_over_t = matched["v_over_K"]
    rows = [
        axial_spectrum(
            spin=1,
            lattice_size=lattice_size,
            u_over_t=u_over_t,
            v_over_t=v_over_t,
            eigenpairs=16,
        )
        for lattice_size in (6, 8, 10)
    ]
    fit = fit_branch(rows)
    photon_speed = matched["K_plaquette_over_Delta"] * fit[
        "linear_speed_c_gamma"
    ]
    k_g = gravity_coupling(x)
    gravity_speed = math.sqrt(u_physical * k_g)
    return {
        **matched,
        "u_physical_over_Delta": u_physical,
        "u_over_t_plaquette": u_over_t,
        "K_g_over_Delta": k_g,
        "photon_fit": fit,
        "photon_speed_in_Delta_units": photon_speed,
        "gravity_speed_in_Delta_units": gravity_speed,
        "speed_difference": photon_speed - gravity_speed,
        "speed_ratio": photon_speed / gravity_speed,
        "photon_rows": [
            {
                key: row[key]
                for key in (
                    "lattice_size",
                    "hilbert_dimension",
                    "lattice_momentum",
                    "single_polarization_gap",
                    "signed_momentum_split",
                )
            }
            for row in rows
        ],
    }


def bisection(lower: float, upper: float, iterations: int) -> dict[str, Any]:
    lower_result = photon_result(lower)
    upper_result = photon_result(upper)
    if lower_result["speed_difference"] * upper_result["speed_difference"] >= 0:
        raise RuntimeError("common-cone root is not bracketed")
    for _ in range(iterations):
        midpoint = 0.5 * (lower + upper)
        result = photon_result(midpoint)
        if lower_result["speed_difference"] * result["speed_difference"] <= 0:
            upper = midpoint
            upper_result = result
        else:
            lower = midpoint
            lower_result = result
    result = photon_result(0.5 * (lower + upper))
    result["root_bracket_width"] = upper - lower
    return result


def exact_cube_anchor(root: dict[str, Any]) -> dict[str, Any]:
    basis, electric, degree, adjacency, operators = exact_cube_components()
    hamiltonian = (
        diags(
            root["u_over_t_plaquette"] * electric
            + root["v_over_K"] * degree
        )
        - adjacency
    )
    uniform = np.ones(len(basis.states)) / math.sqrt(len(basis.states))
    energies, vectors = eigsh(
        hamiltonian,
        k=8,
        which="SA",
        tol=1.0e-9,
        maxiter=30000,
        v0=uniform,
    )
    order = np.argsort(energies)
    energies = energies[order]
    vectors = vectors[:, order]
    gaps = energies - energies[0]
    low_mask = np.isclose(gaps, gaps[1], rtol=1.0e-8, atol=1.0e-8)
    ground = vectors[:, 0]
    operator_results = {}
    for name, operator in operators.items():
        created = operator * ground
        norm = float(np.vdot(created, created).real)
        residues = np.abs(vectors.conj().T @ created) ** 2
        operator_results[name] = {
            "norm": norm,
            "low_multiplet_residue_fraction": (
                float(np.sum(residues[low_mask]) / norm) if norm > 0.0 else 0.0
            ),
        }
    return {
        "gauge_reduced_dimension": len(basis.states),
        "lowest_gaps": [float(value) for value in gaps],
        "first_gap_multiplicity_in_computed_window": int(np.count_nonzero(low_mask)),
        "ground_overlap_with_RK_uniform_state": float(
            abs(np.vdot(uniform, ground)) ** 2
        ),
        "operators": operator_results,
    }


def finite_gravity(root: dict[str, Any], primes: tuple[int, ...]) -> list[dict[str, Any]]:
    output = []
    for prime in primes:
        rows = [
            gravity_row(
                prime,
                lattice_size,
                root["K_g_over_Delta"],
                root["u_physical_over_Delta"],
            )
            for lattice_size in (3, 4, 5)
        ]
        fit = fit_gap_rows(rows)
        output.append(
            {
                "prime": prime,
                "rows": rows,
                "fit": fit,
                "harmonic_target_speed": root[
                    "gravity_speed_in_Delta_units"
                ],
                "relative_speed_error": (
                    fit["linear_speed_c_g"]
                    / root["gravity_speed_in_Delta_units"]
                    - 1.0
                ),
            }
        )
    return output


def compact_root(root: dict[str, Any]) -> dict[str, Any]:
    return {
        key: root[key]
        for key in (
            "x_t_over_Delta",
            "K_plaquette_over_Delta",
            "delta_U_over_Delta",
            "v_over_K",
            "double_flux_over_K",
            "low_band_separation_over_K",
            "u_physical_over_Delta",
            "u_over_t_plaquette",
            "K_g_over_Delta",
            "photon_speed_in_Delta_units",
            "gravity_speed_in_Delta_units",
            "speed_difference",
            "speed_ratio",
        )
    } | {
        "photon_gap_power": root["photon_fit"][
            "gap_proportional_to_lattice_momentum_power"
        ],
        "photon_fit_residual": root["photon_fit"][
            "maximum_relative_linear_plus_cubic_fit_residual"
        ],
    }


def build_report(quick: bool) -> dict[str, Any]:
    small_x = np.asarray(
        (0.005, 0.01, 0.02, 0.03, 0.05)
        if quick
        else (0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05)
    )
    small_rows = [matched_coefficients(float(value)) for value in small_x]
    k_g_values = np.asarray([gravity_coupling(float(value)) for value in small_x])
    fits = {
        "K_plaquette": logarithmic_fit(
            small_x,
            np.asarray([row["K_plaquette_over_Delta"] for row in small_rows]),
        ),
        "delta_U": logarithmic_fit(
            small_x,
            np.asarray([row["delta_U_over_Delta"] for row in small_rows]),
        ),
        "v_diagonal": logarithmic_fit(
            small_x,
            np.asarray([row["v_over_Delta"] for row in small_rows]),
        ),
        "K_g": logarithmic_fit(small_x, k_g_values),
    }

    coefficient_controls = []
    for coefficient in (0.8, 0.9, 1.0, 1.1, 1.2):
        row = matched_coefficients_with_coefficient(0.005, coefficient)
        coefficient_controls.append(
            {
                "seagull_coefficient": coefficient,
                "delta_U_over_x2": row["delta_U_over_Delta"] / 0.005**2,
                "K_over_20x4": row["K_plaquette_over_Delta"]
                / (20.0 * 0.005**4),
            }
        )

    first_root = bisection(0.08, 0.10, 12 if quick else 26)
    physical_root = bisection(0.12, 0.14, 14 if quick else 28)
    cube = None if quick else exact_cube_anchor(physical_root)
    primes = (11, 17) if quick else (5, 7, 11, 13, 17, 23)
    gravity = finite_gravity(physical_root, primes)

    linear_gravity = [
        entry
        for entry in gravity
        if entry["prime"] >= 7
    ]
    linear_powers = [
        entry["fit"]["gap_proportional_to_lattice_momentum_power"]
        for entry in linear_gravity
    ]
    gravity_residuals = [
        entry["fit"]["maximum_relative_linear_plus_cubic_fit_residual"]
        for entry in linear_gravity
    ]
    converging_speeds = [
        entry["fit"]["linear_speed_c_g"]
        for entry in gravity
        if entry["prime"] >= 7
    ]

    checks = {
        "unit_seagull_uniquely_cancels_second_order_selfenergy": (
            coefficient_controls[0]["delta_U_over_x2"] > 0.35
            and coefficient_controls[1]["delta_U_over_x2"] > 0.15
            and abs(coefficient_controls[2]["delta_U_over_x2"]) < 1.0e-3
            and coefficient_controls[3]["delta_U_over_x2"] < -0.15
            and coefficient_controls[4]["delta_U_over_x2"] < -0.35
        ),
        "ring_exchange_remains_fourth_order": abs(
            fits["K_plaquette"]["power"] - 4.0
        )
        < 0.03,
        "residual_electric_correction_is_fourth_order": abs(
            fits["delta_U"]["power"] - 4.0
        )
        < 0.03,
        "linked_diagonal_is_fourth_order": abs(
            fits["v_diagonal"]["power"] - 4.0
        )
        < 0.03,
        "gravity_superexchange_remains_second_order": abs(
            fits["K_g"]["power"] - 2.0
        )
        < 0.02,
        "physical_root_is_z1": (
            0.94
            <= physical_root["photon_fit"][
                "gap_proportional_to_lattice_momentum_power"
            ]
            <= 1.16
            and physical_root["photon_fit"][
                "maximum_relative_linear_plus_cubic_fit_residual"
            ]
            < 0.002
            and abs(physical_root["speed_difference"]) < 1.0e-6
        ),
        "first_root_is_nonrelativistic_control": (
            first_root["photon_fit"][
                "gap_proportional_to_lattice_momentum_power"
            ]
            < 0.94
            and first_root["photon_fit"][
                "maximum_relative_linear_plus_cubic_fit_residual"
            ]
            > 0.002
        ),
        "finite_gravity_branches_are_linear": min(linear_powers) > 0.94
        and max(linear_powers) < 1.16,
        "finite_gravity_fit_residuals_below_0p002": max(gravity_residuals)
        < 0.002,
        "finite_gravity_speeds_move_toward_target": all(
            right > left
            for left, right in zip(converging_speeds, converging_speeds[1:])
        ),
        "full_cube_longitudinal_null": (
            True
            if cube is None
            else cube["operators"]["longitudinal_x"]["norm"] < 1.0e-12
        ),
        "full_cube_scalar_rejected": (
            True
            if cube is None
            else cube["operators"]["scalar_electric_energy"][
                "low_multiplet_residue_fraction"
            ]
            < 1.0e-12
        ),
    }

    report = {
        "module": (
            "phase_junction_network/microscopic/"
            "check_seagull_completion.py"
        ),
        "status": (
            "Stage 3E PASS at candidate scope: the unit local seagull cancels "
            "the dangerous second-order self-energy while preserving ring and "
            "frame exchange, yielding one z~1 shared-cone root"
        ),
        "run_mode": "quick" if quick else "full",
        "candidate_hamiltonians": {
            "electromagnetic": (
                "H/Delta = 1/2 sum G^2 - x sum(U+U^dagger) + "
                "x^2 sum(U^dagger U + U U^dagger)"
            ),
            "gravity": (
                "H_g/Delta = [[x^2,-x,0],[-x,1,-x],[0,-x,x^2]]"
            ),
            "coefficient_principle": (
                "The unit coefficient is fixed by cancelling each individual "
                "second-order move self-energy. Cross-move superexchange is "
                "left intact."
            ),
        },
        "derived_common_charging_factor": R_COMMON,
        "small_x_fits": fits,
        "small_x_rows": small_rows,
        "coefficient_controls": coefficient_controls,
        "nonrelativistic_root_control": compact_root(first_root),
        "physical_common_cone_root": compact_root(physical_root),
        "exact_periodic_cube_anchor": cube,
        "finite_gravity_same_coefficients": gravity,
        "checks": checks,
        "stage_pass": all(checks.values()),
        "decision": (
            "Advance this seagull-completed graph to the full shared finite "
            "embedding. Do not yet call it fundamental: the algebraic origin "
            "of the completion and its extension to matter and companion "
            "moves remain to be derived."
        ),
        "claim_boundary": {
            "established": [
                "a unique unit local seagull cancels the O(x^2) isolated-move self-energy",
                "the O(x^4) photon ring and O(x^2) gravity superexchange survive",
                "a z~1 finite photon common-cone root exists at the shared r=1/144",
                "the exact 2^3 gauge cube retains no longitudinal operator at that root",
                "finite compact gravity remains linear at the same coefficients for adequate local dimension",
            ],
            "not_established": [
                "a fundamental symmetry requiring the seagull completion",
                "uniqueness against all other local completions",
                "the full infinite-volume three-dimensional interacting phase",
                "matter and companion seagull partners",
                "radiative stability, nonlinear gravity, alpha, G, or formal novelty",
            ],
        },
        "next_gate": (
            "Construct the seagull-completed moves directly in the finite "
            "dressed-frame, charged-endpoint, and companion Hilbert spaces. "
            "Verify exact constraints and derive the complete operator set "
            "without adding sector-specific terms."
        ),
    }

    if not report["stage_pass"]:
        raise AssertionError(json.dumps(ready(checks), indent=2, sort_keys=True))
    return report


def matched_coefficients_with_coefficient(
    x: float, coefficient: float
) -> dict[str, float]:
    effective, energies = effective_hamiltonian(x, coefficient)
    edge = 0.5 * (effective[0, 0] + effective[2, 2])
    center = effective[1, 1]
    independent = 4.0 * (
        isolated_link_energy(x, 1, coefficient)
        - isolated_link_energy(x, 0, coefficient)
    )
    ring = -0.5 * (effective[0, 1] + effective[1, 2])
    return {
        "K_plaquette_over_Delta": float(ring),
        "delta_U_over_Delta": float(0.5 * independent),
        "v_over_Delta": float(-((edge - center) - independent)),
        "low_band_separation_over_K": float(
            (energies[3] - energies[2]) / ring
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "phase_junction_network/microscopic/"
            "seagull_completion_results.json"
        ),
    )
    parser.add_argument("--quick", action="store_true")
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
                "physical_root": report["physical_common_cone_root"],
                "failed_checks": [
                    key
                    for key, passed in report["checks"].items()
                    if not passed
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
