#!/usr/bin/env python3
"""Stage 3D: derive all low-band plaquette operators from the shared swap graph.

Stage 3A derived the fourth-order ring exchange but left the diagonal
flippability coefficient open. This script starts from the actual four-link
spin-1 Hilbert space, eliminates the Gauss-violating states exactly, and
matches the resulting three-state gauge-band Hamiltonian to

    2 U E_loop^2 - K (W + W^dagger) + v D

plus the first omitted W^2 operator.

A linked-cluster subtraction removes the independent one-link self-energy from
the plaquette-local diagonal term. The independent self-energy is retained in
the predicted renormalized electric coefficient. The script then tests whether
the complete matched coefficients, together with the Stage-3C r=1/144
charging mode, admit a photon/gravity common cone without a counterterm.
"""
from __future__ import annotations

import argparse
import functools
import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from check_shared_junction_move import gravity_coupling
from finite_em_dynamics_impl import axial_spectrum, fit_branch


R_COMMON = 1.0 / 144.0


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


def spin_one_operators() -> tuple[np.ndarray, np.ndarray]:
    electric = np.diag((-1.0, 0.0, 1.0))
    raising = np.zeros((3, 3), dtype=float)
    raising[1, 0] = 1.0
    raising[2, 1] = 1.0
    return electric, raising


def kron_all(operators: list[np.ndarray]) -> np.ndarray:
    output = operators[0]
    for operator in operators[1:]:
        output = np.kron(output, operator)
    return output


def embed(operator: np.ndarray, location: int, count: int = 4) -> np.ndarray:
    identity = np.eye(operator.shape[0])
    return kron_all(
        [operator if index == location else identity for index in range(count)]
    )


E_LINK, U_LINK = spin_one_operators()
E_LINKS = [embed(E_LINK, index) for index in range(4)]
U_LINKS = [embed(U_LINK, index) for index in range(4)]
GAUSS = [
    E_LINKS[0] + E_LINKS[3],
    E_LINKS[1] - E_LINKS[0],
    -E_LINKS[1] - E_LINKS[2],
    E_LINKS[2] - E_LINKS[3],
]
H0 = 0.5 * sum(generator @ generator for generator in GAUSS)
V_SWAP = -sum((operator + operator.T for operator in U_LINKS))

BASIS_STATES = list(itertools.product((-1, 0, 1), repeat=4))
BASIS_INDEX = {state: index for index, state in enumerate(BASIS_STATES)}
GAUGE_INDICES = [
    BASIS_INDEX[(-1, -1, 1, 1)],
    BASIS_INDEX[(0, 0, 0, 0)],
    BASIS_INDEX[(1, 1, -1, -1)],
]
GAUGE_EMBEDDING = np.zeros((81, 3), dtype=float)
for column, index in enumerate(GAUGE_INDICES):
    GAUGE_EMBEDDING[index, column] = 1.0


def inverse_square_root(matrix: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    return (eigenvectors * (1.0 / np.sqrt(eigenvalues))) @ eigenvectors.T


def effective_gauge_hamiltonian(x: float) -> tuple[np.ndarray, np.ndarray]:
    """Canonical Hermitian low-band Hamiltonian in the three gauge states."""
    energies, eigenvectors = np.linalg.eigh(H0 + x * V_SWAP)
    low_energies = energies[:3]
    overlap = GAUGE_EMBEDDING.T @ eigenvectors[:, :3]
    polar = overlap @ inverse_square_root(overlap.T @ overlap)
    effective = polar @ np.diag(low_energies) @ polar.T
    return effective, energies


def isolated_link_energy(x: float, reference_flux: int) -> float:
    """Exact adiabatic one-link self-energy for fixed endpoint Gauss charges."""
    fluxes = np.asarray((-1.0, 0.0, 1.0))
    hamiltonian = np.diag((fluxes - reference_flux) ** 2)
    hamiltonian -= x * (U_LINK + U_LINK.T)
    energies, eigenvectors = np.linalg.eigh(hamiltonian)
    reference_index = reference_flux + 1
    branch = int(np.argmax(np.abs(eigenvectors[reference_index, :]) ** 2))
    return float(energies[branch])


@functools.lru_cache(maxsize=None)
def matched_coefficients(x_rounded: float) -> dict[str, float]:
    x = float(x_rounded)
    effective, full_energies = effective_gauge_hamiltonian(x)
    edge_diagonal = 0.5 * (effective[0, 0] + effective[2, 2])
    center_diagonal = effective[1, 1]
    full_contrast = edge_diagonal - center_diagonal

    edge_link = isolated_link_energy(x, 1)
    center_link = isolated_link_energy(x, 0)
    independent_link_contrast = 4.0 * (edge_link - center_link)

    # In the three-state gauge band, the electric term gives edge-center 2U,
    # while D=diag(1,2,1) gives edge-center -v.
    delta_u = 0.5 * independent_link_contrast
    linked_diagonal_contrast = full_contrast - independent_link_contrast
    v_diagonal = -linked_diagonal_contrast

    ring_exchange = -0.5 * (effective[0, 1] + effective[1, 2])
    double_flux_exchange = effective[0, 2]
    low_band_separation = full_energies[3] - full_energies[2]

    return {
        "x_t_over_Delta": x,
        "K_plaquette_over_Delta": float(ring_exchange),
        "delta_U_over_Delta": float(delta_u),
        "v_over_Delta": float(v_diagonal),
        "v_over_K": float(v_diagonal / ring_exchange),
        "double_flux_exchange_over_Delta": float(double_flux_exchange),
        "double_flux_exchange_over_K": float(
            double_flux_exchange / ring_exchange
        ),
        "full_edge_center_contrast": float(full_contrast),
        "independent_link_contrast": float(independent_link_contrast),
        "linked_diagonal_contrast": float(linked_diagonal_contrast),
        "low_band_separation_over_K": float(
            low_band_separation / ring_exchange
        ),
    }


def logarithmic_fit(xs: np.ndarray, values: np.ndarray) -> dict[str, float]:
    power, log_prefactor = np.polyfit(np.log(xs), np.log(np.abs(values)), 1)
    return {"power": float(power), "prefactor": float(np.exp(log_prefactor))}


def photon_result(
    x: float,
    include_electric_self_energy: bool,
    mirror_self_energy_into_gravity: bool = False,
) -> dict[str, Any]:
    matched = matched_coefficients(round(float(x), 12))
    k_plaquette = matched["K_plaquette_over_Delta"]
    u_photon = R_COMMON
    if include_electric_self_energy:
        u_photon += matched["delta_U_over_Delta"]
    u_over_t = u_photon / k_plaquette
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
    photon_speed = k_plaquette * fit["linear_speed_c_gamma"]

    u_gravity = R_COMMON
    if mirror_self_energy_into_gravity:
        u_gravity += matched["delta_U_over_Delta"]
    gravity_speed = math.sqrt(u_gravity * gravity_coupling(x)[0])

    return {
        **matched,
        "u_photon_over_Delta": u_photon,
        "u_over_t_plaquette": u_over_t,
        "photon_fit": fit,
        "photon_speed_in_Delta_units": photon_speed,
        "u_gravity_over_Delta": u_gravity,
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


def bisection_root(
    lower: float,
    upper: float,
    include_electric_self_energy: bool,
    mirror_self_energy_into_gravity: bool,
    iterations: int,
) -> dict[str, Any]:
    lower_result = photon_result(
        lower, include_electric_self_energy, mirror_self_energy_into_gravity
    )
    upper_result = photon_result(
        upper, include_electric_self_energy, mirror_self_energy_into_gravity
    )
    if lower_result["speed_difference"] * upper_result["speed_difference"] >= 0:
        raise RuntimeError("requested root is not bracketed")
    for _ in range(iterations):
        midpoint = 0.5 * (lower + upper)
        midpoint_result = photon_result(
            midpoint,
            include_electric_self_energy,
            mirror_self_energy_into_gravity,
        )
        if (
            lower_result["speed_difference"]
            * midpoint_result["speed_difference"]
            <= 0
        ):
            upper = midpoint
            upper_result = midpoint_result
        else:
            lower = midpoint
            lower_result = midpoint_result
    result = photon_result(
        0.5 * (lower + upper),
        include_electric_self_energy,
        mirror_self_energy_into_gravity,
    )
    result["root_bracket_width"] = upper - lower
    return result


def compact_result(result: dict[str, Any]) -> dict[str, Any]:
    return {
        key: result[key]
        for key in (
            "x_t_over_Delta",
            "K_plaquette_over_Delta",
            "delta_U_over_Delta",
            "v_over_K",
            "double_flux_exchange_over_K",
            "low_band_separation_over_K",
            "u_photon_over_Delta",
            "u_over_t_plaquette",
            "u_gravity_over_Delta",
            "photon_speed_in_Delta_units",
            "gravity_speed_in_Delta_units",
            "speed_difference",
            "speed_ratio",
        )
    } | {
        "photon_gap_power": result["photon_fit"][
            "gap_proportional_to_lattice_momentum_power"
        ],
        "photon_fit_residual": result["photon_fit"][
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
    small_fits = {
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
    }

    scan_x = (
        np.asarray((0.05, 0.10, 0.14, 0.20, 0.30, 0.40))
        if quick
        else np.linspace(0.05, 0.40, 36)
    )
    scan = [
        photon_result(float(value), True, False) for value in scan_x
    ]
    z1_scan = [
        result
        for result in scan
        if 0.94
        <= result["photon_fit"]["gap_proportional_to_lattice_momentum_power"]
        <= 1.16
        and result["photon_fit"][
            "maximum_relative_linear_plus_cubic_fit_residual"
        ]
        < 0.002
    ]

    counterterm_root = bisection_root(
        0.13,
        0.15,
        include_electric_self_energy=False,
        mirror_self_energy_into_gravity=False,
        iterations=12 if quick else 28,
    )
    mirrored_root = bisection_root(
        0.18,
        0.22,
        include_electric_self_energy=True,
        mirror_self_energy_into_gravity=True,
        iterations=12 if quick else 28,
    )

    minimum_scan = min(scan, key=lambda result: result["speed_difference"])
    minimum_z1 = min(z1_scan, key=lambda result: result["speed_difference"])
    x_04 = matched_coefficients(0.4)

    checks = {
        "ring_exchange_is_fourth_order": abs(
            small_fits["K_plaquette"]["power"] - 4.0
        )
        < 0.03,
        "electric_self_energy_is_second_order": abs(
            small_fits["delta_U"]["power"] - 2.0
        )
        < 0.02,
        "linked_diagonal_is_fourth_order": abs(
            small_fits["v_diagonal"]["power"] - 4.0
        )
        < 0.10,
        "small_x_v_over_K_approaches_minus_two_over_fifteen": abs(
            small_rows[0]["v_over_K"] + 2.0 / 15.0
        )
        < 2.0e-4,
        "counterterm_free_scan_has_no_common_cone": minimum_scan[
            "speed_difference"
        ]
        > 0.0,
        "z1_subscan_has_no_common_cone": minimum_z1["speed_difference"] > 0.0,
        "omitting_predicted_self_energy_creates_spurious_z1_root": (
            0.94
            <= counterterm_root["photon_fit"][
                "gap_proportional_to_lattice_momentum_power"
            ]
            <= 1.16
            and abs(counterterm_root["speed_difference"]) < 1.0e-6
        ),
        "mirroring_self_energy_into_gravity_still_fails_z1_gate": (
            abs(mirrored_root["speed_difference"]) < 1.0e-6
            and mirrored_root["photon_fit"][
                "gap_proportional_to_lattice_momentum_power"
            ]
            < 0.94
            and mirrored_root["photon_fit"][
                "maximum_relative_linear_plus_cubic_fit_residual"
            ]
            > 0.002
        ),
        "additional_double_flux_operator_becomes_material": (
            x_04["double_flux_exchange_over_K"] > 0.20
        ),
        "low_band_still_separated_at_scan_edge": (
            x_04["low_band_separation_over_K"] > 3.0
        ),
    }

    report = {
        "module": (
            "phase_junction_network/microscopic/"
            "check_plaquette_diagonal_matching.py"
        ),
        "status": (
            "Stage 3D PASS as a rejection gate: exact diagonal matching "
            "rejects the current shared-move graph without a derived "
            "self-energy cancellation or a more complete microscopic frame "
            "charging process"
        ),
        "run_mode": "quick" if quick else "full",
        "matching_scheme": {
            "low_band": (
                "Canonical polar projection of the three lowest eigenstates "
                "onto the exact spin-1 gauge states m=-1,0,+1."
            ),
            "electric_matching": (
                "The exact independent one-link self-energy is assigned to "
                "the extensive electric coefficient delta_U."
            ),
            "plaquette_diagonal_matching": (
                "The linked remainder of the edge-center contrast is matched "
                "to D=diag(1,2,1), so v=-(linked contrast)."
            ),
        },
        "derived_common_charging_factor": R_COMMON,
        "small_x_fits": small_fits,
        "small_x_rows": small_rows,
        "leading_behavior": {
            "K_plaquette": "20*x^4",
            "delta_U": "2*x^2",
            "v_diagonal": "-(8/3)*x^4",
            "v_over_K": "-2/15",
        },
        "counterterm_free_scan": {
            "x_min": float(scan_x[0]),
            "x_max": float(scan_x[-1]),
            "points": len(scan),
            "minimum_speed_difference": minimum_scan["speed_difference"],
            "minimum_speed_difference_x": minimum_scan["x_t_over_Delta"],
            "z1_like_point_count": len(z1_scan),
            "minimum_z1_speed_difference": minimum_z1["speed_difference"],
            "rows": [compact_result(result) for result in scan],
        },
        "unearned_counterterm_control": {
            "description": (
                "Dropping the predicted O(x^2) electric self-energy creates "
                "an apparently successful common cone. This root is a negative "
                "control because no cancellation mechanism has been derived."
            ),
            "root": compact_result(counterterm_root),
            "omitted_delta_U_over_Delta": counterterm_root[
                "delta_U_over_Delta"
            ],
        },
        "mirrored_selfenergy_control": {
            "description": (
                "Copying the photon self-energy into U_g can match speeds, but "
                "the photon branch fails the declared z~1 and fit-residual "
                "gates. Matching one number is not enough."
            ),
            "root": compact_result(mirrored_root),
        },
        "generated_operator_inventory": {
            "double_flux_exchange": (
                "The exact low band also generates a direct m=-1 to m=+1 "
                "operator, equivalent to W^2+W^dagger^2."
            ),
            "double_flux_over_K_at_x_0p4": x_04[
                "double_flux_exchange_over_K"
            ],
            "low_band_gap_over_K_at_x_0p4": x_04[
                "low_band_separation_over_K"
            ],
        },
        "checks": checks,
        "stage_pass": all(checks.values()),
        "decision": (
            "Reject the present pure-swap shared graph as the complete "
            "microscopic unification. It generates an O(t^2/Delta) photon "
            "charging self-energy that dominates the O(t^4/Delta^3) loop "
            "term and prevents a z~1 common cone over the controlled scan."
        ),
        "claim_boundary": {
            "established": [
                "the exact spin-1 low-band ring, electric, diagonal, and double-flux operators",
                "the O(x^2) independent-link electric renormalization",
                "the O(x^4) linked diagonal coefficient and its negative sign",
                "absence of a counterterm-free common cone in the controlled x=0.05..0.40 scan",
                "the apparent common-cone root obtained by deleting delta_U is spurious",
                "copying delta_U into gravity does not pass the photon z~1 gate",
            ],
            "not_established": [
                "a symmetry that cancels the electric self-energy",
                "the complete frame analogue of the diagonal matching",
                "that no more elaborate shared mediator can work",
                "the infinite-volume three-dimensional phase outside the scan",
                "physical alpha, G, nonlinear gravity, or formal novelty",
            ],
        },
        "next_gate": (
            "Enumerate symmetry-complete elementary junction Hamiltonians that "
            "include the diagonal/seagull partners of each swap. A viable "
            "candidate must cancel or share the O(x^2) charging correction by "
            "symmetry, derive the frame analogue, and reproduce all generated "
            "operators without sector-specific counterterms."
        ),
    }

    if not report["stage_pass"]:
        raise AssertionError(json.dumps(ready(checks), indent=2, sort_keys=True))
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "phase_junction_network/microscopic/"
            "plaquette_diagonal_matching_results.json"
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
                "small_x_fits": report["small_x_fits"],
                "minimum_counterterm_free_speed_difference": report[
                    "counterterm_free_scan"
                ]["minimum_speed_difference"],
                "unearned_counterterm_root": report[
                    "unearned_counterterm_control"
                ]["root"],
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
